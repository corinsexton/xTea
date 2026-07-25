// xtea-scan (prototype) — single-pass clip + discordant read collector (Stage 1 of xTea2).
//
// Streaming, bounded-memory design (spec §3): one pass over BAM/CRAM detects clip + disc signal per
// the rules in docs/clipdisc_contract_spec.md. Parallelism = fixed-size bp TILES over a thread pool
// (NOT per-chromosome); a read is owned by the tile containing its breakpoint (clip) / POS (disc).
//
// MEMORY: we aggregate into small per-position counters as we scan — peak RAM is O(distinct candidate
// sites), not O(reads). Clipped sequences are streamed to ONE consolidated FASTQ file (on disk, not
// RAM; not thousands of files). The per-read dump (--out, the oracle gate) is opt-in and region-scale.
//
// Outputs (all optional):
//   --out FILE|-        per-read TSV (opt-in; accumulates — region/testing only):
//                         CLIP chrom clip_pos side read_name is_rev clip_len
//                         DISC chrom pos read_name mate_chrom mate_pos is_rev is_mate_rev
//   --clip-sites FILE   aggregated clip positions: chrom pos n_left n_right          (sorted)
//   --disc-sites FILE   aggregated disc counts:    chrom pos n_disc                  (sorted)
//   --clip-fastq FILE   clipped soft-clip parts, streamed (@name~chrom~pos~side)
//   --disc-pos  FILE    disc records for mate-gathering: chrom pos is_rev mate_chrom mate_pos
//                         read_name is_mate_rev                                       (streamed)
//   --realign-out FILE  run bwa on the clip FASTQ vs --l1/alu/sva-cns; per-site/family L/R counts
// Positions 0-based; region coords 0-based half-open.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cinttypes>
#include <string>
#include <string_view>
#include <sstream>
#include <memory>
#include <optional>
#include <vector>
#include <map>
#include <array>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <getopt.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/regidx.h"
#include "realign.hpp"
#include "cns_emit.hpp"

namespace {

constexpr int  kMinClipMapqDefault   = 12;
constexpr int  kMinDiscMapqDefault   = 20;
constexpr int  kClipPhredCutoff      = 15;
constexpr long kDiscThresholdDefault = 2000;
constexpr hts_pos_t kTileSizeDefault = 20'000'000;
constexpr size_t kFlushBytes = 1u << 18;  // 256 KB per-tile buffer flush

struct SamFileDeleter { void operator()(samFile* p) const { if (p) sam_close(p); } };
struct HdrDeleter     { void operator()(sam_hdr_t* p) const { if (p) sam_hdr_destroy(p); } };
struct BamDeleter     { void operator()(bam1_t* p) const { if (p) bam_destroy1(p); } };
struct IdxDeleter     { void operator()(hts_idx_t* p) const { if (p) hts_idx_destroy(p); } };
struct ItrDeleter     { void operator()(hts_itr_t* p) const { if (p) hts_itr_destroy(p); } };
using SamFile = std::unique_ptr<samFile, SamFileDeleter>;
using Hdr = std::unique_ptr<sam_hdr_t, HdrDeleter>;
using Bam = std::unique_ptr<bam1_t, BamDeleter>;
using Idx = std::unique_ptr<hts_idx_t, IdxDeleter>;
using Itr = std::unique_ptr<hts_itr_t, ItrDeleter>;

struct Options {
    // Multi-file input: --bam takes ONE argument value that may contain several
    // whitespace-separated paths ("a.cram b.cram"), pooling multiple libraries of
    // the same sample. Functionally identical for one path or many. All files must
    // share an identical @SQ dictionary (validate_headers, called from main before
    // any tile is built) -- that precondition is what lets every (tid,pos)-keyed
    // aggregation below stay unchanged, with no per-file tid remapping.
    std::vector<std::string> bams;
    std::string reference, region;
    std::string out, clip_sites, disc_sites, clip_fastq, disc_pos;
    std::string realign_out, bwa{"bwa"}, l1_cns, alu_cns, sva_cns;
    // per-family repeat-annotation BEDs (extended, from .out) for the col5
    // mate-in-repeat clip count; indexed L1=0, Alu=1, SVA=2 (matches fams order).
    std::string l1_bed, alu_bed, sva_bed;
    std::string disc_pos_in, mate_fastq, realign_fastq_in;  // Part 2: mate-gather + realign-any-fastq
    std::string clip_realign_sam;  // move B.1: write per-family clip-realign SAM to <prefix><FAM>.sam
    // move B.1 cns-emit: reproduce cns_filter's temp_clip.sam/temp_disc.sam from the scan pass.
    std::string cns_clip_sam, cns_disc_sam, cns_candidates, cns_ref, cns_tmp_prefix, cns_lowq;
    std::string cns_manifest;  // fold-in: collect once, emit for all types listed in this manifest
    // full fold-in: spill cns collection DURING the detection pass, then fan from the spill later.
    std::string cns_collect_clip, cns_collect_anchor, cns_collect_lowq;   // detection-pass spill outputs
    std::string cns_fan_manifest, cns_fan_clip, cns_fan_anchor, cns_fan_lowq;  // post-candidate fan inputs
    long cns_extnd = 400;
    int cns_sample = 0;
    int min_clip_mapq = kMinClipMapqDefault, min_disc_mapq = kMinDiscMapqDefault;
    long disc_threshold = kDiscThresholdDefault;
    int threads = 1, bwa_threads = 4;
    hts_pos_t tile_size = kTileSizeDefault;
    bool no_baseq = false;
    bool clip_supplementary = true;  // count clips on supplementary (split-read) alignments
};

struct Region { std::string chrom; hts_pos_t start, end; };
struct Tile   { int tid; hts_pos_t start, end; };

// Small per-read record for the opt-in --out dump (NO sequence — bounded even in dump mode).
struct ClipRec { int tid; hts_pos_t pos; char side; int is_rev; int clen; std::string qname; };
struct DiscRec { int tid; hts_pos_t pos; int mtid; hts_pos_t mpos; int is_rev; int is_mrev; std::string qname; };

struct TileResult {
    std::map<std::pair<int, hts_pos_t>, std::pair<int, int>> clip_agg;  // ->(nL,nR)
    std::map<std::pair<int, hts_pos_t>, int> disc_agg;
    // col5: per clip site (tid,cp), # clipped reads whose mate lands in ANY of the
    // MEI repeat annotations (L1 | Alu | SVA), a single union count. Only when rep
    // BEDs are given. (xTea2 design change DC1 vs master's Alu-only col5 — see
    // docs/clipdisc_contract_spec.md §4 Q2/DC1.)
    std::map<std::pair<int, hts_pos_t>, int> clip_matein;
    std::vector<ClipRec> clip_recs;   // only if --out
    std::vector<DiscRec> disc_recs;   // only if --out
    std::string fastq_buf;            // only if streaming clip FASTQ
    std::string discpos_buf;          // only if streaming disc-pos
    // full fold-in: cns collection spilled during the detection pass (candidate-agnostic).
    std::string cns_clip_buf, cns_anchor_buf, cns_lowq_buf;
};

// Thread-safe append-only file (one consolidated file, many writers).
// Thread-safe append-only sink. NOTE the ordering consequence, measured 2026-07-25 and
// deliberately accepted: appends land in tile-completion order, so every file written through
// this class (the cns spill files in particular) has a run-to-run NONDETERMINISTIC line order.
// Content is stable -- two runs of the same binary produce spill files that are identical after
// sorting -- but they are never byte-identical, so do NOT use a spill file as a byte-comparison
// baseline. Verified that this order does not reach the VCF: two runs whose spill order differed
// produced byte-identical per-type VCFs. See HANDOFF_2026-07-25.md "cns spill is not byte-stable".
class ConcFile {
    FILE* f_ = nullptr; std::mutex m_;
public:
    bool open(const std::string& p) { f_ = std::fopen(p.c_str(), "w"); return f_ != nullptr; }
    void append(const std::string& s) { if (s.empty()) return; std::lock_guard<std::mutex> lk(m_); std::fwrite(s.data(), 1, s.size(), f_); }
    void close() { if (f_) { std::fclose(f_); f_ = nullptr; } }
    ~ConcFile() { close(); }
};

bool is_qualified_clip(const uint8_t* qual, int n) {
    if (n <= 0) return false;
    const double half = n / 2.0; int low = 0;
    for (int i = 0; i < n; ++i) if (qual[i] < kClipPhredCutoff && ++low > half) return false;
    return true;
}

std::optional<Region> parse_region(std::string_view s) {
    auto colon = s.rfind(':');
    auto dash = s.find('-', colon == std::string_view::npos ? 0 : colon + 1);
    if (colon == std::string_view::npos || dash == std::string_view::npos) return std::nullopt;
    Region r; r.chrom = std::string(s.substr(0, colon));
    r.start = (hts_pos_t)std::strtoll(std::string(s.substr(colon + 1, dash - colon - 1)).c_str(), nullptr, 10);
    r.end = (hts_pos_t)std::strtoll(std::string(s.substr(dash + 1)).c_str(), nullptr, 10);
    if (r.chrom.empty() || r.end <= r.start) return std::nullopt;
    return r;
}

std::vector<Tile> build_tiles(sam_hdr_t* hdr, const std::optional<Region>& reg, hts_pos_t tsz) {
    std::vector<Tile> tiles;
    auto add = [&](int tid, hts_pos_t lo, hts_pos_t hi) {
        for (hts_pos_t s = lo; s < hi; s += tsz) tiles.push_back({tid, s, std::min(hi, s + tsz)});
    };
    if (reg) { int tid = sam_hdr_name2tid(hdr, reg->chrom.c_str()); if (tid >= 0) add(tid, reg->start, reg->end); }
    else { int n = sam_hdr_nref(hdr); for (int tid = 0; tid < n; ++tid) add(tid, 0, (hts_pos_t)sam_hdr_tid2len(hdr, tid)); }
    return tiles;
}

void append_clip_fastq(std::string& buf, const bam1_t* b, bool left, int clen,
                       const char* chrom, hts_pos_t pos, char side, int is_rev) {
    const uint8_t* s = bam_get_seq(b); const uint8_t* q = bam_get_qual(b);
    const int qlen = b->core.l_qseq; const int off = left ? 0 : (qlen - clen);
    const bool have_q = (q && q[0] != 0xff);
    buf += '@'; buf += bam_get_qname(b); buf += '~'; buf += chrom; buf += '~';
    buf += std::to_string(pos); buf += '~'; buf += side; buf += '\n';
    for (int i = 0; i < clen; ++i) buf += seq_nt16_str[bam_seqi(s, off + i)];
    buf += "\n+\n";
    for (int i = 0; i < clen; ++i) buf += have_q ? (char)(q[off + i] + 33) : 'I';
    buf += '\n';
    (void)is_rev;
}

// file_idx identifies which --bam input this call is scanning. It is used ONLY to tag
// the cns anchor spill (see collect_cns_read): the spill is later read back by
// emit_cns_from_spill, which has no BAM pass of its own to learn file identity from,
// so without the tag it could resolve a mate against the wrong file on a qname
// collision. Everything else here is (tid,pos)-keyed and file-agnostic by construction.
void scan_tile(const Tile& tile, const Options& o, samFile* fp, hts_idx_t* idx, sam_hdr_t* hdr,
               TileResult& res, ConcFile* fastq, ConcFile* discposf,
               regidx_t* const rep_idx[3], ConcFile* cns_clip_f, ConcFile* cns_anchor_f,
               ConcFile* cns_lowq_f, int file_idx) {
    const bool want_fastq = (fastq != nullptr);
    const bool want_out = !o.out.empty();
    const bool want_cns = (cns_clip_f != nullptr);  // fold-in: spill cns collection this pass
    // per-thread regidx iterators (regidx_t is shared read-only; each iterator
    // holds per-query scratch state, so it must be thread-local to this tile).
    const bool rep_any = rep_idx && (rep_idx[0] || rep_idx[1] || rep_idx[2]);
    regitr_t* rep_itr[3] = {nullptr, nullptr, nullptr};
    if (rep_any) for (int f = 0; f < 3; ++f) if (rep_idx[f]) rep_itr[f] = regitr_init(rep_idx[f]);
    const bool want_discpos = (discposf != nullptr);
    Itr it(sam_itr_queryi(idx, tile.tid, tile.start, tile.end));
    if (!it) return;
    Bam bh(bam_init1());
    while (sam_itr_next(fp, it.get(), bh.get()) >= 0) {
        const bam1_t* r = bh.get();
        const uint16_t flag = r->core.flag;
        // Full fold-in: spill the cns collection for this read (owned by ref-start so each read is
        // seen once). Runs BEFORE the detection skip because cns has its own gates (its lowq counts
        // supp/sec; its clip write skips them internally). Candidate-agnostic; fanned later.
        if (want_cns && r->core.pos >= tile.start && r->core.pos < tile.end) {
            xtea::cns::collect_cns_read(r, hdr, res.cns_clip_buf, res.cns_anchor_buf, res.cns_lowq_buf,
                                        file_idx);
            if (res.cns_clip_buf.size()   >= kFlushBytes) { cns_clip_f->append(res.cns_clip_buf);     res.cns_clip_buf.clear(); }
            if (res.cns_anchor_buf.size() >= kFlushBytes) { cns_anchor_f->append(res.cns_anchor_buf); res.cns_anchor_buf.clear(); }
            if (res.cns_lowq_buf.size()   >= kFlushBytes) { cns_lowq_f->append(res.cns_lowq_buf);     res.cns_lowq_buf.clear(); }
        }
        // Split reads spanning an MEI breakpoint carry the clip signal on their
        // supplementary alignment (common for SVA/L1); master's clip collection counts
        // them (verified: golden nL=11 = 4 primary + 7 non-dup supplementary at
        // chr1:168451690). Kept by default; toggle off with --no-clip-supplementary.
        // DISC always excludes supplementary (below), independent of this toggle.
        uint16_t skip = BAM_FDUP | BAM_FSECONDARY | BAM_FUNMAP;
        if (!o.clip_supplementary) skip |= BAM_FSUPPLEMENTARY;
        if (flag & skip) continue;
        const hts_pos_t pos = r->core.pos;
        const int is_rev = bam_is_rev(r) ? 1 : 0;
        const uint32_t* cig = bam_get_cigar(r);
        const int ncig = r->core.n_cigar;
        const char* chrom = sam_hdr_tid2name(hdr, r->core.tid);

        // ---- CLIP ----
        if (r->core.qual >= o.min_clip_mapq && ncig > 0) {
            const uint8_t* qual = bam_get_qual(r);
            const bool have_qual = (qual && qual[0] != 0xff && !o.no_baseq);
            // col5 (xTea2 design change DC1): does this read's mate land in ANY MEI
            // repeat (L1 | Alu | SVA)? Master counted only mates in the Alu annotation
            // -- an accident of the shared-pub_clip race (spec §4 Q2); xTea2 uses the
            // union (each family's own boundary extend baked into its BED) as a
            // deterministic, more accurate repeat-anchored support signal. Computed
            // once per clipped read (lazy), reused for its L and R clip ends. The raw
            // 0-based mate pos vs raw .out coords preserve master's off-by-one.
            int mate_rep = -1;  // -1 = not computed yet; 0/1 = mate-in-any-rep result
            auto mate_in_rep = [&]() -> int {
                if (mate_rep >= 0) return mate_rep;
                mate_rep = 0;
                if (rep_any && !(flag & BAM_FMUNMAP) && r->core.mtid >= 0) {
                    const char* mchrom = sam_hdr_tid2name(hdr, r->core.mtid);
                    const hts_pos_t mp = r->core.mpos;
                    for (int f = 0; f < 3; ++f)
                        if (rep_idx[f] && regidx_overlap(rep_idx[f], mchrom, mp, mp, rep_itr[f])) {
                            mate_rep = 1; break;
                        }
                }
                return mate_rep;
            };
            auto add_matein = [&](int tid, hts_pos_t cp) {
                if (!rep_any) return;
                res.clip_matein[{tid, cp}] += mate_in_rep();
            };
            const int op0 = bam_cigar_op(cig[0]);
            if (op0 == BAM_CSOFT_CLIP || op0 == BAM_CHARD_CLIP) {
                const int clen = bam_cigar_oplen(cig[0]); const hts_pos_t cp = pos;
                bool ok = !(op0 == BAM_CSOFT_CLIP && have_qual) || is_qualified_clip(qual, clen);
                if (ok && cp >= tile.start && cp < tile.end) {
                    res.clip_agg[{r->core.tid, cp}].first += 1;
                    add_matein(r->core.tid, cp);
                    if (want_out) res.clip_recs.push_back({r->core.tid, cp, 'L', is_rev, clen, bam_get_qname(r)});
                    if (want_fastq && op0 == BAM_CSOFT_CLIP) {
                        append_clip_fastq(res.fastq_buf, r, true, clen, chrom, cp, 'L', is_rev);
                        if (res.fastq_buf.size() >= kFlushBytes) { fastq->append(res.fastq_buf); res.fastq_buf.clear(); }
                    }
                }
            }
            const int opN = bam_cigar_op(cig[ncig - 1]);
            if (opN == BAM_CSOFT_CLIP || opN == BAM_CHARD_CLIP) {
                const int clen = bam_cigar_oplen(cig[ncig - 1]); hts_pos_t cp = pos;
                for (int i = 0; i < ncig - 1; ++i) {
                    const int op = bam_cigar_op(cig[i]);
                    if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CINS) continue;
                    cp += bam_cigar_oplen(cig[i]);
                }
                bool ok = true;
                if (opN == BAM_CSOFT_CLIP && have_qual) ok = is_qualified_clip(qual + (r->core.l_qseq - clen), clen);
                if (ok && cp >= tile.start && cp < tile.end) {
                    res.clip_agg[{r->core.tid, cp}].second += 1;
                    add_matein(r->core.tid, cp);
                    if (want_out) res.clip_recs.push_back({r->core.tid, cp, 'R', is_rev, clen, bam_get_qname(r)});
                    if (want_fastq && opN == BAM_CSOFT_CLIP) {
                        append_clip_fastq(res.fastq_buf, r, false, clen, chrom, cp, 'R', is_rev);
                        if (res.fastq_buf.size() >= kFlushBytes) { fastq->append(res.fastq_buf); res.fastq_buf.clear(); }
                    }
                }
            }
        }
        // ---- DISC ---- (exclude supplementary: split-read alignments carry no
        // meaningful insert-size/mate discordance signal; master skips them here too)
        if (!(flag & (BAM_FMUNMAP | BAM_FSUPPLEMENTARY)) && r->core.qual >= o.min_disc_mapq) {
            const int mtid = r->core.mtid;
            if (mtid >= 0 && pos >= tile.start && pos < tile.end) {
                const hts_pos_t mpos = r->core.mpos;
                if (mtid != r->core.tid || llabs((long long)(mpos - pos)) > o.disc_threshold) {
                    const int is_mrev = (flag & BAM_FMREVERSE) ? 1 : 0;
                    res.disc_agg[{r->core.tid, pos}] += 1;
                    if (want_out) res.disc_recs.push_back({r->core.tid, pos, mtid, mpos, is_rev, is_mrev, bam_get_qname(r)});
                    if (want_discpos) {
                        // Trailing field = read reference END (bam_endpos, i.e. one past
                        // the last aligned ref base = pysam reference_end). The disc
                        // adapter needs it to reproduce master's bamfile.fetch overlap
                        // semantics exactly (a read is "in" a window [s,e) iff
                        // ref_start < e AND ref_end > s); a read-start-only test
                        // mis-counts breakpoint-spanning reads by +-1. Appended last so
                        // the mate-gather --disc-pos-in parser (fields 0-6) is unaffected.
                        res.discpos_buf += chrom; res.discpos_buf += '\t';
                        res.discpos_buf += std::to_string(pos); res.discpos_buf += '\t';
                        res.discpos_buf += std::to_string(is_rev); res.discpos_buf += '\t';
                        res.discpos_buf += sam_hdr_tid2name(hdr, mtid); res.discpos_buf += '\t';
                        res.discpos_buf += std::to_string(mpos); res.discpos_buf += '\t';
                        res.discpos_buf += bam_get_qname(r); res.discpos_buf += '\t';
                        res.discpos_buf += std::to_string(is_mrev); res.discpos_buf += '\t';
                        res.discpos_buf += std::to_string((long long)bam_endpos(r)); res.discpos_buf += '\n';
                        if (res.discpos_buf.size() >= kFlushBytes) { discposf->append(res.discpos_buf); res.discpos_buf.clear(); }
                    }
                }
            }
        }
    }
    if (want_fastq && !res.fastq_buf.empty()) { fastq->append(res.fastq_buf); res.fastq_buf.clear(); }
    if (want_discpos && !res.discpos_buf.empty()) { discposf->append(res.discpos_buf); res.discpos_buf.clear(); }
    if (want_cns) {
        if (!res.cns_clip_buf.empty())   { cns_clip_f->append(res.cns_clip_buf);     res.cns_clip_buf.clear(); }
        if (!res.cns_anchor_buf.empty()) { cns_anchor_f->append(res.cns_anchor_buf); res.cns_anchor_buf.clear(); }
        if (!res.cns_lowq_buf.empty())   { cns_lowq_f->append(res.cns_lowq_buf);     res.cns_lowq_buf.clear(); }
    }
    for (int f = 0; f < 3; ++f) if (rep_itr[f]) regitr_destroy(rep_itr[f]);
}

// ---- Part 2: disc mate-gathering ----
// One anchor record needs its mate (same qname) fetched at (mate_tid, mate_pos). Keyed by that
// position so a tiled sweep can pull the mate sequences by name (no per-read random seeks) —
// mirrors master's extract_mate_reads_by_name. The mate's FASTQ header carries the ANCHOR site.
struct MateNeed { std::string qname, anchor_chrom; long anchor_pos; };

// Read the --disc-pos file into {(mate_tid,mate_pos): [MateNeed...]}.
std::map<std::pair<int, hts_pos_t>, std::vector<MateNeed>>
load_mate_needs(const std::string& path, sam_hdr_t* hdr) {
    std::map<std::pair<int, hts_pos_t>, std::vector<MateNeed>> needs;
    FILE* f = std::fopen(path.c_str(), "r");
    if (!f) return needs;
    char* line = nullptr; size_t cap = 0; ssize_t len;
    while ((len = getline(&line, &cap, f)) != -1) {
        // chrom pos is_rev mate_chrom mate_pos read_name is_mate_rev
        char *save = nullptr;
        char* f0 = strtok_r(line, "\t\n", &save);      // anchor chrom
        char* f1 = strtok_r(nullptr, "\t\n", &save);   // anchor pos
        strtok_r(nullptr, "\t\n", &save);              // is_rev
        char* f3 = strtok_r(nullptr, "\t\n", &save);   // mate chrom
        char* f4 = strtok_r(nullptr, "\t\n", &save);   // mate pos
        char* f5 = strtok_r(nullptr, "\t\n", &save);   // read name
        if (!f0 || !f1 || !f3 || !f4 || !f5) continue;
        int mtid = sam_hdr_name2tid(hdr, f3);
        if (mtid < 0) continue;
        needs[{mtid, (hts_pos_t)atoll(f4)}].push_back({f5, f0, atoll(f1)});
    }
    free(line);
    std::fclose(f);
    return needs;
}

// Sweep a tile; for reads whose (pos,qname) is a needed mate, emit the FULL read as FASTQ tagged
// with the anchor site: @qname~anchor_chrom~anchor_pos~D
void gather_tile(const Tile& tile, samFile* fp, hts_idx_t* idx,
                 const std::map<std::pair<int, hts_pos_t>, std::vector<MateNeed>>& needs,
                 ConcFile& mate_fastq) {
    Itr it(sam_itr_queryi(idx, tile.tid, tile.start, tile.end));
    if (!it) return;
    Bam bh(bam_init1());
    std::string buf;
    while (sam_itr_next(fp, it.get(), bh.get()) >= 0) {
        const bam1_t* r = bh.get();
        if (r->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;  // want the real mate
        const hts_pos_t pos = r->core.pos;
        if (pos < tile.start || pos >= tile.end) continue;                    // own by POS (once)
        auto nit = needs.find({r->core.tid, pos});
        if (nit == needs.end()) continue;
        const char* qn = bam_get_qname(r);
        for (const auto& need : nit->second) {
            if (need.qname != qn) continue;
            const uint8_t* s = bam_get_seq(r); const uint8_t* q = bam_get_qual(r);
            const int L = r->core.l_qseq; const bool have_q = (q && q[0] != 0xff);
            buf += '@'; buf += qn; buf += '~'; buf += need.anchor_chrom; buf += '~';
            buf += std::to_string(need.anchor_pos); buf += "~D\n";
            for (int i = 0; i < L; ++i) buf += seq_nt16_str[bam_seqi(s, i)];
            buf += "\n+\n";
            for (int i = 0; i < L; ++i) buf += have_q ? (char)(q[i] + 33) : 'I';
            buf += '\n';
        }
        if (buf.size() >= kFlushBytes) { mate_fastq.append(buf); buf.clear(); }
    }
    if (!buf.empty()) mate_fastq.append(buf);
}

bool open_one(const std::string& path, const std::string& reference,
              SamFile& fp, Idx& idx, Hdr& hdr) {
    fp.reset(sam_open(path.c_str(), "r"));
    if (!fp) return false;
    if (!reference.empty()) hts_set_fai_filename(fp.get(), reference.c_str());
    hdr.reset(sam_hdr_read(fp.get()));
    if (!hdr) return false;
    idx.reset(sam_index_load(fp.get(), path.c_str()));
    return idx != nullptr;
}

// Thin wrapper over the first input file, for the call sites that genuinely only
// need one header (tile construction, contig-name lookup) once headers are known
// to be identical across all inputs.
bool open_private(const Options& o, SamFile& fp, Idx& idx, Hdr& hdr) {
    return open_one(o.bams.at(0), o.reference, fp, idx, hdr);
}

// Split one --bam argument value on whitespace into N paths.
std::vector<std::string> split_ws(const char* s) {
    std::vector<std::string> out;
    std::istringstream is(s);
    std::string tok;
    while (is >> tok) out.push_back(tok);
    return out;
}

// Fail-fast precondition for multi-file input: every input must present the same
// @SQ dictionary -- same reference count, and same name+length at every tid, in the
// same order. This keeps all (tid,pos) aggregation keys valid across files without
// remapping. Real multi-center data aligned to one reference FASTA satisfies this
// naturally; when it doesn't, that's a data problem the user needs told about
// loudly, naming the offending file and the exact discrepancy.
bool validate_headers(const std::vector<std::string>& bams, const std::string& reference) {
    if (bams.size() < 2) return true;
    SamFile fp0; Idx idx0; Hdr hdr0;
    if (!open_one(bams[0], reference, fp0, idx0, hdr0)) {
        std::fprintf(stderr, "error: cannot open %s (need index; reference for CRAM)\n", bams[0].c_str());
        return false;
    }
    const int n0 = sam_hdr_nref(hdr0.get());
    for (size_t b = 1; b < bams.size(); ++b) {
        SamFile fp; Idx idx; Hdr hdr;
        if (!open_one(bams[b], reference, fp, idx, hdr)) {
            std::fprintf(stderr, "error: cannot open %s (need index; reference for CRAM)\n", bams[b].c_str());
            return false;
        }
        const int n = sam_hdr_nref(hdr.get());
        if (n != n0) {
            std::fprintf(stderr,
                "error: header mismatch: %s has %d reference sequences but %s has %d.\n"
                "       All --bam inputs must be aligned/decoded against the same reference.\n",
                bams[b].c_str(), n, bams[0].c_str(), n0);
            return false;
        }
        for (int t = 0; t < n0; ++t) {
            const char* nm0 = sam_hdr_tid2name(hdr0.get(), t);
            const char* nm  = sam_hdr_tid2name(hdr.get(), t);
            if (!nm0 || !nm || std::strcmp(nm0, nm) != 0) {
                std::fprintf(stderr,
                    "error: header mismatch at reference index %d: %s has '%s' but %s has '%s'.\n"
                    "       @SQ order must be identical across all --bam inputs.\n",
                    t, bams[b].c_str(), nm ? nm : "(null)", bams[0].c_str(), nm0 ? nm0 : "(null)");
                return false;
            }
            const hts_pos_t l0 = sam_hdr_tid2len(hdr0.get(), t);
            const hts_pos_t l  = sam_hdr_tid2len(hdr.get(), t);
            if (l0 != l) {
                std::fprintf(stderr,
                    "error: header mismatch for reference '%s': length %" PRIhts_pos " in %s but %"
                    PRIhts_pos " in %s.\n",
                    nm0, l, bams[b].c_str(), l0, bams[0].c_str());
                return false;
            }
        }
    }
    return true;
}

void usage(FILE* f) {
    std::fprintf(f,
      "xtea-scan (prototype) — streaming clip+disc collector (bp-tile parallel)\n"
      "usage: xtea_scan --bam \"FILE [FILE ...]\" [--reference FA] [--region chr:start-end]\n"
      "  --bam takes ONE value that may hold several whitespace-separated paths, pooling\n"
      "        libraries of the same sample; all must share an identical @SQ dictionary.\n"
      "  outputs: [--out FILE|-] [--clip-sites F] [--disc-sites F] [--clip-fastq F] [--disc-pos F]\n"
      "  realign: [--realign-out F --l1-cns FA --alu-cns FA --sva-cns FA [--bwa PATH] [--bwa-threads N]]\n"
      "  col5:    [--l1-bed F --alu-bed F --sva-bed F]  (mate-in-repeat clip count, 6th realign-out col)\n"
      "  tuning:  [--threads N] [--tile-size BP] [--min-clip-mapq N] [--min-disc-mapq N]\n"
      "           [--disc-threshold N] [--no-baseq] [--no-clip-supplementary]\n"
      "region coords 0-based half-open. clip-sites/disc-sites/out are deterministic across --threads.\n");
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    static struct option lopts[] = {
        {"bam",required_argument,0,'b'},{"reference",required_argument,0,'r'},{"region",required_argument,0,'R'},
        {"out",required_argument,0,'o'},{"clip-sites",required_argument,0,10},{"disc-sites",required_argument,0,11},
        {"clip-fastq",required_argument,0,12},{"disc-pos",required_argument,0,14},{"threads",required_argument,0,'t'},
        {"tile-size",required_argument,0,13},{"min-clip-mapq",required_argument,0,1},{"min-disc-mapq",required_argument,0,2},
        {"disc-threshold",required_argument,0,3},{"no-baseq",no_argument,0,4},
        {"realign-out",required_argument,0,20},{"bwa",required_argument,0,21},{"l1-cns",required_argument,0,22},
        {"alu-cns",required_argument,0,23},{"sva-cns",required_argument,0,24},{"bwa-threads",required_argument,0,25},
        {"disc-pos-in",required_argument,0,26},{"mate-fastq",required_argument,0,27},{"realign-fastq-in",required_argument,0,28},
        {"no-clip-supplementary",no_argument,0,29},
        {"l1-bed",required_argument,0,30},{"alu-bed",required_argument,0,31},{"sva-bed",required_argument,0,32},
        {"clip-realign-sam",required_argument,0,33},
        {"cns-clip-sam",required_argument,0,34},{"cns-disc-sam",required_argument,0,35},
        {"cns-candidates",required_argument,0,36},{"cns-ref",required_argument,0,37},
        {"cns-extnd",required_argument,0,38},{"cns-sample",required_argument,0,39},
        {"cns-tmp-prefix",required_argument,0,40},{"cns-lowq",required_argument,0,41},
        {"cns-emit-manifest",required_argument,0,42},
        {"cns-collect-clip",required_argument,0,43},{"cns-collect-anchor",required_argument,0,44},
        {"cns-collect-lowq",required_argument,0,45},{"cns-fan-manifest",required_argument,0,46},
        {"cns-fan-clip",required_argument,0,47},{"cns-fan-anchor",required_argument,0,48},
        {"cns-fan-lowq",required_argument,0,49},
        {"help",no_argument,0,'h'},{0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "b:r:R:o:t:h", lopts, nullptr)) != -1) {
        switch (c) {
            case 'b': o.bams = split_ws(optarg); break; case 'r': o.reference = optarg; break;
            case 'R': o.region = optarg; break; case 'o': o.out = optarg; break;
            case 10: o.clip_sites = optarg; break; case 11: o.disc_sites = optarg; break;
            case 12: o.clip_fastq = optarg; break; case 14: o.disc_pos = optarg; break;
            case 't': o.threads = std::max(1, atoi(optarg)); break;
            case 13: o.tile_size = std::max((hts_pos_t)1, (hts_pos_t)atoll(optarg)); break;
            case 1: o.min_clip_mapq = atoi(optarg); break; case 2: o.min_disc_mapq = atoi(optarg); break;
            case 3: o.disc_threshold = atol(optarg); break; case 4: o.no_baseq = true; break;
            case 20: o.realign_out = optarg; break; case 21: o.bwa = optarg; break;
            case 22: o.l1_cns = optarg; break; case 23: o.alu_cns = optarg; break;
            case 24: o.sva_cns = optarg; break; case 25: o.bwa_threads = std::max(1, atoi(optarg)); break;
            case 26: o.disc_pos_in = optarg; break; case 27: o.mate_fastq = optarg; break;
            case 28: o.realign_fastq_in = optarg; break;
            case 29: o.clip_supplementary = false; break;
            case 30: o.l1_bed = optarg; break; case 31: o.alu_bed = optarg; break;
            case 32: o.sva_bed = optarg; break;
            case 33: o.clip_realign_sam = optarg; break;
            case 34: o.cns_clip_sam = optarg; break; case 35: o.cns_disc_sam = optarg; break;
            case 36: o.cns_candidates = optarg; break; case 37: o.cns_ref = optarg; break;
            case 38: o.cns_extnd = atol(optarg); break; case 39: o.cns_sample = atoi(optarg); break;
            case 40: o.cns_tmp_prefix = optarg; break; case 41: o.cns_lowq = optarg; break;
            case 42: o.cns_manifest = optarg; break;
            case 43: o.cns_collect_clip = optarg; break; case 44: o.cns_collect_anchor = optarg; break;
            case 45: o.cns_collect_lowq = optarg; break; case 46: o.cns_fan_manifest = optarg; break;
            case 47: o.cns_fan_clip = optarg; break; case 48: o.cns_fan_anchor = optarg; break;
            case 49: o.cns_fan_lowq = optarg; break;
            case 'h': usage(stdout); return 0; default: usage(stderr); return 2;
        }
    }
    if (o.bams.empty()) { usage(stderr); return 2; }
    // Multi-file precondition, checked once before anything is built off hdr0.
    if (!validate_headers(o.bams, o.reference)) return 2;

    SamFile fp0; Idx idx0; Hdr hdr0;
    // hdr0/tiles come from bams[0] alone -- safe because validate_headers has just
    // proven every input presents the identical @SQ dictionary.
    if (!open_private(o, fp0, idx0, hdr0)) { std::fprintf(stderr, "error: cannot open %s (need index; reference for CRAM)\n", o.bams[0].c_str()); return 1; }
    std::optional<Region> reg;
    if (!o.region.empty()) { reg = parse_region(o.region); if (!reg) { std::fprintf(stderr, "error: bad --region\n"); return 2; } }
    std::vector<Tile> tiles = build_tiles(hdr0.get(), reg, o.tile_size);

    // ---- fold-in cns-emit: collect ONCE, emit temp SAMs for ALL types in a manifest ----
    // manifest line (whitespace-sep): tag candidate_list consensus out_clip_sam out_disc_sam out_lowq
    if (!o.cns_manifest.empty()) {
        std::vector<xtea::cns::CnsTypeJob> jobs;
        FILE* mf = std::fopen(o.cns_manifest.c_str(), "r");
        if (!mf) { std::fprintf(stderr, "error: cannot open manifest %s\n", o.cns_manifest.c_str()); return 1; }
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, mf)) != -1) {
            char *sv = nullptr;
            char* f0 = strtok_r(line, " \t\n", &sv); char* f1 = strtok_r(nullptr, " \t\n", &sv);
            char* f2 = strtok_r(nullptr, " \t\n", &sv); char* f3 = strtok_r(nullptr, " \t\n", &sv);
            char* f4 = strtok_r(nullptr, " \t\n", &sv); char* f5 = strtok_r(nullptr, " \t\n", &sv);
            if (!f0 || !f1 || !f2 || !f3) continue;
            xtea::cns::CnsTypeJob J;
            J.tag = f0; J.candidate_list = f1; J.consensus = f2; J.out_clip_sam = f3;
            J.out_disc_sam = f4 ? f4 : ""; J.out_lowq = f5 ? f5 : "";
            J.cand.load(J.candidate_list);
            jobs.push_back(std::move(J));
        }
        free(line); std::fclose(mf);
        if (jobs.empty()) { std::fprintf(stderr, "error: empty manifest\n"); return 2; }
        std::string tmp = o.cns_tmp_prefix;
        if (tmp.empty()) { char t[] = "/tmp/xtea_cnsm_XXXXXX"; int fd = mkstemp(t); if (fd>=0){ close(fd); ::unlink(t); } tmp = t; }
        if (!xtea::cns::emit_cns_multi(o.bams, o.reference, jobs, o.cns_extnd, o.cns_sample,
                                       o.bwa, o.bwa_threads, o.region, tmp)) {
            std::fprintf(stderr, "error: cns emit (manifest) failed\n"); return 1;
        }
        return 0;
    }

    // ---- full fold-in FAN: load the detection-pass spill, gather mates (1 BAM pass), emit temp SAMs ----
    if (!o.cns_fan_manifest.empty()) {
        std::vector<xtea::cns::CnsTypeJob> jobs;
        FILE* mf = std::fopen(o.cns_fan_manifest.c_str(), "r");
        if (!mf) { std::fprintf(stderr, "error: cannot open fan manifest %s\n", o.cns_fan_manifest.c_str()); return 1; }
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, mf)) != -1) {
            char *sv = nullptr;
            char* f0 = strtok_r(line, " \t\n", &sv); char* f1 = strtok_r(nullptr, " \t\n", &sv);
            char* f2 = strtok_r(nullptr, " \t\n", &sv); char* f3 = strtok_r(nullptr, " \t\n", &sv);
            char* f4 = strtok_r(nullptr, " \t\n", &sv); char* f5 = strtok_r(nullptr, " \t\n", &sv);
            if (!f0 || !f1 || !f2 || !f3) continue;
            xtea::cns::CnsTypeJob J;
            J.tag = f0; J.candidate_list = f1; J.consensus = f2; J.out_clip_sam = f3;
            J.out_disc_sam = f4 ? f4 : ""; J.out_lowq = f5 ? f5 : "";
            J.cand.load(J.candidate_list);
            jobs.push_back(std::move(J));
        }
        free(line); std::fclose(mf);
        if (jobs.empty()) { std::fprintf(stderr, "error: empty fan manifest\n"); return 2; }
        std::string tmp = o.cns_tmp_prefix;
        if (tmp.empty()) { char t[] = "/tmp/xtea_cnsf_XXXXXX"; int fd = mkstemp(t); if (fd>=0){ close(fd); ::unlink(t); } tmp = t; }
        if (!xtea::cns::emit_cns_from_spill(o.bams, o.reference, jobs, o.cns_fan_clip, o.cns_fan_anchor,
                                            o.cns_fan_lowq, o.cns_extnd, o.cns_sample, o.bwa,
                                            o.bwa_threads, o.region, tmp)) {
            std::fprintf(stderr, "error: cns fan (from spill) failed\n"); return 1;
        }
        return 0;
    }

    // ---- move B.1 cns-emit mode: reproduce cns_filter's temp_clip.sam (and temp_disc.sam) ----
    if (!o.cns_clip_sam.empty() || !o.cns_disc_sam.empty()) {
        if (o.cns_candidates.empty() || o.cns_ref.empty()) {
            std::fprintf(stderr, "error: cns-emit needs --cns-candidates and --cns-ref\n"); return 2;
        }
        std::string tmp = o.cns_tmp_prefix;
        if (tmp.empty()) { char t[] = "/tmp/xtea_cns_XXXXXX"; int fd = mkstemp(t); if (fd>=0){ close(fd); ::unlink(t); } tmp = t; }
        xtea::cns::Candidates cand; cand.load(o.cns_candidates);
        if (!xtea::cns::emit_cns(o.bams, o.reference, cand, o.cns_extnd, o.cns_sample,
                                 o.cns_ref, o.bwa, o.bwa_threads, o.region,
                                 o.cns_clip_sam, o.cns_disc_sam, o.cns_lowq, tmp)) {
            std::fprintf(stderr, "error: cns emit failed\n"); return 1;
        }
        return 0;
    }

    // ---- Part 2 gather-mates mode: fetch disc mate sequences by name into --mate-fastq ----
    if (!o.disc_pos_in.empty()) {
        if (o.mate_fastq.empty()) { std::fprintf(stderr, "error: --disc-pos-in needs --mate-fastq\n"); return 2; }
        auto needs = load_mate_needs(o.disc_pos_in, hdr0.get());
        ConcFile mate_out;
        if (!mate_out.open(o.mate_fastq)) { std::fprintf(stderr, "error: cannot write %s\n", o.mate_fastq.c_str()); return 1; }
        std::atomic<size_t> gn{0};
        const int gt = std::min<int>(o.threads, (int)std::max<size_t>(1, tiles.size()));
        auto gworker = [&]() {
            // Unlike the cns-emit modes, this mode's input (the disc-pos stream) was
            // already pooled across all files by the detection pass, so an anchor here
            // carries no file identity -- every file's tiles must be searched for a
            // matching mate, not just one. Hence one handle per file, gather_tile once
            // per file per tile.
            //
            // Residual, accepted: if two files coincidentally share a qname at the exact
            // position needed, this emits one extra genuinely-matching read rather than
            // corrupting anything. This mode is not invoked by run_xtea.py today
            // (superseded by the pure-Python disc_adapter.py reading the disc-pos stream
            // directly), so the exposure is theoretical.
            std::vector<SamFile> fps(o.bams.size());
            std::vector<Idx>     idxs(o.bams.size());
            std::vector<Hdr>     hdrs(o.bams.size());
            for (size_t b = 0; b < o.bams.size(); ++b)
                if (!open_one(o.bams[b], o.reference, fps[b], idxs[b], hdrs[b])) return;
            for (;;) { size_t i = gn.fetch_add(1); if (i >= tiles.size()) break;
                for (size_t b = 0; b < o.bams.size(); ++b)
                    gather_tile(tiles[i], fps[b].get(), idxs[b].get(), needs, mate_out); }
        };
        if (gt <= 1) gworker();
        else { std::vector<std::thread> pool; for (int i = 0; i < gt; ++i) pool.emplace_back(gworker); for (auto& t : pool) t.join(); }
        mate_out.close();
        // Optionally realign the gathered mates directly.
        if (!o.realign_out.empty()) {
            std::vector<std::pair<std::string, std::string>> fams;
            if (!o.l1_cns.empty()) fams.emplace_back("LINE1", o.l1_cns);
            if (!o.alu_cns.empty()) fams.emplace_back("ALU", o.alu_cns);
            if (!o.sva_cns.empty()) fams.emplace_back("SVA", o.sva_cns);
            if (!fams.empty()) xtea::realign_sites_file(o.mate_fastq, fams, o.bwa, o.bwa_threads, o.realign_out);
        }
        return 0;
    }

    std::vector<TileResult> results(tiles.size());

    // Streaming clip FASTQ: to --clip-fastq if given, else a temp file when realign is requested.
    std::string fastq_path = o.clip_fastq;
    bool fastq_is_temp = false;
    if (fastq_path.empty() && !o.realign_out.empty()) {
        char tmpl[] = "/tmp/xtea_clipfq_XXXXXX";
        int fd = mkstemp(tmpl); if (fd >= 0) { close(fd); fastq_path = tmpl; fastq_is_temp = true; }
    }
    ConcFile fastq_file, discpos_file;
    ConcFile* fastq_ptr = nullptr; ConcFile* discpos_ptr = nullptr;
    if (!fastq_path.empty()) { if (fastq_file.open(fastq_path)) fastq_ptr = &fastq_file; }
    if (!o.disc_pos.empty()) { if (discpos_file.open(o.disc_pos)) discpos_ptr = &discpos_file; }

    // Full fold-in: spill the cns collection during THIS detection pass (candidate-agnostic).
    ConcFile cns_clip_file, cns_anchor_file, cns_lowq_file;
    ConcFile *cns_clip_ptr = nullptr, *cns_anchor_ptr = nullptr, *cns_lowq_ptr = nullptr;
    if (!o.cns_collect_clip.empty() && !o.cns_collect_anchor.empty() && !o.cns_collect_lowq.empty()) {
        if (cns_clip_file.open(o.cns_collect_clip) && cns_anchor_file.open(o.cns_collect_anchor)
            && cns_lowq_file.open(o.cns_collect_lowq)) {
            cns_clip_ptr = &cns_clip_file; cns_anchor_ptr = &cns_anchor_file; cns_lowq_ptr = &cns_lowq_file;
        }
    }

    // Load per-family repeat BEDs into in-memory htslib region indexes (shared
    // read-only across tile threads) for the col5 mate-in-repeat count.
    regidx_t* rep_idx[3] = {nullptr, nullptr, nullptr};
    const std::string* bed_paths[3] = {&o.l1_bed, &o.alu_bed, &o.sva_bed};
    for (int f = 0; f < 3; ++f) {
        if (bed_paths[f]->empty()) continue;
        rep_idx[f] = regidx_init(bed_paths[f]->c_str(), regidx_parse_bed, nullptr, 0, nullptr);
        if (!rep_idx[f]) { std::fprintf(stderr, "error: failed to load repeat BED %s\n", bed_paths[f]->c_str()); return 2; }
    }

    const int nthreads = std::min<int>(o.threads, (int)std::max<size_t>(1, tiles.size()));
    std::atomic<size_t> next{0};
    auto worker = [&]() {
        // One private handle set PER INPUT FILE. Each tile is then scanned once per
        // file, accumulating into the SAME results[i] TileResult -- a flat sum across
        // files, which is exactly master's own pooling semantics for clip/disc counts
        // (x_basic_info/x_clip_disc_filter/x_TEI_locator flat-sum across bam_list
        // lines). scan_tile itself needs no change: clip_agg/disc_agg/clip_matein are
        // (tid,pos)-keyed count maps with no qname component, so pooling cannot collide
        // regardless of file count, and clip_recs/disc_recs (the opt-in --out dump) are
        // plain vectors that concatenate safely.
        std::vector<SamFile> fps(o.bams.size());
        std::vector<Idx>     idxs(o.bams.size());
        std::vector<Hdr>     hdrs(o.bams.size());
        for (size_t b = 0; b < o.bams.size(); ++b)
            if (!open_one(o.bams[b], o.reference, fps[b], idxs[b], hdrs[b])) return;
        for (;;) {
            size_t i = next.fetch_add(1);
            if (i >= tiles.size()) break;
            for (size_t b = 0; b < o.bams.size(); ++b)
                scan_tile(tiles[i], o, fps[b].get(), idxs[b].get(), hdrs[b].get(), results[i],
                          fastq_ptr, discpos_ptr, rep_idx,
                          cns_clip_ptr, cns_anchor_ptr, cns_lowq_ptr, (int)b);
        }
    };
    if (nthreads <= 1) worker();
    else { std::vector<std::thread> pool; for (int i = 0; i < nthreads; ++i) pool.emplace_back(worker); for (auto& t : pool) t.join(); }
    fastq_file.close(); discpos_file.close();
    cns_clip_file.close(); cns_anchor_file.close(); cns_lowq_file.close();

    auto name = [&](int tid){ return sam_hdr_tid2name(hdr0.get(), tid); };

    // Merge col5 union mate-in-repeat counts across tiles, keyed by (chrom name, pos),
    // for realign_sites_file to attach as the 6th column (same value for every family).
    std::map<std::pair<std::string, long>, int> matein_merged;
    if (rep_idx[0] || rep_idx[1] || rep_idx[2]) {
        for (auto& r : results)
            for (auto& kv : r.clip_matein)
                matein_merged[{name(kv.first.first), (long)kv.first.second}] += kv.second;
    }

    // Merge aggregation maps (std::map merges in sorted order) -> emit sorted.
    if (!o.clip_sites.empty()) {
        std::map<std::pair<int, hts_pos_t>, std::pair<int, int>> agg;
        for (auto& r : results) for (auto& kv : r.clip_agg) { auto& a = agg[kv.first]; a.first += kv.second.first; a.second += kv.second.second; }
        if (FILE* f = std::fopen(o.clip_sites.c_str(), "w")) {
            for (auto& kv : agg) std::fprintf(f, "%s\t%" PRIhts_pos "\t%d\t%d\n", name(kv.first.first), kv.first.second, kv.second.first, kv.second.second);
            std::fclose(f);
        }
    }
    if (!o.disc_sites.empty()) {
        std::map<std::pair<int, hts_pos_t>, int> agg;
        for (auto& r : results) for (auto& kv : r.disc_agg) agg[kv.first] += kv.second;
        if (FILE* f = std::fopen(o.disc_sites.c_str(), "w")) {
            for (auto& kv : agg) std::fprintf(f, "%s\t%" PRIhts_pos "\t%d\n", name(kv.first.first), kv.first.second, kv.second);
            std::fclose(f);
        }
    }
    // Per-read dump (opt-in): merge + sort for determinism.
    if (!o.out.empty()) {
        std::vector<ClipRec> clips; std::vector<DiscRec> discs;
        for (auto& r : results) { for (auto& e : r.clip_recs) clips.push_back(std::move(e)); for (auto& e : r.disc_recs) discs.push_back(std::move(e)); }
        std::sort(clips.begin(), clips.end(), [](const ClipRec& a, const ClipRec& b){
            if (a.tid != b.tid) return a.tid < b.tid; if (a.pos != b.pos) return a.pos < b.pos;
            if (a.side != b.side) return a.side < b.side; return a.qname < b.qname; });
        std::sort(discs.begin(), discs.end(), [](const DiscRec& a, const DiscRec& b){
            if (a.tid != b.tid) return a.tid < b.tid; if (a.pos != b.pos) return a.pos < b.pos; return a.qname < b.qname; });
        FILE* f = (o.out == "-") ? stdout : std::fopen(o.out.c_str(), "w");
        if (f) {
            std::string buf; buf.reserve(1 << 20); char line[512];
            for (auto& e : clips) { int n = std::snprintf(line, sizeof line, "CLIP\t%s\t%" PRIhts_pos "\t%c\t%s\t%d\t%d\n", name(e.tid), e.pos, e.side, e.qname.c_str(), e.is_rev, e.clen); buf.append(line, n); if (buf.size() >= (1u<<20)) { std::fwrite(buf.data(),1,buf.size(),f); buf.clear(); } }
            for (auto& e : discs) { int n = std::snprintf(line, sizeof line, "DISC\t%s\t%" PRIhts_pos "\t%s\t%s\t%" PRIhts_pos "\t%d\t%d\n", name(e.tid), e.pos, e.qname.c_str(), name(e.mtid), e.mpos, e.is_rev, e.is_mrev); buf.append(line, n); if (buf.size() >= (1u<<20)) { std::fwrite(buf.data(),1,buf.size(),f); buf.clear(); } }
            if (!buf.empty()) std::fwrite(buf.data(),1,buf.size(),f);
            if (f != stdout) std::fclose(f);
        }
    }

    // ---- Stage 2: realign the streamed clip FASTQ vs the family consensi (concurrent bwa) ----
    if (!o.realign_out.empty()) {
        std::vector<std::pair<std::string, std::string>> fams;
        if (!o.l1_cns.empty()) fams.emplace_back("LINE1", o.l1_cns);
        if (!o.alu_cns.empty()) fams.emplace_back("ALU", o.alu_cns);
        if (!o.sva_cns.empty()) fams.emplace_back("SVA", o.sva_cns);
        if (fams.empty()) { std::fprintf(stderr, "error: --realign-out needs --l1-cns/--alu-cns/--sva-cns\n"); return 2; }
        // realign an explicit FASTQ if given (e.g. a mate FASTQ), else the scanned clip FASTQ.
        std::string rfq = o.realign_fastq_in.empty() ? fastq_path : o.realign_fastq_in;
        if (rfq.empty()) { std::fprintf(stderr, "error: no FASTQ produced for realign\n"); return 1; }
        xtea::realign_sites_file(rfq, fams, o.bwa, o.bwa_threads, o.realign_out, &matein_merged,
                                 o.clip_realign_sam);
    }
    for (int f = 0; f < 3; ++f) if (rep_idx[f]) regidx_destroy(rep_idx[f]);
    if (fastq_is_temp) ::unlink(fastq_path.c_str());
    return 0;
}
