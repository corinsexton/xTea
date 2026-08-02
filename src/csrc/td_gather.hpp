// td_gather.hpp — native per-site read gather for the transduction step (Stage 1).
//
// WHAT IT REPLACES, and why this and nothing else:
//   The transduction step's read gather is 73-87% of the step's wall time (measured, job
//   49437990: clip_gather 577.6s L1 / 622.8s SVA, disc_gather 280.1s / 275.5s, out of 995.6s /
//   1231.9s), which is what made it Stage 1's lead item.  It is two Python functions:
//
//     collect_clipped_reads      (x_clip_disc_filter.py:3309)
//       -> collect_clipped_reads_of_given_list (:336): a fresh multiprocessing Pool PER
//          CHROMOSOME (~25 fork storms), one task per candidate site (~54K for L1 on WGS), each
//          task writing TWO tiny files the parent immediately re-reads and unlinks; produces
//          the merged clip FASTQ and, as a side effect, all_disc_pos.disc_pos.
//     collect_disc_reads         (:3335)
//       -> extract_mate_reads_by_name (x_alignments.py:478): a Pool over 50 Mb genome bins where
//          EVERY bin re-parses the whole all_disc_pos file and fetches its whole 50 Mb range --
//          i.e. a full genome scan plus O(bins x file) text parsing -- then merges per-bin temp
//          FASTAs.
//
//   Both phases run here instead, in one process: the candidate windows are fetched MERGED and
//   threaded (no fork, no per-site temp files), and the mate gather fetches only the mate
//   positions the anchors actually recorded, after ONE parse.
//
// WHAT IT MUST NOT CHANGE:
//   Stage 1 is agreement-gated on byte-identical calls, so this reproduces the legacy outputs
//   BYTE FOR BYTE, not merely record-for-record.  Two ordering facts make that load-bearing:
//
//   1. bwa's tie-break among EQUAL-SCORING placements depends on a read's index in its input
//      (cns_emit.hpp:798-817, measured).  Reordering the clip FASTQ silently moves which repeat
//      copy a multi-mapping clip lands on, and with it borderline calls.  So the FASTQ is emitted
//      CANDIDATE-MAJOR in the SITE FILE's own order: chromosomes in first-appearance order, sites
//      in file order within a chromosome -- which is what the legacy per-chromosome Pool plus
//      per-site file concatenation produces (m_chrm_pos is an insertion-ordered dict; the merge
//      loop walks l_chrm_records).  The site file (sf_raw_slct) is NOT position-sorted -- the cns
//      candidates are written first, then the raw-disc lines -- so "sorted by position" would be
//      the wrong order here, unlike in cns_emit.
//   2. Master's n_cnt (read-name field 5) has GAPS: the right-clip branch builds the name and
//      increments the counter BEFORE applying that side's length/polyA/quality gates
//      (x_clip_disc_filter.py:200-224), so a right clip that fails burns an index.  No parser
//      reads field 5 (parse_transduction_clip_algnmt uses 0,1,2,4), so the values are free -- but
//      reproducing them is what makes `cmp` a valid oracle for this port, so the burns are held
//      as no-emit ClipEvents and counted.
//
// The disc FASTA's order is likewise reproduced: BAM-header chromosome order (decoy/contig names
// skipped, x_reference.py:18), coordinate order within a chromosome, and for one mate read its
// insertion sites in the order they first appear in the disc-pos stream.
//
// DIVERGENCES, both deliberate and both unreachable from the Stage 1 gate:
//   * MULTI-BAM.  Legacy writes all_disc_pos.disc_pos to ONE fixed path, truncated once per input
//     BAM (x_clip_disc_filter.py:349), so with N>1 inputs only the LAST BAM's disc positions
//     survive and every sample's mates are then gathered against them -- a bug, silent, and with
//     no baseline (all germline runs are single-BAM; the multi-CRAM mosaic path is Stage 4).  Here
//     each input BAM is collected and gathered against ITS OWN anchors.  Single-BAM output is
//     unaffected, which is what the gate compares.
//   * An exact duplicate (chrom,pos) line in the site file is collapsed here.  Legacy CRASHES on
//     one (the second Pool task overwrites the first's temp file, which the merge loop then opens
//     after unlinking it, x_clip_disc_filter.py:373 -- no isfile guard, unlike the cns variant),
//     so this can only differ on input legacy rejects.
//   * A site chromosome absent from the BAM header is skipped with a warning; legacy raises out of
//     a Pool worker and fails the step.
//   * A read with NO base qualities skips the clip quality gate here; legacy raises TypeError on
//     `None[:n]`. Same shape as the two above: it differs only where legacy cannot run.
//
// ONE ASSUMPTION, stated because it is not checked: an anchor's recorded mate position equals the
// mate's own reference start (r->core.mpos == the primary mate's pos). That is what lets the mate
// gather fetch by position instead of scanning; a BAM whose mate fields disagree with its mate
// records would gather a different set than the legacy 50 Mb-bin scan. cns_emit.hpp's gather rests
// on the same assumption (:454-456), and the byte-identical disc FASTA on real data is evidence
// for it rather than proof.
//
// Reuses, deliberately, from cns_emit.hpp: the read-level gate ports (contain_enough_A_T,
// is_qualified_clip, is_two_side_clipped, is_decoy_contig), the constants mirrored from
// global_values.py, and the merged-window FetchRegion plan with its exact ownership rule
// (:502-523).  Those are the two patterns the plan names; the divergence in emission ORDER
// (site-file order here vs position order there) is why this is a separate collector rather than
// a parameter on emit_cns_multi.
#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <future>
#include <unistd.h>

#include "htslib/sam.h"
#include "cns_emit.hpp"

namespace xtea {
namespace td {

using xtea::cns::contain_enough_A_T;
using xtea::cns::is_qualified_clip;
using xtea::cns::is_two_side_clipped;
using xtea::cns::is_decoy_contig;
using xtea::cns::kMinClipMapq;
using xtea::cns::kMinDiscMapq;
using xtea::cns::kMinPolyaClip;
using xtea::cns::kBwaRealignCut;
using xtea::cns::kTwoSideClipMin;
using xtea::cns::kDiscThreshold;

// ---------------------------------------------------------------------------------------------
// The site list, kept twice on purpose: once in FILE ORDER (that is the emission order, see
// header note 1) and once sorted for overlap queries.
// ---------------------------------------------------------------------------------------------
struct TdSites {
    std::vector<std::string> chrom_order;                 // first appearance in the file
    std::map<std::string, std::vector<long>> file_order;  // per chrom, file order, exact dups dropped
    cns::Candidates lookup;                               // sorted+unique, for the fetch plan
    size_t n_sites = 0, n_dup = 0;

    bool load(const std::string& path) {
        FILE* f = std::fopen(path.c_str(), "r");
        if (!f) return false;
        std::map<std::string, std::set<long>> seen;
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, f)) != -1) {
            char* sv = nullptr;
            char* c0 = strtok_r(line, " \t\n", &sv);
            char* c1 = strtok_r(nullptr, " \t\n", &sv);
            if (!c0 || !c1) continue;
            std::string chrom(c0);
            const long pos = atol(c1);
            auto ins = seen[chrom].insert(pos);
            if (!ins.second) { ++n_dup; continue; }
            auto it = file_order.find(chrom);
            if (it == file_order.end()) { chrom_order.push_back(chrom); it = file_order.emplace(chrom, std::vector<long>()).first; }
            it->second.push_back(pos);
            lookup.by_chrom[chrom].push_back(pos);
            ++n_sites;
        }
        free(line);
        std::fclose(f);
        for (auto& kv : lookup.by_chrom) {
            std::sort(kv.second.begin(), kv.second.end());
            kv.second.erase(std::unique(kv.second.begin(), kv.second.end()), kv.second.end());
        }
        return true;
    }
};

// The one site set the BAM is decoded against when several entries share a pass.
//
// It carries the SORTED index only, deliberately: `file_order` is what fixes emission order, and
// that is a property of each entry's own site file, not of the union -- fanning from a union file
// order would silently reorder every entry's FASTQ. Only build_fetch_plan and collect_stripe's
// hold test read this, and both want the sorted view.
//
// The union is nearly free on real data: 53,945 of ~54K transduction sites are shared by all three
// families, so one decode serves every entry. Two entries whose site files spell the same contig
// differently ("1" vs "chr1") are NOT merged -- they stay separate keys and that contig is decoded
// twice. Wasteful, still correct (each entry fans from its own key), and it cannot happen within
// one run, where every site file comes from the same pipeline.
inline TdSites union_sites(const std::vector<TdSites>& all) {
    TdSites U;
    for (const TdSites& S : all) {
        for (const auto& chrom : S.chrom_order) {
            auto it = U.lookup.by_chrom.find(chrom);
            if (it == U.lookup.by_chrom.end()) {
                U.chrom_order.push_back(chrom);
                it = U.lookup.by_chrom.emplace(chrom, std::vector<long>()).first;
            }
            const auto& v = S.lookup.by_chrom.at(chrom);
            it->second.insert(it->second.end(), v.begin(), v.end());
        }
    }
    for (auto& kv : U.lookup.by_chrom) {
        std::sort(kv.second.begin(), kv.second.end());
        kv.second.erase(std::unique(kv.second.begin(), kv.second.end()), kv.second.end());
        U.n_sites += kv.second.size();
    }
    return U;
}

// _process_chrm_name (x_clip_disc_filter.py:63): the site file's style is authoritative for the
// NAME we emit; the BAM's style is what we fetch by.
inline std::string to_bam_chrom(const std::string& chrm, bool b_tmplt_with_chr) {
    const bool b_with_chr = (chrm.size() > 3 && chrm.compare(0, 3, "chr") == 0);
    if (b_tmplt_with_chr && !b_with_chr) return "chr" + chrm;
    if (!b_tmplt_with_chr && b_with_chr) return chrm.substr(3);
    return chrm;
}

// ---------------------------------------------------------------------------------------------
// Held, candidate-agnostic collection products.  A read can produce up to two ClipEvents (left
// then right, which is master's write order) and one AnchorEvent.
// ---------------------------------------------------------------------------------------------
struct ClipEvent {
    long ref_start, ref_end;   // read span; decides which candidate windows it belongs to
    long map_pos;              // clip position -> read-name field 1
    char side;                 // 'L' / 'R' (FLAG_LEFT_CLIP / FLAG_RIGHT_CLIP)
    int is_rc;
    bool emit;                 // false = right clip that burned an n_cnt index then failed a gate
    std::string seq, qual;     // empty when !emit
};
struct AnchorEvent {
    long ref_start, ref_end;
    long map_pos;                       // clip-adjusted, exactly as master's disc line carries it
    int is_rc, is_mate_rc, mate_first, b_lclip, b_rclip;
    long mate_pos;
    std::string qname, mate_chrm;       // mate_chrm in BAM style (master writes it that way)
};
struct ChromEvents {
    std::vector<ClipEvent> clips;
    std::vector<AnchorEvent> anchors;
};

// One family's gather. Several of these share ONE decode of the BAM (see run_td_gather_multi);
// alone, it is exactly the pre-existing single-family job.
struct TdGatherEntry {
    std::string tag = "td";      // logging only; the manifest's first column
    std::string sites_file;      // sf_raw_slct: "chrom pos ..." per line; ORDER IS SIGNIFICANT
    long extnd = 400;            // iextnd: window is [max(1,pos-extnd), pos+extnd)
    std::string out_clip_fastq;  // == the legacy raw_candidate_sites_all_clip.fq
    std::string out_disc_fa;     // == the legacy raw_candidate_sites_all_disc.fa
    std::string out_disc_pos;    // optional: the legacy all_disc_pos.disc_pos, for oracle diffs
};

struct TdGatherJob {
    std::string sites_file;
    long extnd = 400;
    std::string out_clip_fastq;
    std::string out_disc_fa;
    std::string out_disc_pos;
    int threads = 1;
};

// ---------------------------------------------------------------------------------------------
// Phase 1 -- collect from the merged candidate windows.
// ---------------------------------------------------------------------------------------------
struct FetchRegion { std::string chrom; int tid; long beg, end, own_from; };

// Build the merged fetch plan.  Identical rule to cns_emit.hpp:502-523: a read belongs to the
// region it is fetched by iff its START is >= own_from, the previous region's end on that contig.
// Regions therefore partition the reads exactly and can be striped across threads; the FAN below
// is region-independent (it queries the whole site index), so a read that reaches into a later
// region's window is still emitted for that window.
//
// UNTRIED IDEA, recorded because it is the obvious next lever if this phase ever needs one.  At
// WGS scale collect is NOT the slower half -- 22.0s against the mate gather's 29.1s for 49,042 L1
// sites at 12 threads -- so this is speculative, not indicated.  (An earlier version of this note
// said collect dominated; that came from the chr1 numbers, where it does, and it does not
// generalise.)  What it costs is seeks: 33,007 regions, each a separate index lookup and BGZF
// seek.  Merging regions that are merely CLOSE (say within 64 KB, not just overlapping)
// would share decompressed blocks between neighbouring windows.  It cannot change the output:
// the ownership rule holds for any merge distance and the fan does not consult regions at all.
// But it reads the gaps too -- with extnd~750 that is ~1.5 KB of window per up-to-64 KB of gap --
// so it is a seek-cost-versus-byte-cost trade that has to be MEASURED, on a real WGS site list,
// before it is believed.  Do not "optimise" it blind.
inline std::vector<FetchRegion> build_fetch_plan(const TdSites& sites, long extnd,
                                                 sam_hdr_t* hdr, bool b_with_chr) {
    std::vector<FetchRegion> fetch;
    for (const auto& chrom : sites.chrom_order) {
        const std::string bamc = to_bam_chrom(chrom, b_with_chr);
        const int tid = sam_hdr_name2tid(hdr, bamc.c_str());
        if (tid < 0) {
            // lookup, not file_order: the UNION site set built by union_sites() populates only
            // the sorted index (file order is per entry and meaningless for a union). The two
            // have identical size for a single loaded entry -- both drop exact duplicates --
            // so this warning's number is unchanged for the one-entry path.
            std::fprintf(stderr, "[td-gather] WARNING: site contig '%s' (as '%s') is not in the "
                                 "BAM header; %zu sites skipped\n",
                         chrom.c_str(), bamc.c_str(), sites.lookup.by_chrom.at(chrom).size());
            continue;
        }
        const auto& v = sites.lookup.by_chrom.at(chrom);   // sorted, unique
        long prev_end = 0;
        for (size_t i = 0; i < v.size(); ) {
            long b = std::max(0L, v[i] - extnd), e = v[i] + extnd;
            while (++i < v.size() && (v[i] - extnd) <= e) e = std::max(e, v[i] + extnd);
            fetch.push_back({chrom, tid, b, e, prev_end});
            prev_end = e;
        }
    }
    return fetch;
}

// One thread's stripe of the fetch plan.  This is collect_clipped_disc_reads_by_region
// (x_clip_disc_filter.py:80) read for read, including its control flow: a gate failure on either
// clip side `continue`s the WHOLE read, so it also suppresses that read's disc record.
inline bool collect_stripe(const std::string& bam, const std::string& reference,
                           const std::vector<FetchRegion>& fetch, size_t lo, size_t hi,
                           long extnd, const TdSites& sites,
                           std::map<std::string, ChromEvents>& out) {
    samFile* fp = sam_open(bam.c_str(), "r");
    if (!fp) return false;
    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); return false; }
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());
    if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); return false; }
    bam1_t* b = bam_init1();
    std::string sseq, squal;
    for (size_t ri = lo; ri < hi; ++ri) {
        const FetchRegion& R = fetch[ri];
        ChromEvents& CE = out[R.chrom];
        // This contig's sorted site positions, resolved once per region rather than per read:
        // the window test below runs on every read in the region.
        const std::vector<long>& spos = sites.lookup.by_chrom.at(R.chrom);
        hts_itr_t* it = sam_itr_queryi(idx, R.tid, R.beg, R.end);
        if (!it) continue;
        while (sam_itr_next(fp, it, b) >= 0) {
            const bam1_t* r = b;
            if (r->core.pos < R.own_from) continue;          // owned by an earlier region
            const uint16_t flag = r->core.flag;
            if (flag & (BAM_FDUP | BAM_FUNMAP)) continue;
            const uint32_t* cig = bam_get_cigar(r);
            const int ncig = r->core.n_cigar;
            if (ncig < 1) continue;
            // Low-mapq reads are dropped here.  Master also tallies them into
            // n_cnt_low_mapq_clip, but collect_clipped_reads throws that return value away
            // (x_clip_disc_filter.py:3320) -- the count is only consumed by the cns path.
            if (r->core.qual < kMinClipMapq) continue;
            if (is_two_side_clipped(cig, ncig, kTwoSideClipMin)) continue;

            const long rs = r->core.pos, re = bam_endpos(r);
            // Only reads that reach at least one candidate window can ever be fanned; dropping
            // the rest here is what keeps peak memory O(window reads) rather than O(region).
            // Membership is pysam.fetch's: C-extnd < re AND C+extnd > rs, i.e. some site in
            // (rs-extnd, re+extnd).  Deliberately WITHOUT the max(1,.) clamp the fan applies, so
            // this stays a superset -- holding is permissive, emission is exact.
            const bool in_window = [&]{
                auto p = std::upper_bound(spos.begin(), spos.end(), rs - extnd);
                return p != spos.end() && *p < re + extnd;
            }();

            const int is_rc = bam_is_rev(r) ? 1 : 0;
            const bool b_first = !(flag & BAM_FREAD2);
            const bool is_supp_sec = (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0;
            const uint8_t* seq = bam_get_seq(r);
            const uint8_t* qual = bam_get_qual(r);
            const int qlen = r->core.l_qseq;
            const bool have_q = (qual && qual[0] != 0xff);
            std::string mate_chrm = "*"; long mate_pos = 0;
            if (!(flag & BAM_FMUNMAP) && r->core.mtid >= 0) {
                mate_chrm = sam_hdr_tid2name(hdr, r->core.mtid);
                mate_pos = r->core.mpos;
            }
            const int is_mate_rc = (flag & BAM_FMREVERSE) ? 1 : 0;

            long map_pos = rs;
            int b_left_clip = 0, b_right_clip = 0;
            if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) {
                if (is_supp_sec) continue;
                const int clen = bam_cigar_oplen(cig[0]);
                sseq.clear();
                for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, i)];
                bool pass = true;
                if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
                else if (clen < kBwaRealignCut) pass = false;
                if (!pass) continue;                          // no n_cnt burn: master gates first
                if (have_q && !is_qualified_clip(qual, clen)) continue;
                if (in_window) {
                    squal.clear();
                    for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[i] + 33) : 'I';
                    CE.clips.push_back({rs, re, rs, 'L', is_rc, true, sseq, squal});
                }
                b_left_clip = 1;
            }
            if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {
                for (int i = 0; i < ncig - 1; ++i) {
                    const int op = bam_cigar_op(cig[i]);
                    if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CINS) continue;
                    map_pos += bam_cigar_oplen(cig[i]);
                }
                if (is_supp_sec) continue;
                const int clen = bam_cigar_oplen(cig[ncig - 1]);
                const int off = qlen - clen;
                sseq.clear();
                for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, off + i)];
                // Master builds the read name -- and increments n_cnt -- HERE, before these
                // three gates.  A right clip that fails still consumed an index, so the event is
                // held with emit=false to keep the counter's gaps.
                bool pass = true;
                if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
                if (clen < kBwaRealignCut) pass = false;
                if (pass && have_q && !is_qualified_clip(qual + off, clen)) pass = false;
                if (in_window) {
                    if (pass) {
                        squal.clear();
                        for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[off + i] + 33) : 'I';
                        CE.clips.push_back({rs, re, map_pos, 'R', is_rc, true, sseq, squal});
                    } else {
                        CE.clips.push_back({rs, re, map_pos, 'R', is_rc, false, std::string(), std::string()});
                    }
                }
                if (!pass) continue;
                b_right_clip = 1;
            }

            // ---- disc anchor ----
            if (mate_chrm == "*") continue;
            if (is_decoy_contig(mate_chrm)) continue;
            if (r->core.qual < kMinDiscMapq) continue;
            const char* rchrom = sam_hdr_tid2name(hdr, r->core.tid);
            const bool disc = (mate_chrm != rchrom) || (llabs(mate_pos - map_pos) > kDiscThreshold);
            if (!disc || !in_window) continue;
            CE.anchors.push_back({rs, re, map_pos, is_rc, is_mate_rc, b_first ? 0 : 1,
                                  b_left_clip, b_right_clip, mate_pos,
                                  bam_get_qname(r), mate_chrm});
        }
        hts_itr_destroy(it);
    }
    bam_destroy1(b);
    hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
    return true;
}

// ---------------------------------------------------------------------------------------------
// Phase 2 -- what the mate gather needs, built while the disc-pos records are emitted so the
// per-read site ORDER is the disc-pos file's order by construction (that is the order
// _parse_read_name_pos_of_region's dict ends up in, and the FASTA follows it).
// ---------------------------------------------------------------------------------------------
struct MateSite {
    std::string ins_chrm; long ins_pos;
    int b_lclip, b_rclip, anchor_rc, mate_rc; long anchor_map_pos;
};
struct MateNeed {
    std::string mate_chrm; long mate_pos = -1;
    std::vector<MateSite> sites;    // first-appearance order; one entry per distinct site
};

inline void note_mate_need(std::unordered_map<std::string, MateNeed>& needs,
                          const std::string& key, const std::string& mate_chrm, long mate_pos,
                          const MateSite& S) {
    MateNeed& N = needs[key];
    N.mate_chrm = mate_chrm; N.mate_pos = mate_pos;
    for (auto& e : N.sites)
        if (e.ins_pos == S.ins_pos && e.ins_chrm == S.ins_chrm) { e = S; return; }  // dict overwrite
    N.sites.push_back(S);
}

// ---------------------------------------------------------------------------------------------
// Phase 3 -- fan ONE entry's windows out of the shared, candidate-agnostic collection.
//
// This is per entry and nothing about it is shared: the emission order is the entry's own site
// file order (header note 1), the window is the entry's own extnd, and the n_cnt index is counted
// within the entry's own window.  `ev` may hold reads no window of this entry covers -- the union
// hold test is deliberately permissive -- and the two explicit `ref_end <= wbeg` tests below are
// what make emission exact anyway.  That is the whole reason a shared decode cannot move bytes.
// ---------------------------------------------------------------------------------------------
struct FanStats { size_t n_fq = 0, n_dpos = 0; };

inline bool fan_entry(const TdGatherEntry& E, const TdSites& S,
                      const std::map<std::string, ChromEvents>& ev,
                      long max_clip_span, long max_anchor_span, size_t bi,
                      std::unordered_map<std::string, MateNeed>& needs, FanStats& st) {
    FILE* fq = E.out_clip_fastq.empty() ? nullptr : std::fopen(E.out_clip_fastq.c_str(), "a");
    FILE* dp = E.out_disc_pos.empty()   ? nullptr : std::fopen(E.out_disc_pos.c_str(), "a");
    if (!E.out_clip_fastq.empty() && !fq) {
        std::fprintf(stderr, "[td-gather] error: cannot append %s\n", E.out_clip_fastq.c_str());
        if (dp) std::fclose(dp);
        return false;
    }
    for (const auto& chrom : S.chrom_order) {
        auto eit = ev.find(chrom);
        if (eit == ev.end()) continue;
        const auto& clips = eit->second.clips;
        const auto& anchors = eit->second.anchors;
        for (long C : S.file_order.at(chrom)) {
            const long wbeg = std::max(1L, C - E.extnd), wend = C + E.extnd;
            // A read is in this window iff ref_start < wend AND ref_end > wbeg (pysam.fetch
            // overlap).  max_span bounds how far back a read can start and still reach wbeg,
            // which is what makes the lower bound a binary search.
            long nc = 0;
            auto clo = std::lower_bound(clips.begin(), clips.end(), wbeg - max_clip_span,
                [](const ClipEvent& a, long x) { return a.ref_start < x; });
            auto chi = std::lower_bound(clips.begin(), clips.end(), wend,
                [](const ClipEvent& a, long x) { return a.ref_start < x; });
            for (auto p = clo; p != chi; ++p) {
                if (p->ref_end <= wbeg) continue;
                const long idx = nc++;                    // burns even when !emit
                if (!p->emit || !fq) continue;
                std::fprintf(fq, "@%s~%ld~%c~%d~%ld~%ld~%zu\n%s\n+\n%s\n",
                             chrom.c_str(), p->map_pos, p->side, p->is_rc, C, idx, bi,
                             p->seq.c_str(), p->qual.c_str());
                ++st.n_fq;
            }
            auto alo = std::lower_bound(anchors.begin(), anchors.end(), wbeg - max_anchor_span,
                [](const AnchorEvent& a, long x) { return a.ref_start < x; });
            auto ahi = std::lower_bound(anchors.begin(), anchors.end(), wend,
                [](const AnchorEvent& a, long x) { return a.ref_start < x; });
            for (auto p = alo; p != ahi; ++p) {
                if (p->ref_end <= wbeg) continue;
                if (dp)
                    std::fprintf(dp, "%s\t%ld\t%d\t%d\t%s\t%ld\t%s\t%d\t%ld\t%d\t%d\n",
                                 chrom.c_str(), p->map_pos, p->is_rc, p->is_mate_rc,
                                 p->mate_chrm.c_str(), p->mate_pos, p->qname.c_str(),
                                 p->mate_first, C, p->b_lclip, p->b_rclip);
                ++st.n_dpos;
                if (E.out_disc_fa.empty()) continue;
                note_mate_need(needs, p->qname + "~" + std::to_string(p->mate_first),
                               p->mate_chrm, p->mate_pos,
                               {chrom, C, p->b_lclip, p->b_rclip, p->is_rc, p->is_mate_rc, p->map_pos});
            }
        }
    }
    if (fq) std::fclose(fq);
    if (dp) std::fclose(dp);
    return true;
}

// ---------------------------------------------------------------------------------------------
// The whole gather, for every input BAM and every entry.  Sample-major, matching the legacy
// merge: all of sample 0's records, then sample 1's, with `~<sample>` appended to each record's
// name.
//
// WHY A VECTOR OF ENTRIES (plan §4 Stage 2 / 2b).  L1's and SVA's transduction gathers ask the
// same BAM for almost exactly the same windows -- 53,945 of ~54K sites are shared by all three
// families -- and today they are two separate processes, each decoding the whole candidate set on
// its own.  Running as one job decodes once.  Note what this does NOT buy: the two processes are
// already concurrent, so wall time is unchanged (~57s either way); what it saves is CPU, and the
// win only becomes wall time once the per-locus loop makes the downstream one job at full width.
//
// WHY IT CANNOT MOVE BYTES, which is the property that matters more than the saving:
//
//   * COLLECTION is a superset.  build_fetch_plan runs over the union site index at the LARGEST
//     entry's extnd, and collect_stripe's `in_window` hold test with it, so `ev` holds every read
//     any single-entry run would have held, plus some.  The header already states the invariant
//     the extras rely on -- "holding is permissive, emission is exact".
//   * ORDER survives the wider plan.  `ev` is accumulated per contig in fetch-plan order and each
//     region owns reads starting at or after the previous region's end, so starts are ascending
//     however the regions happen to be merged; the stable_sort that follows is therefore a no-op
//     and equal-start reads keep BAM order, which is what bwa's tie-break rides on.
//   * EMISSION is per entry and untouched: its own site file order, its own extnd, its own n_cnt.
//   * THE MATE GATHER merges its region plan over the union of the entries' wanted mate positions.
//     Different boundaries, same ownership rule -- a read is emitted by exactly one region either
//     way, and contigs are still walked in BAM-header order with ascending starts inside each --
//     so each entry's FASTA is the byte-identical sequence it would get alone.
//
// That argument is checked rather than trusted: tests/test_td_gather.py runs a two-entry manifest
// against two solo runs and requires `cmp` on every product.
// ---------------------------------------------------------------------------------------------
inline bool run_td_gather_multi(const std::vector<std::string>& bams, const std::string& reference,
                                const std::vector<TdGatherEntry>& entries, int threads) {
    if (entries.empty()) {
        std::fprintf(stderr, "[td-gather] error: no entries\n");
        return false;
    }
    // Two entries sharing an output path would interleave into it and neither would be readable.
    // Cheap to check, impossible to debug from the resulting file.
    {
        std::set<std::string> seen_out;
        for (const auto& E : entries)
            for (const std::string* p : {&E.out_clip_fastq, &E.out_disc_fa, &E.out_disc_pos}) {
                if (p->empty()) continue;
                if (!seen_out.insert(*p).second) {
                    std::fprintf(stderr, "[td-gather] error: two entries write the same output "
                                         "path '%s'\n", p->c_str());
                    return false;
                }
            }
    }

    std::vector<TdSites> sites(entries.size());
    long max_extnd = 1;
    for (size_t ei = 0; ei < entries.size(); ++ei) {
        const TdGatherEntry& E = entries[ei];
        if (!sites[ei].load(E.sites_file)) {
            std::fprintf(stderr, "[td-gather] error: cannot read site list %s\n",
                         E.sites_file.c_str());
            return false;
        }
        if (sites[ei].n_dup)
            std::fprintf(stderr, "[td-gather] %s: NOTE: %zu duplicate (chrom,pos) site lines "
                                 "collapsed (the legacy Python path crashes on these)\n",
                         E.tag.c_str(), sites[ei].n_dup);
        std::fprintf(stderr, "[td-gather] %s: %zu sites over %zu contigs, extnd=%ld, threads=%d\n",
                     E.tag.c_str(), sites[ei].n_sites, sites[ei].chrom_order.size(), E.extnd,
                     threads);
        max_extnd = std::max(max_extnd, E.extnd);
    }
    const TdSites uni = union_sites(sites);
    if (entries.size() > 1)
        std::fprintf(stderr, "[td-gather] union: %zu distinct sites over %zu contigs, "
                             "extnd=%ld (one decode serves %zu entries)\n",
                     uni.n_sites, uni.chrom_order.size(), max_extnd, entries.size());

    // Truncate the outputs once; every per-BAM pass appends.
    for (const auto& E : entries)
        for (const std::string* p : {&E.out_clip_fastq, &E.out_disc_fa, &E.out_disc_pos}) {
            if (p->empty()) continue;
            FILE* f = std::fopen(p->c_str(), "w");
            if (!f) { std::fprintf(stderr, "[td-gather] error: cannot write %s\n", p->c_str()); return false; }
            std::fclose(f);
        }

    // Attribution, three numbers, because "the gather" is now one [TD-PHASE] line on the Python
    // side and the next stage needs to know which third of it to attack: collect = the candidate
    // windows (what the per-site Pool did), fan = pure formatting, mate = the genome-wide mate
    // fetch (what the 50 Mb bins did).
    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point a, clk::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    for (size_t bi = 0; bi < bams.size(); ++bi) {
        // ---- header facts, once per file ----
        samFile* fp0 = sam_open(bams[bi].c_str(), "r");
        if (!fp0) { std::fprintf(stderr, "[td-gather] error: cannot open %s\n", bams[bi].c_str()); return false; }
        if (!reference.empty()) hts_set_fai_filename(fp0, reference.c_str());
        sam_hdr_t* hdr0 = sam_hdr_read(fp0);
        if (!hdr0) { sam_close(fp0); std::fprintf(stderr, "[td-gather] error: no header in %s\n", bams[bi].c_str()); return false; }
        const bool b_with_chr = (sam_hdr_name2tid(hdr0, "chr1") >= 0);   // is_chrm_contain_chr
        const int nref = sam_hdr_nref(hdr0);

        std::vector<FetchRegion> fetch = build_fetch_plan(uni, max_extnd, hdr0, b_with_chr);
        const auto t_collect0 = clk::now();

        // ---- collect, striped ----
        std::map<std::string, ChromEvents> ev;
        {
            const int nthreads = std::max(1, std::min<int>(threads, (int)std::max<size_t>(fetch.size(), 1)));
            std::vector<std::future<std::pair<bool, std::map<std::string, ChromEvents>>>> futs;
            for (int ti = 0; ti < nthreads; ++ti) {
                const size_t lo = fetch.size() * ti / nthreads;
                const size_t hi = fetch.size() * (ti + 1) / nthreads;
                futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
                    std::map<std::string, ChromEvents> local;
                    const bool ok = collect_stripe(bams[bi], reference, fetch, lo, hi, max_extnd, uni, local);
                    return std::make_pair(ok, std::move(local));
                }));
            }
            // Stripes are contiguous slices of a plan that walks each contig's regions in
            // ascending order, and accumulation is per contig, so merging in thread order
            // restores per-contig ascending order -- the order a per-candidate pysam fetch sees.
            // (A stripe may straddle two contigs; that is harmless for the same reason.)
            for (auto& f : futs) {
                auto pr = f.get();
                if (!pr.first) {
                    std::fprintf(stderr, "[td-gather] error: collection failed on %s\n", bams[bi].c_str());
                    sam_hdr_destroy(hdr0); sam_close(fp0); return false;
                }
                for (auto& kv : pr.second) {
                    ChromEvents& dst = ev[kv.first];
                    dst.clips.insert(dst.clips.end(), std::make_move_iterator(kv.second.clips.begin()),
                                     std::make_move_iterator(kv.second.clips.end()));
                    dst.anchors.insert(dst.anchors.end(), std::make_move_iterator(kv.second.anchors.begin()),
                                       std::make_move_iterator(kv.second.anchors.end()));
                }
            }
        }
        // Belt and braces: the merge above is already ascending, and a stable sort cannot
        // reorder ties, so the within-position BAM order that bwa's tie-break rides on survives.
        long max_clip_span = 1, max_anchor_span = 1;
        for (auto& kv : ev) {
            std::stable_sort(kv.second.clips.begin(), kv.second.clips.end(),
                             [](const ClipEvent& a, const ClipEvent& b) { return a.ref_start < b.ref_start; });
            std::stable_sort(kv.second.anchors.begin(), kv.second.anchors.end(),
                             [](const AnchorEvent& a, const AnchorEvent& b) { return a.ref_start < b.ref_start; });
            for (const auto& c : kv.second.clips) max_clip_span = std::max(max_clip_span, c.ref_end - c.ref_start);
            for (const auto& a : kv.second.anchors) max_anchor_span = std::max(max_anchor_span, a.ref_end - a.ref_start);
        }

        // ---- fan to candidates, per entry, CANDIDATE-MAJOR IN ITS OWN SITE-FILE ORDER ----
        const auto t_fan0 = clk::now();
        std::vector<std::unordered_map<std::string, MateNeed>> needs(entries.size());
        std::vector<FanStats> st(entries.size());
        for (size_t ei = 0; ei < entries.size(); ++ei) {
            if (!fan_entry(entries[ei], sites[ei], ev, max_clip_span, max_anchor_span, bi,
                           needs[ei], st[ei])) {
                sam_hdr_destroy(hdr0); sam_close(fp0); return false;
            }
        }
        ev.clear();
        const auto t_mate0 = clk::now();
        for (size_t ei = 0; ei < entries.size(); ++ei)
            std::fprintf(stderr, "[td-gather] bam %zu %s: %zu clip records, %zu disc-pos records, "
                                 "%zu mates to gather\n",
                         bi, entries[ei].tag.c_str(), st[ei].n_fq, st[ei].n_dpos, needs[ei].size());

        // ---- gather the mates, by position, and emit each entry's disc FASTA ----
        bool any_disc_fa = false;
        for (const auto& E : entries) any_disc_fa = any_disc_fa || !E.out_disc_fa.empty();
        if (any_disc_fa) {
            // Targets, merged per contig, over the UNION of what the entries want.  Contigs are
            // walked in BAM-HEADER order and decoy / contig names are skipped, because that is
            // the order and the filter extract_mate_reads_by_name's 50 Mb bin list has
            // (x_alignments.py:487).  Within a contig the legacy bins tile [0,len) and own a read
            // by its start, so per-contig ascending emission is the same sequence; own_from does
            // the same job here, and it does it for any merge boundaries -- which is why widening
            // the plan to the union leaves every entry's own record sequence untouched.
            std::map<std::string, std::vector<long>> want;
            for (size_t ei = 0; ei < entries.size(); ++ei) {
                if (entries[ei].out_disc_fa.empty()) continue;
                for (const auto& kv : needs[ei])
                    if (kv.second.mate_pos >= 0) want[kv.second.mate_chrm].push_back(kv.second.mate_pos);
            }
            constexpr long kMergeSlack = 4096;
            std::vector<FetchRegion> mreg;
            for (int tid = 0; tid < nref; ++tid) {
                const char* nm = sam_hdr_tid2name(hdr0, tid);
                if (!nm) continue;
                const std::string chrom(nm);
                if (is_decoy_contig(chrom)) continue;
                auto wit = want.find(chrom);
                if (wit == want.end()) continue;
                auto& v = wit->second;
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                long prev_end = 0;
                for (size_t i = 0; i < v.size(); ) {
                    long beg = v[i], end = v[i] + 1;
                    while (++i < v.size() && v[i] - end <= kMergeSlack) end = v[i] + 1;
                    mreg.push_back({chrom, tid, beg, end, prev_end});
                    prev_end = end;
                }
            }
            const int nthreads = std::max(1, std::min<int>(threads, (int)std::max<size_t>(mreg.size(), 1)));
            // Each stripe writes its own part file PER ENTRY and the parts are concatenated in
            // stripe order; stripes are contiguous in an already-ordered region list, so that IS
            // the (contig, position) order.  Parts rather than in-memory buffers because the
            // merged FASTA is ~1 GB per family on WGS.
            std::vector<std::vector<std::string>> parts(entries.size(),
                                                        std::vector<std::string>(nthreads));
            std::vector<std::future<bool>> futs;
            for (int ti = 0; ti < nthreads && !mreg.empty(); ++ti) {
                std::vector<std::string> myparts(entries.size());
                for (size_t ei = 0; ei < entries.size(); ++ei) {
                    if (entries[ei].out_disc_fa.empty()) continue;
                    myparts[ei] = entries[ei].out_disc_fa + ".part" + std::to_string(bi) + "_"
                                + std::to_string(ti);
                    parts[ei][ti] = myparts[ei];
                }
                const size_t lo = mreg.size() * ti / nthreads;
                const size_t hi = mreg.size() * (ti + 1) / nthreads;
                futs.push_back(std::async(std::launch::async, [&, lo, hi, myparts]() {
                    std::vector<FILE*> out(entries.size(), nullptr);
                    bool opened = true;
                    for (size_t ei = 0; ei < entries.size(); ++ei) {
                        if (myparts[ei].empty()) continue;
                        out[ei] = std::fopen(myparts[ei].c_str(), "w");
                        if (!out[ei]) opened = false;
                    }
                    auto close_all = [&]() { for (FILE* f : out) if (f) std::fclose(f); };
                    if (!opened) { close_all(); return false; }
                    samFile* fp = sam_open(bams[bi].c_str(), "r");
                    if (!fp) { close_all(); return false; }
                    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
                    sam_hdr_t* hdr = sam_hdr_read(fp);
                    if (!hdr) { sam_close(fp); close_all(); return false; }
                    hts_idx_t* idx = sam_index_load(fp, bams[bi].c_str());
                    if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); close_all(); return false; }
                    bam1_t* b = bam_init1();
                    std::string mseq;
                    for (size_t ri = lo; ri < hi; ++ri) {
                        hts_itr_t* it = sam_itr_queryi(idx, mreg[ri].tid, mreg[ri].beg, mreg[ri].end);
                        if (!it) continue;
                        while (sam_itr_next(fp, it, b) >= 0) {
                            const bam1_t* r = b;
                            // Master excludes secondary/supplementary here for a measured
                            // reason: a soft-clipped supplementary carries the full query
                            // sequence and would emit a byte-identical duplicate record,
                            // inflating that site's disc support (x_alignments.py:439-456).
                            if (r->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
                            if (r->core.pos < mreg[ri].own_from) continue;
                            const int is_first = (r->core.flag & BAM_FREAD1) ? 1 : 0;
                            const std::string key = std::string(bam_get_qname(r)) + "~"
                                                  + std::to_string(is_first);
                            // The sequence is decoded at most once per read, however many
                            // entries claim it -- that, and the single BGZF decode above, are
                            // the whole saving.
                            bool decoded = false;
                            for (size_t ei = 0; ei < entries.size(); ++ei) {
                                if (!out[ei]) continue;
                                auto nit = needs[ei].find(key);
                                if (nit == needs[ei].end()) continue;
                                // The bin parse only admitted records whose mate_chrm was this
                                // contig, so a same-name read elsewhere in the genome is not a hit.
                                if (nit->second.mate_chrm != mreg[ri].chrom) continue;
                                if (!decoded) {
                                    const uint8_t* s = bam_get_seq(r); const int L = r->core.l_qseq;
                                    mseq.clear();
                                    for (int i = 0; i < L; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
                                    decoded = true;
                                }
                                for (const auto& S : nit->second.sites)
                                    std::fprintf(out[ei], ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%zu\n%s\n",
                                                 bam_get_qname(r), is_first, S.b_lclip, S.b_rclip,
                                                 S.anchor_rc, S.mate_rc, S.anchor_map_pos,
                                                 S.ins_chrm.c_str(), S.ins_pos, bi, mseq.c_str());
                            }
                        }
                        hts_itr_destroy(it);
                    }
                    bam_destroy1(b);
                    hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
                    close_all();
                    return true;
                }));
            }
            bool ok = true;
            for (auto& f : futs) ok = f.get() && ok;
            if (!ok) {
                std::fprintf(stderr, "[td-gather] error: mate gather failed on %s\n", bams[bi].c_str());
                sam_hdr_destroy(hdr0); sam_close(fp0); return false;
            }
            std::vector<char> buf(1u << 20);
            for (size_t ei = 0; ei < entries.size(); ++ei) {
                if (entries[ei].out_disc_fa.empty()) continue;
                FILE* fa = std::fopen(entries[ei].out_disc_fa.c_str(), "a");
                if (!fa) {
                    std::fprintf(stderr, "[td-gather] error: cannot append %s\n",
                                 entries[ei].out_disc_fa.c_str());
                    sam_hdr_destroy(hdr0); sam_close(fp0); return false;
                }
                size_t n_fa_bytes = 0;
                for (const auto& part : parts[ei]) {
                    if (part.empty()) continue;
                    FILE* in = std::fopen(part.c_str(), "r");
                    if (!in) continue;
                    size_t got;
                    while ((got = std::fread(buf.data(), 1, buf.size(), in)) > 0) {
                        std::fwrite(buf.data(), 1, got, fa);
                        n_fa_bytes += got;
                    }
                    std::fclose(in);
                    ::unlink(part.c_str());
                }
                std::fclose(fa);
                std::fprintf(stderr, "[td-gather] bam %zu %s: disc FASTA %zu bytes over %zu "
                                     "fetch regions\n",
                             bi, entries[ei].tag.c_str(), n_fa_bytes, mreg.size());
            }
        }

        const auto t_end = clk::now();
        std::fprintf(stderr, "[td-gather] bam %zu: collect %.1fs (%zu regions, %d threads), "
                             "fan %.1fs, mate-gather %.1fs\n",
                     bi, secs(t_collect0, t_fan0), fetch.size(), threads,
                     secs(t_fan0, t_mate0), secs(t_mate0, t_end));
        std::fflush(stderr);

        sam_hdr_destroy(hdr0); sam_close(fp0);
    }
    return true;
}

// The single-family form, kept because that is what XClipDiscFilter.native_td_gather still calls.
// It is a one-entry vector and NOTHING else -- there is no separate single-entry code path to
// drift, which is what makes "the multi-entry version cannot change the one-entry output" a
// structural claim rather than a tested coincidence.
inline bool run_td_gather(const std::vector<std::string>& bams, const std::string& reference,
                          const TdGatherJob& J) {
    TdGatherEntry E;
    E.tag = "td";
    E.sites_file = J.sites_file; E.extnd = J.extnd;
    E.out_clip_fastq = J.out_clip_fastq; E.out_disc_fa = J.out_disc_fa;
    E.out_disc_pos = J.out_disc_pos;
    return run_td_gather_multi(bams, reference, {E}, J.threads);
}

}  // namespace td
}  // namespace xtea
