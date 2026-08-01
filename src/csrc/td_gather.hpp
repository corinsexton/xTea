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
//
// Reuses, deliberately, from cns_emit.hpp: the read-level gate ports (contain_enough_A_T,
// is_qualified_clip, is_two_side_clipped, is_decoy_contig), the constants mirrored from
// global_values.py, and the merged-window FetchRegion plan with its exact ownership rule
// (:502-523).  Those are the two patterns the plan names; the divergence in emission ORDER
// (site-file order here vs position order there) is why this is a separate collector rather than
// a parameter on emit_cns_multi.
#pragma once

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

struct TdGatherJob {
    std::string sites_file;      // sf_raw_slct: "chrom pos ..." per line; ORDER IS SIGNIFICANT
    long extnd = 400;            // iextnd: window is [max(1,pos-extnd), pos+extnd)
    std::string out_clip_fastq;  // == the legacy raw_candidate_sites_all_clip.fq
    std::string out_disc_fa;     // == the legacy raw_candidate_sites_all_disc.fa
    std::string out_disc_pos;    // optional: the legacy all_disc_pos.disc_pos, for oracle diffs
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
inline std::vector<FetchRegion> build_fetch_plan(const TdSites& sites, long extnd,
                                                 sam_hdr_t* hdr, bool b_with_chr) {
    std::vector<FetchRegion> fetch;
    for (const auto& chrom : sites.chrom_order) {
        const std::string bamc = to_bam_chrom(chrom, b_with_chr);
        const int tid = sam_hdr_name2tid(hdr, bamc.c_str());
        if (tid < 0) {
            std::fprintf(stderr, "[td-gather] WARNING: site contig '%s' (as '%s') is not in the "
                                 "BAM header; %zu sites skipped\n",
                         chrom.c_str(), bamc.c_str(), sites.file_order.at(chrom).size());
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
// The whole gather, for every input BAM.  Sample-major, matching the legacy merge: all of
// sample 0's records, then sample 1's, with `~<sample>` appended to each record's name.
// ---------------------------------------------------------------------------------------------
inline bool run_td_gather(const std::vector<std::string>& bams, const std::string& reference,
                          const TdGatherJob& J) {
    TdSites sites;
    if (!sites.load(J.sites_file)) {
        std::fprintf(stderr, "[td-gather] error: cannot read site list %s\n", J.sites_file.c_str());
        return false;
    }
    if (sites.n_dup)
        std::fprintf(stderr, "[td-gather] NOTE: %zu duplicate (chrom,pos) site lines collapsed "
                             "(the legacy Python path crashes on these)\n", sites.n_dup);
    std::fprintf(stderr, "[td-gather] %zu sites over %zu contigs, extnd=%ld, threads=%d\n",
                 sites.n_sites, sites.chrom_order.size(), J.extnd, J.threads);

    // Truncate the outputs once; every per-BAM pass appends.
    for (const std::string* p : {&J.out_clip_fastq, &J.out_disc_fa, &J.out_disc_pos}) {
        if (p->empty()) continue;
        FILE* f = std::fopen(p->c_str(), "w");
        if (!f) { std::fprintf(stderr, "[td-gather] error: cannot write %s\n", p->c_str()); return false; }
        std::fclose(f);
    }

    for (size_t bi = 0; bi < bams.size(); ++bi) {
        // ---- header facts, once per file ----
        samFile* fp0 = sam_open(bams[bi].c_str(), "r");
        if (!fp0) { std::fprintf(stderr, "[td-gather] error: cannot open %s\n", bams[bi].c_str()); return false; }
        if (!reference.empty()) hts_set_fai_filename(fp0, reference.c_str());
        sam_hdr_t* hdr0 = sam_hdr_read(fp0);
        if (!hdr0) { sam_close(fp0); std::fprintf(stderr, "[td-gather] error: no header in %s\n", bams[bi].c_str()); return false; }
        const bool b_with_chr = (sam_hdr_name2tid(hdr0, "chr1") >= 0);   // is_chrm_contain_chr
        const int nref = sam_hdr_nref(hdr0);

        std::vector<FetchRegion> fetch = build_fetch_plan(sites, J.extnd, hdr0, b_with_chr);

        // ---- collect, striped ----
        std::map<std::string, ChromEvents> ev;
        {
            const int nthreads = std::max(1, std::min<int>(J.threads, (int)std::max<size_t>(fetch.size(), 1)));
            std::vector<std::future<std::pair<bool, std::map<std::string, ChromEvents>>>> futs;
            for (int ti = 0; ti < nthreads; ++ti) {
                const size_t lo = fetch.size() * ti / nthreads;
                const size_t hi = fetch.size() * (ti + 1) / nthreads;
                futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
                    std::map<std::string, ChromEvents> local;
                    const bool ok = collect_stripe(bams[bi], reference, fetch, lo, hi, J.extnd, sites, local);
                    return std::make_pair(ok, std::move(local));
                }));
            }
            // Stripes are contiguous slices of a plan that is already in (contig, position)
            // order, so merging in thread order restores per-contig ascending order -- the order
            // a per-candidate pysam fetch sees.
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

        // ---- fan to candidates, CANDIDATE-MAJOR IN SITE-FILE ORDER ----
        FILE* fq = J.out_clip_fastq.empty() ? nullptr : std::fopen(J.out_clip_fastq.c_str(), "a");
        FILE* dp = J.out_disc_pos.empty() ? nullptr : std::fopen(J.out_disc_pos.c_str(), "a");
        if (!J.out_clip_fastq.empty() && !fq) {
            std::fprintf(stderr, "[td-gather] error: cannot append %s\n", J.out_clip_fastq.c_str());
            sam_hdr_destroy(hdr0); sam_close(fp0); return false;
        }
        std::unordered_map<std::string, MateNeed> needs;
        size_t n_fq = 0, n_dpos = 0;
        for (const auto& chrom : sites.chrom_order) {
            auto eit = ev.find(chrom);
            if (eit == ev.end()) continue;
            const auto& clips = eit->second.clips;
            const auto& anchors = eit->second.anchors;
            for (long C : sites.file_order.at(chrom)) {
                const long wbeg = std::max(1L, C - J.extnd), wend = C + J.extnd;
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
                    ++n_fq;
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
                    ++n_dpos;
                    if (J.out_disc_fa.empty()) continue;
                    note_mate_need(needs, p->qname + "~" + std::to_string(p->mate_first),
                                   p->mate_chrm, p->mate_pos,
                                   {chrom, C, p->b_lclip, p->b_rclip, p->is_rc, p->is_mate_rc, p->map_pos});
                }
            }
        }
        if (fq) std::fclose(fq);
        if (dp) std::fclose(dp);
        ev.clear();
        std::fprintf(stderr, "[td-gather] bam %zu: %zu clip records, %zu disc-pos records, "
                             "%zu mates to gather\n", bi, n_fq, n_dpos, needs.size());

        // ---- gather the mates, by position, and emit the disc FASTA ----
        if (!J.out_disc_fa.empty()) {
            // Targets, merged per contig.  Contigs are walked in BAM-HEADER order and decoy /
            // contig names are skipped, because that is the order and the filter
            // extract_mate_reads_by_name's 50 Mb bin list has (x_alignments.py:487).  Within a
            // contig the legacy bins tile [0,len) and own a read by its start, so per-contig
            // ascending emission is the same sequence; own_from does the same job here.
            std::map<std::string, std::vector<long>> want;
            for (const auto& kv : needs)
                if (kv.second.mate_pos >= 0) want[kv.second.mate_chrm].push_back(kv.second.mate_pos);
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
            const int nthreads = std::max(1, std::min<int>(J.threads, (int)std::max<size_t>(mreg.size(), 1)));
            // Each stripe writes its own part file and the parts are concatenated in stripe
            // order; stripes are contiguous in an already-ordered region list, so that IS the
            // (contig, position) order.  Parts rather than in-memory buffers because the merged
            // FASTA is ~1 GB per family on WGS.
            std::vector<std::string> parts(nthreads);
            std::vector<std::future<bool>> futs;
            for (int ti = 0; ti < nthreads && !mreg.empty(); ++ti) {
                parts[ti] = J.out_disc_fa + ".part" + std::to_string(bi) + "_" + std::to_string(ti);
                const size_t lo = mreg.size() * ti / nthreads;
                const size_t hi = mreg.size() * (ti + 1) / nthreads;
                const std::string part = parts[ti];
                futs.push_back(std::async(std::launch::async, [&, lo, hi, part]() {
                    FILE* out = std::fopen(part.c_str(), "w");
                    if (!out) return false;
                    samFile* fp = sam_open(bams[bi].c_str(), "r");
                    if (!fp) { std::fclose(out); return false; }
                    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
                    sam_hdr_t* hdr = sam_hdr_read(fp);
                    if (!hdr) { sam_close(fp); std::fclose(out); return false; }
                    hts_idx_t* idx = sam_index_load(fp, bams[bi].c_str());
                    if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); std::fclose(out); return false; }
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
                            auto nit = needs.find(std::string(bam_get_qname(r)) + "~" + std::to_string(is_first));
                            if (nit == needs.end()) continue;
                            // The bin parse only admitted records whose mate_chrm was this
                            // contig, so a same-name read elsewhere in the genome is not a hit.
                            if (nit->second.mate_chrm != mreg[ri].chrom) continue;
                            const uint8_t* s = bam_get_seq(r); const int L = r->core.l_qseq;
                            mseq.clear();
                            for (int i = 0; i < L; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
                            for (const auto& S : nit->second.sites)
                                std::fprintf(out, ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%zu\n%s\n",
                                             bam_get_qname(r), is_first, S.b_lclip, S.b_rclip,
                                             S.anchor_rc, S.mate_rc, S.anchor_map_pos,
                                             S.ins_chrm.c_str(), S.ins_pos, bi, mseq.c_str());
                        }
                        hts_itr_destroy(it);
                    }
                    bam_destroy1(b);
                    hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
                    std::fclose(out);
                    return true;
                }));
            }
            bool ok = true;
            for (auto& f : futs) ok = f.get() && ok;
            if (!ok) {
                std::fprintf(stderr, "[td-gather] error: mate gather failed on %s\n", bams[bi].c_str());
                sam_hdr_destroy(hdr0); sam_close(fp0); return false;
            }
            FILE* fa = std::fopen(J.out_disc_fa.c_str(), "a");
            if (!fa) {
                std::fprintf(stderr, "[td-gather] error: cannot append %s\n", J.out_disc_fa.c_str());
                sam_hdr_destroy(hdr0); sam_close(fp0); return false;
            }
            size_t n_fa_bytes = 0;
            std::vector<char> buf(1u << 20);
            for (const auto& part : parts) {
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
            std::fprintf(stderr, "[td-gather] bam %zu: disc FASTA %zu bytes over %zu fetch regions\n",
                         bi, n_fa_bytes, mreg.size());
        }

        sam_hdr_destroy(hdr0); sam_close(fp0);
    }
    return true;
}

}  // namespace td
}  // namespace xtea
