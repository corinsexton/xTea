// cns_emit.hpp — move B.1: emit the cns_filter temp SAMs (temp_clip.sam / temp_disc.sam) from
// the scan pass, so x_clip_disc_filter's collect+realign (a per-candidate BAM re-read + bwa re-run)
// can be bypassed. This is a faithful port of collect_clipped_disc_reads_by_region
// (x_clip_disc_filter.py:81) + realign_clipped_read_with_polyA / realign_disc_reads (bwa_align.py),
// verified against ACTUAL code/constants (docs/phase_c_unified_downstream.md §7.1) and the golden's
// preserved temp SAMs (docs/cns_realign_sam_contract.md).
//
// Relationship to master: same ACCESS PATTERN — fetch each candidate's window with the index —
// because that is the part of master's design that is right. ~2,800 candidates means ~4 Mb of
// reference touched instead of all 3.1 Gb, and no C++ makes a whole-genome stream competitive with
// that. What C++ replaces is the per-read Python overhead, the process Pool, and the per-candidate
// temp files, not the algorithm. Two refinements over master, both cost-only: candidate windows are
// MERGED so overlapping ones are read once (with an exact ownership rule so no read is collected
// twice — see emit_cns_multi), and the disc mate gather fetches recorded mate positions instead of
// scanning 50 Mb bins. bwa realign params + read-name formats are byte-identical to master, so the
// emitted SAMs feed straight into parse_clip_realignment_consensus / parse_disc_algnmt_consensus.
//
// A candidate-AGNOSTIC variant used to live here (collect during the detection scan, spill to disk,
// fan afterwards). It is deleted; see the note where it was, near the bottom of this file.

#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <chrono>
#include <future>
#include <unistd.h>

#include "htslib/sam.h"
#include "realign.hpp"

namespace xtea {
namespace cns {

// --- constants mirrored from global_values.py (verified 2026-07-18) ---
constexpr int  kMinClipMapq   = 12;   // MINIMUM_CLIP_MAPQ
constexpr int  kMinDiscMapq    = 20;   // MINIMUM_DISC_MAPQ
constexpr int  kMinPolyaClip   = 7;    // MINIMUM_POLYA_CLIP
constexpr int  kBwaRealignCut  = 9;    // BWA_REALIGN_CUTOFF
constexpr int  kTwoSideClipMin = 8;    // TWO_SIDE_CLIP_MIN_LEN
constexpr int  kClipPhredCut   = 15;   // CLIP_PHRED_SCORE_CUTOFF
constexpr long kDiscThreshold  = 2000; // DISC_THRESHOLD
constexpr int  kClipLowMapq    = 5;    // CLIP_LOW_MAPQ  (mapq < 5 => background low-q clip)
constexpr int  kNearbyLowQClip = 35;   // NEARBY_LOW_Q_CLIP
// READ_LENGTH, used only to place a low-MAPQ right clip (pos + READ_LENGTH - clip length).
// master's cns step sets READ_LENGTH to the measured read length (x_TEA_main --cns calls
// global_values.set_read_length(rlth) before collecting), so the driver passes that value
// with --cns-read-length; 100 is only the fallback when it is not given.
inline int g_read_length = 100;

// Which MAPQ<5 clipped reads count toward low_mapq_clip_cnt.txt: primary alignments only.
// Master counted supplementary/secondary too (its count ran before the supp/sec skip). That only
// mattered for `bwa mem -Y` input, where supplementaries are soft-clipped: a MAPQ-0 supplementary
// clipped at a junction is usually a read whose primary alignment is in the inserted element --
// evidence for the site -- and on the 1KG CRAMs (-Y) such reads were ~40% of the count and made
// the 65% low-MAPQ ratio filter drop true insertions. Without -Y they are hard-clipped and never
// counted, so this makes the count identical with and without -Y. Mirrors
// x_clip_disc_filter.low_mapq_clip_anchor.
inline bool count_as_lowq_clip(uint16_t flag) {
    return (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) == 0;
}

// --- small ports of the master helpers ---

// contain_enough_A_T(seq, n_min): true iff #A >= n_min or #T >= n_min (x_polyA.py:97).
inline bool contain_enough_A_T(const std::string& seq, int n_min) {
    int nA = 0, nT = 0;
    for (char ch : seq) { if (ch=='A'||ch=='a') ++nA; if (ch=='T'||ch=='t') ++nT; }
    return nA >= n_min || nT >= n_min;
}

// _is_qualified_clip(l_score): false iff #(phred < 15) > len/2 (x_clip_disc_filter.py:436).
inline bool is_qualified_clip(const uint8_t* qual, int n) {
    const double half = n / 2.0; int low = 0;
    for (int i = 0; i < n; ++i) if (qual[i] < kClipPhredCut && ++low > half) return false;
    return true;
}

// is_two_side_clipped: both ends soft-clip AND both lens > i_min (x_clip_disc_filter.py:270).
inline bool is_two_side_clipped(const uint32_t* cig, int ncig, int i_min) {
    if (ncig < 1) return false;
    if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP && bam_cigar_op(cig[ncig-1]) == BAM_CSOFT_CLIP)
        return (int)bam_cigar_oplen(cig[ncig-1]) > i_min && (int)bam_cigar_oplen(cig[0]) > i_min;
    return false;
}

// Strip the "<file_idx>:" tag that multi-BAM collection prepends to anchor qnames so two inputs
// sharing a qname cannot cross-pair. The index is a decimal integer, so the FIRST colon is always
// the boundary even though Illumina qnames contain colons of their own.
inline std::string strip_file_tag(const std::string& tagged) {
    const size_t c = tagged.find(':');
    return (c == std::string::npos) ? tagged : tagged.substr(c + 1);
}

// --- candidate list ---
struct Candidates {
    std::map<std::string, std::vector<long>> by_chrom;  // sorted ascending

    void load(const std::string& path) {
        FILE* f = std::fopen(path.c_str(), "r");
        if (!f) return;
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, f)) != -1) {
            char* save = nullptr;
            char* c0 = strtok_r(line, " \t\n", &save);
            char* c1 = strtok_r(nullptr, " \t\n", &save);
            if (!c0 || !c1) continue;
            by_chrom[c0].push_back(atol(c1));
        }
        free(line);
        std::fclose(f);
        for (auto& kv : by_chrom) {
            std::sort(kv.second.begin(), kv.second.end());
            kv.second.erase(std::unique(kv.second.begin(), kv.second.end()), kv.second.end());
        }
    }

    // Candidate positions C whose window [max(1,C-extnd), C+extnd) overlaps read span [rs, re).
    // pysam.fetch membership: rs < C+extnd AND re > C-extnd  <=>  C in (rs-extnd, re+extnd).
    void overlapping(const std::string& chrom, long rs, long re, long extnd,
                     std::vector<long>& out) const {
        auto it = by_chrom.find(chrom);
        if (it == by_chrom.end()) return;
        const auto& v = it->second;
        auto lo = std::upper_bound(v.begin(), v.end(), rs - extnd);      // C > rs-extnd
        auto hi = std::lower_bound(v.begin(), v.end(), re + extnd);      // C < re+extnd
        for (auto p = lo; p != hi; ++p) out.push_back(*p);
    }

    // Candidate positions C with |P - C| < dist (lowq clip proximity, NEARBY_LOW_Q_CLIP).
    void near(const std::string& chrom, long P, long dist, std::vector<long>& out) const {
        auto it = by_chrom.find(chrom);
        if (it == by_chrom.end()) return;
        const auto& v = it->second;
        auto lo = std::lower_bound(v.begin(), v.end(), P - dist + 1);    // C > P-dist
        auto hi = std::lower_bound(v.begin(), v.end(), P + dist);        // C < P+dist
        for (auto p = lo; p != hi; ++p) out.push_back(*p);
    }
};

// is_decoy_contig_chrms port (x_reference.py:19) — decides whether a mate chrom is a decoy/contig
// (excluded from disc collection). Mitochondria treated as decoy (default GLOBAL_MITCHONDRION_SWITCH=OFF).
inline bool is_decoy_contig(const std::string& c) {
    if (c.find('_') != std::string::npos) return true;
    if (c == "hs37d5") return true;
    if (c.size() > 3) {
        auto p3 = c.substr(0, 3);
        if (p3 == "HLA" || p3 == "HPV" || p3 == "HIV" || p3 == "CMV" || p3 == "MCV") return true;
        if (c.substr(0, 2) == "SV") return true;
        if (c.substr(0, 5) == "decoy") return true;
        if (c.substr(0, 6) == "random") return true;
        if (c.substr(0, 4) == "KSHV") return true;
    }
    if (c.find('-') != std::string::npos) return true;
    if (c.find("hs38") != std::string::npos || c.find("hs37") != std::string::npos) return true;
    if (c == "MT" || c == "chrMT" || c == "chrM") return true;  // mito off by default
    if (c.find('.') != std::string::npos) return true;
    return false;
}

// One disc mate to gather + realign: keyed for lookup by (anchor read name, mate's is_read1). Carries
// the per-candidate metadata that becomes the mate FASTA header (temp_disc.sam qname).
struct DiscNeed {
    std::string ins_chrm; long ins_pos;
    long anchor_map_pos; int b_lclip, b_rclip, anchor_rc, mate_rc;
    std::string mate_chrm; long mate_pos;   // where the anchor's mate maps (gather target)
};

// ============================ CLIP + DISC emission ============================
//
// One BAM pass reproduces collect_clipped_disc_reads_by_region (x_clip_disc_filter.py:81): per read it
// runs the shared clip/disc control flow (a gate failure on either clip side `continue`s the whole read;
// map_pos is advanced by the right-clip block and reused as the disc anchor pos; b_left/b_right flag the
// disc record). Clip parts fan out to overlapping candidates -> two-cohort bwa -> temp_clip.sam. Disc
// anchors fan out -> a second pass gathers each mate's full seq by (name,is_read1) -> plain bwa ->
// temp_disc.sam. out_* left empty skips that product (clip logic still runs, needed for disc state).
//
// Returns false on a fatal I/O/bwa error.
// Realign a clip cohort pair (non-polyA -T9 -k9; polyA -T7 -k7 -c70) and merge -> out_sam.
inline bool realign_clip_cohorts(const std::string& bwa_bin, const std::string& consensus,
                                 const std::string& fq_non, const std::string& fq_pa,
                                 const std::string& tstr, const std::string& out_sam,
                                 const std::string& tmp_prefix) {
    const std::string sam_non = tmp_prefix + ".non_polyA.sam";
    const std::string sam_pa  = tmp_prefix + ".polyA.sam";
    if (!run_bwa_sam(bwa_bin, consensus, fq_non, {"-t",tstr,"-K",kBwaBatch,"-T","9","-k","9"}, sam_non)) return false;
    if (!run_bwa_sam(bwa_bin, consensus, fq_pa,  {"-t",tstr,"-K",kBwaBatch,"-T","7","-k","7","-c","70"}, sam_pa)) return false;
    FILE* out = std::fopen(out_sam.c_str(), "w");
    if (!out) return false;
    if (FILE* in = std::fopen(sam_non.c_str(), "r")) {
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, in)) != -1) std::fwrite(line, 1, len, out);
        free(line); std::fclose(in);
    }
    if (FILE* in = std::fopen(sam_pa.c_str(), "r")) {
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, in)) != -1) { if (len > 0 && line[0] == '@') continue; std::fwrite(line, 1, len, out); }
        free(line); std::fclose(in);
    }
    std::fclose(out);
    ::unlink(sam_non.c_str()); ::unlink(sam_pa.c_str());
    return true;
}

// low_mapq_clip_cnt.txt — one line per candidate ("chrm pos count"), matching master's
// m_lowq_clip_cnt (built with a 0 default for every candidate). The bypass needs this file.
inline bool emit_cns_lowq_file(const std::string& out_lowq, const Candidates& cand,
                               const std::map<std::pair<std::string,long>, long>& lowq) {
    if (out_lowq.empty()) return true;
    FILE* fl = std::fopen(out_lowq.c_str(), "w");
    if (!fl) return false;
    for (const auto& kv : cand.by_chrom)
        for (long pos : kv.second) {
            auto f = lowq.find({kv.first, pos});
            std::fprintf(fl, "%s %ld %ld\n", kv.first.c_str(), pos, f == lowq.end() ? 0L : f->second);
        }
    std::fclose(fl);
    return true;
}

inline bool emit_cns(const std::vector<std::string>& bams, const std::string& reference,
                     const Candidates& cand, long extnd, int sample_id,
                     const std::string& consensus, const std::string& bwa_bin,
                     int bwa_threads, const std::string& region,
                     const std::string& out_clip_sam, const std::string& out_disc_sam,
                     const std::string& out_lowq, const std::string& tmp_prefix) {
    const bool want_clip = !out_clip_sam.empty();
    const bool want_disc = !out_disc_sam.empty();
    // Two cohort FASTQs (bwa params differ): non-polyA (clip len >= 9), polyA (len < 9). Only if clip.
    // These, the disc FASTA, and ncnt/lowq are all POOLED across input files -- bwa runs once
    // at the end over the fully pooled read set, exactly as in the single-file case.
    const std::string fq_non = tmp_prefix + ".clip.non_polyA.fq";
    const std::string fq_pa  = tmp_prefix + ".clip.polyA.fq";
    FILE* f_non = want_clip ? std::fopen(fq_non.c_str(), "w") : nullptr;
    FILE* f_pa  = want_clip ? std::fopen(fq_pa.c_str(), "w") : nullptr;
    if (want_clip && (!f_non || !f_pa)) { if (f_non) std::fclose(f_non); if (f_pa) std::fclose(f_pa);
                           return false; }
    const std::string fa_disc = tmp_prefix + ".disc.fa";
    FILE* fa = want_disc ? std::fopen(fa_disc.c_str(), "w") : nullptr;
    if (want_disc && !fa) { if (f_non) std::fclose(f_non); if (f_pa) std::fclose(f_pa); return false; }

    // per-candidate running clip index (n_cnt); value is ignored by the parser (field 5) but kept
    // for format fidelity / uniqueness.
    std::map<std::pair<std::string,long>, long> ncnt;
    std::map<std::pair<std::string,long>, long> lowq;  // per-candidate background low-mapq clip count
    std::string sseq, squal;  // scratch

    // Multi-file: both passes run per input file, back to back. disc_needs is declared INSIDE
    // this loop, so it is fresh for every file -- a read's mate is always in the same physical
    // file as its anchor, so pass 2 for file N only ever needs the map pass 1 for file N built.
    // Two files sharing a qname therefore cannot cross-pair: their maps never coexist. (The
    // emit_cns_multi pools across files, so it cannot use this argument and tags qnames instead.)
    for (size_t bi = 0; bi < bams.size(); ++bi) {
    samFile* fp = sam_open(bams[bi].c_str(), "r");
    if (!fp) return false;
    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); return false; }
    hts_idx_t* idx = sam_index_load(fp, bams[bi].c_str());

    // disc mate needs, keyed by (anchor read name, mate's is_read1 = s_mate_first). File-local.
    std::map<std::pair<std::string,int>, std::vector<DiscNeed>> disc_needs;

    hts_itr_t* it = nullptr;
    if (!region.empty() && idx) it = sam_itr_querys(idx, hdr, region.c_str());
    else if (idx)               it = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    bam1_t* b = bam_init1();

    auto write_clip = [&](const char* chrom, long map_pos, char side, int is_rc,
                          const std::string& clip_seq, const std::string& clip_qual, long rs, long re) {
        std::vector<long> cands;
        cand.overlapping(chrom, rs, re, extnd, cands);
        if (cands.empty()) return;
        FILE* dst = ((int)clip_seq.size() < kBwaRealignCut) ? f_pa : f_non;  // cohort by length
        for (long ins_pos : cands) {
            long idx_cnt = ncnt[{std::string(chrom), ins_pos}]++;
            // @chrm~map_pos~L|R~is_rc~ins_pos~n_cnt~sample_id
            std::fprintf(dst, "@%s~%ld~%c~%d~%ld~%ld~%d\n%s\n+\n%s\n",
                         chrom, map_pos, side, is_rc, ins_pos, idx_cnt, sample_id,
                         clip_seq.c_str(), clip_qual.c_str());
        }
    };

    while (it && sam_itr_next(fp, it, b) >= 0) {
        const bam1_t* r = b;
        const uint16_t flag = r->core.flag;
        if (flag & (BAM_FDUP | BAM_FUNMAP)) continue;
        const uint32_t* cig = bam_get_cigar(r);
        const int ncig = r->core.n_cigar;
        if (ncig < 1) continue;
        const char* chrom = sam_hdr_tid2name(hdr, r->core.tid);
        const long rs = r->core.pos;

        // ---- background low-mapq clip count (mapq<12 gate; only mapq<5 counted) ----
        // Primary alignments only (2026-11-13; master also counted supp/sec, see
        // count_as_lowq_clip).
        if (r->core.qual < kMinClipMapq) {
            if (r->core.qual < kClipLowMapq && count_as_lowq_clip(flag)) {
                std::vector<long> nc;
                if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) {           // left (priority, master elif)
                    cand.near(chrom, rs, kNearbyLowQClip, nc);
                } else if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {  // right
                    long tmp = g_read_length - (long)bam_cigar_oplen(cig[ncig - 1]);
                    cand.near(chrom, rs + tmp, kNearbyLowQClip, nc);
                }
                for (long C : nc) lowq[{std::string(chrom), C}] += 1;
            }
            continue;                                                   // clip mapq gate (>=12)
        }
        if (is_two_side_clipped(cig, ncig, kTwoSideClipMin)) continue;

        const long re = bam_endpos(r);
        const int is_rc = bam_is_rev(r) ? 1 : 0;
        const bool b_first = !(flag & BAM_FREAD2);      // master: b_first = not is_read2
        const bool is_supp_sec = (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0;
        const uint8_t* seq = bam_get_seq(r);
        const uint8_t* qual = bam_get_qual(r);
        const int qlen = r->core.l_qseq;
        const bool have_q = (qual && qual[0] != 0xff);

        // mate info (only set when the mate is mapped with a valid ref) -- master lines 130-139.
        std::string mate_chrm = "*"; long mate_pos = 0;
        if (!(flag & BAM_FMUNMAP) && r->core.mtid >= 0) {
            mate_chrm = sam_hdr_tid2name(hdr, r->core.mtid);
            mate_pos  = r->core.mpos;
        }
        const int is_mate_rc = (flag & BAM_FMREVERSE) ? 1 : 0;

        long map_pos = rs;         // anchor pos; the right-clip block advances it (reused for disc)
        int b_left_clip = 0, b_right_clip = 0;

        // ---- LEFT clip: cigar[0] == S ----
        // NB: master processes left THEN right within one read; a gate failure on EITHER side is a
        // `continue` that skips the WHOLE read (so a read whose left clip fails never contributes its
        // right clip OR its disc record). We mirror that with `continue`.
        if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) {
            if (is_supp_sec) continue;                              // clip block skips supp/sec (whole read)
            const int clen = bam_cigar_oplen(cig[0]);
            sseq.clear();
            for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, i)];
            bool pass = true;                                      // master left: 7<=len<=9 need polyA; elif <9 skip
            if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
            else if (clen < kBwaRealignCut) pass = false;
            if (!pass) continue;
            if (have_q && !is_qualified_clip(qual, clen)) continue;
            if (want_clip) {
                squal.clear();
                for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[i] + 33) : 'I';
                write_clip(chrom, rs, 'L', is_rc, sseq, squal, rs, re);
            }
            b_left_clip = 1;
        }
        // ---- RIGHT clip: cigar[-1] == S ----
        if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {
            for (int i = 0; i < ncig - 1; ++i) {                   // advance map_pos to the right breakpoint
                const int op = bam_cigar_op(cig[i]);
                if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CINS) continue;
                map_pos += bam_cigar_oplen(cig[i]);
            }
            if (is_supp_sec) continue;                             // (master checks after advancing map_pos)
            const int clen = bam_cigar_oplen(cig[ncig - 1]);
            const int off = qlen - clen;
            sseq.clear();
            for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, off + i)];
            bool pass = true;                                      // master right: 7<=len<=9 need polyA; then <9 skip
            if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
            if (clen < kBwaRealignCut) pass = false;
            if (!pass) continue;
            if (have_q && !is_qualified_clip(qual + off, clen)) continue;
            if (want_clip) {
                squal.clear();
                for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[off + i] + 33) : 'I';
                write_clip(chrom, map_pos, 'R', is_rc, sseq, squal, rs, re);
            }
            b_right_clip = 1;
        }

        // ---- DISC anchor ---- (master lines 235-263)
        if (!want_disc) continue;
        if (mate_chrm == "*") continue;                            // mate unmapped
        if (is_decoy_contig(mate_chrm)) continue;
        if (r->core.qual < kMinDiscMapq) continue;                 // anchor mapq >= 20
        // is_discordant: different chrom, OR same chrom with |mate_pos - map_pos| > threshold.
        bool disc = (std::string(chrom) != mate_chrm) || (llabs(mate_pos - map_pos) > kDiscThreshold);
        if (!disc) continue;
        const int s_mate_first = b_first ? 0 : 1;                  // is_read1 of the MATE
        std::vector<long> cands;
        cand.overlapping(chrom, rs, re, extnd, cands);
        if (cands.empty()) continue;
        auto& vec = disc_needs[{std::string(bam_get_qname(r)), s_mate_first}];
        for (long ins_pos : cands)
            vec.push_back({std::string(chrom), ins_pos, map_pos, b_left_clip, b_right_clip,
                           is_rc, is_mate_rc, mate_chrm, mate_pos});
    }
    bam_destroy1(b);
    if (it) hts_itr_destroy(it);

    // ---- DISC pass 2 (this file): gather each mate's full seq by (name,is_read1) -> pooled FASTA --
    if (want_disc) {
        hts_itr_t* it2 = nullptr;
        if (!region.empty() && idx) it2 = sam_itr_querys(idx, hdr, region.c_str());
        else if (idx)               it2 = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
        bam1_t* b2 = bam_init1();
        std::string mseq;
        while (it2 && sam_itr_next(fp, it2, b2) >= 0) {
            const bam1_t* r = b2;
            const int is_r1 = (r->core.flag & BAM_FREAD1) ? 1 : 0;
            auto nit = disc_needs.find({std::string(bam_get_qname(r)), is_r1});
            if (nit == disc_needs.end()) continue;
            const char* rchrom = sam_hdr_tid2name(hdr, r->core.tid);
            const uint8_t* s = bam_get_seq(r); const int L = r->core.l_qseq;
            mseq.clear();
            for (int i = 0; i < L; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
            for (const auto& d : nit->second) {
                // The mate must map on the anchor's mate chrom (master fetches mate at (mate_chrm,
                // mate_pos)); this rejects same-name supplementary alignments that land elsewhere
                // when the true mate is off this contig (e.g. a cross-chrom mate absent from a
                // single-chrom BAM).
                if (d.mate_chrm != rchrom) continue;
                // >rname~is_first~b_lclip~b_rclip~anchor_rc~mate_rc~anchor_map_pos~ins_chrm~ins_pos~sample
                std::fprintf(fa, ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%d\n%s\n",
                             bam_get_qname(r), is_r1, d.b_lclip, d.b_rclip, d.anchor_rc, d.mate_rc,
                             d.anchor_map_pos, d.ins_chrm.c_str(), d.ins_pos, sample_id, mseq.c_str());
            }
        }
        bam_destroy1(b2);
        if (it2) hts_itr_destroy(it2);
    }

    if (idx) hts_idx_destroy(idx);
    sam_hdr_destroy(hdr); sam_close(fp);
    }   // end per-input-file loop

    if (f_non) std::fclose(f_non);
    if (f_pa)  std::fclose(f_pa);
    if (fa)    std::fclose(fa);

    const std::string tstr = std::to_string(bwa_threads);
    if (!emit_cns_lowq_file(out_lowq, cand, lowq)) return false;
    if (want_clip && !realign_clip_cohorts(bwa_bin, consensus, fq_non, fq_pa, tstr,
                                           out_clip_sam, tmp_prefix + ".clip")) return false;
    if (want_clip) { ::unlink(fq_non.c_str()); ::unlink(fq_pa.c_str()); }
    if (want_disc) {
        // realign_disc_reads: plain `bwa mem` (bwa_align.py:120) -> temp_disc.sam
        if (!run_bwa_sam(bwa_bin, consensus, fa_disc, {"-t",tstr,"-K",kBwaBatch}, out_disc_sam))
            return false;
        ::unlink(fa_disc.c_str());
    }
    return true;
}

// ==================== fold-in: collect ONCE, emit for ALL types ====================
//
// The per-type emit_cns above re-reads the BAM twice PER TYPE (6 whole-genome passes for L1/Alu/SVA).
// But the clip/disc collection gates are type- and candidate-agnostic: only the fan-out target
// (per-type candidate list) and the realign reference (per-type consensus) differ. So we collect the
// qualifying reads ONCE (2 BAM passes total, against the UNION of all types' candidate windows), hold
// them, then per type fan the held reads to that type's candidates and realign to its consensus.
// Net: 8 whole-genome passes (1 detect + 1 detect-mate + 6 cns) -> 4 (1+1 detect + 1+1 cns).

// A collected clip part, held candidate-agnostic (fanned to each type's candidates later).
struct HeldClip {
    std::string chrm; long map_pos; char side; int is_rc;
    long ref_start, ref_end; std::string seq, qual;
};
// One disc anchor occurrence, held candidate-agnostic. mate seq is filled by the gather pass.
// mate_pos is the gather TARGET (the mate's reference start, from the anchor's mpos). It is what
// lets the spill-driven gather fetch mates by position instead of scanning the genome for names;
// emit_cns/emit_cns_multi never needed it because they gather in the same pass that built the map.
struct HeldAnchor {
    std::string ins_chrm; long anchor_map_pos; int b_lclip, b_rclip, anchor_rc, mate_rc;
    long ref_start, ref_end; std::string mate_chrm; long mate_pos = -1;
};
// One MEI type's fan-out target + outputs.
struct CnsTypeJob {
    std::string tag, candidate_list, consensus, out_clip_sam, out_disc_sam, out_lowq;
    Candidates cand;
};

inline bool emit_cns_multi(const std::vector<std::string>& bams, const std::string& reference,
                           std::vector<CnsTypeJob>& jobs, long extnd, int sample_id,
                           const std::string& bwa_bin, int bwa_threads,
                           const std::string& region, const std::string& tmp_prefix) {
    // `region` is unused now and that is the intended consequence of the port: the CANDIDATES
    // define what gets read. Restricting a whole-genome stream was only ever a way to bound its
    // cost, and a candidate window that fell outside the region would have been silently
    // truncated. Kept in the signature so it still matches emit_cns.
    (void)region;

    // union of all candidate positions, so overlapping windows are read once, not once per type.
    Candidates uni;
    for (auto& j : jobs)
        for (auto& kv : j.cand.by_chrom) {
            auto& v = uni.by_chrom[kv.first];
            v.insert(v.end(), kv.second.begin(), kv.second.end());
        }
    for (auto& kv : uni.by_chrom) {
        std::sort(kv.second.begin(), kv.second.end());
        kv.second.erase(std::unique(kv.second.begin(), kv.second.end()), kv.second.end());
    }

    // ---- Fetch plan: the candidate windows, merged, with an ownership boundary each ----
    //
    // This is the whole point of the design, so it is worth stating plainly. Master collects by
    // looping over candidates and doing ONE indexed fetch per candidate window
    // ([pos-extnd, pos+extnd), x_clip_disc_filter.py:108). That touches ~4 Mb of reference for
    // ~2,800 candidates instead of all 3.1 Gb, and no amount of C++ makes a whole-genome stream
    // competitive with that -- it is the access pattern, not the language.
    //
    // We fetch the same windows, but MERGED, so overlapping candidate windows are read once
    // rather than once per candidate. Merging needs an ownership rule, because htslib returns
    // every read OVERLAPPING a region: without one, a read spanning the gap between two merged
    // regions is returned by both and collected twice.
    //
    // The rule: a read belongs to the region it is fetched by iff its START is at or after
    // `own_from`, the previous region's end on that contig (0 for the first). That is exact, not
    // approximate. Regions only fail to merge when the next window starts beyond the current
    // end, so beg > own_from always; a read with pos < own_from that reaches into this region
    // must span own_from, hence also overlapped the previous region and was collected there.
    // Reads starting in the gap (own_from <= pos < beg) are collected here. So every qualifying
    // read is collected exactly once, and each region's boundary depends only on its
    // predecessor -- which makes the regions independent, so they parallelise.
    //
    // Fan-out is against the FULL union candidate list, not the region's own candidates, so
    // membership is unaffected by where the merge boundaries happen to fall.
    struct FetchRegion { std::string chrom; long beg, end, own_from; };
    std::vector<FetchRegion> fetch;
    for (const auto& kv : uni.by_chrom) {
        const auto& v = kv.second;              // sorted, unique
        long prev_end = 0;
        for (size_t i = 0; i < v.size(); ) {
            long b = std::max(0L, v[i] - extnd), e = v[i] + extnd;
            while (++i < v.size() && (v[i] - extnd) <= e) e = std::max(e, v[i] + extnd);
            fetch.push_back({kv.first, b, e, prev_end});
            prev_end = e;
        }
    }

    std::vector<HeldClip> clips;
    std::map<std::pair<std::string,int>, std::vector<HeldAnchor>> anchors;
    std::map<std::pair<std::string,long>, long> lowq;   // union candidate -> low-mapq clip count

    // Multi-file: both phases run per input file. `clips` and `lowq` are positional and pool
    // trivially. `anchors`/`mate_seq` are qname-keyed, so their keys carry a "<file_idx>:" tag --
    // two inputs sharing a qname therefore cannot cross-pair, and the mate gather for file N only
    // resolves anchors that file N's collection created. The tag is stripped when the disc FASTA
    // is written, so emitted records are byte-identical to a single-file run.
    std::map<std::pair<std::string,int>, std::string> mate_seq;
    for (size_t bi = 0; bi < bams.size(); ++bi) {
    const std::string ftag = std::to_string(bi) + ":";

    // ---- Phase 1: collect clips + disc anchors from the candidate windows, in parallel ----
    //
    // Per-thread accumulators merged after, so nothing is locked on the hot path. Each thread
    // takes a contiguous stripe of the coordinate-sorted region list so its fetches stay in
    // increasing order and reuse the BGZF block it just decompressed. One samFile + index handle
    // per thread: htslib file handles are not thread-safe.
    {
        const int nthreads = std::max(1, std::min<int>(bwa_threads, (int)std::max<size_t>(fetch.size(), 1)));
        struct Part {
            std::vector<HeldClip> clips;
            std::map<std::pair<std::string,int>, std::vector<HeldAnchor>> anchors;
            std::map<std::pair<std::string,long>, long> lowq;
            bool ok = true;
        };
        std::vector<std::future<Part>> futs;
        for (int ti = 0; ti < nthreads; ++ti) {
            futs.push_back(std::async(std::launch::async, [&, ti]() {
                Part P;
                samFile* fp = sam_open(bams[bi].c_str(), "r");
                if (!fp) { P.ok = false; return P; }
                if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
                sam_hdr_t* hdr = sam_hdr_read(fp);
                if (!hdr) { sam_close(fp); P.ok = false; return P; }
                hts_idx_t* idx = sam_index_load(fp, bams[bi].c_str());
                if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); P.ok = false; return P; }
                bam1_t* b = bam_init1();
                std::string sseq, squal;
                const size_t lo = fetch.size() * ti / nthreads;
                const size_t hi = fetch.size() * (ti + 1) / nthreads;
                for (size_t ri = lo; ri < hi; ++ri) {
                    const int tid = sam_hdr_name2tid(hdr, fetch[ri].chrom.c_str());
                    if (tid < 0) continue;
                    hts_itr_t* it = sam_itr_queryi(idx, tid, fetch[ri].beg, fetch[ri].end);
                    if (!it) continue;
                    while (sam_itr_next(fp, it, b) >= 0) {
                        const bam1_t* r = b;
                        if (r->core.pos < fetch[ri].own_from) continue;   // owned by an earlier region
                        const uint16_t flag = r->core.flag;
                        if (flag & (BAM_FDUP | BAM_FUNMAP)) continue;
                        const uint32_t* cig = bam_get_cigar(r);
                        const int ncig = r->core.n_cigar;
                        if (ncig < 1) continue;
                        const char* chrom = sam_hdr_tid2name(hdr, r->core.tid);
                        const long rs = r->core.pos;

                        if (r->core.qual < kMinClipMapq) {          // lowq background, primary only
                            if (r->core.qual < kClipLowMapq && count_as_lowq_clip(flag)) {
                                std::vector<long> nc;
                                if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP)
                                    uni.near(chrom, rs, kNearbyLowQClip, nc);
                                else if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {
                                    long tmp = g_read_length - (long)bam_cigar_oplen(cig[ncig - 1]);
                                    uni.near(chrom, rs + tmp, kNearbyLowQClip, nc);
                                }
                                for (long C : nc) P.lowq[{std::string(chrom), C}] += 1;
                            }
                            continue;
                        }
                        if (is_two_side_clipped(cig, ncig, kTwoSideClipMin)) continue;

                        const long re = bam_endpos(r);
                        const int is_rc = bam_is_rev(r) ? 1 : 0;
                        const bool b_first = !(flag & BAM_FREAD2);
                        const bool is_supp_sec = (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0;
                        const uint8_t* seq = bam_get_seq(r);
                        const uint8_t* qual = bam_get_qual(r);
                        const int qlen = r->core.l_qseq;
                        const bool have_q = (qual && qual[0] != 0xff);
                        std::string mate_chrm = "*"; long mate_pos = 0;
                        if (!(flag & BAM_FMUNMAP) && r->core.mtid >= 0) {
                            mate_chrm = sam_hdr_tid2name(hdr, r->core.mtid); mate_pos = r->core.mpos;
                        }
                        const int is_mate_rc = (flag & BAM_FMREVERSE) ? 1 : 0;
                        const bool overlaps_union = [&]{
                            std::vector<long> c; uni.overlapping(chrom, rs, re, extnd, c); return !c.empty();
                        }();

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
                            if (!pass) continue;
                            if (have_q && !is_qualified_clip(qual, clen)) continue;
                            if (overlaps_union) {
                                squal.clear();
                                for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[i] + 33) : 'I';
                                P.clips.push_back({chrom, rs, 'L', is_rc, rs, re, sseq, squal});
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
                            bool pass = true;
                            if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
                            if (clen < kBwaRealignCut) pass = false;
                            if (!pass) continue;
                            if (have_q && !is_qualified_clip(qual + off, clen)) continue;
                            if (overlaps_union) {
                                squal.clear();
                                for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[off + i] + 33) : 'I';
                                P.clips.push_back({chrom, map_pos, 'R', is_rc, rs, re, sseq, squal});
                            }
                            b_right_clip = 1;
                        }

                        // disc anchor
                        if (mate_chrm == "*") continue;
                        if (is_decoy_contig(mate_chrm)) continue;
                        if (r->core.qual < kMinDiscMapq) continue;
                        bool disc = (std::string(chrom) != mate_chrm) || (llabs(mate_pos - map_pos) > kDiscThreshold);
                        if (!disc || !overlaps_union) continue;
                        const int s_mate_first = b_first ? 0 : 1;
                        P.anchors[{ftag + bam_get_qname(r), s_mate_first}].push_back(
                            {std::string(chrom), map_pos, b_left_clip, b_right_clip, is_rc, is_mate_rc,
                             rs, re, mate_chrm, mate_pos});
                    }
                    hts_itr_destroy(it);
                }
                bam_destroy1(b);
                hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
                return P;
            }));
        }
        // Merge. Clip ORDER matters: the per-type fan writes the FASTQ in `clips` order, and a
        // candidate's n_cnt indices -- and therefore bwa's input order, whose tie-break among
        // equal-scoring consensus placements is input-index dependent -- follow it. Stripes are
        // contiguous and taken in ascending region order, so concatenating them in thread order
        // reproduces coordinate order, which is the order a per-candidate fetch sees.
        for (auto& f : futs) {
            Part P = f.get();
            if (!P.ok) return false;
            clips.insert(clips.end(), std::make_move_iterator(P.clips.begin()),
                         std::make_move_iterator(P.clips.end()));
            for (auto& kv : P.anchors) {
                auto& dst = anchors[kv.first];
                dst.insert(dst.end(), std::make_move_iterator(kv.second.begin()),
                           std::make_move_iterator(kv.second.end()));
            }
            for (auto& kv : P.lowq) lowq[kv.first] += kv.second;
        }
    }

    // ---- Phase 2: gather each needed mate's PRIMARY seq, by POSITION, in parallel ----
    //
    // Disc mates are scattered genome-wide, so this half does NOT become sparse just because
    // collection did -- it is the one part that has to touch the whole genome. Master handles it
    // by splitting the genome into 50 Mb bins under a process Pool (x_alignments.py:407); we fetch
    // the mate positions the anchors recorded, merged into coordinate-sorted regions and striped
    // across threads. Same shape, without the fork or the per-bin temp files.
    {
        std::map<std::string, std::vector<long>> want;
        for (const auto& kv : anchors) {
            if (kv.first.first.compare(0, ftag.size(), ftag) != 0) continue;
            for (const auto& o : kv.second)
                if (o.mate_pos >= 0) want[o.mate_chrm].push_back(o.mate_pos);
        }
        constexpr long kMergeSlack = 4096;   // coalesce targets closer than this into one fetch
        std::vector<std::pair<std::string, std::pair<long,long>>> mregions;
        for (auto& kv : want) {
            auto& v = kv.second;
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            for (size_t i = 0; i < v.size(); ) {
                long beg = v[i], end = v[i] + 1;
                while (++i < v.size() && v[i] - end <= kMergeSlack) end = v[i] + 1;
                mregions.push_back({kv.first, {beg, end}});
            }
        }
        const int nthreads = std::max(1, std::min<int>(bwa_threads, (int)std::max<size_t>(mregions.size(), 1)));
        std::vector<std::future<std::map<std::pair<std::string,int>, std::string>>> futs;
        for (int ti = 0; ti < nthreads && !mregions.empty(); ++ti) {
            futs.push_back(std::async(std::launch::async, [&, ti]() {
                std::map<std::pair<std::string,int>, std::string> local;
                samFile* fp = sam_open(bams[bi].c_str(), "r");
                if (!fp) return local;
                if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
                sam_hdr_t* hdr = sam_hdr_read(fp);
                if (!hdr) { sam_close(fp); return local; }
                hts_idx_t* idx = sam_index_load(fp, bams[bi].c_str());
                if (!idx) { sam_hdr_destroy(hdr); sam_close(fp); return local; }
                bam1_t* b = bam_init1();
                std::string mseq;
                const size_t lo = mregions.size() * ti / nthreads;
                const size_t hi = mregions.size() * (ti + 1) / nthreads;
                for (size_t ri = lo; ri < hi; ++ri) {
                    const int tid = sam_hdr_name2tid(hdr, mregions[ri].first.c_str());
                    if (tid < 0) continue;
                    hts_itr_t* it = sam_itr_queryi(idx, tid, mregions[ri].second.first,
                                                  mregions[ri].second.second);
                    if (!it) continue;
                    while (sam_itr_next(fp, it, b) >= 0) {
                        const bam1_t* r = b;
                        if (r->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
                        const int is_r1 = (r->core.flag & BAM_FREAD1) ? 1 : 0;
                        std::pair<std::string,int> key{ftag + bam_get_qname(r), is_r1};
                        auto ait = anchors.find(key);
                        if (ait == anchors.end() || local.count(key)) continue;
                        const char* rchrom = sam_hdr_tid2name(hdr, r->core.tid);
                        bool chrom_ok = false;
                        for (auto& occ : ait->second) if (occ.mate_chrm == rchrom) { chrom_ok = true; break; }
                        if (!chrom_ok) continue;
                        const uint8_t* s = bam_get_seq(r); const int L = r->core.l_qseq;
                        mseq.clear();
                        for (int i = 0; i < L; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
                        local[key] = mseq;
                    }
                    hts_itr_destroy(it);
                }
                bam_destroy1(b);
                hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp);
                return local;
            }));
        }
        // First writer wins, matching the old single-pass "skip if already present". Order-
        // independent regardless: there is exactly one primary alignment per (qname, is_read1).
        for (auto& f : futs) for (auto& kv : f.get()) mate_seq.emplace(kv.first, std::move(kv.second));
    }
    }   // end per-input-file loop

    const std::string tstr = std::to_string(bwa_threads);

    // Index the held clips by contig, ordered by reference start, so the emit below can walk a
    // single candidate's window without rescanning every clip. max_span bounds how far back a
    // clip can start and still reach a given position, which is what makes the lower bound of
    // that walk a binary search rather than a scan from the beginning.
    std::map<std::string, std::vector<const HeldClip*>> clip_ix;
    long max_span = 1;
    for (const auto& c : clips) {
        clip_ix[c.chrm].push_back(&c);
        max_span = std::max(max_span, c.ref_end - c.ref_start);
    }
    for (auto& kv : clip_ix)
        std::stable_sort(kv.second.begin(), kv.second.end(),
                         [](const HeldClip* a, const HeldClip* b) { return a->ref_start < b->ref_start; });

    // ---- Per type: fan held reads to that type's candidates, realign, emit ----
    for (size_t ti = 0; ti < jobs.size(); ++ti) {
        CnsTypeJob& J = jobs[ti];
        const std::string pfx = tmp_prefix + "." + std::to_string(ti);
        // CLIP
        if (!J.out_clip_sam.empty()) {
            const std::string fq_non = pfx + ".clip.non.fq", fq_pa = pfx + ".clip.pa.fq";
            FILE* f_non = std::fopen(fq_non.c_str(), "w");
            FILE* f_pa  = std::fopen(fq_pa.c_str(), "w");
            // CANDIDATE-MAJOR, and that nesting is load-bearing rather than incidental.
            //
            // Master emits one candidate's window at a time (a fetch per candidate, then the next),
            // so its FASTQ is grouped BY CANDIDATE. Fanning read-major instead groups it BY READ,
            // which is a different order for the same records -- and order is not cosmetic here:
            // bwa's tie-break among EQUAL-SCORING consensus placements depends on a read's index in
            // its input, so a reordered FASTQ silently moves which copy multi-mapping clip parts
            // land on, and with it a few borderline calls. Measured read-major on chr1: identical
            // read sets but 570/2021/1751 clip records placed differently, and one extra SVA call.
            //
            // Candidates are iterated in Candidates order, which is position-sorted. That matches
            // master because master's per-candidate files are concatenated in candidate-FILE order
            // and the candidate lists are already position-sorted (verified on all three types).
            //
            // n_cnt (field 5) is NOT reproduced exactly: master builds the right-clip read name and
            // increments the counter BEFORE applying that side's polyA/length/quality gates
            // (x_clip_disc_filter.py:205-218), so a right clip that fails burns an index and leaves
            // a gap. parse_clip_realignment_consensus reads fields 0,1,2,4,6 and never field 5, so
            // the values are free -- only their ORDER matters, and that now matches. Reproducing
            // the gaps would mean collecting reads solely to consume indices.
            std::vector<const HeldClip*> hits;
            for (const auto& kv : J.cand.by_chrom) {
                auto ixit = clip_ix.find(kv.first);
                if (ixit == clip_ix.end()) continue;
                const auto& v = ixit->second;
                for (long C : kv.second) {
                    // clip overlaps candidate C iff ref_start < C+extnd AND ref_end > C-extnd
                    auto lo = std::lower_bound(v.begin(), v.end(), C - extnd - max_span,
                        [](const HeldClip* a, long x) { return a->ref_start < x; });
                    auto hi = std::lower_bound(v.begin(), v.end(), C + extnd,
                        [](const HeldClip* a, long x) { return a->ref_start < x; });
                    hits.clear();
                    for (auto p = lo; p != hi; ++p) if ((*p)->ref_end > C - extnd) hits.push_back(*p);
                    long nc = 0;               // per candidate, as master's n_cnt_clip is
                    for (const HeldClip* c : hits) {
                        FILE* dst = ((int)c->seq.size() < kBwaRealignCut) ? f_pa : f_non;
                        std::fprintf(dst, "@%s~%ld~%c~%d~%ld~%ld~%d\n%s\n+\n%s\n",
                                     c->chrm.c_str(), c->map_pos, c->side, c->is_rc, C, nc++, sample_id,
                                     c->seq.c_str(), c->qual.c_str());
                    }
                }
            }
            std::fclose(f_non); std::fclose(f_pa);
            if (!realign_clip_cohorts(bwa_bin, J.consensus, fq_non, fq_pa, tstr, J.out_clip_sam, pfx + ".clip"))
                return false;
            ::unlink(fq_non.c_str()); ::unlink(fq_pa.c_str());
        }
        // DISC
        if (!J.out_disc_sam.empty()) {
            const std::string fa = pfx + ".disc.fa";
            FILE* faf = std::fopen(fa.c_str(), "w");
            std::vector<long> cs;
            for (const auto& kv : anchors) {
                auto sit = mate_seq.find(kv.first);
                if (sit == mate_seq.end()) continue;
                const int s_mate_first = kv.first.second;
                // Strip the file-index tag -- emitted qname must match master byte for byte.
                const std::string qn = strip_file_tag(kv.first.first);
                for (const auto& occ : kv.second) {
                    cs.clear();
                    J.cand.overlapping(occ.ins_chrm, occ.ref_start, occ.ref_end, extnd, cs);
                    for (long ins_pos : cs)
                        std::fprintf(faf, ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%d\n%s\n",
                                     qn.c_str(), s_mate_first, occ.b_lclip, occ.b_rclip,
                                     occ.anchor_rc, occ.mate_rc, occ.anchor_map_pos, occ.ins_chrm.c_str(),
                                     ins_pos, sample_id, sit->second.c_str());
                }
            }
            std::fclose(faf);
            if (!run_bwa_sam(bwa_bin, J.consensus, fa, {"-t",tstr,"-K",kBwaBatch}, J.out_disc_sam)) return false;
            ::unlink(fa.c_str());
        }
        // LOWQ (per type: this type's candidates, counts from the union scan)
        if (!J.out_lowq.empty()) {
            if (FILE* fl = std::fopen(J.out_lowq.c_str(), "w")) {
                for (const auto& kv : J.cand.by_chrom)
                    for (long pos : kv.second) {
                        auto it2 = lowq.find({kv.first, pos});
                        std::fprintf(fl, "%s %ld %ld\n", kv.first.c_str(), pos, it2 == lowq.end() ? 0L : it2->second);
                    }
                std::fclose(fl);
            }
        }
    }
    return true;
}

// The "full fold-in" spill path lived here and is DELETED (2026-07-28).
//
// It collected the cns reads during the detection scan -- before candidates existed -- spilled
// them to three text files (1.10 GB on NA12878 WGS), and fanned them to candidates afterwards.
// The point was to reuse the detection pass and cut whole-genome BAM passes from 8 to 4.
//
// It was removed because emit_cns_multi above now does the same job by fetching the candidate
// windows directly, which is strictly better on all three axes the fold-in traded against:
// no spill (vs 1.10 GB), no genome-wide read set held in memory (vs ~4.15M anchors, of which
// ~0.5M were ever used), and ~4 Mb of reference touched instead of 3.1 Gb. Collecting before
// candidates are known is what forced all three; knowing them first removes the need.
//
// Worth keeping in mind if a future step wants candidate-agnostic collection again: the reason
// this was slow was not the missing index. It was that a candidate-agnostic collector cannot
// know what to discard, so it pays for everything. Git history has the implementation.

}  // namespace cns
}  // namespace xtea
