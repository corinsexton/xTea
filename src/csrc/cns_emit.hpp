// cns_emit.hpp — move B.1: emit the cns_filter temp SAMs (temp_clip.sam / temp_disc.sam) from
// the scan pass, so x_clip_disc_filter's collect+realign (a per-candidate BAM re-read + bwa re-run)
// can be bypassed. This is a faithful port of collect_clipped_disc_reads_by_region
// (x_clip_disc_filter.py:81) + realign_clipped_read_with_polyA / realign_disc_reads (bwa_align.py),
// verified against ACTUAL code/constants (docs/phase_c_unified_downstream.md §7.1) and the golden's
// preserved temp SAMs (docs/cns_realign_sam_contract.md).
//
// Difference from master's per-candidate collection: we stream the BAM ONCE (candidate-agnostic),
// then fan each qualifying read out to every candidate whose window [ins_pos-extnd, ins_pos+extnd]
// overlaps the read's reference span — reproducing pysam.fetch's per-candidate membership without the
// N-fetch BAM re-read. bwa realign params + read-name formats are byte-identical to master, so the
// emitted SAMs feed straight into parse_clip_realignment_consensus / parse_disc_algnmt_consensus.

#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
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
constexpr int  kReadLength     = 100;  // READ_LENGTH (default; used only in lowq right-clip pos)

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
inline bool emit_cns(const std::string& bam, const std::string& reference,
                     const Candidates& cand, long extnd, int sample_id,
                     const std::string& consensus, const std::string& bwa_bin,
                     int bwa_threads, const std::string& region,
                     const std::string& out_clip_sam, const std::string& out_disc_sam,
                     const std::string& out_lowq, const std::string& tmp_prefix) {
    const bool want_clip = !out_clip_sam.empty();
    const bool want_disc = !out_disc_sam.empty();
    samFile* fp = sam_open(bam.c_str(), "r");
    if (!fp) return false;
    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); return false; }
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());

    // Two cohort FASTQs (bwa params differ): non-polyA (clip len >= 9), polyA (len < 9). Only if clip.
    const std::string fq_non = tmp_prefix + ".clip.non_polyA.fq";
    const std::string fq_pa  = tmp_prefix + ".clip.polyA.fq";
    FILE* f_non = want_clip ? std::fopen(fq_non.c_str(), "w") : nullptr;
    FILE* f_pa  = want_clip ? std::fopen(fq_pa.c_str(), "w") : nullptr;
    if (want_clip && (!f_non || !f_pa)) { if (f_non) std::fclose(f_non); if (f_pa) std::fclose(f_pa);
                           if (idx) hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp); return false; }

    // per-candidate running clip index (n_cnt); value is ignored by the parser (field 5) but kept
    // for format fidelity / uniqueness.
    std::map<std::pair<std::string,long>, long> ncnt;
    std::map<std::pair<std::string,long>, long> lowq;  // per-candidate background low-mapq clip count
    // disc mate needs, keyed by (anchor read name, mate's is_read1 = s_mate_first).
    std::map<std::pair<std::string,int>, std::vector<DiscNeed>> disc_needs;
    std::string sseq, squal;  // scratch

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
        // master counts these BEFORE the supp/sec skip, so supp/sec low-q clips ARE counted.
        if (r->core.qual < kMinClipMapq) {
            if (r->core.qual < kClipLowMapq) {
                std::vector<long> nc;
                if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) {           // left (priority, master elif)
                    cand.near(chrom, rs, kNearbyLowQClip, nc);
                } else if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {  // right
                    long tmp = kReadLength - (long)bam_cigar_oplen(cig[ncig - 1]);
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
    if (f_non) std::fclose(f_non);
    if (f_pa)  std::fclose(f_pa);

    const std::string tstr = std::to_string(bwa_threads);

    // low_mapq_clip_cnt.txt — one line per candidate ("chrm pos count"), matching master's
    // m_lowq_clip_cnt (built with a 0 default for every candidate). The bypass needs this file.
    if (!out_lowq.empty()) {
        if (FILE* fl = std::fopen(out_lowq.c_str(), "w")) {
            for (const auto& kv : cand.by_chrom)
                for (long pos : kv.second) {
                    auto it2 = lowq.find({kv.first, pos});
                    std::fprintf(fl, "%s %ld %ld\n", kv.first.c_str(), pos, it2 == lowq.end() ? 0L : it2->second);
                }
            std::fclose(fl);
        }
    }

    // ---- CLIP realign: two cohorts (non-polyA -T9 -k9; polyA -T7 -k7 -c70), merge -> temp_clip.sam ----
    if (want_clip) {
        const std::string sam_non = tmp_prefix + ".clip.non_polyA.sam";
        const std::string sam_pa  = tmp_prefix + ".clip.polyA.sam";
        if (!run_bwa_sam(bwa_bin, consensus, fq_non, {"-t",tstr,"-T","9","-k","9"}, sam_non)) return false;
        if (!run_bwa_sam(bwa_bin, consensus, fq_pa,  {"-t",tstr,"-T","7","-k","7","-c","70"}, sam_pa)) return false;
        FILE* out = std::fopen(out_clip_sam.c_str(), "w");
        if (!out) return false;
        if (FILE* in = std::fopen(sam_non.c_str(), "r")) {         // full (headers + records)
            char* line = nullptr; size_t cap = 0; ssize_t len;
            while ((len = getline(&line, &cap, in)) != -1) std::fwrite(line, 1, len, out);
            free(line); std::fclose(in);
        }
        if (FILE* in = std::fopen(sam_pa.c_str(), "r")) {          // records only (drop @ headers)
            char* line = nullptr; size_t cap = 0; ssize_t len;
            while ((len = getline(&line, &cap, in)) != -1) { if (len > 0 && line[0] == '@') continue; std::fwrite(line, 1, len, out); }
            free(line); std::fclose(in);
        }
        std::fclose(out);
        ::unlink(fq_non.c_str()); ::unlink(fq_pa.c_str());
        ::unlink(sam_non.c_str()); ::unlink(sam_pa.c_str());
    }

    // ---- DISC: gather each mate's full seq by (name,is_read1), realign to consensus -> temp_disc.sam --
    if (want_disc) {
        const std::string fa_disc = tmp_prefix + ".disc.fa";
        FILE* fa = std::fopen(fa_disc.c_str(), "w");
        if (!fa) { if (idx) hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp); return false; }
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
        std::fclose(fa);
        // realign_disc_reads: plain `bwa mem` (bwa_align.py:120) -> temp_disc.sam
        if (!run_bwa_sam(bwa_bin, consensus, fa_disc, {"-t",tstr}, out_disc_sam)) {
            if (idx) hts_idx_destroy(idx); sam_hdr_destroy(hdr); sam_close(fp); return false;
        }
        ::unlink(fa_disc.c_str());
    }

    if (idx) hts_idx_destroy(idx);
    sam_hdr_destroy(hdr); sam_close(fp);
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
struct HeldAnchor {
    std::string ins_chrm; long anchor_map_pos; int b_lclip, b_rclip, anchor_rc, mate_rc;
    long ref_start, ref_end; std::string mate_chrm;
};
// One MEI type's fan-out target + outputs.
struct CnsTypeJob {
    std::string tag, candidate_list, consensus, out_clip_sam, out_disc_sam, out_lowq;
    Candidates cand;
};

// Realign a clip cohort pair (non-polyA -T9 -k9; polyA -T7 -k7 -c70) and merge -> out_sam.
inline bool realign_clip_cohorts(const std::string& bwa_bin, const std::string& consensus,
                                 const std::string& fq_non, const std::string& fq_pa,
                                 const std::string& tstr, const std::string& out_sam,
                                 const std::string& tmp_prefix) {
    const std::string sam_non = tmp_prefix + ".non_polyA.sam";
    const std::string sam_pa  = tmp_prefix + ".polyA.sam";
    if (!run_bwa_sam(bwa_bin, consensus, fq_non, {"-t",tstr,"-T","9","-k","9"}, sam_non)) return false;
    if (!run_bwa_sam(bwa_bin, consensus, fq_pa,  {"-t",tstr,"-T","7","-k","7","-c","70"}, sam_pa)) return false;
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

inline bool emit_cns_multi(const std::string& bam, const std::string& reference,
                           std::vector<CnsTypeJob>& jobs, long extnd, int sample_id,
                           const std::string& bwa_bin, int bwa_threads,
                           const std::string& region, const std::string& tmp_prefix) {
    // union of all candidate positions (for candidate-agnostic collection).
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

    samFile* fp = sam_open(bam.c_str(), "r");
    if (!fp) return false;
    if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
    sam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) { sam_close(fp); return false; }
    hts_idx_t* idx = sam_index_load(fp, bam.c_str());

    std::vector<HeldClip> clips;
    std::map<std::pair<std::string,int>, std::vector<HeldAnchor>> anchors;
    std::map<std::pair<std::string,long>, long> lowq;   // union candidate -> low-mapq clip count
    std::string sseq, squal;

    // ---- Pass 1: collect clips + disc anchors (candidate-agnostic vs the union) ----
    {
        hts_itr_t* it = nullptr;
        if (!region.empty() && idx) it = sam_itr_querys(idx, hdr, region.c_str());
        else if (idx)               it = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
        bam1_t* b = bam_init1();
        while (it && sam_itr_next(fp, it, b) >= 0) {
            const bam1_t* r = b;
            const uint16_t flag = r->core.flag;
            if (flag & (BAM_FDUP | BAM_FUNMAP)) continue;
            const uint32_t* cig = bam_get_cigar(r);
            const int ncig = r->core.n_cigar;
            if (ncig < 1) continue;
            const char* chrom = sam_hdr_tid2name(hdr, r->core.tid);
            const long rs = r->core.pos;

            if (r->core.qual < kMinClipMapq) {                       // lowq background (mapq<5), incl supp/sec
                if (r->core.qual < kClipLowMapq) {
                    std::vector<long> nc;
                    if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) uni.near(chrom, rs, kNearbyLowQClip, nc);
                    else if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {
                        long tmp = kReadLength - (long)bam_cigar_oplen(cig[ncig - 1]);
                        uni.near(chrom, rs + tmp, kNearbyLowQClip, nc);
                    }
                    for (long C : nc) lowq[{std::string(chrom), C}] += 1;
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
            const bool overlaps_union = [&]{ std::vector<long> c; uni.overlapping(chrom, rs, re, extnd, c); return !c.empty(); }();

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
                    clips.push_back({chrom, rs, 'L', is_rc, rs, re, sseq, squal});
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
                    clips.push_back({chrom, map_pos, 'R', is_rc, rs, re, sseq, squal});
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
            anchors[{std::string(bam_get_qname(r)), s_mate_first}].push_back(
                {std::string(chrom), map_pos, b_left_clip, b_right_clip, is_rc, is_mate_rc, rs, re, mate_chrm});
        }
        bam_destroy1(b);
        if (it) hts_itr_destroy(it);
    }

    // ---- Pass 2: gather each needed mate's PRIMARY seq by (name, is_read1) on its mate chrom ----
    std::map<std::pair<std::string,int>, std::string> mate_seq;
    {
        hts_itr_t* it = nullptr;
        if (!region.empty() && idx) it = sam_itr_querys(idx, hdr, region.c_str());
        else if (idx)               it = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
        bam1_t* b = bam_init1();
        std::string mseq;
        while (it && sam_itr_next(fp, it, b) >= 0) {
            const bam1_t* r = b;
            if (r->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;   // primary mate only
            const int is_r1 = (r->core.flag & BAM_FREAD1) ? 1 : 0;
            std::pair<std::string,int> key{std::string(bam_get_qname(r)), is_r1};
            auto ait = anchors.find(key);
            if (ait == anchors.end() || mate_seq.count(key)) continue;
            const char* rchrom = sam_hdr_tid2name(hdr, r->core.tid);
            bool chrom_ok = false;
            for (auto& occ : ait->second) if (occ.mate_chrm == rchrom) { chrom_ok = true; break; }
            if (!chrom_ok) continue;
            const uint8_t* s = bam_get_seq(r); const int L = r->core.l_qseq;
            mseq.clear();
            for (int i = 0; i < L; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
            mate_seq[key] = mseq;
        }
        bam_destroy1(b);
        if (it) hts_itr_destroy(it);
    }
    if (idx) hts_idx_destroy(idx);
    sam_hdr_destroy(hdr); sam_close(fp);

    const std::string tstr = std::to_string(bwa_threads);

    // ---- Per type: fan held reads to that type's candidates, realign, emit ----
    for (size_t ti = 0; ti < jobs.size(); ++ti) {
        CnsTypeJob& J = jobs[ti];
        const std::string pfx = tmp_prefix + "." + std::to_string(ti);
        // CLIP
        if (!J.out_clip_sam.empty()) {
            const std::string fq_non = pfx + ".clip.non.fq", fq_pa = pfx + ".clip.pa.fq";
            FILE* f_non = std::fopen(fq_non.c_str(), "w");
            FILE* f_pa  = std::fopen(fq_pa.c_str(), "w");
            std::map<std::pair<std::string,long>, long> ncnt;
            std::vector<long> cs;
            for (const auto& c : clips) {
                cs.clear();
                J.cand.overlapping(c.chrm, c.ref_start, c.ref_end, extnd, cs);
                if (cs.empty()) continue;
                FILE* dst = ((int)c.seq.size() < kBwaRealignCut) ? f_pa : f_non;
                for (long ins_pos : cs) {
                    long nc = ncnt[{c.chrm, ins_pos}]++;
                    std::fprintf(dst, "@%s~%ld~%c~%d~%ld~%ld~%d\n%s\n+\n%s\n",
                                 c.chrm.c_str(), c.map_pos, c.side, c.is_rc, ins_pos, nc, sample_id,
                                 c.seq.c_str(), c.qual.c_str());
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
                for (const auto& occ : kv.second) {
                    cs.clear();
                    J.cand.overlapping(occ.ins_chrm, occ.ref_start, occ.ref_end, extnd, cs);
                    for (long ins_pos : cs)
                        std::fprintf(faf, ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%d\n%s\n",
                                     kv.first.first.c_str(), s_mate_first, occ.b_lclip, occ.b_rclip,
                                     occ.anchor_rc, occ.mate_rc, occ.anchor_map_pos, occ.ins_chrm.c_str(),
                                     ins_pos, sample_id, sit->second.c_str());
                }
            }
            std::fclose(faf);
            if (!run_bwa_sam(bwa_bin, J.consensus, fa, {"-t",tstr}, J.out_disc_sam)) return false;
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

// ==================== full fold-in: collect during the detection scan pass ====================
//
// The detection scan (scan_tile) already streams every read once. collect_cns_read runs the cns
// collection gates on that same read and SPILLS the qualifying clip parts + disc anchors + low-mapq
// clip positions to text buffers (candidate-agnostic -- candidates are not known until after the
// scan). A later fan step (emit_cns_from_spill) loads the spill, does the ONE mate-gather pass, and
// fans per type. Net BAM passes for detection+cns: 2 (detection pass + mate-gather), down from 3.
//
// Spill line formats (tab-sep):
//   clip:   chrm  map_pos  L|R  is_rc  ref_start  ref_end  seq  qual
//   anchor: rname  s_mate_first  anchor_chrm  mate_chrm  anchor_map_pos  b_lclip  b_rclip  anchor_rc  mate_rc  ref_start  ref_end
//   lowq:   chrm  pos          (one per mapq<5 soft-clip near-position; fan counts within 35 bp)
inline void collect_cns_read(const bam1_t* r, sam_hdr_t* hdr,
                             std::string& clip_buf, std::string& anchor_buf, std::string& lowq_buf) {
    const uint16_t flag = r->core.flag;
    if (flag & (BAM_FDUP | BAM_FUNMAP)) return;
    const uint32_t* cig = bam_get_cigar(r);
    const int ncig = r->core.n_cigar;
    if (ncig < 1) return;
    const char* chrom = sam_hdr_tid2name(hdr, r->core.tid);
    const long rs = r->core.pos;
    char tmp[64];

    if (r->core.qual < kMinClipMapq) {                        // lowq background (mapq<5), incl supp/sec
        if (r->core.qual < kClipLowMapq) {
            long P = -1;
            if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) P = rs;
            else if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP)
                P = rs + (kReadLength - (long)bam_cigar_oplen(cig[ncig - 1]));
            if (P >= 0) { lowq_buf += chrom; int n = std::snprintf(tmp, sizeof tmp, "\t%ld\n", P); lowq_buf.append(tmp, n); }
        }
        return;
    }
    if (is_two_side_clipped(cig, ncig, kTwoSideClipMin)) return;

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
    std::string sseq, squal;
    long map_pos = rs;
    int b_left_clip = 0, b_right_clip = 0;

    if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) {
        if (is_supp_sec) return;
        const int clen = bam_cigar_oplen(cig[0]);
        sseq.clear();
        for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, i)];
        bool pass = true;
        if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
        else if (clen < kBwaRealignCut) pass = false;
        if (!pass) return;
        if (have_q && !is_qualified_clip(qual, clen)) return;
        squal.clear();
        for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[i] + 33) : 'I';
        int n = std::snprintf(tmp, sizeof tmp, "\t%ld\tL\t%d\t%ld\t%ld\t", rs, is_rc, rs, re);
        clip_buf += chrom; clip_buf.append(tmp, n); clip_buf += sseq; clip_buf += '\t'; clip_buf += squal; clip_buf += '\n';
        b_left_clip = 1;
    }
    if (bam_cigar_op(cig[ncig - 1]) == BAM_CSOFT_CLIP) {
        for (int i = 0; i < ncig - 1; ++i) {
            const int op = bam_cigar_op(cig[i]);
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP || op == BAM_CINS) continue;
            map_pos += bam_cigar_oplen(cig[i]);
        }
        if (is_supp_sec) return;
        const int clen = bam_cigar_oplen(cig[ncig - 1]);
        const int off = qlen - clen;
        sseq.clear();
        for (int i = 0; i < clen; ++i) sseq += seq_nt16_str[bam_seqi(seq, off + i)];
        bool pass = true;
        if (clen >= kMinPolyaClip && clen <= kBwaRealignCut) pass = contain_enough_A_T(sseq, clen - 1);
        if (clen < kBwaRealignCut) pass = false;
        if (!pass) return;
        if (have_q && !is_qualified_clip(qual + off, clen)) return;
        squal.clear();
        for (int i = 0; i < clen; ++i) squal += have_q ? (char)(qual[off + i] + 33) : 'I';
        int n = std::snprintf(tmp, sizeof tmp, "\t%ld\tR\t%d\t%ld\t%ld\t", map_pos, is_rc, rs, re);
        clip_buf += chrom; clip_buf.append(tmp, n); clip_buf += sseq; clip_buf += '\t'; clip_buf += squal; clip_buf += '\n';
        b_right_clip = 1;
    }

    if (mate_chrm == "*") return;
    if (is_decoy_contig(mate_chrm)) return;
    if (r->core.qual < kMinDiscMapq) return;
    bool disc = (std::string(chrom) != mate_chrm) || (llabs(mate_pos - map_pos) > kDiscThreshold);
    if (!disc) return;
    const int s_mate_first = b_first ? 0 : 1;
    // rname  s_mate_first  anchor_chrm  mate_chrm  anchor_map_pos  b_lclip b_rclip anchor_rc mate_rc ref_start ref_end
    anchor_buf += bam_get_qname(r);
    int n = std::snprintf(tmp, sizeof tmp, "\t%d\t", s_mate_first);
    anchor_buf.append(tmp, n); anchor_buf += chrom; anchor_buf += '\t'; anchor_buf += mate_chrm;
    n = std::snprintf(tmp, sizeof tmp, "\t%ld\t%d\t%d\t%d\t%d\t%ld\t%ld\n",
                      map_pos, b_left_clip, b_right_clip, is_rc, is_mate_rc, rs, re);
    anchor_buf.append(tmp, n);
}

// Load the three spill files back into the in-memory collection used by the per-type fan.
inline void load_cns_spill(const std::string& clip_file, const std::string& anchor_file,
                           const std::string& lowq_file, std::vector<HeldClip>& clips,
                           std::map<std::pair<std::string,int>, std::vector<HeldAnchor>>& anchors,
                           std::map<std::string, std::vector<long>>& lowq_pos) {
    if (FILE* f = std::fopen(clip_file.c_str(), "r")) {
        char* ln = nullptr; size_t cap = 0; ssize_t L;
        while ((L = getline(&ln, &cap, f)) != -1) {
            char* sv = nullptr;
            char* a[8]; int k = 0;
            for (char* t = strtok_r(ln, "\t\n", &sv); t && k < 8; t = strtok_r(nullptr, "\t\n", &sv)) a[k++] = t;
            if (k < 8) continue;
            clips.push_back({a[0], atol(a[1]), a[2][0], atoi(a[3]), atol(a[4]), atol(a[5]), a[6], a[7]});
        }
        free(ln); std::fclose(f);
    }
    if (FILE* f = std::fopen(anchor_file.c_str(), "r")) {
        char* ln = nullptr; size_t cap = 0; ssize_t L;
        while ((L = getline(&ln, &cap, f)) != -1) {
            char* sv = nullptr;
            char* a[11]; int k = 0;
            for (char* t = strtok_r(ln, "\t\n", &sv); t && k < 11; t = strtok_r(nullptr, "\t\n", &sv)) a[k++] = t;
            if (k < 11) continue;
            // rname s_mate_first anchor_chrm mate_chrm anchor_map_pos b_lclip b_rclip anchor_rc mate_rc ref_start ref_end
            anchors[{std::string(a[0]), atoi(a[1])}].push_back(
                {std::string(a[2]),   /* ins_chrm = anchor chrom */
                 atol(a[4]), atoi(a[5]), atoi(a[6]), atoi(a[7]), atoi(a[8]),
                 atol(a[9]), atol(a[10]), std::string(a[3]) /* mate_chrm */});
        }
        free(ln); std::fclose(f);
    }
    if (FILE* f = std::fopen(lowq_file.c_str(), "r")) {
        char* ln = nullptr; size_t cap = 0; ssize_t L;
        while ((L = getline(&ln, &cap, f)) != -1) {
            char* sv = nullptr;
            char* c = strtok_r(ln, "\t\n", &sv); char* p = strtok_r(nullptr, "\t\n", &sv);
            if (c && p) lowq_pos[c].push_back(atol(p));
        }
        free(ln); std::fclose(f);
    }
    for (auto& kv : lowq_pos) std::sort(kv.second.begin(), kv.second.end());
}

// Count spilled low-mapq clip positions within kNearbyLowQClip of candidate C (== emit_cns_multi's
// union-near count, computed at fan time from raw positions).
inline long lowq_near_count(const std::map<std::string, std::vector<long>>& lowq_pos,
                            const std::string& chrom, long C) {
    auto it = lowq_pos.find(chrom);
    if (it == lowq_pos.end()) return 0;
    const auto& v = it->second;
    auto lo = std::lower_bound(v.begin(), v.end(), C - kNearbyLowQClip + 1);  // P > C-dist
    auto hi = std::lower_bound(v.begin(), v.end(), C + kNearbyLowQClip);      // P < C+dist
    return (long)(hi - lo);
}

// Fold-in fan: given the spilled collection, gather mates (the ONE remaining BAM pass) and emit each
// type's temp SAMs. anchors[].ins_chrm is empty (loaded); we use the candidate list's chrm as the
// output ins_chrm and occ.mate_chrm as the gather target.
inline bool emit_cns_from_spill(const std::string& bam, const std::string& reference,
                                std::vector<CnsTypeJob>& jobs,
                                const std::string& clip_file, const std::string& anchor_file,
                                const std::string& lowq_file, long extnd, int sample_id,
                                const std::string& bwa_bin, int bwa_threads,
                                const std::string& region, const std::string& tmp_prefix) {
    std::vector<HeldClip> clips;
    std::map<std::pair<std::string,int>, std::vector<HeldAnchor>> anchors;
    std::map<std::string, std::vector<long>> lowq_pos;
    load_cns_spill(clip_file, anchor_file, lowq_file, clips, anchors, lowq_pos);

    // ---- the ONE remaining BAM pass: gather each mate's PRIMARY seq by (name,is_read1) ----
    std::map<std::pair<std::string,int>, std::string> mate_seq;
    {
        samFile* fp = sam_open(bam.c_str(), "r");
        if (!fp) return false;
        if (!reference.empty()) hts_set_fai_filename(fp, reference.c_str());
        sam_hdr_t* hdr = sam_hdr_read(fp);
        if (!hdr) { sam_close(fp); return false; }
        hts_idx_t* idx = sam_index_load(fp, bam.c_str());
        hts_itr_t* it = nullptr;
        if (!region.empty() && idx) it = sam_itr_querys(idx, hdr, region.c_str());
        else if (idx)               it = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
        bam1_t* b = bam_init1();
        std::string mseq;
        while (it && sam_itr_next(fp, it, b) >= 0) {
            const bam1_t* r = b;
            if (r->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
            const int is_r1 = (r->core.flag & BAM_FREAD1) ? 1 : 0;
            std::pair<std::string,int> key{std::string(bam_get_qname(r)), is_r1};
            auto ait = anchors.find(key);
            if (ait == anchors.end() || mate_seq.count(key)) continue;
            const char* rchrom = sam_hdr_tid2name(hdr, r->core.tid);
            bool ok = false;
            for (auto& occ : ait->second) if (occ.mate_chrm == rchrom) { ok = true; break; }
            if (!ok) continue;
            const uint8_t* s = bam_get_seq(r); const int Lq = r->core.l_qseq;
            mseq.clear();
            for (int i = 0; i < Lq; ++i) mseq += seq_nt16_str[bam_seqi(s, i)];
            mate_seq[key] = mseq;
        }
        bam_destroy1(b);
        if (it) hts_itr_destroy(it);
        if (idx) hts_idx_destroy(idx);
        sam_hdr_destroy(hdr); sam_close(fp);
    }

    const std::string tstr = std::to_string(bwa_threads);
    for (size_t ti = 0; ti < jobs.size(); ++ti) {
        CnsTypeJob& J = jobs[ti];
        const std::string pfx = tmp_prefix + "." + std::to_string(ti);
        if (!J.out_clip_sam.empty()) {
            const std::string fq_non = pfx + ".clip.non.fq", fq_pa = pfx + ".clip.pa.fq";
            FILE* f_non = std::fopen(fq_non.c_str(), "w");
            FILE* f_pa  = std::fopen(fq_pa.c_str(), "w");
            std::map<std::pair<std::string,long>, long> ncnt;
            std::vector<long> cs;
            for (const auto& c : clips) {
                cs.clear();
                J.cand.overlapping(c.chrm, c.ref_start, c.ref_end, extnd, cs);
                if (cs.empty()) continue;
                FILE* dst = ((int)c.seq.size() < kBwaRealignCut) ? f_pa : f_non;
                for (long ins_pos : cs) {
                    long nc = ncnt[{c.chrm, ins_pos}]++;
                    std::fprintf(dst, "@%s~%ld~%c~%d~%ld~%ld~%d\n%s\n+\n%s\n",
                                 c.chrm.c_str(), c.map_pos, c.side, c.is_rc, ins_pos, nc, sample_id,
                                 c.seq.c_str(), c.qual.c_str());
                }
            }
            std::fclose(f_non); std::fclose(f_pa);
            if (!realign_clip_cohorts(bwa_bin, J.consensus, fq_non, fq_pa, tstr, J.out_clip_sam, pfx + ".clip"))
                return false;
            ::unlink(fq_non.c_str()); ::unlink(fq_pa.c_str());
        }
        if (!J.out_disc_sam.empty()) {
            const std::string fa = pfx + ".disc.fa";
            FILE* faf = std::fopen(fa.c_str(), "w");
            std::vector<long> cs;
            for (const auto& kv : anchors) {
                auto sit = mate_seq.find(kv.first);
                if (sit == mate_seq.end()) continue;
                const int s_mate_first = kv.first.second;
                for (const auto& occ : kv.second) {
                    cs.clear();
                    J.cand.overlapping(occ.ins_chrm, occ.ref_start, occ.ref_end, extnd, cs);
                    for (long ins_pos : cs)
                        std::fprintf(faf, ">%s~%d~%d~%d~%d~%d~%ld~%s~%ld~%d\n%s\n",
                                     kv.first.first.c_str(), s_mate_first, occ.b_lclip, occ.b_rclip,
                                     occ.anchor_rc, occ.mate_rc, occ.anchor_map_pos, occ.ins_chrm.c_str(),
                                     ins_pos, sample_id, sit->second.c_str());
                }
            }
            std::fclose(faf);
            if (!run_bwa_sam(bwa_bin, J.consensus, fa, {"-t",tstr}, J.out_disc_sam)) return false;
            ::unlink(fa.c_str());
        }
        if (!J.out_lowq.empty()) {
            if (FILE* fl = std::fopen(J.out_lowq.c_str(), "w")) {
                for (const auto& kv : J.cand.by_chrom)
                    for (long pos : kv.second)
                        std::fprintf(fl, "%s %ld %ld\n", kv.first.c_str(), pos, lowq_near_count(lowq_pos, kv.first, pos));
                std::fclose(fl);
            }
        }
    }
    return true;
}

}  // namespace cns
}  // namespace xtea
