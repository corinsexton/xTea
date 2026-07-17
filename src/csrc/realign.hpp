// realign.hpp — Stage 2 of xtea-scan: realign clipped/mate parts to the TE consensi with bwa.
//
// Memory-lean design: the scan streams candidate sequences to a consolidated FASTQ *file* (one file,
// not thousands; on disk, not RAM). Here we run one `bwa mem` per MEI family (L1/Alu/SVA)
// concurrently against that file, reusing pre-built consensus indexes (never rebuilt). Each candidate
// read's FASTQ header encodes its site (name~chrom~pos~side), so we attribute consensus hits back to
// sites straight from the qname — no in-memory read list. Per site we count how many left/right reads
// aligned to each family's consensus ("aligned to copies" evidence).
//
// bwa params mirror master's consensus-realign step (-T 7 -k 7 -c 70); -K fixes the batch size so
// output is independent of bwa's thread count.
//
// Scope: validates mechanism (concurrent bwa, cached index, per-site/family aggregation) + biological
// sanity. Exact master LCLIP/RCLIP/LDISC/RDISC integer counts are pinned later by VCF concordance.

#pragma once
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_set>
#include <map>
#include <future>

#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>

#include "htslib/sam.h"

namespace xtea {

inline const char* kBwaMinScore = "7";        // -T
inline const char* kBwaSeedLen  = "7";        // -k
inline const char* kBwaClipPen  = "70";       // -c
inline const char* kBwaBatch    = "10000000"; // -K (deterministic across threads)

// Run one bwa mem over a FASTQ file; return qnames whose PRIMARY alignment mapped.
inline std::unordered_set<std::string>
realign_one_file(const std::string& bwa_bin, const std::string& consensus,
                 const std::string& fastq_path, int threads) {
    std::unordered_set<std::string> mapped;
    int from_bwa[2];
    if (pipe(from_bwa) != 0) return mapped;
    pid_t pid = fork();
    if (pid == 0) {
        dup2(from_bwa[1], STDOUT_FILENO);
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) dup2(devnull, STDERR_FILENO);
        close(from_bwa[0]); close(from_bwa[1]);
        std::string tstr = std::to_string(threads);
        execlp(bwa_bin.c_str(), "bwa", "mem", "-t", tstr.c_str(),
               "-K", kBwaBatch, "-T", kBwaMinScore, "-k", kBwaSeedLen, "-c", kBwaClipPen,
               consensus.c_str(), fastq_path.c_str(), (char*)nullptr);
        _exit(127);
    }
    close(from_bwa[1]);
    FILE* in = fdopen(from_bwa[0], "r");
    if (in) {
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, in)) != -1) {
            if (len == 0 || line[0] == '@') continue;
            char* tab1 = (char*)memchr(line, '\t', len);
            if (!tab1) continue;
            int flag = atoi(tab1 + 1);
            if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
            mapped.emplace(line, tab1 - line);
        }
        free(line);
        fclose(in);
    } else {
        close(from_bwa[0]);
    }
    int status = 0; waitpid(pid, &status, 0);
    return mapped;
}

struct FamilyCounts { int nL = 0; int nR = 0; };

// qname format: "<readname>~<chrom>~<pos>~<side>"  (readname uses ':' not '~', so ~ is a safe delim)
inline bool parse_site_qname(const std::string& qn, std::string& chrom, long& pos, char& side) {
    auto p3 = qn.rfind('~'); if (p3 == std::string::npos) return false;
    auto p2 = qn.rfind('~', p3 - 1); if (p2 == std::string::npos) return false;
    auto p1 = qn.rfind('~', p2 - 1); if (p1 == std::string::npos) return false;
    side  = qn[p3 + 1];
    pos   = std::strtol(qn.substr(p2 + 1, p3 - p2 - 1).c_str(), nullptr, 10);
    chrom = qn.substr(p1 + 1, p2 - p1 - 1);
    return true;
}

// Realign a candidate FASTQ against all families concurrently; write per-site/family L/R counts.
// Output (sorted by chrom,pos then family): chrom \t pos \t family \t nL \t nR
inline void realign_sites_file(const std::string& fastq_path,
                               const std::vector<std::pair<std::string, std::string>>& families,
                               const std::string& bwa_bin, int threads, const std::string& out_path) {
    std::vector<std::future<std::unordered_set<std::string>>> futs;
    for (const auto& fam : families)
        futs.push_back(std::async(std::launch::async, realign_one_file, bwa_bin, fam.second,
                                  fastq_path, threads));
    // family -> site -> counts
    std::vector<std::map<std::pair<std::string, long>, FamilyCounts>> per_family(families.size());
    for (size_t fi = 0; fi < families.size(); ++fi) {
        auto mapped = futs[fi].get();
        for (const auto& qn : mapped) {
            std::string chrom; long pos; char side;
            if (!parse_site_qname(qn, chrom, pos, side)) continue;
            auto& fc = per_family[fi][{chrom, pos}];
            if (side == 'L') ++fc.nL; else ++fc.nR;
        }
    }
    FILE* out = std::fopen(out_path.c_str(), "w");
    if (!out) return;
    // Merge to a site-major ordering for stable output.
    std::map<std::pair<std::string, long>, std::vector<std::pair<size_t, FamilyCounts>>> bysite;
    for (size_t fi = 0; fi < families.size(); ++fi)
        for (const auto& kv : per_family[fi]) bysite[kv.first].push_back({fi, kv.second});
    for (const auto& s : bysite)
        for (const auto& fc : s.second)
            std::fprintf(out, "%s\t%ld\t%s\t%d\t%d\n", s.first.first.c_str(), s.first.second,
                         families[fc.first].first.c_str(), fc.second.nL, fc.second.nR);
    std::fclose(out);
}

}  // namespace xtea
