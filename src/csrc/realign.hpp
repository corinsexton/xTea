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
#include <array>
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

// Split a total CPU-core budget across `n` bwa processes that run CONCURRENTLY, so the
// sum of their `-t` never exceeds the allocation (no oversubscription). Job `i` (0-based)
// gets floor(total/n); the first (total % n) jobs get one extra core, so every core is
// used even when total isn't divisible by n (e.g. total=8,n=3 -> 3,3,2). Every job gets
// at least 1 thread, so total < n (e.g. a single core, total=1) degrades gracefully to
// 1 thread each rather than 0 -- minimally oversubscribed only when unavoidable.
inline int bwa_threads_for(int total, int n, int i) {
    if (n < 1) n = 1;
    if (total < 1) total = 1;
    if (n == 1) return total;
    int t = total / n + (i < (total % n) ? 1 : 0);
    return t < 1 ? 1 : t;
}

// Run one bwa mem over a FASTQ file; return qnames whose PRIMARY alignment mapped.
// If sam_out is non-empty, ALSO tee the complete bwa SAM (headers + every record,
// incl. supplementary/secondary) to that path -- the raw material for the cns
// temp_clip.sam / temp_disc.sam that move B.1 emits from the scan pass.
inline std::unordered_set<std::string>
realign_one_file(const std::string& bwa_bin, const std::string& consensus,
                 const std::string& fastq_path, int threads,
                 const std::string& sam_out = "") {
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
    FILE* samf = sam_out.empty() ? nullptr : std::fopen(sam_out.c_str(), "w");
    FILE* in = fdopen(from_bwa[0], "r");
    if (in) {
        char* line = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&line, &cap, in)) != -1) {
            if (samf) std::fwrite(line, 1, len, samf);   // tee full SAM before filtering
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
    if (samf) std::fclose(samf);
    int status = 0; waitpid(pid, &status, 0);
    return mapped;
}

// Run `bwa mem <args...> <ref> <fastq>` writing the FULL SAM (headers + every record)
// to out_sam. args are the flags between "mem" and the ref (e.g. {"-t","1","-T","9","-k","9"}).
// Returns true on a clean bwa exit. Used by the cns-emit path (move B.1) where downstream
// needs full alignment records (RNAME/POS/CIGAR/SEQ), not just the mapped-qname set.
inline bool run_bwa_sam(const std::string& bwa_bin, const std::string& ref,
                        const std::string& fastq, const std::vector<std::string>& args,
                        const std::string& out_sam) {
    pid_t pid = fork();
    if (pid < 0) return false;
    if (pid == 0) {
        int outfd = open(out_sam.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (outfd < 0) _exit(126);
        dup2(outfd, STDOUT_FILENO); close(outfd);
        int devnull = open("/dev/null", O_WRONLY);
        if (devnull >= 0) { dup2(devnull, STDERR_FILENO); close(devnull); }
        std::vector<const char*> argv;
        argv.push_back("bwa"); argv.push_back("mem");
        for (const auto& a : args) argv.push_back(a.c_str());
        argv.push_back(ref.c_str()); argv.push_back(fastq.c_str());
        argv.push_back(nullptr);
        execvp(bwa_bin.c_str(), const_cast<char* const*>(argv.data()));
        _exit(127);
    }
    int status = 0; waitpid(pid, &status, 0);
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
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
// matein (optional): (chrom,pos) -> union mate-in-repeat clip count (mate in L1|Alu|SVA),
// emitted as col5 -- the SAME value for every family at a site. Null => col5 is 0.
// sam_out_prefix (optional): if non-empty, each family's full bwa SAM is written to
// "<prefix><FAMILY>.sam" (move B.1: raw material for the scan-emitted cns temp SAMs).
inline void realign_sites_file(const std::string& fastq_path,
                               const std::vector<std::pair<std::string, std::string>>& families,
                               const std::string& bwa_bin, int threads, const std::string& out_path,
                               const std::map<std::pair<std::string, long>, int>* matein = nullptr,
                               const std::string& sam_out_prefix = "") {
    // These families realign CONCURRENTLY (one bwa each), so split the core budget
    // across them -- otherwise every family gets `threads` and we oversubscribe the
    // node by families.size()x. bwa here uses -K (kBwaBatch), so the alignment is
    // identical regardless of the per-family thread count.
    const int nf = (int)families.size();
    std::vector<std::future<std::unordered_set<std::string>>> futs;
    for (int fi = 0; fi < nf; ++fi)
        futs.push_back(std::async(std::launch::async, realign_one_file, bwa_bin, families[fi].second,
                                  fastq_path, bwa_threads_for(threads, nf, fi),
                                  sam_out_prefix.empty() ? std::string()
                                                         : sam_out_prefix + families[fi].first + ".sam"));
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
        for (const auto& fc : s.second) {
            const std::string& fam = families[fc.first].first;
            int col5 = 0;
            if (matein) {
                auto it = matein->find(s.first);
                if (it != matein->end()) col5 = it->second;
            }
            std::fprintf(out, "%s\t%ld\t%s\t%d\t%d\t%d\n", s.first.first.c_str(), s.first.second,
                         fam.c_str(), fc.second.nL, fc.second.nR, col5);
        }
    std::fclose(out);
}

}  // namespace xtea
