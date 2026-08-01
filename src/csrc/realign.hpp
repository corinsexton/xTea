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
#include <mutex>
#include <chrono>

#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>

#include "htslib/sam.h"

namespace xtea {

// ---------------------------------------------------------------------------------------------
// Phase instrumentation. The realign phase was opaque: bwa's stderr went to /dev/null and there
// were no internal timers, so "20-22 min of a ~60 min run, of which only ~2 min is detection"
// could be measured from outside but nothing inside it could. These make it observable.
//
// The number that actually settles the standing hypothesis -- that the single-threaded PARENT
// parsing bwa's SAM through a pipe is the bottleneck, one reader per family against N bwa
// threads -- is the getline/body split:
//
//   * blocked in getline >> body  => the parent is waiting on bwa. bwa is the bottleneck; the
//     fix is more bwa throughput (more processes, or a faster consensus index), not a faster
//     reader.
//   * body >> blocked in getline  => the pipe is staying full because the parent cannot drain
//     it; bwa is blocking on write(). The parent IS the bottleneck, and the fix is bwa -> file
//     plus a parallel parse, or more bwa processes with fewer threads each.
//
// Both are measured with one clock read per line, which is cheap relative to the getline and
// the unordered_set insert it brackets, and is what makes the two mutually exclusive readings
// distinguishable rather than inferred.
// ---------------------------------------------------------------------------------------------
struct RealignStats {
    std::string family;
    double t_spawn = 0;        // fork+exec to the parent's first read attempt
    double t_first_record = 0; // spawn to the first non-header SAM line (bwa startup + index load)
    double t_total = 0;        // spawn to bwa reaped
    double t_getline = 0;      // wall spent blocked in getline (waiting on bwa)
    double t_body = 0;         // wall spent parsing/teeing/inserting (the parent's own work)
    long   lines = 0;
    long   records = 0;        // non-header SAM lines
    long   mapped = 0;         // distinct qnames with a mapped primary
    long   bytes = 0;
    int    threads = 0;
    int    exit_status = 0;
};

// One line per family to stderr, mutex-guarded so the concurrent per-family threads do not
// interleave. run_xtea captures the scan binary's stderr, so this lands in the run log next to
// the phase timings it explains.
inline void log_realign_stats(const RealignStats& s) {
    static std::mutex m;
    std::lock_guard<std::mutex> lk(m);
    const double pct_wait = s.t_total > 0 ? 100.0 * s.t_getline / s.t_total : 0.0;
    std::fprintf(stderr,
        "[realign] %-6s t=%d total=%.1fs spawn=%.3fs first_rec=%.1fs "
        "getline=%.1fs(%.0f%%) body=%.1fs lines=%ld recs=%ld mapped=%ld MB=%.1f exit=%d %s\n",
        s.family.empty() ? "?" : s.family.c_str(), s.threads, s.t_total, s.t_spawn,
        s.t_first_record, s.t_getline, pct_wait, s.t_body, s.lines, s.records, s.mapped,
        s.bytes / 1048576.0, s.exit_status,
        // The verdict, stated so a reader does not have to recompute it. 60% is not a magic
        // number: it just has to be far enough from 50 that the two readings are not a
        // coin flip, and the WGS gap is expected to be lopsided either way.
        pct_wait >= 60.0 ? "=> BWA-BOUND (parent waits on bwa)"
                         : (pct_wait <= 40.0 ? "=> PARENT-BOUND (parse cannot drain the pipe)"
                                             : "=> BALANCED (neither side clearly dominates)"));
}

inline const char* kBwaMinScore = "7";        // -T
inline const char* kBwaSeedLen  = "7";        // -k
inline const char* kBwaClipPen  = "70";       // -c
// -K fixes bwa's batch size in bases so the output does not depend on the thread count. ANY
// fixed value gives that, and a smaller one holds far less resident: measured on a real clip
// FASTQ (293,758 reads, 24 cores) 10000000 -> 2.59 GB, 1000000 -> 0.33 GB, 200000 -> 0.07 GB,
// with no time cost (10.44 s / 9.87 s / 9.80 s). Three of these run concurrently, so the old
// value was holding ~7.8 GB for nothing on a run whose peak forced --mem=64G.
//
// The alignment records are byte-identical across all three, verified for EVERY param set bwa
// is invoked with here -- detection realign (-T7 -k7 -c70), cns clip non-polyA (-T9 -k9), cns
// clip polyA (-T7 -k7 -c70) and cns disc (no extra flags). The only differing line is bwa's own
// @PG CL: header, which records the -K value; the cns parsers read the SAM through pysam, so
// they never see it. That check mattered here rather than being pro forma: the cns temp SAMs are
// sensitive to alignment PLACEMENT (see B10's note on bwa's read-index-dependent tie-break), so
// "same records" is the property that had to hold, not "same counts".
//
// 1000000 rather than 200000: it captures essentially all of the memory win (2.59 -> 0.33 GB)
// while leaving chunk headroom for longer reads than this workload's 13-70 bp clip fragments.
inline const char* kBwaBatch    = "1000000";  // -K (deterministic across threads)

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

// Serializes the "create fds + mark close-on-exec + fork" sequence across the concurrent
// per-family bwa spawns (realign_sites_file launches one realign_one_file per family via
// std::async). WITHOUT this, family B's fork can snapshot family A's freshly-pipe()'d fds
// during the window before they're FD_CLOEXEC-tagged, carrying A's stdout WRITE end into B's
// bwa process across exec. A's reader (this parent) then never sees EOF on that pipe -> it
// blocks reading -> never drains B's stdout -> B blocks writing -> permanent deadlock.
// (Observed: WGS run 48091342 wedged here with all bwa+scan threads in futex_do_wait.)
// The mutex closes the race; FD_CLOEXEC ensures any inherited pipe fd is dropped at exec.
// (pipe2(O_CLOEXEC) would be atomic, but it's not declared under strict -std=c++17.)
inline std::mutex& bwa_spawn_mutex() { static std::mutex m; return m; }

// Run one bwa mem over a FASTQ file; return qnames whose PRIMARY alignment mapped.
// If sam_out is non-empty, ALSO tee the complete bwa SAM (headers + every record,
// incl. supplementary/secondary) to that path -- the raw material for the cns
// temp_clip.sam / temp_disc.sam that move B.1 emits from the scan pass.
//
// bwa_log: where bwa's own stderr goes. bwa mem reports its per-chunk progress there
// ("[M::mem_process_seqs] Processed N reads in X CPU sec, Y real sec"), which is the only
// view of what bwa itself is doing and is worth keeping; it used to be sent to /dev/null.
// Empty => /dev/null, preserving the old behaviour for callers that do not want the file.
// stats: filled in if non-null (see RealignStats).
inline std::unordered_set<std::string>
realign_one_file(const std::string& bwa_bin, const std::string& consensus,
                 const std::string& fastq_path, int threads,
                 const std::string& sam_out = "", const std::string& bwa_log = "",
                 RealignStats* stats = nullptr) {
    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    auto secs = [](clk::duration d) {
        return std::chrono::duration_cast<std::chrono::duration<double>>(d).count();
    };
    std::unordered_set<std::string> mapped;
    int from_bwa[2];
    pid_t pid;
    {
        // pipe + CLOEXEC + fork as one critical section (see bwa_spawn_mutex above).
        std::lock_guard<std::mutex> lk(bwa_spawn_mutex());
        if (pipe(from_bwa) != 0) return mapped;
        fcntl(from_bwa[0], F_SETFD, FD_CLOEXEC);   // read end: never inherited past exec
        fcntl(from_bwa[1], F_SETFD, FD_CLOEXEC);   // write end: closed in every child at exec
        pid = fork();
        if (pid == 0) {
            // dup2 target fds are NOT close-on-exec, so this child's own stdout->pipe and
            // stderr->/dev/null survive exec; all other (CLOEXEC) fds are dropped.
            dup2(from_bwa[1], STDOUT_FILENO);
            // bwa's own progress goes to bwa_log when asked for, /dev/null otherwise. Opened
            // WITHOUT O_CLOEXEC on purpose: dup2 is what carries it past exec.
            int errfd = bwa_log.empty()
                          ? open("/dev/null", O_WRONLY)
                          : open(bwa_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (errfd < 0) errfd = open("/dev/null", O_WRONLY);
            if (errfd >= 0) dup2(errfd, STDERR_FILENO);
            close(from_bwa[0]); close(from_bwa[1]);
            std::string tstr = std::to_string(threads);
            execlp(bwa_bin.c_str(), "bwa", "mem", "-t", tstr.c_str(),
                   "-K", kBwaBatch, "-T", kBwaMinScore, "-k", kBwaSeedLen, "-c", kBwaClipPen,
                   consensus.c_str(), fastq_path.c_str(), (char*)nullptr);
            _exit(127);
        }
    }
    if (pid < 0) { close(from_bwa[0]); close(from_bwa[1]); return mapped; }
    close(from_bwa[1]);
    const auto t_spawned = clk::now();
    FILE* samf = sam_out.empty() ? nullptr : std::fopen(sam_out.c_str(), "w");
    FILE* in = fdopen(from_bwa[0], "r");
    clk::duration d_getline{0}, d_body{0};
    clk::time_point t_first{};
    long n_lines = 0, n_records = 0, n_bytes = 0;
    if (in) {
        char* line = nullptr; size_t cap = 0; ssize_t len;
        for (;;) {
            // Bracket ONLY the getline: whatever it blocks on is time the parent spent waiting
            // for bwa. Everything after it, up to the next getline, is the parent's own work.
            const auto a = clk::now();
            len = getline(&line, &cap, in);
            const auto b = clk::now();
            d_getline += b - a;
            if (len == -1) break;
            ++n_lines; n_bytes += len;
            if (samf) std::fwrite(line, 1, len, samf);   // tee full SAM before filtering
            if (!(len == 0 || line[0] == '@')) {
                if (n_records == 0) t_first = b;
                ++n_records;
                char* tab1 = (char*)memchr(line, '\t', len);
                if (tab1) {
                    int flag = atoi(tab1 + 1);
                    if (!(flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP)))
                        mapped.emplace(line, tab1 - line);
                }
            }
            d_body += clk::now() - b;
        }
        free(line);
        fclose(in);
    } else {
        close(from_bwa[0]);
    }
    if (samf) std::fclose(samf);
    int status = 0; waitpid(pid, &status, 0);
    if (stats) {
        stats->t_spawn = secs(t_spawned - t0);
        stats->t_first_record = n_records ? secs(t_first - t0) : 0.0;
        stats->t_total = secs(clk::now() - t0);
        stats->t_getline = secs(d_getline);
        stats->t_body = secs(d_body);
        stats->lines = n_lines;
        stats->records = n_records;
        stats->mapped = (long)mapped.size();
        stats->bytes = n_bytes;
        stats->threads = threads;
        stats->exit_status = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    }
    return mapped;
}

// Run `bwa mem <args...> <ref> <fastq>` writing the FULL SAM (headers + every record)
// to out_sam. args are the flags between "mem" and the ref (e.g. {"-t","1","-T","9","-k","9"}).
// Returns true on a clean bwa exit. Used by the cns-emit path (move B.1) where downstream
// needs full alignment records (RNAME/POS/CIGAR/SEQ), not just the mapped-qname set.
//
// This path has no pipe to read, so there is no getline/body split to measure -- bwa writes
// the SAM itself and the only cost is bwa. What was missing is visibility: it discarded bwa's
// stderr too, so a cns realign that took 20 minutes reported nothing about why. bwa's log now
// goes to "<out_sam>.bwa.log" and the wall time is reported per call.
inline bool run_bwa_sam(const std::string& bwa_bin, const std::string& ref,
                        const std::string& fastq, const std::vector<std::string>& args,
                        const std::string& out_sam) {
    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    const std::string bwa_log = out_sam.empty() ? std::string() : out_sam + ".bwa.log";
    pid_t pid;
    {
        // Serialize the fork with all other bwa spawns (the cns-emit fan runs several of these
        // concurrently via std::async) so no child inherits another thread's in-flight fds; the
        // child's own outfd/devnull are O_CLOEXEC (dropped at exec) and survive only via dup2.
        std::lock_guard<std::mutex> lk(bwa_spawn_mutex());
        pid = fork();
        if (pid == 0) {
            int outfd = open(out_sam.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
            if (outfd < 0) _exit(126);
            dup2(outfd, STDOUT_FILENO); close(outfd);
            int errfd = bwa_log.empty()
                          ? open("/dev/null", O_WRONLY | O_CLOEXEC)
                          : open(bwa_log.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
            if (errfd < 0) errfd = open("/dev/null", O_WRONLY | O_CLOEXEC);
            if (errfd >= 0) { dup2(errfd, STDERR_FILENO); close(errfd); }
            std::vector<const char*> argv;
            argv.push_back("bwa"); argv.push_back("mem");
            for (const auto& a : args) argv.push_back(a.c_str());
            argv.push_back(ref.c_str()); argv.push_back(fastq.c_str());
            argv.push_back(nullptr);
            execvp(bwa_bin.c_str(), const_cast<char* const*>(argv.data()));
            _exit(127);
        }
    }
    if (pid < 0) return false;
    int status = 0; waitpid(pid, &status, 0);
    const int rc = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    {
        static std::mutex m;   // several of these run concurrently (the per-type cns fan)
        std::lock_guard<std::mutex> lk(m);
        const double t = std::chrono::duration_cast<std::chrono::duration<double>>(
                             clk::now() - t0).count();
        // Basename only: the full temp paths are long and the family/round is the informative part.
        const size_t slash = out_sam.rfind('/');
        std::fprintf(stderr, "[bwa_sam] %-28s %.1fs exit=%d\n",
                     (slash == std::string::npos ? out_sam : out_sam.substr(slash + 1)).c_str(),
                     t, rc);
    }
    return rc == 0;
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
// matein (optional): (chrom,pos) -> mate_in_rep_clip_cnt, the union mate-in-repeat clip count
// (mate in L1|Alu|SVA), emitted as the 6th column -- the SAME value for every family at a
// site. Null => the count is 0. ("col5" was its old name: a column position, not a concept.)
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
    // Per-family stats live here (not in the lambda) so they outlive the async call and can be
    // logged after the join; bwa's own stderr goes next to out_path, one log per family.
    std::vector<RealignStats> stats(nf);
    std::vector<std::future<std::unordered_set<std::string>>> futs;
    for (int fi = 0; fi < nf; ++fi) {
        stats[fi].family = families[fi].first;
        futs.push_back(std::async(std::launch::async, realign_one_file, bwa_bin, families[fi].second,
                                  fastq_path, bwa_threads_for(threads, nf, fi),
                                  sam_out_prefix.empty() ? std::string()
                                                         : sam_out_prefix + families[fi].first + ".sam",
                                  out_path + "." + families[fi].first + ".bwa.log",
                                  &stats[fi]));
    }
    // family -> site -> counts
    std::vector<std::map<std::pair<std::string, long>, FamilyCounts>> per_family(families.size());
    for (size_t fi = 0; fi < families.size(); ++fi) {
        auto mapped = futs[fi].get();
        log_realign_stats(stats[fi]);
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
            int mate_in_rep_clip_cnt = 0;
            if (matein) {
                auto it = matein->find(s.first);
                if (it != matein->end()) mate_in_rep_clip_cnt = it->second;
            }
            std::fprintf(out, "%s\t%ld\t%s\t%d\t%d\t%d\n", s.first.first.c_str(), s.first.second,
                         fam.c_str(), fc.second.nL, fc.second.nR, mate_in_rep_clip_cnt);
        }
    std::fclose(out);
}

}  // namespace xtea
