#!/usr/bin/env python3
"""Whole-genome three-way comparison: master, xTea2 with Python cns, xTea2 with C++ cns.

  master  job 48825977  master 46520ab, three types serially, -n 24. The reference
                        implementation AS PUBLISHED -- so it still carries the B10 disc
                        double-count, which is the point of keeping it as an arm.
  python  job 49065573  xTea2, Python cns collection, B10-fixed.
  cpp     job 49070030  xTea2, C++ cns collection (--cns-emit-scan).

python and cpp differ ONLY in which code collects the cns reads. master differs in
everything the fork has changed, so master-vs-either is a cumulative comparison, not an
attribution -- read the python-vs-cpp column when you want the port's own effect.

Three levels, because they answer different questions:

  A. CALLS         the deliverable. Per type and merged, with gained/lost positions.
  B. READ SETS     did the collectors gather the same reads? Primary records only, keyed on
                   what the collector decides, SEQ canonicalised for strand -- so bwa's
                   order-dependent placement does not masquerade as a collection difference.
  C. PLACEMENT     read set + RNAME/POS/CIGAR. Everything downstream actually parses. Where
                   B matches and C does not, the difference is bwa tie-breaking among
                   equal-scoring consensus copies, which follows from FASTQ order.

Usage: python compare_wgs_three_way.py [--quick]   (--quick skips B and C, which read ~5 GB)
"""
import collections
import os
import subprocess
import sys

SCRATCH = "/n/scratch/users/c/cos689/xtea2_testdata"

# (label, sample root, job id, note). master's tree uses NA12878/, xTea2's NA12878.wgs/.
ARMS = [
    ("master", f"{SCRATCH}/master_serial_out/NA12878",   "48825977", "46520ab, as published (pre-B10)"),
    ("python", f"{SCRATCH}/b10_out/NA12878.wgs",         "49065573", "xTea2, Python cns, B10-fixed"),
    ("cpp",    f"{SCRATCH}/port_wgs_out/NA12878.wgs",    "49070030", "xTea2, C++ cns collection"),
]
VCF_PAT = {"L1": "_LINE1.vcf", "Alu": "_ALU.vcf", "SVA": "_SVA.vcf"}
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def find_vcf(root, t):
    """master keeps per-type VCFs under <type>/; xTea2 under merged/."""
    for sub in ("merged", t):
        d = os.path.join(root, sub)
        if not os.path.isdir(d):
            continue
        hits = sorted(f for f in os.listdir(d) if f.endswith(VCF_PAT[t]))
        if hits:
            return os.path.join(d, hits[0])
    return None


def positions(vcf):
    if not vcf or not os.path.exists(vcf):
        return None
    out = set()
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 2)
            if len(f) >= 2:
                out.add((f[0], f[1]))
    return out


def canon(s):
    return min(s, s.translate(COMP)[::-1])


def cns_keys(path, kind):
    """-> (readset Counter, placement Counter). Primary records only."""
    rs, pl = collections.Counter(), collections.Counter()
    if not os.path.exists(path):
        return None, None
    with open(path) as fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or (int(f[1]) & 0x900):
                continue
            q = f[0].split("~")
            n = 5 if kind == "clip" else 9
            if len(q) <= n:
                continue
            base = tuple(q[:n])
            seq = canon(f[9]) if f[9] != "*" else "*"
            rs[base + (seq,)] += 1
            pl[base + (seq, f[2], f[3], f[5])] += 1
    return rs, pl


def sacct(job, field):
    try:
        out = subprocess.run(["sacct", "-j", job, "-n", "-P", "--format=" + field],
                             capture_output=True, text=True, timeout=30).stdout
        for line in out.splitlines():
            v = line.strip()
            if v and v != "":
                return v
    except Exception:
        pass
    return "?"


def phases(root):
    p = os.path.join(root, "step_timings.tsv")
    if not os.path.exists(p):
        return {}
    out = {}
    with open(p) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == "phase":
                out[f[1]] = float(f[2])
    return out


def main():
    quick = "--quick" in sys.argv
    present = [(lab, root, job, note) for lab, root, job, note in ARMS if os.path.isdir(root)]
    for lab, root, job, note in ARMS:
        if not os.path.isdir(root):
            print(f"!! {lab} missing: {root}", file=sys.stderr)
    if len(present) < 2:
        print("need at least two arms present", file=sys.stderr)
        return 1

    print("=" * 78)
    print("WGS three-way comparison")
    for lab, root, job, note in present:
        print(f"  {lab:<7} job {job}  {note}")
        print(f"          {root}")
    print("=" * 78)

    # ---------------- A. calls ----------------
    print("\n### A. CALLS ###\n")
    print("%-6s" % "type" + "".join("%10s" % l for l, _, _, _ in present))
    pos = {}
    for t in ("L1", "Alu", "SVA"):
        row = "%-6s" % t
        for lab, root, _, _ in present:
            p = positions(find_vcf(root, t))
            pos[(lab, t)] = p
            row += "%10s" % ("--" if p is None else len(p))
        print(row)
    # merged VCF (xTea2 only writes one; master has none)
    row = "%-6s" % "merged"
    for lab, root, _, _ in present:
        hits = [f for f in os.listdir(root) if f.endswith(".vcf")] if os.path.isdir(root) else []
        row += "%10s" % (len(positions(os.path.join(root, hits[0]))) if hits else "--")
    print(row)

    print("\n-- pairwise position differences (gained by the second / lost from the first) --")
    for i in range(len(present)):
        for j in range(i + 1, len(present)):
            a, b = present[i][0], present[j][0]
            print(f"\n  {a} -> {b}")
            for t in ("L1", "Alu", "SVA"):
                pa, pb = pos.get((a, t)), pos.get((b, t))
                if pa is None or pb is None:
                    print(f"    {t:<4} n/a")
                    continue
                print(f"    {t:<4} +{len(pb - pa):<5} -{len(pa - pb):<5} shared {len(pa & pb)}")

    # ---------------- timing / memory / disk ----------------
    print("\n### TIMING / MEMORY / DISK ###\n")
    print("%-7s %10s %10s %12s" % ("arm", "elapsed", "MaxRSS", "output"))
    for lab, root, job, _ in present:
        du = subprocess.run(["du", "-sh", os.path.dirname(root)], capture_output=True,
                            text=True).stdout.split("\t")[0] or "?"
        print("%-7s %10s %10s %12s" % (lab, sacct(job, "Elapsed"), sacct(job, "MaxRSS"), du))
    print("\n-- phases (s) --")
    allp = {lab: phases(root) for lab, root, _, _ in present}
    keys = []
    for d in allp.values():
        for k in d:
            if k not in keys:
                keys.append(k)
    print("  %-22s" % "phase" + "".join("%12s" % l for l in allp))
    for k in keys:
        print("  %-22s" % k + "".join("%12s" % (f"{allp[l].get(k, 0):.1f}" if k in allp[l] else "--")
                                      for l in allp))
    print("  %-22s" % "TOTAL" + "".join("%12.1f" % sum(allp[l].values()) for l in allp))

    if quick:
        print("\n(--quick: skipped read-set and placement comparison)")
        return 0

    # ---------------- B/C. read sets and placement ----------------
    print("\n### B. READ SETS  and  C. PLACEMENT (cns temp SAMs, primary records) ###\n")
    print("B identical + C differing  =>  same reads, different equal-scoring bwa placement.\n")
    print("%-6s %-5s %-16s %12s %12s" % ("type", "kind", "pair", "read set", "placement"))
    for t in ("L1", "Alu", "SVA"):
        for kind in ("clip", "disc"):
            loaded = {}
            for lab, root, _, _ in present:
                p = os.path.join(root, t, "tmp", "cns", f"temp_{kind}.sam")
                loaded[lab] = cns_keys(p, kind)
            for i in range(len(present)):
                for j in range(i + 1, len(present)):
                    a, b = present[i][0], present[j][0]
                    ra, pa = loaded[a]
                    rb, pb = loaded[b]
                    if ra is None or rb is None:
                        print("%-6s %-5s %-16s %12s %12s" % (t, kind, f"{a}/{b}", "n/a", "n/a"))
                        continue
                    d1 = (sum((ra - rb).values()), sum((rb - ra).values()))
                    d2 = (sum((pa - pb).values()), sum((pb - pa).values()))
                    f = lambda d: "IDENTICAL" if d == (0, 0) else "+%d/-%d" % (d[1], d[0])
                    print("%-6s %-5s %-16s %12s %12s" % (t, kind, f"{a}/{b}", f(d1), f(d2)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
