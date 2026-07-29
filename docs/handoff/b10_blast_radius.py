#!/usr/bin/env python3
"""How much did B10's double-count inflate per-site disc support, and at how many sites?

Re-baselining needs this before it starts. The fix is strictly more conservative, so goldens sit
above the fixed pipeline by construction -- but "how far above, at how many sites" decides whether
the losses are a handful of borderline calls or a broad shift.

Counts, per candidate site, the MAPPED disc records in the pre-fix run's temp_disc.sam, both raw
(what the filter actually saw) and deduplicated (what it would have seen). Only mapped records
count: an unmapped mate contributes no consensus position, so duplicating it changes nothing.

Run against the pre-fix indexed run (48842713).
"""
import collections
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    "/n/scratch/users/c/cos689/xtea2_testdata/fix_nofan_out/NA12878.wgs"

print("site-level disc support inflation from B10 (pre-fix run: %s)\n" % ROOT)
print("%-5s %9s %9s %9s %9s %9s" % ("type", "sites", "inflated", "pct", "median_x", "worst_x"))

for t in ("L1", "Alu", "SVA"):
    path = "%s/%s/tmp/cns/temp_disc.sam" % (ROOT, t)
    raw = collections.Counter()          # site -> mapped records seen by the filter
    seen = collections.defaultdict(set)  # site -> distinct record keys
    try:
        fh = open(path)
    except FileNotFoundError:
        print("%-5s MISSING" % t)
        continue
    with fh:
        for line in fh:
            if line[0] == "@":
                continue
            f = line.split("\t", 10)
            flag = int(f[1])
            if flag & 0x904:             # secondary/supplementary/unmapped
                continue
            q = f[0].split("~")
            if len(q) < 10:
                continue
            site = (q[7], q[8])
            raw[site] += 1
            seen[site].add("~".join(q[:9]))

    ratios = []
    inflated = 0
    for site, n in raw.items():
        d = len(seen[site])
        if n > d:
            inflated += 1
            ratios.append(n / d)
    ratios.sort()
    med = ratios[len(ratios) // 2] if ratios else 1.0
    worst = ratios[-1] if ratios else 1.0
    print("%-5s %9d %9d %8.1f%% %9.2f %9.2f"
          % (t, len(raw), inflated, 100.0 * inflated / max(1, len(raw)), med, worst))

    # The band that actually decides calls: sites whose support crossed a small integer boundary.
    if t == "Alu":
        drops = collections.Counter()
        for site, n in raw.items():
            d = len(seen[site])
            if n > d:
                drops[n - d] += 1
        print("      support lost per inflated site: %s"
              % ", ".join("%d reads: %d sites" % (k, v) for k, v in sorted(drops.items())[:8]))
        low = sum(1 for s, n in raw.items() if n > len(seen[s]) and len(seen[s]) <= 6)
        print("      of those, %d end up with <=6 deduplicated mapped disc reads -- the band where"
              % low)
        print("      a 1-2 read change can cross a cutoff, so these are the sites to inspect")
