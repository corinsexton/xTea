#!/usr/bin/env python
"""Move-A prep: characterize the multi-type conflict loci on the WGS golden.
Replicates x_vcf_merger.merge_single_sample_ins: load L1, then SVA, then Alu;
conflict = incoming site within +-50bp of an already-placed site; resolution =
keep existing UNLESS incoming has two_side_tprt_both and existing does not.
Reports every locus where >=2 types call within 50bp = move A's divergence surface."""
import glob, collections

ROOT = "/home/cos689/scratch/xtea2_testdata/master_golden_wgs/NA12878"
WIN = 50
TPRT = "two_side_tprt_both"
FILES = [("L1", f"{ROOT}/L1/NA12878.wgs_LINE1.vcf"),
         ("SVA", f"{ROOT}/SVA/NA12878.wgs_SVA.vcf"),
         ("Alu", f"{ROOT}/Alu/NA12878.wgs_ALU.vcf")]

def load(path):
    out = []
    for ln in open(path):
        if ln.startswith("#") or len(ln.split("\t")) < 3:
            continue
        f = ln.split("\t")
        out.append((f[0], int(f[1]), TPRT in ln))
    return out

# per-type call counts
per_type = {t: load(p) for t, p in FILES}
for t, _ in FILES:
    n = len(per_type[t]); ntprt = sum(1 for _, _, b in per_type[t] if b)
    print(f"{t}: {n} calls ({ntprt} two_side_tprt_both)")

# Replicate merger: placed[chrm] = list of (pos, type, tprt); conflict within WIN.
placed = collections.defaultdict(list)
conflicts = []   # (chrm, existing(pos,type,tprt), incoming(pos,type,tprt), winner)
for t, _ in FILES:  # L1, then SVA, then Alu  (load order = priority)
    for chrm, pos, tprt in per_type[t]:
        hit = None
        for i, (p2, t2, tprt2) in enumerate(placed[chrm]):
            if abs(p2 - pos) <= WIN:
                hit = i; break
        if hit is not None:
            p2, t2, tprt2 = placed[chrm][hit]
            # keep existing unless incoming is tprt-both and existing is not
            keep_first = not (tprt and not tprt2)
            winner = t2 if keep_first else t
            conflicts.append((chrm, (p2, t2, tprt2), (pos, t, tprt), winner))
            if not keep_first:
                placed[chrm][hit] = (pos, t, tprt)  # incoming replaces
        else:
            placed[chrm].append((pos, t, tprt))

total_placed = sum(len(v) for v in placed.values())
total_calls = sum(len(per_type[t]) for t, _ in FILES)
print(f"\nper-type total = {total_calls}; merged (deduped) = {total_placed}; "
      f"conflicts collapsed = {total_calls - total_placed}")
print(f"conflict events = {len(conflicts)}")

# Break down conflicts by type pair + who wins + whether TPRT flipped the default priority
pair_ct = collections.Counter()
tprt_override = 0
for chrm, (p2, t2, b2), (pos, t, b), winner in conflicts:
    pair_ct[tuple(sorted((t2, t)))] += 1
    default_winner = t2  # load-order default = existing (L1>SVA>Alu)
    if winner != default_winner:
        tprt_override += 1
print("\nconflict type-pairs:", dict(pair_ct))
print(f"TPRT overrode load-order priority in {tprt_override}/{len(conflicts)} conflicts")

print("\n--- every conflict locus (existing | incoming -> winner) ---")
for chrm, (p2, t2, b2), (pos, t, b), winner in conflicts:
    print(f"  {chrm}:{p2}({t2}{'*' if b2 else ''}) vs {pos}({t}{'*' if b else ''}) -> {winner}"
          f"{'  [TPRT-override]' if winner != t2 else ''}")
