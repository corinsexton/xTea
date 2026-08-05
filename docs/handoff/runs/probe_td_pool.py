#!/usr/bin/env python3
"""The Stage 3 transduction gate's POOL: every (insertion locus, source) pair some caller asserts.

Offline, seconds, no BAM. Two things, and the second is the one that needs the arms to exist --

    # one arm: the shape table the gate spec calls sec 2.1
    python docs/handoff/runs/probe_td_pool.py --run-dir <out>/<sample> [--ledger master_na12878]

    # the pool: N arms merged at +-200 bp, one-to-one
    python docs/handoff/runs/probe_td_pool.py \
        --arm master=<master_out>/HG002 \
        --arm BASELINE-0=<base0>/BASELINE-0/HG002.wgs.vcf \
        --arm decoy=<gate_out>/HG002.wgs \
        --pool-tsv /path/to/td_pool.tsv

WHAT THIS IS FOR. `HANDOFF_2026-08-04_stage3_td_gate_spec.md` sec 2 defines the gate's pool as the
union of `(insertion locus X, source Y)` pairs from `TD_SRC` over three arms, stratified by sec 4's
two axes, and sec 7 step 2 lists "measure HG002's pool size" as unknown -- only NA12878's 115 was
ever counted. This produces both: the per-arm table, and the merged pool as a TSV the verifier
consumes. It decides nothing; a pool entry is a HYPOTHESIS, and the reads return the verdict.

IT DOES NOT PARSE `TD_SRC` ITSELF. Every classification is imported from
`benchmark/benchmark_truth_set.py` -- `parse_td_src`, `td_src_call_class`, `td_src_source_class`,
`match_one_to_one`, `DEFAULT_POS_WINDOW` -- for the same reason `probe_td_source_classes.py` imports
`transduction_evidence`: a probe with its own copy of a rule can agree with a scorer that disagrees
with it. The +-200 bp one-to-one merge below IS the harness's matcher, not a second implementation
of it, so a pool built here and a score computed there cannot disagree about which loci are the
same locus.

WHAT WRITING THIS MEASURED, and both corrections are in the spec's own sec 2.1 table:

  1. **`TD_SRC` has FIVE shapes, not four, and two of them used to parse to `None`** -- which is
     indistinguishable from `not_transduction`. `{source}~orphan` and
     `{source}~{cluster_l}-{cluster_r}` both did. That silently dropped 7 of master's 115 NA12878
     assertions, 3 of BASELINE-0's, and 5 in the `--td-decoy-structural` arm, and reached the
     scorer as `TP_POS_SOURCE_MISSED` for calls that named their source perfectly well. Fixed in
     `parse_td_src` by stripping the provenance tail; the one number it moves on BASELINE-0 is L1's
     `FP_TRANSDUCTION`, 0 -> 1.
  2. **"orphan" does not imply "polymorphic".** The spec's table read master's SVA orphans
     (`chr6:32546845-32546845~0~orphan`, start == end) as the shape of an orphan. Master's three L1
     orphans are reference copies spanning ~6 kb (`chr4:79966908-79972933~0~orphan`). Call class and
     source class are INDEPENDENT axes, which is why they are two columns here and not one.

     Consequence for the spec's L1 row: reference-copy sources are **36**, not 33 -- the three
     tailed orphans were invisible.

WHAT IT CANNOT SAY, which is sec 6.1 and belongs next to every number below: the pool is bounded by
what some caller found. It cannot say what all the arms missed, so it is not a denominator for
sensitivity and does not state how many transductions HG002 has.
"""
import argparse
import os
import sys
from collections import Counter, OrderedDict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "benchmark"))
import benchmark_truth_set as bts  # noqa: E402

FAMILIES = ("L1", "Alu", "SVA")
# The per-family VCF suffixes xTea writes, in this order, so a discovered file set is ordered.
ALT_SUFFIX = OrderedDict((("L1", "_LINE1.vcf"), ("Alu", "_ALU.vcf"), ("SVA", "_SVA.vcf")))

# EVERY source class has a reference window the verifier can align to, including polymorphic --
# which is not obvious and is worth stating, because the obvious reading is wrong. What a
# polymorphic source lacks in the reference is the ELEMENT BODY; the transduced sequence comes
# from DOWNSTREAM of that element's 3' end, which is ordinary genomic sequence and is present.
# `[ins_pos, ins_pos + flklen]` is exactly the window `gnrt_flank_regions_of_polymerphic_
# insertions` builds its flank from, which is why it writes start == end in the first place.
# So the source class does not decide verifiability; it decides where the window STARTS -- at the
# insertion point rather than after a body. Kept as a named tuple of all three so a future reader
# does not re-derive "polymorphic means unverifiable".
SRC_CLASSES_WITH_A_REFERENCE_WINDOW = (bts.SRC_CLASS_REFERENCE, bts.SRC_CLASS_UNKNOWN_START,
                                       bts.SRC_CLASS_POLYMORPHIC)

# The transduction step's own picked-source file, per family. Read ONLY for the reconciliation
# below -- it is a different population from the pool and the whole point is to say so.
PICKED_SOURCES = "candidate_disc_filtered_cns2.txt.all_non_sibling_td.txt.unique_trsdct_disc_only"


# ---------------------------------------------------------------------------------------
# the ledgers: what these two runs measured, so a re-run that disagrees says so out loud
# ---------------------------------------------------------------------------------------
# Keyed (family, call_class, source_class); source_class is None for a sibling, which has no
# interval. `calls` is the family's whole VCF. Measured 2026-08-04 from the runs named.
LEDGERS = {
    # master xTea's own NA12878 WGS run, job 48825977. The gate spec's sec 2.1 table, corrected
    # for the two tailed shapes: L1's reference copies are 36, not 33.
    "master_na12878": {
        "calls": {"L1": 203, "Alu": 1262, "SVA": 73},
        "strata": {
            ("L1", "regular", bts.SRC_CLASS_REFERENCE): 33,
            ("L1", "regular", bts.SRC_CLASS_UNKNOWN_START): 23,
            ("L1", "sibling", None): 27,
            ("L1", "orphan", bts.SRC_CLASS_REFERENCE): 3,
            ("SVA", "sibling", None): 25,
            ("SVA", "orphan", bts.SRC_CLASS_POLYMORPHIC): 4,
        },
    },
    # xTea2 BASELINE-0 on HG002, jobs 49441918-22, the frozen reference point.
    "baseline0_hg002": {
        "calls": {"L1": 188, "Alu": 1223, "SVA": 82},
        "strata": {
            ("L1", "regular", bts.SRC_CLASS_REFERENCE): 24,
            ("L1", "regular", bts.SRC_CLASS_UNKNOWN_START): 14,
            ("L1", "regular", bts.SRC_CLASS_POLYMORPHIC): 1,
            ("L1", "sibling", None): 26,
            ("SVA", "regular", bts.SRC_CLASS_POLYMORPHIC): 1,
            ("SVA", "sibling", None): 24,
        },
    },
}


# ---------------------------------------------------------------------------------------
# resolving an arm to its VCFs
# ---------------------------------------------------------------------------------------

def resolve_vcfs(path):
    """An arm -> the VCF path(s) to read, and the layout name, PRINTED rather than inferred.

    A `.vcf` is read as-is: `read_calls_vcf` takes the family from the ALT tag, so a merged VCF
    and a per-family one read identically. A directory is tried against the three layouts that
    exist on disk, in order, and the FIRST one that yields a file wins:

        <dir>/{L1,Alu,SVA}/*_{LINE1,ALU,SVA}.vcf     master's own pipeline
        <dir>/merged/*_{LINE1,ALU,SVA}.vcf           xTea2's run tree
        <dir>/*_{LINE1,ALU,SVA}.vcf                  a frozen BASELINE-0 directory

    First-match rather than union on purpose: an xTea2 tree contains BOTH `merged/` and a
    top-level merged VCF, and a union would count every call twice.
    """
    if path.endswith(".vcf"):
        return [path], "explicit vcf"
    if not os.path.isdir(path):
        raise SystemExit("not a .vcf and not a directory: {0}".format(path))
    layouts = [
        ("master per-type", lambda fam, suf: os.path.join(path, fam, "*" + suf)),
        ("xtea2 merged/", lambda fam, suf: os.path.join(path, "merged", "*" + suf)),
        ("flat", lambda fam, suf: os.path.join(path, "*" + suf)),
    ]
    import glob
    for name, mk in layouts:
        found = []
        for fam, suf in ALT_SUFFIX.items():
            found.extend(sorted(glob.glob(mk(fam, suf))))
        if found:
            return found, name
    raise SystemExit("no per-family VCFs under {0} in any known layout".format(path))


def read_arm(path):
    """(calls, vcf_paths, layout). Calls are deduplicated on (family, chrom, pos, TD_SRC),
    because the master layout's three VCFs are disjoint by family but a mis-resolved layout
    would not be -- and a pool built from a double-counted arm is silently wrong."""
    paths, layout = resolve_vcfs(path)
    seen = set()
    calls = []
    for p in paths:
        for c in bts.read_calls_vcf(p):
            key = (c.family, c.chrom, c.pos, c.td_src)
            if key in seen:
                continue
            seen.add(key)
            calls.append(c)
    return calls, paths, layout


# ---------------------------------------------------------------------------------------
# the per-arm table: sec 2.1
# ---------------------------------------------------------------------------------------

class Entry(object):
    """One pool candidate: a call that asserts a transduction, with both strata resolved.

    `chrom`/`pos` are the INSERTION locus and are named that way because `match_one_to_one`
    reads those two attributes -- this object stands in for a CallRecord in the merge.
    """
    __slots__ = ("arm", "family", "chrom", "pos", "call_class", "source_class",
                 "src_chrom", "src_start", "src_end", "raw")

    def __init__(self, arm, call):
        parsed = bts.parse_td_src(call.td_src)
        self.arm = arm
        self.family = call.family
        self.chrom = call.chrom
        self.pos = call.pos
        self.call_class = bts.td_src_call_class(call.td_src)
        self.source_class = bts.td_src_source_class(parsed)
        self.raw = call.td_src
        if parsed[0] == "sibling":
            # sec 2.1: use xTea's OWN convention for a sibling point, which is what
            # `x_orphan_transduction.py:1190-1191` does when it consumes the field --
            # [POS-100, POS+100]. Do not invent a window.
            self.src_chrom, self.src_start, self.src_end = parsed[1], parsed[2] - 100, parsed[2] + 100
        else:
            _, self.src_chrom, self.src_start, self.src_end = parsed

    @property
    def source_point(self):
        """The one comparable source coordinate, for cross-arm identity."""
        return (self.src_chrom, self.src_end if self.src_start < 0 else self.src_start)

    @property
    def stratum(self):
        return (self.family, self.call_class, self.source_class)


def entries_of(arm, calls):
    """Every call that asserts a transduction, as an Entry. `not_transduction` and the
    unparseable are dropped -- and the unparseable are COUNTED, because that count going
    above zero means a sixth TD_SRC shape has appeared and this probe is now blind to it."""
    out, n_unparsed, unparsed_examples = [], 0, []
    for c in calls:
        if not c.td_src or c.td_src == "not_transduction":
            continue
        if bts.parse_td_src(c.td_src) is None:
            n_unparsed += 1
            if len(unparsed_examples) < 3:
                unparsed_examples.append(c.td_src)
            continue
        out.append(Entry(arm, c))
    return out, n_unparsed, unparsed_examples


def report_arm(arm, calls, entries, n_unparsed, unparsed_examples, paths, layout):
    per_fam_calls = Counter(c.family for c in calls)
    strata = Counter(e.stratum for e in entries)

    print("--- arm {0!r} [{1}] ---".format(arm, layout))
    for p in paths:
        print("      {0}".format(p))
    print("      {0:<5} {1:>7} {2:>9}   {3}".format("fam", "calls", "assert", "by stratum"))
    for fam in FAMILIES:
        n_assert = sum(v for (f, _cc, _sc), v in strata.items() if f == fam)
        print("      {0:<5} {1:>7} {2:>9}".format(fam, per_fam_calls.get(fam, 0), n_assert))
        rows = sorted((k, v) for k, v in strata.items() if k[0] == fam)
        if not rows:
            if fam == "Alu":
                print("            none, BY CONSTRUCTION -- Alu has no transduction flank, so the"
                      " step never runs. Biology, not a gap.")
            continue
        for (_f, cc, sc), v in rows:
            note = ""
            if sc is None:
                note = "  (a sibling has no interval, so no source class -- not 'unknown')"
            elif sc == bts.SRC_CLASS_POLYMORPHIC:
                note = "  (the ELEMENT is absent from the reference; its 3' flank is not)"
            print("            {0:<8} {1:<15} {2:>5}{3}".format(cc, sc or "-", v, note))
    n_ref_window = sum(v for (_f, _cc, sc), v in strata.items()
                       if sc in SRC_CLASSES_WITH_A_REFERENCE_WINDOW)
    print("      TOTAL assert {0}   of which {1} name a source INTERVAL (a sibling names a point,"
          " so it gets xTea's own +-100)".format(len(entries), n_ref_window))
    if n_unparsed:
        print("      *** {0} TD_SRC value(s) DID NOT PARSE -- e.g. {1}".format(
            n_unparsed, ", ".join(repr(x) for x in unparsed_examples)))
        print("          A sixth TD_SRC shape exists and every count above is a floor. Fix")
        print("          `parse_td_src` before using this pool for anything.")
    return per_fam_calls, strata


def reconcile_with_picked_sources(path, entries):
    """The step's picked sources vs the pool, PER FAMILY -- two populations, not one number.

    THIS EXISTS BECAUSE ONE NUMBER WAS READ AT THE WRONG LEVEL. The gate spec's sec 2 puts the
    decoy-fixed arm in the pool on the strength of "241 reference-copy SVA sources", which is
    measured in `unique_trsdct_disc_only` -- the transduction step's own output, BEFORE the
    downstream filtering that produces the call set. The pool is defined on `TD_SRC` in the VCF,
    and there the same arm contributes **4** reference-copy SVA sources (1 regular + 3 orphan),
    against BASELINE-0's 0. Measured on job 49724056 vs BASELINE-0:

        family  picked (locus,source) pairs   non-sibling TD_SRC in the VCF
        L1                              159                            39     4.1x
        SVA, BASELINE-0                 105                             1   105.0x
        SVA, decoy-structural           346                             5    69.2x

    Neither number is wrong; they answer different questions. 241 sizes the STEP's recovered
    hypotheses, which is the right figure for Stage 3 threshold work because that is where the
    thresholds act. 4 sizes what the CALL SET asserts, which is the right figure for a precision
    gate -- sec 1 defines the question as "whether the transductions it asserts are real". So
    sec 5's GROWTH arm for SVA reference copies has **N = 4**, and sec 6.2's warning that SVA's
    regular class has N ~ 0 survives the fix rather than being retired by it.

    Printed for every arm that has the files, so this cannot be misread twice.
    """
    if not os.path.isdir(path):
        return
    rows = []
    for fam in ("L1", "SVA"):
        f = os.path.join(path, fam, PICKED_SOURCES)
        if not os.path.isfile(f):
            continue
        with open(f) as fin:
            n_picked = sum(1 for line in fin if line.strip())
        n_pool = sum(1 for e in entries
                     if e.family == fam and e.call_class != bts.CALL_CLASS_SIBLING)
        rows.append((fam, n_picked, n_pool))
    if not rows:
        return
    print("      --- reconciliation: the STEP's picked sources vs the POOL ---")
    print("      Different populations. The pool is `TD_SRC` in the call set; picked sources are")
    print("      the step's own hypotheses before downstream filtering. Do not quote one as the")
    print("      other -- the gate spec's '241 reference-copy SVA sources' is the picked level.")
    for fam, n_picked, n_pool in rows:
        ratio = ("{0:.1f}x".format(float(n_picked) / n_pool)) if n_pool else "all lost"
        print("      {0:<5} picked {1:>5}   non-sibling TD_SRC {2:>4}   attrition {3}".format(
            fam, n_picked, n_pool, ratio))


def check_ledger(name, per_fam_calls, strata):
    """PASS/FAIL against a measured run. Reports every disagreement, not the first."""
    led = LEDGERS[name]
    bad = []
    for fam, exp in sorted(led["calls"].items()):
        got = per_fam_calls.get(fam, 0)
        if got != exp:
            bad.append("calls {0}: expected {1}, got {2}".format(fam, exp, got))
    for k in sorted(set(led["strata"]) | set(strata), key=lambda t: tuple(str(x) for x in t)):
        exp, got = led["strata"].get(k, 0), strata.get(k, 0)
        if exp != got:
            bad.append("stratum {0}: expected {1}, got {2}".format(k, exp, got))
    print()
    if bad:
        print("  LEDGER {0}: FAIL ({1} disagreement(s))".format(name, len(bad)))
        for b in bad:
            print("      {0}".format(b))
        return False
    n = sum(led["strata"].values())
    print("  LEDGER {0}: PASS -- {1} calls in {2} strata reproduced exactly".format(
        name, n, len(led["strata"])))
    return True


# ---------------------------------------------------------------------------------------
# the pool: N arms merged
# ---------------------------------------------------------------------------------------

class PoolEntry(object):
    __slots__ = ("family", "chrom", "pos", "call_class", "source_class",
                 "src_chrom", "src_start", "src_end", "arms", "raws")

    def __init__(self, e):
        self.family, self.chrom, self.pos = e.family, e.chrom, e.pos
        self.call_class, self.source_class = e.call_class, e.source_class
        self.src_chrom, self.src_start, self.src_end = e.src_chrom, e.src_start, e.src_end
        self.arms = [e.arm]
        self.raws = {e.arm: e.raw}

    def absorb(self, e):
        self.arms.append(e.arm)
        self.raws[e.arm] = e.raw


def merge_pool(arm_entries, window, src_window):
    """The union of (locus, source) pairs over the arms, merged at +-`window` bp ONE-TO-ONE.

    TWO decisions, both of which change the count:

      * **locus identity is the harness's matcher**, `match_one_to_one` at
        `DEFAULT_POS_WINDOW`, so one call cannot satisfy two pooled entries (spec sec 2). Plan
        sec 7.3 records that greedy-vs-nearest and one-to-one-vs-not "moves sensitivity/FDR by
        percentage points"; the harness's choice is the one to keep, and importing it is how that
        stays true without anyone re-deciding it.
      * **matching is per (family, source), not per locus alone.** Two arms calling the same
        locus with DIFFERENT sources are asserting two different things, each independently
        verifiable, so they stay two pool entries. Folding them would silently pick a winner.
        Sources are the same when their comparable points are within `src_window`, the same
        window `benchmark_truth_set score` compares a called source to a derived one with.

    Arms are absorbed in the order given, so the first arm seeds the pool. That is a reporting
    convention only -- the union does not depend on it.
    """
    pool = []
    for arm, entries in arm_entries:
        # bucket the existing pool and this arm's entries the same way, so the matcher only
        # ever sees one family and one source group at a time
        def group(x):
            return (x.family, x.src_chrom)
        by_group = {}
        for p in pool:
            by_group.setdefault(group(p), []).append(p)
        for key, incoming in _grouped(entries, group):
            existing = by_group.get(key, [])
            if not existing:
                pool.extend(PoolEntry(e) for e in incoming)
                continue
            matched = bts.match_one_to_one(existing, incoming, window)
            taken = set()
            for ei, (ii, _d) in matched.items():
                # a locus match is only a POOL match if the source agrees too
                p, e = existing[ei], incoming[ii]
                if abs(_point(p) - _point(e)) <= src_window:
                    p.absorb(e)
                    taken.add(ii)
            for ii, e in enumerate(incoming):
                if ii not in taken:
                    pool.append(PoolEntry(e))
    return pool


def _point(x):
    return x.src_end if x.src_start < 0 else x.src_start


def _grouped(items, keyfn):
    out = {}
    for it in items:
        out.setdefault(keyfn(it), []).append(it)
    return sorted(out.items(), key=lambda kv: tuple(str(x) for x in kv[0]))


def report_pool(pool, arm_names, window, src_window):
    print("=" * 78)
    print("THE POOL -- {0} entries over {1} arm(s), merged at +-{2} bp one-to-one, sources "
          "matched at +-{3} bp".format(len(pool), len(arm_names), window, src_window))
    print("=" * 78)
    print("A pool entry is a HYPOTHESIS. Nothing here is verified, and the pool is bounded by")
    print("what these arms found -- it is not a sensitivity denominator (spec sec 6.1).")
    print()
    print("  --- by family x call class x source class ---")
    strata = Counter((p.family, p.call_class, p.source_class) for p in pool)
    for k in sorted(strata, key=lambda t: tuple(str(x) for x in t)):
        print("      {0:<5} {1:<8} {2:<15} {3:>5}".format(k[0], k[1], k[2] or "-", strata[k]))
    print()
    print("  --- arm agreement: how many arms assert each entry ---")
    for n in sorted(set(len(p.arms) for p in pool)):
        subset = [p for p in pool if len(p.arms) == n]
        print("      {0} arm(s): {1:>5} entr(ies)".format(n, len(subset)))
    print()
    print("  --- per arm: contributed / unique to it ---")
    for arm in arm_names:
        mine = [p for p in pool if arm in p.arms]
        only = [p for p in mine if len(p.arms) == 1]
        print("      {0:<14} {1:>5} contributed   {2:>5} unique".format(arm, len(mine), len(only)))
    print()
    print("      Read 'unique' as scope, NOT as error. An entry only one arm asserts is exactly")
    print("      what the gate exists to adjudicate from the reads -- and for the decoy arm it")
    print("      is the expected shape of the fix, which sec 5 scores as GROWTH in VERIFIED.")


def write_pool_tsv(pool, path, arm_names):
    with open(path, "w") as fout:
        fout.write("# the Stage 3 transduction gate's pool. One row = one (insertion locus,\n")
        fout.write("# source) HYPOTHESIS. `arms` is which callers assert it; that is provenance,\n")
        fout.write("# not evidence. src_start == src_end means a POLYMORPHIC source (absent from\n")
        fout.write("# the reference), src_start == -1 an unknown start whose END is the usable\n")
        fout.write("# coordinate, and a sibling's interval is xTea's own [POS-100, POS+100].\n")
        fout.write("# arms in this pool: {0}\n".format(",".join(arm_names)))
        fout.write("\t".join(("family", "ins_chrom", "ins_pos", "call_class", "source_class",
                             "src_chrom", "src_start", "src_end", "n_arms", "arms",
                             "td_src_raw")) + "\n")
        for p in sorted(pool, key=lambda x: (x.family, x.chrom, x.pos, x.src_chrom, x.src_start)):
            fout.write("\t".join(str(v) for v in (
                p.family, p.chrom, p.pos, p.call_class, p.source_class or "-",
                p.src_chrom, p.src_start, p.src_end, len(p.arms), ",".join(p.arms),
                p.raws[p.arms[0]])) + "\n")
    print()
    print("  wrote {0}".format(path))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", help="a single arm, named 'run-dir' in the report")
    p.add_argument("--arm", action="append", default=[], metavar="NAME=PATH",
                   help="an arm; repeatable. PATH is a .vcf or a run/sample directory")
    p.add_argument("--vcf-prefix", help="accepted and ignored -- layouts are discovered, not "
                                       "templated. Kept so callers can pass it harmlessly")
    p.add_argument("--pool-tsv", help="write the merged pool here")
    p.add_argument("--ledger", choices=sorted(LEDGERS),
                   help="check a single arm against a measured run; exits nonzero on FAIL")
    p.add_argument("--window", type=int, default=bts.DEFAULT_POS_WINDOW,
                   help="insertion-locus merge window (default %(default)s, the harness's)")
    p.add_argument("--src-window", type=int, default=1000,
                   help="source-agreement window (default %(default)s, `score`'s)")
    a = p.parse_args()

    arms = []
    if a.run_dir:
        arms.append(("run-dir", a.run_dir))
    for spec in a.arm:
        if "=" not in spec:
            raise SystemExit("--arm wants NAME=PATH, got {0!r}".format(spec))
        name, path = spec.split("=", 1)
        arms.append((name, path))
    if not arms:
        raise SystemExit("give at least one --run-dir or --arm NAME=PATH")

    print("=" * 78)
    print("TD_SRC -> THE STAGE 3 TRANSDUCTION GATE'S POOL")
    print("  strata: family x call class (regular/sibling/orphan) x source class "
          "(reference_copy/polymorphic/unknown_start)")
    print("=" * 78)
    print()

    arm_entries = []
    ok = True
    for name, path in arms:
        calls, paths, layout = read_arm(path)
        entries, n_unparsed, examples = entries_of(name, calls)
        per_fam, strata = report_arm(name, calls, entries, n_unparsed, examples, paths, layout)
        reconcile_with_picked_sources(path, entries)
        if n_unparsed:
            ok = False
        if a.ledger and len(arms) == 1:
            ok = check_ledger(a.ledger, per_fam, strata) and ok
        arm_entries.append((name, entries))
        print()

    if len(arm_entries) > 1 or a.pool_tsv:
        pool = merge_pool(arm_entries, a.window, a.src_window)
        report_pool(pool, [n for n, _ in arm_entries], a.window, a.src_window)
        if a.pool_tsv:
            write_pool_tsv(pool, a.pool_tsv, [n for n, _ in arm_entries])

    if a.ledger and not ok:
        return 1
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
