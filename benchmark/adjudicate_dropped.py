#!/usr/bin/env python3
"""Stage 3's check 2b: did the ON arm drop a hypothesis the OFF arm's own reads VERIFIED?

    python benchmark/adjudicate_dropped.py <lost_keys.txt> <off_arm_verdicts.tsv>

    exit 0  no dropped hypothesis was VERIFIED  (the count is still printed)
    exit 1  at least one was, or at least one could not be adjudicated at all
    exit 2  the inputs are not what this program can read

THE RULE, decided by Corinne 2026-08-07 after the gate FAILED on a truth-set false positive
(job 49970296): **a dropped hypothesis fails the gate only if the OFF arm's own reads VERIFIED
it.** Dropping one the reads called UNDECIDABLE, REFUTED or SUB_THRESHOLD is REPORTED and does
not fail.

WHY THE OLD RULE WAS THE WRONG SHAPE. Check 2a counts hypotheses, and a change whose whole
mechanism is REMOVING spurious sources reduces that count by construction -- so the test could
never be satisfied by the class of change it was being applied to, however good the change was.
`--td-decoy-concat` dropped exactly one entry, `chr17:45585847`, which the truth set counts a
FALSE POSITIVE and which the verifier (after the self-window guard, 2026-08-06) calls
UNDECIDABLE. Failing a gate for deleting a false positive inverts what the gate is for.

WHAT THE OLD RULE WAS ACTUALLY PROTECTING, and this keeps it exactly: a per-stratum VERIFIED
count cannot see a dropped entry, because the entry is absent from BOTH sides of it. So evidence
can leave silently. That hazard is real and is entirely about VERIFIED entries -- which is what
this tests, and nothing more.

**THE COUNT IS PRINTED ON PASS AS WELL AS FAIL.** The residual risk in the new rule is a change
that drops many un-verified but nonetheless real hypotheses, and no automatic rule catches that.
A number on the page is what makes it noticeable.

WHY THIS IS A MODULE AND NOT THE HEREDOC IT USED TO BE. It lived inline in
`docs/handoff/runs/td_gate_no_regression.sbatch` as a `python - <<EOF` block -- load-bearing
pass/fail logic in the one surface of this project with no tests and no mutation coverage. Both
harness defects of 2026-08-06 (a byte-`cmp` against a file whose column count had changed, and a
reader taking columns by INDEX) lived in that same surface while the unit suite stayed green, and
handoff sec 5.2 asks its reader to *check* that this check would still fail a drop of a VERIFIED
entry rather than take the claim. That check is only possible once the logic can be imported.
`tests/test_adjudicate_dropped.py` is that check.

COLUMNS ARE READ BY NAME, NEVER BY INDEX. `verdict` is field 13 in the frozen `verdicts.tsv` and
field 15 in the current writer's, because the restratification added `shape` and `sidedness`
after the file was frozen. An index-based reader silently compares `n_clip_support` against
`verdict` and can report any answer at all. A missing column is a loud exit 2 rather than a
KeyError traceback, because "the harness is reading the wrong file" and "the gate failed" must
not look alike.
"""
import sys

# The verdict string that, on a DROPPED entry, fails the gate. Everything else -- UNDECIDABLE,
# REFUTED, SUB_THRESHOLD -- is reported and waived. One name, because a second copy of this
# string somewhere else is the private-copy-of-a-shipping-rule mistake this project has now paid
# for four times (plan sec 10).
FATAL_VERDICT = "VERIFIED"

# The columns this program needs from verdicts.tsv, by name.
REQUIRED = ("family", "ins_chrom", "ins_pos", "td_src_raw", "verdict", "why")

_UNADJUDICATED = "<not in the OFF arm's verdicts>"

# Exit 2, distinct from the gate's own 1. `SystemExit("...")` would print the message and exit
# **1**, which is precisely the confusion this code is trying to prevent: the caller could not
# tell "the harness is reading the wrong file" from "the flip lost verified evidence".
EXIT_BAD_INPUT = 2


def _bad_input(msg):
    sys.stderr.write(msg + "\n")
    raise SystemExit(EXIT_BAD_INPUT)


def read_lost_keys(path):
    """The dropped keys, as written by check 2a's `cut -f1-8,11` of the pool.

    Fields are family, ins_chrom, ins_pos, ..., td_src_raw -- so the key is (0, 1, 2, 8), the
    same 4-tuple the verdicts file is unique on. Short or blank lines are skipped: `comm` writes
    a trailing newline and the caller may `printf` an empty result into the file, so "no drops"
    arrives as a file containing one empty line rather than as an empty file.
    """
    lost = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[0]:
            lost.append((f[0], f[1], f[2], f[8]))
    return lost


def read_verdicts(path):
    """`(family, ins_chrom, ins_pos, td_src_raw) -> (verdict, why)`, parsed BY HEADER NAME."""
    hdr, verd = None, {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            missing = [c for c in REQUIRED if c not in hdr]
            if missing:
                _bad_input(
                    "FATAL: {0} has no column(s) {1}.\n"
                    "  Its header is: {2}\n"
                    "  This is a wrong-file or wrong-writer-version error, NOT a gate failure.".
                    format(path, ", ".join(missing), "\t".join(hdr)))
            continue
        d = dict(zip(hdr, f))
        verd[(d["family"], d["ins_chrom"], d["ins_pos"], d["td_src_raw"])] = (d["verdict"],
                                                                             d["why"])
    if hdr is None:
        _bad_input("FATAL: {0} is empty -- no header row. The OFF arm did not "
                   "adjudicate, so nothing here can be waived.".format(path))
    return verd


def adjudicate(lost, verd, emit=print):
    """Report every dropped entry and return the ones that fail the gate.

    Two distinct fatal cases, and they are different in kind:
      * the OFF arm's reads VERIFIED it -- verified evidence is leaving where a per-stratum
        count cannot see it, which is the whole hazard this check exists for;
      * it is absent from the OFF arm's verdicts entirely -- it cannot be waived on evidence
        nobody has. Silence is not a pass.
    """
    if not lost:
        emit("  nothing was dropped -- check 2a was purely additive, so there is nothing "
             "to adjudicate.")
        return []
    emit("  {0} dropped hypothes(es), adjudicated against the OFF arm's reads:".format(len(lost)))
    fatal = []
    for k in lost:
        state, why = verd.get(k, (_UNADJUDICATED, ""))
        flag = ""
        if state == FATAL_VERDICT:
            flag = "   <-- FAILS: the OFF arm's reads VERIFIED this"
            fatal.append(k)
        elif state == _UNADJUDICATED:
            flag = "   <-- FAILS: cannot adjudicate, so it cannot be waived"
            fatal.append(k)
        emit("    {0} {1}:{2}  {3} ({4}){5}".format(k[0], k[1], k[2], state, why, flag))
    return fatal


def main(argv):
    if len(argv) != 3:
        raise SystemExit("usage: adjudicate_dropped.py <lost_keys.txt> <off_verdicts.tsv>")
    lost = read_lost_keys(argv[1])
    # Read the verdicts even when nothing was dropped: a malformed or empty OFF-arm verdicts
    # file is a broken gate either way, and finding that out only on the runs that happen to
    # drop something is how a harness defect survives (sec 5.3).
    verd = read_verdicts(argv[2])
    fatal = adjudicate(lost, verd)
    if fatal:
        print("  FAIL: {0} dropped hypothes(es) carried VERIFIED evidence, or could not be "
              "adjudicated at all.".format(len(fatal)))
        print("  A per-stratum count cannot see this, because the entry is absent from both "
              "sides of it.")
        return 1
    if lost:
        print("  PASS: no dropped hypothesis was VERIFIED, so no verified evidence left silently.")
        print("  (The COUNT above is reported and not gated -- Corinne 2026-08-07. A change that "
              "drops")
        print("  many un-verified but real hypotheses would pass here, which is why the number "
              "is shown.)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
