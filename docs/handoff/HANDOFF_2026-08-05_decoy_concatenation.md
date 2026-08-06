# xTea2 handoff — 2026-08-05 (third session), the decoy-consensus concatenation is built and gated

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-05_gate_a_measured.md`. **Its §6 NEXT item 1 — the decoy-consensus
> concatenation, the top of that list — is BUILT.** Items 2-6 are unchanged and are re-listed in §6
> below in the same order.
>
> **§10's four decisions (2026-08-04) remain the standing agreement. Nothing here reopens one.**

**The one-sentence version: a comment that has been false on hg38's SVA path since master can now
be made true, and the fix needed no new data — only a recipe, verified as bytes.**

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | `ba83d7e` the module + wiring + tests + census probe; `7f313e6` the chr1 gate; `10a247a` the WGS gate + the generalised no-regression harness; `05f23ae` the plan's record; `RESULTS_COMMIT` the gate results |
| **`xtea/` DID change this session** | `xtea/td_decoy_library.py` (new), `xtea/x_transduction.py`, `xtea/run_xtea.py`, `xtea/global_values.py`. The last `xtea/` change before this was `53ca400`. **BASELINE-0 and every frozen harness baseline are still reproducible** — the new flag is OFF by default and absent-means-OFF on both entry points |
| tests | **254 passed / 30 skipped** (`python3 -m pytest -m "not slow"`, system 3.9.25), was 234/30. 20 new in `tests/test_td_decoy_library.py` |
| jobs | `49879284` chr1 three-arm gate, **PASS**, 3 arms at ~6 min each. `49879677` WGS one-arm gate — RESULT_WGS |
| worktrees | `/n/scratch/users/c/cos689/xtea2_wt/td_concat_457eb32` (BEFORE) and `td_concat_ba83d7e` (AFTER). Neither has a built `xtea_scan`; both gates borrow `td_decoy_55055ba`'s, and `git diff --name-only 11b0e61 HEAD -- src/` is empty |
| what is NOT done | the default is still OFF, deliberately — see §5 |

## 1. What the change is, and what was actually wrong

`construct_novel_flanks` builds the bwa reference **both** transduction realignment passes use —
the clip pass (`x_transduction.py:272`) and the disc pass (`:299`). The comment above both of its
call sites has always said:

    ####Note, we add L1, Alu and SVA consensus to the flank regions to filter out those mapped to consensus

**It never did.** It cats the packaged flank library with the polymorphic flanks it just built and
nothing else. So whether the step has a decoy channel at all is a property of how that build's
library was **packaged**, and the packaging is not uniform
(`docs/handoff/runs/probe_decoy_library.py`, offline, seconds, no job):

| family | build | records | decoy records |
|---|---|---|---|
| LINE | hg38 | 707 | **45** |
| LINE | hg19 | 705 | **45** |
| LINE | chm13 | 818 | **0** |
| SVA | hg38 | 5454 | **0** |
| SVA | hg19 | 5158 | **44** (all but HERV-K) |
| SVA | chm13 | 6072 | **0** |

**hg38 SVA is the path every measurement in this project is taken on**, and it has none. That is
the same hole two independent instruments already showed: the disc side's `decoy_cross_family`
firing 144× for L1 and 0 for SVA, and the clip side's **3,992 of 99,245** qualified L1 clip records
resting on decoy-consensus support against SVA's **0 of 128,262**
(`HANDOFF_2026-08-05_gate_a_measured.md` §2.4). Plan §5's SVA/Alu mapping hazard — SVA is
SINE-VNTR-**Alu** — is untested rather than measured-absent because of it.

## 2. The recipe is BYTES, not an inference

`xtea/td_decoy_library.py` rebuilds the 45 records from `<rep_lib>/consensus/`:
`transduction_decoy.fa`'s 3 (`>LINE1`, `>SVA`, `>ALU`), `HERV.fa`'s 1, `SVA.fa`'s 6, and `ALU.fa`'s
35 under an `ALU_` prefix.

**The emitted text is byte-identical to the TAIL of both packaged L1 libraries** — each source
file's own line wrapping (60 / 120 / 70 / 71 columns) and the header comment fields
(`>SVA_A  SINE    Homo sapiens`) included. hg19's SVA 44 are the same set less HERV-K, in a
different order. That is `tests/test_td_decoy_library.py`'s central test, and it runs against the
shipped libraries rather than a fixture; it skips when rep_lib is absent, it does not weaken.

**Records are taken by NAME, not by a positional slice.** Both source files continue past the cut —
`ALU.fa` has 47 records and `SVA.fa` 88 — so a positional cut would be silently wrong if a
consensus file were ever re-ordered or extended, and this block ends up in a bwa reference nothing
downstream re-checks. `SVA` and `ALU` are taken from `transduction_decoy.fa` and not from
`SVA.fa`/`ALU.fa`, which also carry records of those names; that is the packaged libraries' choice,
and the byte-identity test is what proves it is theirs and not a guess.

## 3. Three implementation decisions, each of which a gate tests

1. **It appends by NAME, so it is a structural no-op on both L1 libraries** and adds exactly
   HERV-K on hg19's SVA. A duplicate name in a bwa reference makes `reference_name` ambiguous, and
   that string is what the entire decoy classification tests.

   **This is what turns "L1 must not move" from an argument about the code into a prediction about
   the DATA** — and it is a different argument from the structural fix's, which rested on L1 flank
   ids colliding with no decoy tag. Both gates check it **on the reference file**, upstream of
   every gate, because a VCF comparison cannot tell "the reference was identical" from "the
   reference differed and nothing survived to the call set".
2. **The block goes between the library and the polymorphic flanks**, which is where hg19's SVA
   library carries its own 44. An appended block reproduces a packaged layout rather than inventing
   a third.
3. **The step reports the COUNT it appended**, not only that the policy arrived
   (`[td-decoy-concat] N of 45 decoy records appended to ...`). On an L1 run "armed" and "changed
   anything" are different facts, and the decoy *default* flip had to be rescued from exactly that
   confusion.

`run_xtea --td-decoy-concat` / `XTEA_TD_DECOY_CONCAT`, **OFF by default and absent-means-OFF on
both entry points** — deliberately unlike `--td-decoy-structural`, whose driver default was flipped
only after a measured no-regression arm.

## 4. The gates

### 4.1 chr1, three arms — **PASS**, job `49879284`, 6:03 / 5:49 / ~6 min

`docs/handoff/runs/gate_decoy_concat_chr1.sbatch`, BEFORE (`457eb32`) / AFTER_OFF / AFTER_ON.

**Containment:** AFTER_OFF is byte-identical to BEFORE over the whole output tree — **100 identical,
0 differing, 0 missing**, plus the two characterised exceptions (`cns_manifest.txt`'s embedded arm
paths, `scan_disc_pos.tsv`'s completion-order permutation over 302,912 rows). No unexpected result
file in either after arm.

**The data-level prediction, confirmed at pipeline scale and at the bwa reference itself:**

| | L1 | SVA |
|---|---|---|
| `[td-decoy-concat]` appended | **0 of 45** | **45 of 45** |
| reference records OFF → ON | 709 → **709**, and **byte-identical** | 5464 → **5509**, exactly **+45** |
| block position | n/a | verbatim, followed by `>chr1~64384507~64384507~polymerphic~0L` |

So "L1 does not move" is here a statement about the **data** — its packaged library already carries
the block, so `write_missing_decoys` selects nothing — and it is checked *upstream of every gate*
rather than inferred from an unchanged VCF.

**Nothing moved at chr1, exactly as the header predicted.** 0 files differ in AFTER_ON (L1, Alu,
SVA and cross-family all zero), both ledgers are identical across all three arms (L1 1,695 rows,
SVA 2,315, same reason split), and call counts are L1 1 / Alu 14 / SVA 2 in every arm. chr1's gate
B SVA survivors are 0, so this is the population limitation the structural gate's header already
records — **not** evidence that the concatenation does nothing.

### 4.2 THE ONE NUMBER THAT DID MOVE, and what it says about which side can move calls

`probe_clip_gate.py` on each arm's own surviving clip SAM, at the shipping operating point (both
arms reproduce their run's ledger `A_clip` split *and* its `.unique_trsdct_half_clip` file as a set,
so this is a measurement):

| | gate [A] survivors | with ANY decoy-consensus clip support | with ONLY that |
|---|---|---|---|
| L1, BEFORE and AFTER_ON | 217 (both) | **33** (both) | 0 |
| SVA, BEFORE | 270 | **0** | 0 |
| SVA, AFTER_ON | **270** | **19** | 0 |

**SVA has a decoy channel on the clip side for the first time — and gate [A]'s verdicts did not
move at all.** 270 survivors both ways, byte-identical ledgers. That is not a coincidence and it
is worth carrying forward:

> **On the clip side the concatenation changes what a read maps TO, not whether it counts.** The
> clip path does not skip decoy hits (the skip is commented out, `x_transduction.py:923-925`), so a
> record that used to hit a reference flank and now hits `>SVA` is still one qualified clip record
> at the same locus — it just carries a degenerate source. The one asymmetry: a decoy hit has no
> `~`-fields, so it **bypasses the `MIN_POLYMORPHIC_SOURCE_DIST` self-source exclusion** at
> `x_transduction.py:947-948`, which can only *add* records. At chr1 it added none that changed a
> verdict.
>
> **So the side that can move SVA calls is the DISC side**, where `parse_transduction_disc_algnmt`
> skips a decoy hit **before the source is extracted** and tallies it into `gate_decoy_cross_family`.
> Read the WGS arm's SVA delta as a disc-side result and check §4's ledger `decoy_cross_family`
> counts first.

### 4.3 WGS, one arm — job `49879677`

RESULTS_WGS_SECTION

## 5. What is NOT decided

**The default.** It stays OFF. RESULTS_DEFAULT_NOTE

**Direction was not predicted, and that was written down before the run.** The concatenation adds
competition, so its expected effect is to REMOVE spurious sources — a read that came from SVA or
Alu consensus body now maps to the consensus record instead of to some reference copy's flank.
That is precision-shaped and can show up as lost calls. Plan §4's attribution (306 of BASELINE-0's
322 false negatives never reach this step) is why sensitivity is not what is on the table.

## 6. NEXT

**In this order.** The previous handoff's items 2-6, unchanged, with the concatenation removed
from the top because it is built.

1. **Gate [A]'s two semantic decisions, in this order.** They are decisions, not measurements, and
   `HANDOFF_2026-08-05_gate_a_measured.md` §2 is the measurement they need:
   - **Separate `n_polya_cutoff` from `n_clip_cutoff`.** One dial for two thresholds is the §5
     defect; separating them moves calls, so it is a flag and an arm.
   - **Decide whether polyA-map membership should be a veto at all**, given that it drops 68% / 78%
     of the clip-qualified population including 594 / 810 loci with `n_clip >= 10`. Expect
     *precision*: 306 of BASELINE-0's 322 FNs never reach this gate.
2. **`NEARBY_CLIP` at gate [A] is an open question**, and it is the largest lever measured in the
   step. It needs a truth-set-scored two-arm run, and the arm is a one-line env change. **Do not
   raise it on the 4× alone** — read that handoff's §3 for the two reasons the 4× is not 4× more
   evidence.
3. **§11's bare-shape writer fix** — still scoped, not built, per §10.3. The disjointness test
   fails when it lands, deliberately.
4. **§3.4's repeat cutoff by family** — still explicitly NOT the route to justifying anything about
   the decoy fix (§10.1). Argue it on biology, in isolation, in its own change.
5. **A third arm that is neither master nor xTea2** — still the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4).

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the standing decisions** | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four decisions, unchanged |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-05_gate_a_measured.md` |
| **the decoy block** | `xtea/td_decoy_library.py`; the census probe is `docs/handoff/runs/probe_decoy_library.py`; the tests are `tests/test_td_decoy_library.py` |
| the seam | `xtea/x_transduction.py` `construct_novel_flanks` — one `cat`, three files when the flag is on |
| **the gates** | `docs/handoff/runs/gate_decoy_concat_chr1.sbatch` (three arms) and `gate_decoy_concat_wgs.sbatch` (one arm vs job `49792012`) |
| **"two arms, one flag apart"** | `docs/handoff/runs/td_gate_no_regression.sbatch`, now general: `XTEA2_OFF_ARM` / `XTEA2_ON_ARM` / `XTEA2_OFF_LEDGER=none` / `XTEA2_ARMS_LABEL`. Its defaults are still the structural flip's pair |
| the frozen OFF arm for this change | `/n/scratch/users/c/cos689/xtea2_testdata/gate_decoy_structural_wgs_49792012/HG002.wgs` — **not BASELINE-0**, which predates the 2026-08-05 structural flip |
| gate [A]'s instrument | `docs/handoff/runs/probe_clip_gate.py` — also the before/after on the clip-side decoy channel |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 — **254 passed / 30 skipped** |
