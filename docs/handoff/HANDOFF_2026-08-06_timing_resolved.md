# xTea2 handoff — 2026-08-06, the unexplained timing is explained, and it cost no compute

**This is the current entry point for state.** The agenda entry point is still
`docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md`; the gate's design entry point is still
`HANDOFF_2026-08-04_stage3_td_gate_spec.md`.

> Supersedes `HANDOFF_2026-08-05_decoy_concatenation.md`. **Its §5 items 2 and 3 — two of the three
> things a default flip would have to resolve — are CLOSED.** Item 1, the instrument conflict at
> `chr17:45585847`, is untouched and is still Corinne's decision. §6 NEXT items 1 and 2 are
> unchanged and re-listed in §5 below.
>
> **§10's four decisions (2026-08-04) remain the standing agreement. Nothing here reopens one.**

**The one-sentence version: the gate said "one arm on one node has no control" and filed a rerun;
the control was already inside the job, and it convicts the change of 55s rather than the 132s the
phase number showed.**

---

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| commits this session | this handoff + `docs/handoff/runs/probe_step_timing.py` |
| **`xtea/` did NOT change this session** | the last `xtea/` change is still `ba83d7e`. **BASELINE-0 and every frozen harness baseline are untouched, and no gate needs re-running** |
| tests | unchanged — no `xtea/`, `benchmark/` or `tests/` file was touched |
| jobs | **none. This session ran no compute at all**, which is the point of §1 |
| what is NOT done | **the default still stays OFF**, on ONE blocker now instead of three (§4) |

## 1. The rule that made this free, for the third time

`HANDOFF_2026-08-05_decoy_concatenation.md` §5 item 2 filed the timing as *"One rerun settles it"* —
24 cores, ~45 minutes of priopark, against a number the gate itself called un-actionable because
**"one arm on one node has no control."**

That disclaimer is true of the *phase* number and false of the *job*. Three facts already in the
logs turn one arm into a controlled experiment:

1. **`--parallel-downstream` runs the three families concurrently in one process tree on one
   node.** L1, Alu and SVA share a machine, a filesystem and a minute-by-minute load.
2. **A family whose bwa reference is byte-identical between the arms is a NULL ARM** — same node,
   same job, same phase, provably identical work. The concatenation gate's own check 1a already
   proved L1's reference byte-identical. *The gate established the control and then reported the
   timing as uncontrolled.*
3. **`x_transduction.td_phase` already writes a 13-line decomposition of the transduction step**
   to each family's `.steps.log` (`[TD-PHASE]`), so "the step got slower" resolves to *which*
   of `flank_bwa_index` / `native_gather` / `clip_bwa_realign` / `clip_sam_parse` / … got slower.

**This is the third time an intermediate has outlived a handoff's belief that a measurement needed
a run** — `probe_nearby_clip` found the disc-to-consensus SAM, gate [A] found
`temp_transduction_clip.sam`, and now the timing decomposition was in the step logs. Plan §10's
"check the intermediates before scoping a run" now has a companion: **check the logs too, and ask
which arm of the job is already a control.**

## 2. The measurement

`docs/handoff/runs/probe_step_timing.py --off-arm … --on-arm …`, offline, seconds, no job. Run on
job `49792012` (OFF) vs `49879677` (ON).

### 2.1 `downstream_parallel` is SVA's serial chain, so only SVA's steps can move it

| arm | SVA chain | reported phase |
|---|---|---|
| OFF | 63.0 + 354.3 + 51.0 + 2.2 + 0.9 = **471.4s** | 471.5s |
| ON | 66.8 + 433.8 + 58.7 + 30.0 + 14.2 = **603.5s** | 603.6s |

The phase is a **max over families, not a sum over steps**. Worth stating because the per-step
delta table sums to far more than 132.1s, and summing it over-attributes.

### 2.2 The node is exonerated by L1, and it is exonerated *hard*

L1's reference is byte-identical, so L1 is doing identical work in the same phase on the same node:

| L1, the null arm | OFF | ON | delta |
|---|---|---|---|
| `TOTAL_call_candidate_transduction_v3` | 169.95 | 165.60 | **−4.35s** |
| `clip_bwa_realign` | 75.97 | 72.49 | **−3.48s** |

**The null arm got FASTER while the arm under test got slower.** A generic slow-node explanation
predicts both slow. It also predicts `scan_pass` slow, and `scan_pass` — 1,600s, the biggest phase
in the run, untouched by the flag — moved **+2.8%**, while `cns_emit` moved **−5.4%**.

### 2.3 …but 58% of the phase delta is still not the flag

| SVA step | OFF | ON | delta | can the flag reach it? |
|---|---|---|---|---|
| `cns_filter` | 63.0 | 66.8 | +3.8 | no — upstream of transduction |
| **`transduction`** | 354.3 | 433.8 | **+79.5** | **the only step it touches** |
| `sibling` | 51.0 | 58.7 | +7.7 | no — **ran uncontended in BOTH arms** |
| `postF_all` | 2.2 | 30.0 | +27.8 | no — **ran uncontended in BOTH arms** |
| `postF_hc` | 0.9 | 14.2 | +13.3 | no — **ran uncontended in BOTH arms** |

Those last three are downstream of a call set that differs **by one call** between these arms, and
the probe's timeline reconstruction shows nothing else was running during them in *either* arm. A
13× and 16× slowdown on an uncontended, flag-independent step is the machine, full stop.

Corroborating, from the same table: **`Alu:transduction` went 2.2s → 41.9s (+1800%)** on a family
whose VCF is byte-identical — and `HG002.wgs_Alu.steps.log` contains **no `[TD-PHASE]` lines at
all**, so that step does no transduction work whatever. 41.9 seconds of nothing.

### 2.4 What the change actually costs: +55.4s, and it is the two bwa passes

| SVA `td_phase` | OFF | ON | delta |
|---|---|---|---|
| `flank_build` (the `cat` itself) | 0.08 | 0.13 | +0.05 |
| `flank_bwa_index` | 6.58 | 5.73 | **−0.85** |
| `native_gather` | 53.15 | 59.45 | +6.30 *(flag-independent — reads the BAM)* |
| **`clip_bwa_realign`** | 125.95 | 164.15 | **+38.20** |
| `clip_sam_parse` | 14.76 | 12.35 | **−2.41** |
| **`disc_bwa_realign`** | 101.33 | 118.51 | **+17.18** |
| `disc_sam_parse` | 7.35 | 6.57 | **−0.78** |
| `TOTAL_call_candidate_transduction_v3` | 317.83 | 375.23 | +57.40 |

```
downstream_parallel delta                      : +132.1s (28.0%)
  attributable to the flag (the 2 bwa passes)  :  +55.4s (42%)
  NOT reachable by the flag                    :  +76.7s (58%)
    of which, uncontended in both arms         :  +48.8s
```

**Quote it as ~55s on a 472s phase — 12% — not as 28%.**

## 3. The mechanism, and why the gate's own guess missed it

The gate reasoned *"the bwa index grows by 45 records out of 5,778, which cannot plausibly account
for it."* **That reasoning is correct and the conclusion was still wrong**, because record count is
the wrong variable twice over:

- **By bases the reference grows +0.20%** (17,331,894 → 17,367,290) — even smaller than by records.
  And `flank_bwa_index` got *faster*, so index construction is not the cost.
- **The 45 records are repeat CONSENSUS sequences** — full-length L1, HERV-K, six SVA, 35 Alu. They
  are seed magnets. A read that used to fail fast now finds many seed hits to extend.

The output volume proves it is work-per-read and not more reads:

| | OFF | ON | |
|---|---|---|---|
| `temp_transduction_clip.sam` | 1,979,331,242 B | 1,978,591,405 B | **−0.04%** |
| `temp_transduction_disc.sam` | 772,545,001 B | 772,187,258 B | **−0.05%** |

**Same reads, same output, +30% alignment time.**

And note the gate's stated prediction — *"if anything expect SVA's disc parse to get slightly
cheaper"* — **was right**: `clip_sam_parse` −2.41s and `disc_sam_parse` −0.78s, real savings from
skipping decoy hits before source extraction. It was swamped by an alignment cost 20× larger that
nobody predicted. A correct prediction about the right mechanism can still be the wrong sign on the
total.

## 4. So §5's three blockers are down to one

| | status |
|---|---|
| **1. the instrument conflict at `chr17:45585847`** | **OPEN — and it is Corinne's, not a session's.** Unchanged from the previous handoff §4.5 |
| 2. the unexplained timing | **CLOSED.** ~55s / 12%, mechanism named, node exonerated by a null arm (§2, §3) |
| 3. the 241 → 227 | **CLOSED.** §4.1 below |

### 4.1 The 241 → 227, closed, and L1 is the control again

From the gate's own `probe_src_{OFF,ON}.txt`, already on disk:

| picked sources | OFF | ON |
|---|---|---|
| **L1** total / reference copies / polymorphic | 159 / 114 / 45 | **159 / 114 / 45 — identical** |
| **SVA** total / reference copies / polymorphic | 346 / 241 / 105 | 330 / **227** / 103 |

SVA loses 16 picked sources = 14 reference-copy + 2 polymorphic. **The alarming outcome the WGS
harness named — `241 → 0`, the block outcompeting the whole library — did not happen; this is
−5.8%.** L1's source classes are byte-identical, which is the same null-arm control as §2.2.

The interpretation the gate wrote in advance stands: *a hit that prefers the SVA consensus over a
specific reference copy was never evidence for that copy.* The net call-level effect of all 16 is
the single false positive already reported. **The residual obligation is documentation, not
measurement:** anything quoting the structural flip's "0 → 241" must now say which configuration it
means — 241 is concat-OFF, 227 is concat-ON.

## 5. NEXT

**In this order.** Items 1 and 2 are unchanged and are decisions; item 3 has been re-scoped by a
correction (§5.1).

1. **Decide the default.** Now blocked on **one** thing, not three: the instrument conflict at
   `chr17:45585847` (previous handoff §4.5). Either `verify_transduction` grows a minimum
   source-distance guard — the caller has `MIN_POLYMORPHIC_SOURCE_DIST = 1000` and the verifier has
   nothing — or the no-regression rule gains a case for a change that *removes* evidence. **Both
   are §10-kind decisions.** Note a guard is not free: it moves the frozen 3-arm baseline, which
   contains that entry. The timing objection is gone (§2) and the 241 objection is gone (§4.1).
2. **Gate [A]'s two semantic decisions, in this order** — decisions, not measurements;
   `HANDOFF_2026-08-05_gate_a_measured.md` §2 is the measurement they need:
   - **Separate `n_polya_cutoff` from `n_clip_cutoff`.** At the shipping operating point the polyA
     threshold is *unreachable*, so this should move zero calls and come back byte-identical —
     which makes it a cheap change with a strong gate.
   - **Decide whether polyA-map membership should be a veto at all**, given it drops 68% / 78% of
     the clip-qualified population including 594 / 810 loci with `n_clip >= 10`.
3. **`NEARBY_CLIP` at gate [A]** — still the largest lever measured in the step, but **it is not a
   one-line env change; see §5.1.** Scope it before running anything.
4. **§11's bare-shape writer fix** — still scoped, not built, per §10.3.
5. **§3.4's repeat cutoff by family** — still explicitly NOT the route to justifying anything about
   the decoy fix (§10.1).
6. **A third arm that is neither master nor xTea2** — still the only thing that can take the
   no-regression margin out of PROVISIONAL (§10.4).

### 5.1 A correction to item 3: `NEARBY_CLIP` has no env hook, and it is not scoped to gate [A]

The previous two handoffs say *"the arm is a one-line env change rather than a code change."*
**There is no env hook.** `NEARBY_CLIP = 50` is a plain module constant at
`xtea/global_values.py:114`, and it has **13 read sites across 4 modules**:

| module | sites | what they gate |
|---|---|---|
| `x_orphan_transduction.py` | **6** | the orphan-transduction path — **neither gate [A] nor [B2]** |
| `x_clip_disc_filter.py` | 3 | includes gate [B2] (`:2961`) |
| `x_transduction.py` | 2 | `:916` is gate [A]'s clip parse; `:1103` |
| `transduction_evidence.py` | 1 | a comment reference only |

So changing the constant globally is **not** a scoped gate-[A] arm — it silently moves the orphan
path too, and §4 of the gate-[A] handoff established that the orphan writers are exactly where the
two-"orphan" confusion lives. **An honest gate-[A] arm needs the constant threaded as a parameter
to the clip parse, which is a code change and a Stage 3 decision**, or it needs to accept and
report that it is measuring three call paths at once. That changes what this item costs; it does
not change its priority.

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| **the standing decisions** | `HANDOFF_2026-08-04_td_gate_built.md` §10 — four decisions, unchanged |
| the previous state handoff | `docs/handoff/HANDOFF_2026-08-05_decoy_concatenation.md` |
| **the timing instrument (new)** | `docs/handoff/runs/probe_step_timing.py` — two arms, offline, seconds. Reports which families are **null arms** and refuses to attribute anything when there is none |
| the decoy block | `xtea/td_decoy_library.py`; census probe `docs/handoff/runs/probe_decoy_library.py`; tests `tests/test_td_decoy_library.py` |
| the gates | `docs/handoff/runs/gate_decoy_concat_chr1.sbatch`, `gate_decoy_concat_wgs.sbatch` |
| **"two arms, one flag apart"** | `docs/handoff/runs/td_gate_no_regression.sbatch` — `XTEA2_OFF_ARM` / `XTEA2_ON_ARM` / `XTEA2_OFF_LEDGER=none` / `XTEA2_ARMS_LABEL` |
| the two arms this handoff measures | `gate_decoy_structural_wgs_49792012/HG002.wgs` (OFF) and `gate_decoy_concat_wgs_49879677/HG002.wgs` (ON), both under `/n/scratch/users/c/cos689/xtea2_testdata/` |
| **a trap in the neighbouring harness** | `docs/handoff/runs/master_serial_wgs.sbatch` opens `rm -rf "$OUT"` and its default `OUT` is the NA12878 tree the gate spec cites. Never invoke it without `MASTER_OUT` |
| running the tests | `python3 -m pytest -m "not slow"` on the system 3.9.25 — **254 passed / 30 skipped** |
