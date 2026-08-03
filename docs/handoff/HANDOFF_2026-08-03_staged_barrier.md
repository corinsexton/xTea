# xTea2 handoff — 2026-08-03 late (the barrier stopped being a join)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-03_locus_universe_driver.md`, whose §7 NEXT item **2** is what this
> session did. Its item 1 (the per-locus decision loop) is untouched and is still the front edge.

**One-line status: the `cns_filter → transduction` barrier no longer makes every family wait for the
slowest one.** `run_xtea.StagedBarrier` gives each family only the waits it actually has, and the
three things that run at that point turn out to have three different dependency sets — per-family
enumeration needs nobody, the driver gather needs (and is needed by) only the flanked families, and
only the union table and `--family-assign-early` need all three. The join was measured at **+52.9s at
WGS** (job `49550925`) and none of it was work. **Gated PASS at both scales on one commit
(`d565445`): chr1 6/6 (job `49599837`) and WGS 6/6 (job `49599838`, including `OFF == BASELINE-0`
4/4 and 188/1223/82/1493 calls in all three arms).** The flankless family's wait at the barrier went
**45.8s → 0.1s** and `downstream_parallel` **recovered 41.4s of the join's 46.2s**.
`--global-barrier` restores the old form so the numbers stay a one-flag A/B.

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| the feature commits | `2f757b8` (the staged barrier) and `d3cb7f1` (the assignment's timing row). **No HEAD hash recorded on purpose** — an earlier handoff named one and an `--amend` invalidated it the same day. |
| working tree | clean |
| tests | **198 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py`. Was 187 — **+11** in the new `tests/test_downstream_barrier.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`). |
| in flight | see §3 |
| the agenda | unchanged: the loop's **third** box, the per-locus decision. The previous handoff's §4 still says what shape it can have. |

## 1. What was built

### 1.1 `run_xtea.StagedBarrier` — the barrier, decomposed

The split itself is unchanged (`downstream_waves` still returns two waves). What changed is that the
split is no longer *run* as a global join. Each family's worker now owns its whole chain and calls
`barrier.arrive(t)` between the waves, where it waits for **only what it reads**:

| runs at the barrier | needs | waited on by |
|---|---|---|
| per-family enumeration (`--locus-universe`) | that family's own cns list + raw_disc | **nobody** |
| the driver gather (`--td-gather-driver`) | the **flanked** families' `.slct` | the **flanked** families only |
| the union table | all three | **nobody** — nothing in the run reads it |
| `--family-assign-early` | all three | everybody, correctly |

So **`--locus-universe` alone now waits for nothing at all**, which is what makes flipping it on cost
its ~3s of enumeration instead of the slowest family's `cns_filter`. `--family-assign-early` still
joins, and should: its rule compares every family's cns rows against every other's.

**Who runs the shared work:** whichever thread completes the dependency, not a fourth thread and not
the main one. The gather goes to the last **flanked** family to enumerate (it would otherwise be
sitting in `wait()`); the union table and the assignment to the last family overall. Election is
under the lock and latched, so each runs exactly once.

### 1.2 Why it is answer-preserving by construction, not by luck

**The staged interleaving is the one production already runs.** With the flags off, the one-wave
chain already has L1 and SVA in transduction while Alu is still in `cns_filter` — that is documented
in `main()` as the reason the transduction core split can briefly oversubscribe. Staging reproduces
exactly that pattern. The global join was the anomaly, not the staging.

**The one ordering that could lie is preserved for free.** `--family-assign-early` *rewrites* each
family's cns list to its own winners, so enumerating after it would silently produce a
post-assignment universe — a smaller number that still looks like a universe. Every family enumerates
before it registers, and the assignment cannot be elected until every family has registered, so all
three enumerations still precede it. Same guarantee the two-wave order gave, without the join.

### 1.3 The two things staging introduces, both handled

1. **`XTEA_TD_GATHER_PRECOMPUTED` is set BEFORE any worker starts**, not on the gather thread as the
   global form does. Under a join no subprocess exists when the driver mutates `os.environ`; under
   staging families are spawning steps throughout, and mutating environ concurrently with
   `fork`/`exec` is the one race staging would otherwise add. Setting it early is **free**, because
   the step's test is `env AND both files present` (`x_transduction.py`'s `b_td_precomputed`) — a
   family the gather skipped for having zero sites, or a gather that produced no entries at all,
   still falls back to gathering for itself.
2. **A death before the barrier must release the waiters.** A join has one thing going for it:
   nothing waits on a family that died. Here it does, so every failure path calls `abort()`, which
   latches the failure and sets both events. Without it a non-zero `cns_filter` leaves the other
   families blocked inside the executor's shutdown and the run **hangs instead of failing** — at WGS,
   forty minutes in. Two tests cover it (a death in the step and a death inside the enumeration).

### 1.4 `--global-barrier`

Restores the join. Default off. It exists so the measured numbers in `run_td_gather_driver` stay
reproducible as a one-flag A/B rather than a checkout, and so the gate has a control arm. No effect
at all on a run with no barrier, which is the default one.

### 1.5 Driver log lines now carry elapsed seconds

`[xtea2 +  123.4s] …`. `step_timings.tsv` records *durations*, so two rows that overlap and two that
abut are indistinguishable in it — and under `--parallel-downstream` overlap is the whole question.
Nothing greps the `[xtea2]` prefix (checked across `docs/handoff/runs/` and `tests/`); every reader
matches on the message, so existing gate scripts are unaffected.

## 2. The unit gate, and the mutations that prove it bites

`tests/test_downstream_barrier.py`, 11 tests, all driving real threads with a deadline — a
coordination bug is a hang, and a hang with no deadline is a suite that never finishes.

Three mutations were applied to confirm the tests fail for the right reason:

| mutation | test that caught it |
|---|---|
| every family waits for the gather (the bug being removed) | `test_flankless_family_never_waits_for_the_driver_gather` |
| flanked families skip the gather wait (job `49498881`'s failure mode) | `test_every_flanked_family_transduces_after_the_gather` |
| `abort()` no longer releases the events | `test_a_failing_enumeration_aborts_from_inside_arrive` |

## 3. The end-to-end gates

### 3.0 The first chr1 run FAILED, and it failed on the one thing only it could see

Job `49596361`, at `2f757b8`. Five checks of six passed: every call identical across all three arms
(9/9 artifacts, `L1=1 Alu=14 SVA=2 merged=17`), both `.slct` replicas byte-identical to the steps'
own in **both** barrier arms, the bypass firing 6/6, and check 5 clean.

**Check 2 failed**: `locus_universe.tsv` differed between GLOBAL and STAGED over an identical 3,272
loci with identical claims —

```
GLOBAL   chr1  180998  2  L1,SVA  new,new
STAGED   chr1  180998  2  SVA,L1  new,new
```

`LocusUniverse.rows()` writes each locus's `families`/`claims` pairs in the iteration order of the
map it was built from. The global join built that map by walking `type_ctx`; the staged barrier built
it as each family landed, so **the file's bytes depended on whose `cns_filter` finished first**.
Fixed at the barrier (`d565445`) by ordering the map canonically — *not* by sorting inside
`LocusUniverse`, which would silently rewrite every table already gated, including job `49550925`'s.

**Two things worth keeping from this.** First, the reason check 2 exists is exactly the reason it
fired: the universe pass is diagnostic, it cannot move a call, so a staged enumeration producing a
different — or empty — table sails through a call diff in silence. **A concurrency change's failure
mode is usually a byte, not a call.** Second, it should never have needed a cluster job:
`test_the_union_sees_families_in_driver_order_not_completion_order` now forces SVA to enumerate first
and catches it in two seconds.

### 3.1 `docs/handoff/runs/gate_staged_barrier_chr1.sbatch` — three arms, one commit

OFF / GLOBAL (`--td-gather-driver --global-barrier`) / STAGED (`--td-gather-driver`), 18 min.

**RESULT at `d565445`, job `49599837`: PASS, 6/6.**

| check | result |
|---|---|
| 0. arms ran the topologies they claim | OFF 0/0, GLOBAL universe+gather+**global-line**, STAGED universe+gather+**staged-line** |
| 1. GLOBAL vs OFF, STAGED vs OFF — no call may move | **9/9 identical each**, 0 new |
| 2. STAGED vs GLOBAL, the universe artifacts | **3/3 identical** — table 99,665 B / 3,272 rows, L1 2,815 rows, SVA 3,205 |
| 3. driver replica == step `.slct`, in **both** arms | **4/4 identical** |
| 4. bypass fired, per flanked family, all three arms | **6/6** |
| 5. the flankless family stopped waiting | **GLOBAL 20.2s → STAGED 0.1s**, against a 20.06s gather |
| 6. call counts agree | all three **1 / 14 / 2 / 17** |

**Check 5 is the point of the change**: Alu waited the driver gather's entire duration under the join
and essentially nothing under staging. The order-dependent observation agrees — in the STAGED arm the
first transduction starts at `+247.4s`, *before* the last enumeration at `+247.8s`; in GLOBAL it
starts at `+282.8s`, 20s after the last enumeration at `+262.7s`.

**chr1 cannot show the magnitude, and did not:** `downstream_parallel` went 109.6 / 104.8 / 100.6s
across OFF / GLOBAL / STAGED. chr1's three `cns_filter`s are 17.2 / 18.1 / 17.3s — near-identical, so
there is almost no imbalance to recover, and the 4s spread is not worth reading. The whole cost being
removed is Alu's WGS `cns_filter` being ~2× the others', which only exists at WGS. That asymmetry is
why §3.2 carries the magnitude check and this one carries the structure.

### 3.2 `docs/handoff/runs/gate_staged_barrier_wgs.sbatch` — three arms, ~2h

Same three arms at WGS, plus the `OFF == frozen BASELINE-0` leg (without it a clean STAGED-vs-GLOBAL
diff could be two equally-wrong arms) and the **magnitude** check. Submitted with
`--dependency=afterok` on the chr1 job, so two hours of compute cannot start behind a failing chr1.

**RESULT at `d565445`, job `49599838`: PASS, 6/6.**

| check | result |
|---|---|
| 0. arms ran what they claim | OFF 0/0; GLOBAL and STAGED each one universe + one gather, right barrier line |
| 1. **arm OFF == frozen BASELINE-0** | **4/4 identical** |
| 2. GLOBAL vs OFF, STAGED vs OFF | **9/9 identical each**, 0 unexpected new |
| 3. STAGED vs GLOBAL, universe artifacts | **3/3 identical** — table 1,264,244 B / **41,604 loci**, L1 793,811 B / 34,197 rows, SVA 915,815 B / 39,452 rows |
| 4. replica == step `.slct` (both arms) + bypass fired | **4/4** and **6/6** |
| 5. the flankless family stopped waiting | **GLOBAL 45.8s → STAGED 0.1s**, against a 45.49s gather |
| 6. magnitude — `downstream_parallel` | OFF **511.75s**, GLOBAL **557.98s**, STAGED **516.59s** → **recovered 41.4s** |
| calls | all three **188 / 1223 / 82 / 1493** |

Leg 3's byte counts are **identical to the offline gate's and to job `49550925`'s**, and the universe
reproduces exactly: 41,604 loci, 7,910 / 33,139 / 555 by one / two / three families, **38,815
created-only**.

#### The magnitude, read honestly

**This node is ~10% slower than job `49550925`'s** — OFF's `downstream_parallel` is 511.75s where
that job measured 462.7s, and its `cns_filter`s are 65.5 / 155.1 / 70.6s against 57.3 / 124.1 /
62.5s. So **only the within-job arm deltas are comparable**, which is why all three arms run in one
job on one node. Do not compare these absolutes to the previous handoff's.

- **The join cost 46.2s here** (GLOBAL 557.98 − OFF 511.75), consistent with the 52.9s measured on
  the faster node.
- **Staging recovers 41.4s of that — about 90%.**
- **The residual is +4.8s against OFF, not the −13.7s the hoist alone would predict.** The likely
  cause is the one §1.3 flagged: the driver gather now runs while Alu is still in `cns_filter` and
  contends for cores the join used to hand it exclusively. STAGED's Alu `cns_filter` is 171.4s
  against OFF's 155.1s, which is the shape that would produce.
  **State this as the likely cause, not a measured one.** The same table shows GLOBAL's L1 and SVA
  `cns_filter`s at 83.4 / 88.9s against OFF's 65.5 / 70.6s — +18s in an arm where nothing touches
  `cns_filter` at all. So run-to-run variance on these per-family numbers is ±~18s, and the +16.3s
  attributed to contention sits inside it. **Only check 5's 45.8 → 0.1s is a direct measurement of
  the wait**; everything else here is a critical-path inference.
- **Total run wall (2398 / 2375 / 2264s) is not the place to read any of this** — `scan_pass` varies
  by ~66s between arms and no barrier change can touch it.

### 3.3 The one check worth copying: measure the flankless family's own gap

The tempting concurrency assertion — "did some family cross the barrier before the slowest
`cns_filter` landed" — is real but **order-dependent**: whether it holds depends on which family
happens to be slowest and on whether the gather outruns the gap behind it, neither of which is a
property of the code. *A gate that passes for a reason the change does not control is worse than no
gate.*

The structural version is the flankless family's own gap, from its enumeration line to its wave-2
announcement. Under staging that contains only the union-table write; under the join it must contain
the entire driver gather — **on any ordering**. That is check 5 in both scripts.

## 4. What NOT to redo

Plan §10 is the authority and gained the entries below. Everything in the previous handoff's §6 still
stands unchanged.

- **The cost of the barrier was lost OVERLAP, and it has been taken back.** Do not re-derive the
  +52.9s as an open problem — it is `run_xtea.StagedBarrier`, default on whenever a flag opens the
  split. The old numbers still describe `--global-barrier`.
- **The staged interleaving is not new behaviour to be feared** — it is the interleaving the one-wave
  production chain already has. The join was the anomaly.
- **Do not "simplify" the env-var set back onto the gather thread.** It is early on purpose (§1.3).
- **Do not remove `abort()`** — without it a failed `cns_filter` hangs the run instead of failing it.
- **A concurrency check must not be order-dependent** (§3.3).
- **`--family-assign-early` still joins and should.** Staging saves nothing there and does not
  pretend to.

## 5. NEXT

1. **The per-locus decision loop** — unchanged, and now the only thing left in the loop's front half.
   The previous handoff's §4 has its record type (`collect_ins_info`'s 11 fields plus the 10-field
   genotype tuple — do not invent one) and the two things still unknown about its tail.
2. **Reconsider flipping `--td-gather-driver` on.** The previous handoff's item 3 said "item 2 first,
   or the flip ships a 53s regression". Item 2 is done, so the arithmetic has changed sign: the hoist
   itself is **−13.7s and ~66s of CPU returned**, and the barrier that made it a net +42.6s is gone.
   The flip is now defensible on its own merits rather than only as a prerequisite for the loop.
   **Left as a decision, not taken** — it changes the default topology, and the original reasoning
   (flip when the loop consumes it) is about more than wall clock.
3. **The last ~18s: give the driver gather fewer threads while a family is still in `cns_filter`.**
   Staging recovered 41.4s of the join's 46.2s and landed at **+4.8s vs OFF** where the hoist alone
   predicts −13.7s. The likely cause is that the gather now runs concurrently with Alu's `cns_filter`
   at `--threads opt.cores`, contending for cores the join used to hand it exclusively (§3.2).
   **Measure before building**: per-family `cns_filter` varies by ±~18s between arms of the same job,
   so the +16.3s this is inferred from is inside the noise. The cheap experiment is one extra WGS arm
   with the gather at `opt.cores − dcores`, not a redesign.
4. **`--emit-per-type`**, still deliberately last: flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
5. **Stage 3**, whose gate is fully specified (`HANDOFF_2026-08-02_stage2c_barrier.md` §1).

## 6. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§10 gained the staged-barrier block) |
| the barrier | `xtea/run_xtea.py` — `StagedBarrier`, `enumerate_family_universe`, `write_locus_universe_table`, `timed`, `_worker_staged` |
| its unit gate | `tests/test_downstream_barrier.py` |
| its end-to-end gates | `docs/handoff/runs/gate_staged_barrier_chr1.sbatch`, `…_wgs.sbatch` |
| the loop's front half | `xtea/downstream_unified.py` (unchanged this session) |
| everything else | see `HANDOFF_2026-08-03_locus_universe_driver.md` §8 |
