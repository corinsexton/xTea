# xTea2 handoff — 2026-08-03 (site enumeration is solved; §3's diagram order is not)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it; this is the index.
> Supersedes `HANDOFF_2026-08-02_td_multi_and_site_creation.md`, whose §4 NEXT item 1 is the per-locus
> loop. This session did **not** build the loop — it discharged the first of the two dependencies that
> item names, and in doing so found that the architecture diagram the loop would have been built from
> is wrong.

**One-line status: the 52 of 1,515 final records that do not exist when the family barrier runs are
`created`, not moved — every postF_hc position is exact in `cns2` and the nearest cns row is ≥3,293 bp
away, median 2.4 Mb — and every one of them is already a row in
`candidate_list_from_disc.txt.clip_sites_raw_disc.txt`, which the disc adapter writes **before**
`cns_filter`. The selection that turns that file into transduction's candidate list needs no BAM, no
bwa and no cns evidence, so it is now `xtea/site_universe.py`, gated **offline at WGS scale, 4/4
byte-identical** against a real run's own `.slct`/`.slct2` for both production callers. Two
consequences: the per-locus loop's iteration domain is **41,604 union loci** (vs 75,853 per-family
visits, 1.82×), and **plan §3's order — assign family, then loop — cannot work**, because for those 52
records the evidence *is* the transduction evaluation's output. The 7-of-22 ceiling is unchanged.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `03711ea` |
| working tree | clean; all work committed |
| tests | **173 passed, 1 skipped** with `XTEA2_SCAN_BIN` set, excluding `test_downstream_from_cns.py` (the slow chr1 golden). Was 162 — +11 in the new `test_site_universe.py`. The skip is `test_e2e_concordance` (needs `XTEA2_OUTPUT_DIR`) |
| in flight | **nothing.** The gate is offline and already run — `docs/handoff/runs/gate_site_universe.py`, PASS 4/4. No sbatch this session, and none needed: this change touches no BAM path |
| the agenda | the per-locus loop, now with one fewer unknown and one corrected premise (§4) |

## 1. What was measured, and why the first measurement matters most

`probe_site_creation.py` (previous session) established that **52 of 1,515 per-family final records
have no `candidate_disc_filtered_cns.txt` row** — 50 L1, 2 SVA, all from the transduction step. It
matched positions with a ±50 bp window, justified in its own docstring by "the chain shifts a call's
coordinate". That justification admits an alternative reading of the whole result: *maybe those 52 are
ordinary cns rows that the chain moved more than 50 bp.* If so, "unified site creation" would be the
wrong name for the problem and a wider window would be the fix.

**It is not that.** Measured on the same run's kept intermediates:

| L1 postF_hc positions matched within … | cns | cns.hc | cns2 | cns2.hc |
|---|---|---|---|---|
| 0 bp (exact) | 150 | 150 | **200** | **200** |
| ≤ 50 bp | 150 | 150 | 200 | 200 |
| ≤ 10 kb | 151 | 151 | 200 | 200 |

Every final position is **exact** in `cns2` — the post-cns chain does not re-anchor coordinates at
all, so on this hop the ±50 window does no work — and the 50 unmatched records' nearest cns row is
**≥ 3,293 bp away, median 2.4 Mb, max 20.6 Mb**. They are created. Alu: 0 of 1230. SVA: 2 of 85, same
shape.

**Then: where does the creating step get them?** `x_TEA_main`'s `--transduction` branch calls
`re_slct_with_clip_raw_disc_sites(sf_raw_disc, cns_out, …)` → `<raw_disc>.slct`, and evaluates that
list. Checking each of the 52 against every candidate file in the type directory:

| present in | L1 (of 50) | SVA (of 2) |
|---|---|---|
| `<disc>.clip_sites_raw_disc.txt` | **50** | **2** |
| `<disc>.clip_sites_raw_disc.txt.slct` | **50** | **2** |
| `candidate_list_from_disc.txt` | 5 | 0 |
| `…cns2.txt.…new_sites` | 2 | 2 |

**100% in raw_disc, which `disc_adapter.run_disc_from_scan` writes before `cns_filter` runs**
(`xtea/disc_adapter.py:225`). Site creation is a selection over a file that already exists at the
barrier. Note the `new_sites` column is 2 and 2 — that file is a much narrower thing (orphan sites,
`TD_SRC` ending `~orphan`) and is **not** where the 52 come from. Anyone sizing site creation off
`new_sites` will conclude it is 4 records and be wrong by an order of magnitude.

## 2. What was built: `xtea/site_universe.py`

The rule, lifted verbatim out of `XTransduction.re_slct_with_clip_raw_disc_sites`, which is now a
**delegation and nothing else** — the same "one copy, so it cannot drift" discipline
`locus_table.assign_type` uses. `is_matched_rep_type` moved with it, so
`fall_in_low_div_same_type_rep` shares the one copy of the subfamily-name test.

Its inputs are the whole point: **raw_disc, the cns candidate list, the family's rmsk annotation, the
disc cutoff, the rtype.** No BAM. No bwa. No cns *evidence* — only the cns positions, for the
proximity exclusion. So a driver can run it for every family before any per-family chain starts.

### 2.1 How it is gated

**Unit — `tests/test_site_universe.py`, +11.** Its centre is
`test_library_reproduces_the_replaced_implementation`: `_legacy_re_slct` in that file is a **verbatim
copy of the deleted loop**, including its use of the real `XIntemediateSites.are_sites_close`, and the
two are compared file-for-file over randomized inputs — positions drawn densely on the proximity
window's edges, both rtypes, both `b_tumor` values, three cutoffs. That is what makes the half-open
window claim (§5) a *verified property of the original* rather than a reading of it. The other ten
tests pin one rejection each, so a failure says which.

**WGS scale — `docs/handoff/runs/gate_site_universe.py`, PASS 4/4 byte-identical.** Offline: reads the
2c gate's kept intermediates, no BAM, no sbatch, ~2 min. The unit tests are synthetic — a fake
annotation, a few dozen sites — and every remaining risk is scale- or data-shaped: 978,710 real L1
rmsk intervals loaded at extnd 450 / min copy 225, ~39,600 raw-disc rows against 300 candidates over
every contig including decoys.

| | rows | vs the run's own file |
|---|---|---|
| L1 `.slct` — transduction caller, 300 candidates + 33,897 new | 34,197 | **identical**, 793,811 B |
| L1 `.slct2` — sibling caller, `b_with_candidate=False`, extend = mean IS | 33,927 | **identical**, 790,127 B |
| SVA `.slct` — 292 + 39,160 | 39,452 | **identical**, 915,815 B |
| SVA `.slct2` | 39,158 | **identical**, 911,479 B |

**Both production callers, deliberately** — they differ in every argument that could hide a bug: the
sibling one (`x_TEA_main.py:622`) passes `b_with_candidate=False`, a different candidate file
(`…cns2.txt.all_non_sibling_td.txt`) and a different rmsk boundary extend (**mean insert size, a
float**, instead of 450). Alu has no `.slct` and that is a pass, not a skip: no transduction flank, so
the step copies its input through and the rule never runs.

**The gate parses the run's own `*.steps.log` for the argv rather than rebuilding it**, because the
thing being reproduced is *that run*; a rebuilt argv would silently track later driver changes instead
of failing. And it carries the `products compared: N (expected N)` counter — job `49516776`'s lesson.

### 2.2 The trap that would have made this gate lie

**`--nd 5` on the transduction command line is not the cutoff the step used.** `b_automatic` is True
unless `--user_specific` (`x_TEA_main.py:313`), so the step overwrites it with
`Parameters.get_par_by_cov(coverage)[1]` — **4** at HG002's 31.816×. Hand-typing 5 produces a
different `.slct` and reads as "the lift broke something". The gate derives it from the run's own
cached `tmp_basic/basic_cov_is_rlth_info.txt`.

## 3. The number the per-locus loop needed: 41,604

Union over the three families of (cns candidates ∪ new raw-disc sites), exact position:

| | sites |
|---|---|
| L1 (`.slct`) | 34,197 |
| SVA (`.slct`) | 39,452 |
| Alu (cns list — no re-selection) | 2,204 |
| **union, exact position** | **41,604** — 7,910 claimed by one family, 33,139 by two, 555 by all three |

So a loop over the union visits **41,604** loci where the three per-family chains between them visit
**75,853** — **1.82×**. This is the iteration domain plan §3's "one pass over loci" asks for, and it
did not exist before. **It is not the 53,945-of-54K figure** in plan §4 2b: that one is about the td
gather's fetch *windows*, which are per-site read-collection intervals, not candidate loci.

## 4. The correction: §3's diagram order cannot work

Plan §3 draws *locus table → family assignment → per-locus loop*. **Assignment cannot go there.** For
52 of 1,515 records the only per-family evidence is the transduction evaluation's own output, so at
the point the diagram assigns, there is nothing to assign on — which is the same fact, reached a third
time, that makes the barrier's ceiling 7 of 22.

The order the data supports:

```
enumerate loci (site_universe: 41,604)  →  gather per-family evidence  →  assign family  →  filter
```

Only the assignment box moves; the "one pass" shape survives, because the locus set genuinely *is*
knowable up front. **Designing the loop off the picture as drawn rebuilds the barrier that Stage 2c
already measured as unable to close.** Recorded as a correction block at the diagram itself, not only
here.

## 5. What NOT to redo

Plan §10 is the authority and gained four entries. The ones most likely to cost a day:

- **The 52 are created, not moved, and no wider matching window recovers them** (§1). Distance 0 in
  `cns2`, ≥3.3 kb from any cns row.
- **They are all already in raw_disc.** Do not go looking for them in `new_sites` (that file has 4 and
  5 rows) or reason about them as orphan transductions.
- **Site enumeration being solved does NOT raise the 7-of-22 ceiling**, and `--family-assign-early`
  still stays OFF. `.slct` is a 34-40K-row superset; knowing a site exists is not knowing which family
  claims it.
- **The transduction and sibling steps ignore `--nd`** (§2.2). Any offline reproduction of a step that
  reads `n_disc_cutoff` must derive it from the basic-info file.
- **`XIntemediateSites.are_sites_close` implements a HALF-OPEN window**, because it walks
  `range(pos1-iwin, pos1+iwin)`. In this rule that makes the exclusion window on the raw-disc site's
  offset from a candidate **`(-100, +100]`**: exactly 100 bp to the *left* of a candidate is kept,
  exactly 100 bp to the *right* is suppressed. I got the direction backwards on the first write and
  the differential test caught it. Symmetrising it changes which sites transduction re-selects.
- **`is_matched_rep_type` has no Alu arm** and never did, despite the comment above it that said "only
  consider SVA and Alu transduction" (now corrected in the tree). Alu has no transduction flank.
- Everything in the previous handoff's §6 still stands, in particular that holding is permissive /
  emission is exact in the td gather, that sharing the gather buys CPU not wall clock, that
  `scan_clip_realign.tsv` carries no per-family signal, and that L1 precision must always be quoted
  with its stratum (74.5% overall, 92.7% excluding orphan sites).

## 6. NEXT

1. **The per-locus loop**, still item 1, now with one dependency discharged and one premise corrected.
   Its remaining unknowns are the ones §4 implies: what a per-family evidence record looks like at a
   locus that only `.slct` knows about, and where in the loop's tail assignment sits. Note the loop
   has **two callers waiting**: `--td-manifest` (whose "throwaway caller work" objection is now half
   answered — what the driver still lacks is the step's `iextnd` from `automatic_gnrt_parameters`) and
   the collapse of the two decision points.
2. **`--emit-per-type`**, still deliberately last: flipping per-type VCFs to opt-in breaks every gate
   script that `cmp`s them.
3. **Stage 3**, whose gate is fully specified (`HANDOFF_2026-08-02_stage2c_barrier.md` §1):
   position-level per family plus the orphan / sibling / no-`TD_SRC` strata, no transduction pass/fail.

## 7. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§3 correction block, §4 Stage 2c's last block, §8, §10 updated) |
| the lifted rule | `xtea/site_universe.py` — `write_slct_file`, `select_new_sites`, `is_matched_rep_type`, `read_positions` |
| its unit gate | `tests/test_site_universe.py` (11 tests; the differential one is the load-bearing one) |
| its WGS gate | `docs/handoff/runs/gate_site_universe.py` — offline, takes a finished run's sample dir |
| the site-creation probe | `docs/handoff/runs/probe_site_creation.py` (docstring amended: its ±50 window does no work on this hop) |
| the run both gates read | `/n/scratch/users/c/cos689/xtea2_testdata/stage2c_family_assign_49502115/OFF/HG002.wgs` |
| the multi-entry gather | `src/csrc/td_gather.hpp`; CLI `xtea_scan --td-manifest` |
| the family rule, one copy | `xtea/locus_table.py` |
| the barrier + waves | `xtea/run_xtea.py` — `downstream_waves`, `run_family_assign_barrier` |
| BASELINE-0 | `/n/scratch/users/c/cos689/xtea2_testdata/baseline0_hg002_L1-Alu-SVA/` |
| truth sets | `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/truth/` |
