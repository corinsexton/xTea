# xTea2 handoff — 2026-08-04, latest (reference SVA sources are reachable; and the decoy set turns out to be reconstructible)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it (§5's decoy blockquote,
> §5's "alignment passes" claim, §10's `decoy_cross_family` entry); this is the index.
> Supersedes `HANDOFF_2026-08-04_cns_cluster_ledger.md`, whose §7 NEXT item **1** — the SVA decoy fix —
> is what this session did. Its item 1's two open questions are answered separately (§3: one is fixed,
> the other is scoped to the line of code it needs); its items 2, 5 and 6 stand; **its items 3 and 4 are
> CLOSED by measurement** (§4, §4.1) — `MIN_DISC_MAPPED_RATIO` moves 10 loci genome-wide and
> `NEARBY_CLIP` moves 1, both measured offline from artifacts already on disk, which that handoff said
> needed instrumented runs.

**Three things happened. Only the first was on the agenda, and the third corrects the agenda.**

1. **The SVA decoy fix exists, behind `--td-decoy-structural`, off by default.** The classification is
   one shared function; the default path is the pre-fix expression verbatim; and the flip is gated at
   both scales. See §1, §2, §5.
2. **The agenda's second half of that item — "whether the SVA library should also get the decoy
   consensus concatenated in" — is answered offline and needs no library repackaging.** The 45-record
   decoy set is reconstructible **byte-identically from `rep_lib/consensus/` on every build**, and
   `consensus/transduction_decoy.fa` — a 3-record file that is exactly `>LINE1`/`>SVA`/`>ALU` — has
   **zero references in xTea2 or in master**. Not built: it is call-moving and needs its own arm. §3.
3. **NEXT items 3 and 4 are both CLOSED, and the reason they went unmeasured is a file everyone
   thought was deleted.** The disc-to-consensus SAM is **not deleted at the end of the step** — it is
   in every gate's output tree — so both took seconds, offline. `MIN_DISC_MAPPED_RATIO = 0.7` does not
   produce the 146 `no_cns_disc_reads` loci: ceiling **34**, moves **10**. And `NEARBY_CLIP`, swept
   **0 → ∞ through the shipping gate**, moves **one locus in a whole genome** — which also answers
   §5's request for a measured clip-to-insertion distance distribution: it is bimodal and 50 sits in
   the gap. §4 and §4.1.

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| the commits that touch code | `55055ba` — `transduction_evidence.is_flank_record`/`decoy_tags`/`DECOY_TAGS`, the two call sites in `x_transduction.py` (`:1160` primary hit, `:1285` XA tag), `XTransduction.b_decoy_structural`, `global_values.ENV_TD_DECOY_STRUCTURAL`, `run_xtea --td-decoy-structural`. `11b0e61` is the two gate scripts. **No HEAD hash recorded on purpose** — an earlier handoff named one and an `--amend` invalidated it the same day. |
| working tree | clean |
| tests | **43 in `tests/test_transduction_evidence.py`** (was 32 — eleven new, four of them through the shipping seam); whole fast suite green, 247 collected, 213 passed / 31 skipped with the C++ gates skipped. |
| in flight | **nothing — all four gate jobs have reported PASS** (§5). |
| the agenda | plan §7 / the previous handoff's §7, with item 1 discharged and items 3 and 4 **closed by measurement** (§4, §4.1). **1.** decide the default (§7 item 1). **2.** plan §5's semantic half on gate **[A]** — 97% of the rejecting; expect *precision*. **3.** the decoy-consensus concatenation (§3), which is what makes §5's SVA/Alu hazard testable at all. **4.** `ileftmost_cns` 500-vs-5000, the only one of §5's four "interrogate" thresholds left — **`NEARBY_CLIP` and `MIN_DISC_MAPPED_RATIO` are both measured and CLOSED** (§4, §4.1), moving 1 and 10 loci genome-wide. |

## 1. What was built, and the one thing it deliberately did not do

`parse_transduction_disc_algnmt` decides "this disc read hit a decoy consensus, not a source" with a
substring match on the bwa reference name. Every SVA flank record is named
`chr13~51077120~51078923~SVA_D~0L`, so `"SVA" in hit_flank_id` is True for **100% of hg38's 5,454 SVA
flank records** (hg19 99.1%, chm13 100%), the branch counts the hit into the cross-family tally and
`continue`s **before the source is extracted**. Measured consequence at WGS (job `49634290`): of SVA's
105 picked sources, **0 are reference copies and 105 polymorphic**, against L1's 114 of 159.

**And master's own output now says the same thing, which turns "master shares this bug" from a reading
of its source into a measurement.** The probe's run half over master's NA12878 WGS run
(`master_golden_wgs_seeded/NA12878`, job `48825977`): **SVA 92 picked sources, 0 reference copies**,
against **L1 135 with 84 reference copies**. Different sample, different pipeline, same signature —
which is also why the Stage 3 transduction gate cannot use master's calls as transduction truth.

- **The classification is now ONE function, `transduction_evidence.decoy_tags(name, b_structural)`,**
  asked by both the primary-hit test and `check_algnmt_tag_decoy`'s XA-tag copy of it. Two copies of a
  predicate is how the two forms diverged in the first place — the file also contains a *third*, exact
  equality, in the dead `parse_clip_algnmt_rescued` (`:1051`). That one is left alone on purpose: it is
  one of §8's zero-caller functions, and editing dead code to match live code makes the §8 sweep
  harder to reason about, not easier. **If it is ever revived, it must be revived through
  `decoy_tags`.**
- **`b_structural=False` is the pre-fix expression verbatim**, including the HERV asymmetry: there are
  four tags and three counters, so a HERV-only hit enters the skip branch, increments nothing and is
  dropped. `test_the_default_decoy_test_is_the_substring_match_verbatim` pins it against a copy of the
  pre-fix code. **Anything that ships off has to be provably inert, or its default is not a rollback.**
- **`b_structural=True` asks `is_flank_record` first.** A flank id has at least three `~`-delimited
  fields with integer coordinates in fields 2 and 3 — `5~77881127~77881127~polymerphic~0L` included,
  since `gnrt_flank_regions_of_polymerphic_insertions` writes `start == end` — and a consensus id
  (`LINE1`, `ALU_AluYb8`, `SVA_D`) has none of that. It **fails closed**: an id with enough tildes but
  unparseable coordinates is still classified by the tags.
- **The switch crosses a process boundary as an env var, because the step is its own `x_TEA_main`
  subprocess** — `XTEA_TD_DECOY_STRUCTURAL`, resolved once in `XTransduction.__init__` so one run
  cannot answer it two ways, and set by `run_xtea --td-decoy-structural`. Same mechanism as
  `XTEA_SCAN_TD_GATHER` and `XTEA_CNS_PRECOMPUTED_ALIGN`. The driver prints a one-line banner when and
  only when it is on, which is what both gates' check 0 reads.

### 1.1 Why it ships OFF, when it is a fix

Because the default is the **frozen** behaviour, not the correct one. BASELINE-0 and every gate in
`docs/handoff/` is frozen against the substring form; flipping the default silently would invalidate
all of them at once and make the flip's own measurement unattributable. It is the `--td-gather-scan`
shape: ship off, gate both arms, flip on the strength of the numbers. **The flip is a decision this
handoff does not make** — §5 has the evidence for it and §6 says what is still missing.

## 2. Eleven tests, four of them through the seam

The differential-on-the-library habit has a known blind spot in this repo — the return-arity bug (job
`49633492`) that the whole fast suite missed — so four of the eleven drive
`parse_transduction_disc_algnmt` itself on a small SAM whose one reference is a **real** hg38 SVA flank
record:

| test | what it pins |
|---|---|
| `…reaches_an_sva_reference_source_only_with_the_structural_test` | five disc reads of chr1:1000 at MAPQ 60 → **no source at all** with the default test, and exactly `chr13:51077120-51078923~0` with the fix on. The source is the flank record's own interval, i.e. a REFERENCE copy |
| `…a_real_consensus_hit_is_still_dropped_…_with_the_fix_on` | turning the switch on must not turn the decoy skip off (`SVA_D`, a record hg19's library really ships) |
| `…the_l1_flank_path_is_byte_identical_either_way` | the claim that makes an L1-identical WGS arm meaningful |
| `…the_xa_tag_decoy_test_uses_the_same_classification` | one run cannot hold two definitions of "decoy" |
| `…the_switch_is_off_unless_the_env_var_says_otherwise` | `""` and `"0"` are both off, so an inherited variable cannot look like an opt-in |

The remaining six are the classification itself, on **real library record names** rather than invented
ones, including that a malformed id is not promoted to a flank record and that `is_flank_record` is
strict about what "structure" means.

`probe_td_source_classes.py` now **imports** both shipping predicates instead of writing its own copy,
and reproduces every documented number against them (hg19 5,114 of 5,158 collisions, hg38 5,454,
chm13 6,072; 44 / 45 / 0 genuine decoy records). A probe with its own copy of the rule can agree with
itself while the step disagrees.

## 3. THE FINDING THAT WAS NOT ON THE AGENDA: the decoy set is reconstructible, and its packaged form is unreferenced

The previous handoff's item 1 had a second half — *"whether the SVA library should also get the decoy
consensus concatenated in, which is what would make §5's SVA/Alu hazard testable at all"* — recorded as
an open question. It is not open. Measured offline from the libraries and the consensus directory, in
seconds:

**The decoy set is 45 records, and every one of them is byte-identical to something in
`rep_lib/consensus/`:**

| the 45 | source | sequences match |
|---|---|---|
| `>LINE1`, `>SVA`, `>ALU` | `consensus/transduction_decoy.fa` (which has exactly these 3) | 3/3 |
| 35 × `>ALU_Alu*` | `consensus/ALU.fa`, under an `ALU_` prefix | 35/35 |
| 6 × `>SVA_[A-F]` | `consensus/SVA.fa` | 6/6 |
| `>HERV-K` | `consensus/HERV.fa` (which has exactly this 1) | 1/1 |

**Who carries them:** both L1 libraries (hg19 and hg38) carry all 45 and their name sets are
*identical*; hg19's SVA library carries 44 — the same set minus `HERV-K`; hg38's SVA, chm13's SVA and
**chm13's L1** carry none.

Three consequences, and the third is the one to act on:

1. **`TD_DECOY_HERV` is not speculative.** `>HERV-K` is the 45th record, which is why the tag exists
   with no counter beside it. That closes a loop the previous handoff left as an oddity.
2. **`consensus/transduction_decoy.fa` has ZERO references** — not in xTea2, not in master
   (`git grep transduction_decoy 46520ab` is empty). A file packaged to be exactly this step's decoy
   set, that nothing reads. And `construct_novel_flanks`'s own comment (`x_transduction.py:207`) says
   *"we add L1, Alu and SVA consensus to the flank regions to filter out those mapped to consensus"*
   while the function concatenates the family library with the polymorphic flanks and nothing else.
   **Plan §5's claim that "the flank library already has consensus concatenated in" is that comment,
   not behaviour** — corrected in the plan.
3. **So making §5's SVA/Alu hazard testable needs no library repackaging and no build-specific data.**
   It is a concatenation in `construct_novel_flanks`, sourcing the 45 from `rep_lib/consensus/`, which
   works on every build including chm13 (where even L1's tally is dead today). **It is CALL-MOVING and
   is NOT built here** — adding records to a bwa reference changes which hits win, so it moves calls in
   both families and needs its own arm, with `probe_td_source_classes.py` on both sides of it. It is
   also independent of the fix in §1: that one makes flank records reachable, this one makes the decoy
   arm live.

## 4. AND THE NEXT AGENDA ITEM IS SMALLER THAN ITS OWN ENTRY SAYS — measured offline, in seconds

The previous handoff's NEXT item **3** reads: *"`MIN_DISC_MAPPED_RATIO = 0.7` — the knob the split just
exposed … it decides the **146** `no_cns_disc_reads` loci — four times everything [B2]'s own thresholds
decide … Interrogate it before `NEARBY_CLIP`."* **That is the same mistake as the claim it corrected**:
a population attributed to the threshold nearest to it.

**THE ENABLING CORRECTION FIRST, because it is worth more than the finding.** Items 3 and 4 were both
filed as needing an instrumented run, on the strength of "the disc-to-consensus SAM … is deleted at the
end of the step". **It is not.** `run_xtea.cleanup_successful_run` deletes it at the end of the **run**,
and only when `--keep-intermediates` was not passed — and **every gate under `docs/handoff/runs/` passes
it**. So `<fam>/tmp/transduction/temp_transduction_disc_cns.sam` is sitting in every one of their output
trees at **1.4 MB per family at WGS**, next to the 393 MB disc-to-flank SAM and the 2.46 GB
clip-to-flank SAM. Nothing here needed a new run; job `49703322`'s tree already had it.

`docs/handoff/runs/probe_disc_cns_mapped_ratio.py` reads it, **imports** the shipping predicate
(`is_clipped_part_qualified_algnmt` on a real `XClipDiscFilter`) rather than copying the cutoff, and
**reconciles against the ledger** — 113 (L1) and 33 (SVA) absences by two derivations from different
files, a SAM and a TSV. It exits nonzero on a mismatch, so a wrong model of the filter fails loudly
instead of publishing a plausible table. On job `49703322`:

| | L1 | SVA | total |
|---|---|---|---|
| `no_cns_disc_reads` loci | 113 | 33 | **146** |
| …**every** record came back unmapped from bwa — **unreachable at any cutoff** | 89 | 23 | **112** |
| …≥1 record mapped and failed the ratio — **the cutoff's whole ceiling** | 24 | 10 | **34** |
| loci recovered by 0.7 → 0.5 | +4 | +6 | **+10** |
| loci lost by 0.7 → 0.9 | −4 | −8 | −12 |

**So it reaches 34, against `NEARBY_CLIP`'s 118, and moves 10 genome-wide.** It is the *smaller* lever
of the two, and the ranking in the previous handoff's NEXT is inverted. The largest fact in this region
is not a threshold at all: **90.2% (L1) and 84.9% (SVA) of disc-to-consensus alignment records come back
unmapped.**

**And a perf observation that is NOT a sensitivity one.** Of L1's 4,969 distinct focal reads, **2,784
were already unmapped against the flank library** — `re_select_disc_reads` applies no mapped filter — and
they contribute **3 of the 456 records that survive to gate [B2]**. So the third bwa pass and its parse
handle ~2.2× the reads that can matter. Filtering them would move 3 records, so it is **call-moving at
the margin, not a free win**, and it belongs with the other call-moving items rather than in a cleanup.

### 4.1 …and so is the item I promoted above it. `NEARBY_CLIP` is INERT, measured the same way.

Having said "`NEARBY_CLIP` is the larger lever and it is now offline", the honest next step was to
measure it rather than re-rank it. `docs/handoff/runs/probe_nearby_clip.py` sweeps the constant
**0 → ∞ through the shipping `check_seprt_disc_consistency`** over job `49703322`'s own inputs:

| NEARBY_CLIP | L1 survivors | SVA survivors |
|---|---|---|
| 0, 10, 25, **50 (shipping)**, 100, 200, 500 | **24** at every value | **56** at every value |
| ∞ | 24 | 57 |

**One locus in a whole genome, and only past the largest observed distance.** At 50 the probe
reproduces the run's own ledger B2 split exactly — L1 `no_cns_disc_reads` 113 + `no_cns_disc_cluster`
22, SVA 33 + 16 — which is what makes this a measurement and not a model of one.

**Plan §5 asked for "a measured clip-to-insertion distance distribution" to replace the arbitrary
proxy. Here it is, and it says the proxy is fine.** Over the disc records clipped on exactly one side
— 64 of L1's 296 in the domain, 178 of SVA's 505 — `|anchor_pos - ins_pos|` is sharply **bimodal**: L1
min 0, median **3**, p75 **411**, max 762; SVA median 129, p75 358, p90 787, max 998. **50 sits in the
empty gap between the modes**, which is why moving it anywhere in 0-500 changes nothing.

Two facts make this computable from read names alone, and both are worth knowing independently:
`prep_representative_pos` is the **identity** (`x_transduction.py:710` — `clip_pos` IS the insertion
position, despite its name), and `parse_disc_algnmt_consensus` stores the read name's anchor field
verbatim.

**So three of plan §5's four "interrogate" thresholds are now measured near-inert** — the two-sided
polyA/polyT veto (1 and 2 loci), the large-indel check (1 each), and this. The one left is
`ileftmost_cns` 500-vs-5000, and it is not a small lever: `family_params`' docstring measured that
waking SVA's 5000 admits essentially the whole candidate list.

## 5. The gates

| gate | job | result |
|---|---|---|
| HG002 chr1, three arms (BEFORE / AFTER_OFF / AFTER_ON) | `49721312` | **PASS** — 100 identical + 2 characterised with the flag off; with it on, **one** file in the tree moves and it is SVA's ledger. §5.1 |
| NA12878 golden downstream, default path | `49721733` | **PASS** — 3/3 families, exactly the documented notes: 4 L1 carry-forward, Alu genotype flips |
| HG002 WGS, one arm with the flag ON vs frozen BASELINE-0 | `49724055` (20 cores, `short`) and `49724056` (24 cores, `priopark`) | **PASS on both, byte-identical to each other.** L1+Alu identical to BASELINE-0 with the fix ON, SVA +3 calls, **241 reference-copy SVA sources from 0** — every number as predicted in §5.3. §5.2 |

### 5.1 `gate_decoy_structural_chr1.sbatch` — three arms, because "off is inert" and "on works" are different claims

- **AFTER_OFF must be byte-identical to BEFORE over the whole output tree.** The unit differential
  pins the predicate; only a run pins the seam. Both decision ledgers are compared **outright** here,
  unlike the gate [B2] ledger gates — nothing in this change is *supposed* to alter a ledger.
- **AFTER_ON is the wiring proof, not the measurement.** chr1's gate B has 5 L1 survivors and **0
  SVA**, so its SVA delta may be zero: the fix's population is at WGS. A zero delta at chr1 is not
  evidence the fix does nothing, and the script says so where a reader will hit it.
- **Check 0 reads the driver banner in all three arms**, because with two arms sharing one worktree the
  thing most likely to go wrong is the *flag*, not the commit.
- Check 3 is the same prediction the WGS gate makes: L1 and Alu must not move with the fix on.

**RESULT at `11b0e61`, job `49721312`: PASS.** Walls 7:51 / 5:59 / 5:57.

| check | result |
|---|---|
| 0. the arms are the code AND the configuration | `decoy_tags` BEFORE=0 AFTER=1; the `[decoy]` banner 0 / 0 / **1** |
| 1. AFTER_OFF vs BEFORE, whole tree | **100 identical, 2 characterised, 0 DIFFERING, 0 missing** |
| 2. new or removed files, either after arm | **0** |
| 3. what moved with the fix ON | **L1 0, Alu 0, SVA 1, other 0** — and the one SVA file is its **decision ledger** (the other two diffs are the characterised `cns_manifest.txt` and `scan_disc_pos.tsv`) |
| 4. picked sources per arm | L1 5 / 5 reference copies in all three; SVA **0** in all three |
| 5. the ledger populations | L1 1,695 rows identical in all three arms. **SVA 2,281 → 2,315 with the fix on** |
| 6. call counts | **1 / 14 / 2 / 17 in all three arms** |

**AND CHECK 5 IS THE ONE THAT SHOWS THE FIX ARRIVED.** Those 34 new SVA ledger rows are new gate-**B**
decisions — `low_disc_support` **2 → 35** plus one `no_dominant_source` — on loci whose disc hits used
to be discarded at the decoy branch *before any gate saw them*, so they had no row at all. At chr1 the
fix reaches **34 loci and every one is rejected for thin evidence**, which is exactly what a gate-B
domain of 5 L1 loci and 0 SVA loci predicts. **Not one call moves.** That is the wiring proven and the
population deferred, which is what this leg was for.

### 5.2 `gate_decoy_structural_wgs.sbatch` — one arm, and the frozen baseline is the other one

BASELINE-0 was reproduced byte-for-byte by the last WGS gate at this commit's behavioural parent (job
`49703322`, 4/4 VCFs), so it **is** the off-arm at no compute. That makes check 1 two-sided rather than
a diff:

> **L1 and Alu must be byte-identical to BASELINE-0 WITH THE FIX ON, and SVA must not be.**

Both directions are failures. A moved L1 or Alu file falsifies the classification's premise (L1
collides with no tag — `TD_DECOY_LINE` is `"LINE1"`, not `"L1"`; Alu never runs this step). An unmoved
SVA file means the flag reached nothing. **L1's 21,158 ledger rows are gated too**, which is finer than
its VCF: a locus can move between gates without changing whether it is called.

**Check 5 is the payoff and it is the one number the fix exists for** — the probe's SVA reference-copy
picked sources, 0 of 105 today on every build and on master, against a 5,454-record library built to be
exactly that source set.

**Checks 6 and 7 report rather than gate, and §7.1 says why.** HG002's derived transduction denominator
is **0 (L1) and 2 (SVA)**, so the truth set cannot score a transduction call as a transduction. What it
can say is whether SVA *position* sensitivity/precision moved from 85.6% / 93.9% over 90 truth records
— and with only 5 SVA false positives in the whole genome, a handful of new calls moves precision
several points either way. Read it as a direction, not a verdict.

**RESULT at `11b0e61`: PASS on both arms** — `49724055` (20 cores, wall 47:54) and `49724056` (24
cores, 41:08) — **and the two arms are byte-identical to each other**, which retires §5.4's
core-count caveat by measurement rather than argument.

| check | result |
|---|---|
| 0. the flag reached the driver | banner 1 |
| 1. **THE PREDICTION** | **L1 and Alu IDENTICAL to BASELINE-0 with the fix ON; SVA and merged differ.** Three positions gained, none lost: chr12:39235233, chr20:45799810, chrY:10120435 |
| 2. call counts | **L1 188, Alu 1223** (gated, exact); SVA **82 → 85**; merged 1493 → 1496 |
| 3. the ledgers | **L1 21,158 rows, and its by-reason line is identical to job `49703322`'s to the last count**, `decoy_cross_family=144` included. SVA **27,787 → 29,001** |
| 4. `downstream_parallel` | 550.3s at 20 cores, **494.9s at 24** — inside the 462-486s band the no-flag runs sat in, so the extra evidence costs ~10-25s |
| 5. **THE PAYOFF** | **SVA reference-copy picked sources 241 of 346, from 0 of 105.** PASS |
| 6. the truth-set score | SVA sensitivity **unchanged at 85.6%**; precision **93.9% → 90.6%** (FP 5 → 8). All families: sensitivity unchanged 80.4%, precision 88.4% → 88.2% |
| 7. the FN attribution | 322 FNs, unchanged — but **one moves from `no_ledger_row` to `low_disc_support`** (306 → 305), a locus that previously had no row because its hits were discarded before any gate saw them |
| 8. the funnel reconciliation | **OK for both families** |

**AND THE OFFLINE PREDICTION HIT EVERY NUMBER.** SVA picked **346** (predicted 346), reference copies
**241** (241), gate [B] drops **1,086** (1,086), gate [B]'s SVA domain **1,432** (1,432), orphan loci
**87** (87), and L1 identical to the last ledger count. Two derivations from different sources — a
frozen SAM replayed offline, and a 47-minute run — agreeing exactly.

### 5.2.1 WHAT THE THREE NEW CALLS ARE, and why "precision 93.9% → 90.6%" overstates the cost

All three carry `SUBTYPE=one_side_and_half_transduction` and a `TD_SRC` naming a **reference SVA
copy** — `chr2:197618052-197619758~0`, `chr8:100891476-100893344~0`, `chr20:30846336-30848045~1` —
which is precisely the population the fix unlocked. They are not incidental FPs; they are the fix
working.

**And two of the three are ORPHAN transduction sites** (`~orphan` in `TD_SRC`, `SVLEN=0`): chr20 and
chrY. Plan §4 measured what that means — an orphan site is where transduced sequence landed *without*
the element body, so it is **not an MEI insertion**, and this truth set enumerates MEI insertions.
All 37 of L1's orphan calls are off-truth for exactly that reason, and the harness reports L1
precision both ways because of it.

**The harness does not catch these, and that is a real gap.** `benchmark_truth_set.py:675` keys the
orphan stratum on `subtype == "orphan_or_sibling_transduction"` — the string L1's 37 calls carry —
and these two carry `one_side_and_half_transduction` with the orphan marker in `TD_SRC` instead. So
the reported 90.6% counts two category-error calls as ordinary false positives. **The measured
precision cost of the fix is one call, not three**, and fixing the stratum's key is a small,
non-call-moving harness change that should happen before this number is quoted anywhere.

### 5.3 A PREDICTION, WRITTEN BEFORE THE WGS ARMS LANDED — and it found a parameter this repo's tests describe wrongly

`--td-decoy-structural` changes exactly one thing: how `parse_transduction_disc_algnmt` classifies a
hit's reference name. It does not change the bwa reference (`construct_novel_flanks` depends on
neither the flag nor anything downstream of it), it does not change which reads reach that pass, and
it runs strictly **after** the alignment. **So a pre-fix run's `temp_transduction_disc.sam` is the
same file the fixed run will produce** — and per §4 it is still on disk. Replaying the shipping
method over job `49703322`'s SAM with the flag both ways is therefore the fix's exact effect on gate
[B], derived from a different source than the run. `docs/handoff/runs/probe_decoy_fix_delta.py`.

**IT ALSO CAUGHT A MIS-STATED PARAMETER, and the catch is what makes the rest trustworthy.** The
first replay reproduced **94** of the run's 105 picked SVA sources. The reason: `x_TEA_main` sets
`b_automatic = True` unless `--user`, which nothing passes, so it **overrides `--nfclip/--nfdisc`**
with `Parameters.get_par_by_cov(coverage)`. HG002 WGS is **31.816×** → `m_cov_par[30] = (3, 4, 1)`,
so the step runs at **n_clip_cutoff 3 and n_disc_cutoff 4** — `n_half_disc` **2.0**, not the 2.5 that
`--nfdisc 5` implies. `run_xtea.py`'s locus-universe docstring already says this ("wrong by one at
HG002 (4, not 5)"); `tests/test_transduction_evidence.py`'s header comment claiming `--nclip 4
--nd 5` gives "the production ones" does not, and is corrected. **No test result changes** — the
differentials are input-agnostic — but any threshold reasoning done at 2.5 is reasoning about a run
that never happened.

At the corrected arguments the OFF arm reproduces the run **exactly**, on both families and on the
full gate-[B] ledger split, which is what licenses reading the ON arm as a prediction:

| | OFF (= the run) | ON (predicted) |
|---|---|---|
| **SVA** picked sources | **105** (0 reference, 105 polymorphic) | **346** (**241** reference, 105 polymorphic) |
| **SVA** orphan loci | 26 | **87** |
| **SVA** gate [B] ledger | `low_disc_support` 92, `no_dominant_source` 13 = **105** | 898 + 188 = **1,086** |
| **SVA** gate [B] domain (rows + picked) | **210** | **1,432** |
| **L1** everything | 159 picked (114 ref / 45 poly), 30 orphan, ledger 144 + 424 + 108 = **676** | **identical, every number** |

Both OFF columns match what is already recorded: 105 SVA picked with 0 reference copies, 159 L1 with
114, and plan §10's L1 `B_source` split of 676 = 424 + 144 + 108. So:

> **PREDICTED, before the runs reported: SVA reference-copy picked sources 0 → 241, total picked
> 105 → 346, gate [B]'s SVA domain 210 → 1,432 loci, and L1 unchanged to the last row.**

**And note what the L1 row is:** the "L1 collides with no decoy tag" claim, which the whole two-sided
WGS check rests on, measured over 393 MB of real alignments rather than argued from `"LINE1" != "L1"`.
It includes `decoy_cross_family` staying at exactly 144 — the fix does not cost L1 its tally.

**What this CANNOT say is what the run is for.** Everything after gate [B] — [B2]'s cluster
consistency, [C]'s conjunction, [D], [E] — is untouched by the replay, so 346 picked sources is an
*input* to the cascade, not a call count. 3.3× the SVA candidates is also why check 6's precision
number is the one to read carefully.

### 5.4 The WGS leg's scheduling, and the one variable it introduces

`priopark` was saturated when this ran — 18 idle CPUs of 1,152, with slurm's estimated start **10 hours
out** — so the leg went out **twice**: `49724056` at 24 cores on `priopark` as submitted, and `49724055`
at **20 cores on `short`**, which is that partition's per-job cap (`medium` has the same cap, and
neither accepts a 24-core job at all).

**THE CAVEAT BELOW IS NOW RETIRED BY MEASUREMENT — both arms ran, and their SVA VCF is
byte-identical.** Keep the reasoning: it is why the second arm was submitted, and it is the argument
to repeat the next time a gate has to be scheduled at a different width.

**READ THE 20-CORE ARM WITH ONE CAVEAT AND IT IS NOT A LOOSE ONE.** Every byte-identity claim in this
repo's gates was measured between arms at *equal* core counts, and `-n` reaches bwa as `-t`. For these
realignments the input is single-end FASTA, so bwa's per-batch insert-size estimation — the usual reason
`-t` changes paired-end output — does not apply, and the native gather is separately proven
deterministic. But that is an argument, not a measurement. So:

- **check 5 (the payoff) and check 3 (the ledgers) are core-count-independent** and can be read off
  either arm;
- **check 1's L1/Alu byte-identity is the claim at risk.** If the 20-core arm reports an L1 or Alu
  difference, the 24-core arm is what says whether the fix or the thread count caused it. If the
  24-core arm has not landed, do not resolve it by argument.
- **check 4's `downstream_parallel` bound is not comparable across the two** — 485.8s was measured at 24
  cores, and the 900s bound in this script is loose enough to cover 20.

**MEASURED: the two arms' SVA VCFs are BYTE-IDENTICAL** (and every ledger count and call count agrees),
so at this width the thread count is not a variable for this pipeline. `downstream_parallel` is the one
thing that did move: **550.3s at 20 cores against 494.9s at 24**. Treat the identity as one data point,
not a licence — it says nothing about a *paired-end* realignment path.

## 6. What NOT to redo

Plan §10 is the authority. Everything in the previous handoff's §6 still stands except where noted.

- **The decoy classification is `transduction_evidence.decoy_tags`, and there must not be a second
  copy.** There were two live ones and there is still a third in dead code (`:1051`, exact equality).
- **`b_structural=False` is a frozen artifact, not a preference.** Do not "clean up" the substring
  branch: it is what makes every existing gate and BASELINE-0 comparable, and a test pins it against a
  verbatim copy of the pre-fix expression.
- **Do not flip the default on the strength of the fix being correct.** It is correct, and it is
  measured contained; that is not the question. **The question is what the 241 newly reachable
  reference sources are WORTH, and HG002 cannot say** — its derived transduction denominator is 2-3
  (plan §4 Stage 3), so the truth set sees the fix only through 3 extra SVA calls, two of which it
  cannot match by construction (§5.2.1). The instrument that could score the 241 is the Stage 3
  transduction gate, which is **spec'd and unbuilt** — and which needs this flag ON to build its
  candidate pool. So: build `benchmark/verify_transduction.py`, score the 241, then flip. See §7.
- **Do not quote "SVA precision 93.9% → 90.6%" without §5.2.1.** Two of the three new FPs are orphan
  transduction SITES, which are not MEI insertions and which this truth set cannot contain; the harness
  misses them because its orphan stratum is keyed on a `SUBTYPE` string they do not carry. The measured
  cost is **one call**.
- **L1 is untouched BY CONSTRUCTION, and that is a load-bearing claim, not a convenience.** It is what
  lets an L1-identical WGS arm be evidence about the fix rather than a coincidence. If an L1 file ever
  moves under this flag, the premise is wrong, not the gate.
- **The hg38 SVA library's zero decoy records is still a live dead-arm**, so §5's SVA/Alu-competition
  hazard is **still untested** even with the fix on. The fix makes flank records reachable; it does not
  give the cross-family tally anything to count. Do not read a passing WGS gate as testing that hazard.
- **`consensus/transduction_decoy.fa` is unreferenced and the "consensus is concatenated in" comment is
  false** (§3). Anything sized off that comment is sized off nothing.
- **HERV has a tag and no counter on purpose** — `>HERV-K` is a real record in both L1 libraries, and
  giving it a counter is a call-moving change with no evidence behind it.
- **Everything the previous handoff's §6 lists** — six gate groups behind three batched passes, [B2]'s
  body staying in `x_clip_disc_filter.py`, never sorting ledger phase names, `no_cns_disc_reads` being
  an absence, [B2] not being the largest threshold rejection, `NEARBY_CLIP`'s small reach, master's
  calls not being transduction truth, the two preserved defects, the two coverage windows, the ledger
  reasons as public contract, the `.slct` two-independent-screens rule — unchanged.

## 7. NEXT

1. **The default: RECOMMENDED SEQUENCE IS "build the transduction gate, then flip" — and this is
   Corinne's call, not mine.** All the evidence is in and it does not point at flipping *now*:

   | | |
   |---|---|
   | contained? | **yes, measured.** L1 and Alu byte-identical to BASELINE-0 at WGS *with the fix on*, L1's 21,158 ledger rows identical to the count, 100/100 files at chr1 with it off |
   | does what it claims? | **yes.** 241 reference-copy SVA sources from 0, exactly as predicted offline |
   | costs what? | **one off-truth SVA call.** Sensitivity unchanged everywhere; the reported 93.9% → 90.6% counts two orphan SITES the truth set cannot contain (§5.2.1) |
   | worth what? | **unknown, and HG002 cannot say.** Its derived transduction denominator is 2-3, so the 241 are invisible to the only scorer that exists |

   The missing number is the last row, and the instrument for it is the **Stage 3 transduction gate**
   (`docs/handoff/HANDOFF_2026-08-04_stage3_td_gate_spec.md`, deliverable
   `benchmark/verify_transduction.py`) — which the plan already says **needs this flag ON to build its
   third pool arm**, because nothing else produces reference-copy SVA sources. So the flag is exactly
   what that gate was waiting for, and the gate is exactly what the default is waiting for. Build it,
   score the 241, flip on the result.

   **And when the flip happens, its cost is not the one-line change.** It is **re-freezing BASELINE-0**:
   every gate under `docs/handoff/runs/` that compares against it has to restate its expectations in
   the same commit, and this handoff's own numbers (188 / 1223 / 82 / 1493, L1 21,158 rows) become
   stale. Do that in one commit or not at all.

   **Before quoting any precision number, fix the harness's orphan stratum** (§5.2.1): a one-line key
   change in `benchmark_truth_set.py:675`, not call-moving, and it is what makes the cost row above
   readable by someone who did not run this.
2. **Plan §5's semantic half, on gate [A]** (97% of the rejecting). Expect *precision*: 95% of
   BASELINE-0's misses never reach this decision. Unchanged from the previous handoff's item 2, and its
   `decoy_cross_family` caveat stays answered.
3. **The decoy-consensus concatenation (§3).** Its own arm; the recipe is exact and build-independent;
   it is what makes §5's SVA/Alu hazard testable at all. Sequence it *after* the default decision in
   item 1, or the two changes' effects on SVA are not separable.
4. **`ileftmost_cns` 500-vs-5000 — the ONLY one of plan §5's four "interrogate" thresholds still
   worth an item** (§4.1). The other three are measured near-inert, and this one is the opposite of
   small: `family_params`' docstring measured that waking SVA's 5000 admits essentially the whole
   candidate list into the polymorphic flank library, and every realignment with it. Truth-set gated,
   like the decoy default.
5. **NOT `NEARBY_CLIP` and NOT `MIN_DISC_MAPPED_RATIO` — both are measured and closed** (§4, §4.1).
   `NEARBY_CLIP` moves **one** locus swept 0 → ∞; the mapped-ratio's ceiling is 34 and it moves 10 at
   0.5. The question left in that region is not a threshold: it is **why 85-90% of disc-to-consensus
   records come back unmapped**, and whether gate [B2] wanting element-body evidence from reads
   selected for carrying *flank* sequence is the design or an accident.
6. **The per-locus decision loop's remaining box**, and the staged-barrier handoff's items 2-5,
   unchanged.

## 8. Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` (§4 Stage 3, §5, §8, §10) |
| the decision | `xtea/transduction_evidence.py` — `is_flank_record`, `decoy_tags`, `DECOY_TAGS`, `TransductionParams`, `SiteEvidence`, `gate_*`, `PHASE_ORDER`/`phase_rank`, `DecisionLedger` |
| the two call sites | `xtea/x_transduction.py:1160` (primary hit) and `:1285` (`check_algnmt_tag_decoy`, XA tags); the switch is `XTransduction.b_decoy_structural` |
| the switch | `run_xtea --td-decoy-structural` → `global_values.ENV_TD_DECOY_STRUCTURAL` (`XTEA_TD_DECOY_STRUCTURAL`) |
| its unit + seam tests | `tests/test_transduction_evidence.py` (43 tests; four drive the shipping method) |
| its end-to-end gates | `docs/handoff/runs/gate_decoy_structural_{chr1,wgs}.sbatch`, `gate_golden_downstream.sbatch` |
| the instrument | `docs/handoff/runs/probe_td_source_classes.py` — libraries + a finished run, seconds, no BAM; `--genome hg19,hg38,chm13` sweeps builds; now imports the shipping predicates |
| the decoy set's sources | `rep_lib_annotation/consensus/{transduction_decoy,ALU,SVA,HERV}.fa` (§3) |
| the mapped-ratio probe | `docs/handoff/runs/probe_disc_cns_mapped_ratio.py` — sweeps cutoffs, separates unmapped from ratio-failed, reconciles against the ledger, exits nonzero on a mismatch |
| the NEARBY_CLIP probe | `docs/handoff/runs/probe_nearby_clip.py` — sweeps the constant 0 → ∞ through the **shipping** `check_seprt_disc_consistency`, prints the clip-to-insertion distance distribution §5 asked for, and checks itself against the run's ledger B2 split |
| the decoy-fix delta probe | `docs/handoff/runs/probe_decoy_fix_delta.py` — replays the shipping disc parse over a pre-fix run's SAM with the flag both ways; the OFF arm reproducing the run is its own check (§5.3) |
| the SAM everyone thought was deleted | `<run>/<fam>/tmp/transduction/temp_transduction_disc_cns.sam` (1.4 MB/family at WGS), next to `temp_transduction_disc.sam` (393 MB) and `temp_transduction_clip.sam` (2.46 GB). Present whenever `--keep-intermediates` was passed, which every gate does |
| the FN attribution | `benchmark/ledger_attribution.py` |
| the truth-set scorer | `benchmark/benchmark_truth_set.py score --calls <vcf> --truth-dir <dir>` |
| the funnel reconciliation | `docs/handoff/runs/probe_decision_funnel.py` |
| worktrees for this session | `/n/scratch/users/c/cos689/xtea2_wt/td_decoy_55055ba` (AFTER, at `11b0e61`; `xtea_scan` **copied** from the BEFORE worktree, legitimately — `src/` is byte-identical across the range), BEFORE is `…/td_b2_365f352` |
| everything else | see `HANDOFF_2026-08-04_cns_cluster_ledger.md` §8 |
