# xTea2 handoff — 2026-08-04, latest (reference SVA sources are reachable; and the decoy set turns out to be reconstructible)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It still outranks this chain on agenda.
> Everything factual below is written into the plan at the section that owns it (§5's decoy blockquote,
> §5's "alignment passes" claim, §10's `decoy_cross_family` entry); this is the index.
> Supersedes `HANDOFF_2026-08-04_cns_cluster_ledger.md`, whose §7 NEXT item **1** — the SVA decoy fix —
> is what this session did. Its items 2-6 stand, with item 1's two open questions now answered
> separately (see §3 below: one is fixed, the other is scoped to the line of code it needs).

**Two things happened, and the second one was not on the agenda either.**

1. **The SVA decoy fix exists, behind `--td-decoy-structural`, off by default.** The classification is
   one shared function; the default path is the pre-fix expression verbatim; and the flip is gated at
   both scales. See §1, §2, §4.
2. **The agenda's second half of that item — "whether the SVA library should also get the decoy
   consensus concatenated in" — is answered offline and needs no library repackaging.** The 45-record
   decoy set is reconstructible **byte-identically from `rep_lib/consensus/` on every build**, and
   `consensus/transduction_decoy.fa` — a 3-record file that is exactly `>LINE1`/`>SVA`/`>ALU` — has
   **zero references in xTea2 or in master**. Not built: it is call-moving and needs its own arm. §3.

## 0. State at handoff

| | |
|---|---|
| branch | `xTea2` |
| the commits that touch code | `55055ba` — `transduction_evidence.is_flank_record`/`decoy_tags`/`DECOY_TAGS`, the two call sites in `x_transduction.py` (`:1160` primary hit, `:1285` XA tag), `XTransduction.b_decoy_structural`, `global_values.ENV_TD_DECOY_STRUCTURAL`, `run_xtea --td-decoy-structural`. `11b0e61` is the two gate scripts. **No HEAD hash recorded on purpose** — an earlier handoff named one and an `--amend` invalidated it the same day. |
| working tree | clean |
| tests | **43 in `tests/test_transduction_evidence.py`** (was 32 — eleven new, four of them through the shipping seam); whole fast suite green, 247 collected, 213 passed / 31 skipped with the C++ gates skipped. |
| in flight | see §4 for the three gates and their results. |
| the agenda | plan §7 / the previous handoff's §7, with item 1 discharged. **1.** plan §5's semantic half on gate **[A]** (97% of the rejecting; expect *precision*). **2.** `MIN_DISC_MAPPED_RATIO = 0.7`, the largest knob between gates [B] and [C] (146 loci) and still not in §5's list. **3.** the decoy-consensus concatenation from §3 below, which is what makes §5's SVA/Alu hazard testable at all. **4.** `NEARBY_CLIP` — still 4th, still a small lever. |

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
handoff does not make** — §4 has the evidence for it and §5 says what is still missing.

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

## 4. The gates

| gate | job | result |
|---|---|---|
| HG002 chr1, three arms (BEFORE / AFTER_OFF / AFTER_ON) | `49721312` | TBD |
| HG002 WGS, one arm with the flag ON vs frozen BASELINE-0 | `49721731` | TBD |
| NA12878 golden downstream, default path | `49721733` | TBD |

### 4.1 `gate_decoy_structural_chr1.sbatch` — three arms, because "off is inert" and "on works" are different claims

- **AFTER_OFF must be byte-identical to BEFORE over the whole output tree.** The unit differential
  pins the predicate; only a run pins the seam. Both decision ledgers are compared **outright** here,
  unlike the gate [B2] ledger gates — nothing in this change is *supposed* to alter a ledger.
- **AFTER_ON is the wiring proof, not the measurement.** chr1's gate B has 5 L1 survivors and **0
  SVA**, so its SVA delta may be zero: the fix's population is at WGS. A zero delta at chr1 is not
  evidence the fix does nothing, and the script says so where a reader will hit it.
- **Check 0 reads the driver banner in all three arms**, because with two arms sharing one worktree the
  thing most likely to go wrong is the *flag*, not the commit.
- Check 3 is the same prediction the WGS gate makes: L1 and Alu must not move with the fix on.

### 4.2 `gate_decoy_structural_wgs.sbatch` — one arm, and the frozen baseline is the other one

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

## 5. What NOT to redo

Plan §10 is the authority. Everything in the previous handoff's §6 still stands except where noted.

- **The decoy classification is `transduction_evidence.decoy_tags`, and there must not be a second
  copy.** There were two live ones and there is still a third in dead code (`:1051`, exact equality).
- **`b_structural=False` is a frozen artifact, not a preference.** Do not "clean up" the substring
  branch: it is what makes every existing gate and BASELINE-0 comparable, and a test pins it against a
  verbatim copy of the pre-fix expression.
- **Do not flip the default on the strength of the fix being correct.** It is correct; that is not the
  question. The question is what it does to SVA precision on a truth set whose transduction
  denominator is 2, and §4.2 check 6 is a direction rather than an answer.
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

## 6. NEXT

1. **Decide the default.** Read §4.2's checks 1, 5 and 6 together: 1 says the fix is contained, 5 says
   it does the thing, 6 says what it costs on the only truth signal HG002 can give. If precision holds,
   flipping `--td-decoy-structural` on by default is a one-line change plus a re-freeze of BASELINE-0 —
   and re-freezing BASELINE-0 is the expensive part, because every gate in `docs/handoff/runs/` that
   compares against it needs its expectations restated in the same commit.
2. **Plan §5's semantic half, on gate [A]** (97% of the rejecting). Expect *precision*: 95% of
   BASELINE-0's misses never reach this decision. Unchanged from the previous handoff's item 2, and its
   `decoy_cross_family` caveat stays answered.
3. **`MIN_DISC_MAPPED_RATIO = 0.7`** (`global_values.py:122`) — the largest knob between gates [B] and
   [C], deciding all **146** `no_cns_disc_reads` loci, four times everything [B2]'s own thresholds
   decide, and still not in §5's interrogate list.
4. **The decoy-consensus concatenation (§3).** Its own arm; the recipe is exact and build-independent;
   it is what makes §5's SVA/Alu hazard testable at all. Sequence it *after* the default decision in
   item 1, or the two changes' effects on SVA are not separable.
5. **`NEARBY_CLIP=50`** — still 4th, still bounded at 38 rejections plus 80 survivors.
6. **The per-locus decision loop's remaining box**, and the staged-barrier handoff's items 2-5,
   unchanged.

## 7. Where things live

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
| the FN attribution | `benchmark/ledger_attribution.py` |
| the truth-set scorer | `benchmark/benchmark_truth_set.py score --calls <vcf> --truth-dir <dir>` |
| the funnel reconciliation | `docs/handoff/runs/probe_decision_funnel.py` |
| worktrees for this session | `/n/scratch/users/c/cos689/xtea2_wt/td_decoy_55055ba` (AFTER, at `11b0e61`; `xtea_scan` **copied** from the BEFORE worktree, legitimately — `src/` is byte-identical across the range), BEFORE is `…/td_b2_365f352` |
| everything else | see `HANDOFF_2026-08-04_cns_cluster_ledger.md` §8 |
