# xTea2 handoff — 2026-07-30 (the agenda moved to the rewrite plan; HG002 germline inputs acquired)

**This is the current entry point for state. It is NOT the entry point for agenda — the plan is.**

> **Read `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` first.** It was written 2026-07-29 11:52,
> an hour after the previous handoff, and it redirects the whole project: a **full downstream
> overhaul** (unify the three per-type pipelines into one locus-major pass, rewrite transduction around
> an explicit evidence model, judge on truth sets instead of agreement with master), sequenced in five
> stages with **Stage 0 blocking everything**. This handoff records what state that plan is in and what
> was done to Stage 0 — it does not restate the plan.

**One-line status: nothing is in flight, no pipeline code has been touched, the tree is unchanged at
`ca1390d`. Stage 0's *acquisition* item is complete — the HG002 30x alignment and all five truth sets
are staged and verified at `test_data/HG002/` — and everything else in Stage 0 (phase timers,
BASELINE-0, chr1 slices, the `col5` rename, the type-set-dependence fix) is not started. The one
germline blocker that acquisition could not clear is what `_tprt` means, which sets the transduction
sensitivity denominator and needs the paper, not the code.**

## 0. State at handoff

| | |
|---|---|
| branch / HEAD | `xTea2` @ `ca1390d` |
| working tree | clean except the two untracked handoff docs (this one and the 07-29 one) |
| pipeline code | **untouched this session** |
| tests | 82 passed, 18 skipped as of `ca1390d` (skips are the unbuilt C scan bin, pre-existing) |
| in flight | **nothing** |
| germline benchmark inputs | **staged and verified** — see §2 |
| the agenda | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |

## 1. Why the previous handoff's "next steps" are no longer the next steps

`HANDOFF_2026-07-29_transduction_inputs.md` is now bannered SUPERSEDED. Its §5 agenda was: measure the
`fb1d766` BAM cache, then build a shared cross-type collection pass, then re-baseline the NA12878
goldens. Against the plan:

| previous agenda item | status under the plan |
|---|---|
| measure fix 1's cache (resubmit `perf_validate.sbatch`) | **drop it.** A constant-factor fix (~387s ceiling) to code the plan rewrites, measured against a comparator the plan **retires**. Not worth a WGS job. |
| shared collection pass, 1.89x (§4 there) | **survives, as plan Stage 1** — and explicitly conditional. The plan forbids pre-committing the C++ boundary until the phase breakdown exists, because the gather is somewhere between 26% and 80% of the step. |
| re-baseline the NA12878 goldens | **replaced** by BASELINE-0 on **HG002**, for `L1,Alu,SVA` *and* `L1` alone. NA12878's frozen comparators (master 50:18, port 54:39, chr1 goldens) become historical context, not gates. |
| the family-agnostic superset question (§3.3 there) | still a real open design question, now folded into the plan's Stage 2/3 family-assignment work rather than standing alone. |

**Everything factual in that handoff still holds** — the `fb1d766` fixes are correct and output-neutral,
fix 2's core split is verified live at WGS scale, transduction's raw-disc input is family-agnostic by
construction and must not be filtered, `.slct` is not a subset of `raw_disc`, and its whole "What NOT
to redo" list. Do not re-derive any of it.

## 2. What was acquired, and the one deviation that needs honoring

Staged at `/n/data1/hms/dbmi/park/corinne/claude_sandbox/test_data/HG002/` — outside the git repo, per
`.gitignore`'s "test data lives outside the repo". **`test_data/HG002/PROVENANCE.md` is the authority**
on all of this; it records the header dumps, the reasoning, and the dead ends. Do not re-derive from
this summary.

| path | what | verified |
|---|---|---|
| `alignment/HG002.novaseq.pcr-free.30x.dedup.grch38.bam` (+`.bai`) | 43,667,905,004 B | exact size, `md5 7378e27ba9a623ccb5eb1a4353aaebd9`, `samtools quickcheck` |
| `truth/{L1,Alu}_{all,tprt}.vcf`, `truth/SVA_tprt_all.vcf` | the paper's truth sets | pinned to `parklab/xTea_paper@2338ece` |
| `truth/x_cmp.py` | the paper's comparison script, kept as a **calibration reference, not a gate** | |
| `fetch_hg002.sbatch`, `logs/` | the resumable download job (ran in 10:43 as `49171775`) | |

**The deviation: the BAM is aligned to hs38DH, not `no_alt`.** Header reads
`bwa mem -t 16 references/hs38dh_bwa_index/hs38DH.fa`, 3,366 `@SQ`, MarkDuplicates 4.1.2.0,
`SM:HG002`. The plan's §7.1 wanted the paper's `no_alt` CRAM specifically because alt-aware placement
depresses MAPQ in repeat-homologous regions and xTea filters on MAPQ (clip ≥ 12, disc ≥ 20).

That file cannot be had — this is settled, do not re-chase:
- the paper's CRAM was built in-house under `simon_chu/projects/XTEA/`, **which no longer exists**, and
  is not hosted publicly under that or any other name;
- **GIAB hosts no bwa-mem GRCh38 Illumina alignment for HG002 at all** (only novoalign, per
  `data_indexes/AshkenazimTrio/`) — a different aligner is a *worse* confounder;
- `vinay/datasets/for_acb/HG002_30x_2_hg38_bwamem.sorted.bam`, which looks exactly right by name, is a
  **dangling symlink** into a deleted `park/DATA/GIAB/` tree. Every other local HG002 BAM is the
  `1mb_bin` subset.

**Why the substitute was accepted:** the faithful route is a FASTQ download plus a bwa-mem run, buying
comparability with the *published* numbers that the plan explicitly does not gate on (§2 and §7.3 —
`x_cmp.py` is a calibration reference; the real gates are BASELINE-0, an internal before/after diff,
and truth-set sensitivity/PPV, both of which use the same BAM on both sides). And the pipeline is
already proven on this reference shape: **`NA12878.wgs.bam` is bwa-mem 0.7.15, GRCh38 *with*
alt/decoy, also 3,366 `@SQ`.**

**What to honor:** our absolute sensitivity/FDR may sit below the paper's table for alignment reasons
alone. **Never present them side by side without saying so.** A `no_alt` rebuild is cheap to start if
that ever becomes the point — the FASTQs are in the same bucket, and the paper's exact reference is
still local and already bwa-indexed at
`simon_chu/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`.

**Reference to run against:** `/n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa`, the
same one the NA12878 runs use — a 3,366-contig alignment with a 195-contig primary reference is the
already-proven pairing. *Unverified:* whether anything in the pipeline trips on this BAM specifically.
The NA12878 precedent makes it likely fine; the first HG002 run is the actual test.

## 3. Two findings worth not rediscovering

- **The truth sets were already local**, at `corinne/dev_xTea/truth_set_hg002/` (Dec 2024),
  **byte-identical** to the fresh clone for all five VCFs. The staged copies add only the pinned
  commit. That older directory also holds **`all_tes.bed` (39 KB), which the paper repo does not
  ship** — worth identifying before writing `benchmark/benchmark_truth_set.py`, since a combined
  all-family annotation may already be what the harness needs.
- **The truth VCFs really are headerless** (no `##`, no `#CHROM`) — the plan flagged this as a parser
  hazard and it is now confirmed on the staged files rather than assumed.

## 4. NEXT — the rest of Stage 0, in dependency order

Straight from the plan; nothing here is new invention.

1. **Instrument transduction.** Phase timers around the blocks in
   `x_transduction.py:call_candidate_transduction_v3` (line 173) — flank-library build + bwa index
   (`:179-191`), the ~49K-site clip gather (`:196`), clip realign + the 2.4 GB SAM parse (`:199-207`),
   disc collect/realign/parse (`:221-233`), cns realign and the rest (`:250+`). The existing
   `XTEA_TD_PROFILE` gate (`x_orphan_transduction.py:18`) instruments only `parse_td_sibling_by_site`
   (~305 calls), **not** the gather. Consider emitting unconditionally — five clock calls — so this
   never needs re-running instrumented.
2. **Get the breakdown** via `docs/handoff/runs/transduction_attribution.sbatch` (staged cns inputs;
   its `impl_spec` loop is trivially extended). **This is what unblocks the Stage 1 scope decision** —
   2x off 26% of the step is a different call from 2x off 80%.
3. **The two standalone fixes**, before any baseline is frozen so every later baseline is against
   correct behaviour: the `col5` → `mate_in_rep_clip_cnt` rename, and the **type-set-dependence bug**
   (`run_xtea.py:837` intersects requested types with `--col5-families`, so `--types SVA` can emit
   *fewer* SVA calls than `--types L1,Alu,SVA`). Signature of a correct fix: three-type output does not
   move at all; SVA-only output does.
4. **Freeze BASELINE-0 on HG002** for both `L1,Alu,SVA` and `L1` alone, plus **chr1 slices** of HG002
   and SMHT005-3G — without a fast loop every gate is a ~50 min job.
5. **Resolve `_tprt`** (needs the paper) and decide how transduction is scored at all — the truth sets
   carry no source-element annotation, and a right-position/wrong-source call is a different error from
   a miss.

## What NOT to redo

Carried forward from `HANDOFF_2026-07-29_transduction_inputs.md` — **its entire "What NOT to redo"
list still applies**, including: no algorithmic port-vs-master difference exists in transduction; don't
re-derive per-candidate cost from the ~300-line candidate list; don't treat chr1 timings as evidence
about the cache; don't filter `raw_disc` by family or down to the candidate list; compare raw-disc
files across types with `cut -f1,2` only; don't put comma-bearing values in `sbatch --export`; don't
run downstream types serially at `-n 24`. Plus the plan's §10.

New this session:

- **Do not look for the paper's HG002 CRAM, and do not look for a bwa-mem GRCh38 HG002 alignment on
  GIAB.** Neither exists (§2).
- **Do not use `vinay/datasets/for_acb/HG002_30x_2_hg38_bwamem.sorted.bam`** — dangling symlink,
  however right the filename looks.
- **Do not re-download the truth sets** — already local twice over, byte-identical, pinned.
- **Do not measure the `fb1d766` cache.** Superseded; see §1.
- **Do not treat this handoff as the agenda.** The plan is.

## Where things live

| | |
|---|---|
| **the agenda** | `docs/handoff/DOWNSTREAM_OVERHAUL_PLAN.md` |
| germline benchmark inputs + full provenance | `test_data/HG002/`, esp. `PROVENANCE.md` |
| the pre-existing local truth sets (+ `all_tes.bed`) | `/n/data1/hms/dbmi/park/corinne/dev_xTea/truth_set_hg002/` |
| the paper's `no_alt` reference, bwa-indexed | `simon_chu/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` |
| mosaic truth set + inputs | `park-smaht_dac/DATA/…/SMHT005-SMaHT-sMEI-v1.0.0.vcf`; plan §7.2 |
| transduction findings that feed Stage 1 | `HANDOFF_2026-07-29_transduction_inputs.md` §3-§4 |
| the `fb1d766` fixes and their tests | `xtea/x_alignments.py`, `xtea/x_clip_disc_filter.py`, `xtea/run_xtea.py`; `tests/test_transduction_core_split.py` |
| the attribution harness | `docs/handoff/runs/transduction_attribution.sbatch` |
| frozen worktrees | `perf_wt` (`fb1d766`), `port_wt` (`b44a85f`), `parent_wt` (`6fb6247`), `master_cmp_wt` (`46520ab`) |
| the approved architecture this plan re-sequences | `docs/phase_c_unified_downstream.md` |
| mosaic design, already done — do not re-plan | `docs/handoff/MOSAIC_PLAN.md` |
