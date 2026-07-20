# Case/control (somatic) mode: prepared plan

**Prereq:** WGS move-A.2 timing run (job `48228181`) comes back green — check that first if picking
this up later than the same session it was drafted in.

**Decisions already made (user, this session — do not re-litigate):**
- **Single VCF out, somatic calls only — no merged germline+somatic-with-labels VCF.** I raised the
  "run germline on both BAMs and merge with a GERMLINE/SOMATIC label" idea; user said no, keep the
  original design. This mirrors master exactly (see below): the control BAM is never run through a
  full germline call, only mined for targeted evidence at the tumor's own candidate sites.
- **First gate = synthetic self-vs-self smoke test** (control BAM = the tumor BAM itself, expect the
  final somatic count to collapse to ~0), not the real PAUTSB tumor/normal pair. Real-data validation
  deferred (see "Validation" below for a *free* real-data option found this session).

## What I pinned during prep (so implementation is de-risked)

### 1. Master's real case-control pipeline order (confirmed from an ACTUAL prior run, not just code)
Found a real, already-executed case-control run at
`/home/cos689/data1/DATA/LFS_light_et_al/analysis/xtea_somatic/results/2921/L1/run_xTEA_pipeline.sh`
(tumor `3762A.recal.cram` vs control `3763.recal.cram`, LFS project, L1 only, master's own
`gnrt_pipeline_local_v38.py`-generated script, already run — outputs on disk). This is the ground
truth for exact command order/args, not a re-derivation from `x_TEA_main.py` reading alone:

1. **Per-type chain runs on the tumor BAM only**, exactly like today's germline chain (clip → disc →
   cns_filter → transduction → sibling → postF_all → postF_hc), but with **`--tumor --purity <p>`**
   appended to clip/disc/cns_filter/transduction/sibling, and **`--tumor`** (no `--purity`) on both
   postF steps. **xTea2's `build_type_steps` currently has no `--tumor`/`--purity` plumbing at all** —
   this needs adding as a prerequisite, independent of case-control (it changes the automatic
   clip/disc cutoff derivation via `f_cov*=f_purity` in `x_TEA_main.py:automatic_gnrt_parameters`).
2. `--case_control -b <control_bam_list> -i <postF_all output, NOT postF_hc>` → filters the RAW
   (non-HC) post-filter list against the control BAM via `CaseControlMode.call_somatic_TE_insertion`
   → a somatic-filtered raw list (`..._post_filtering_somatic.txt`).
3. `--case_control --somatic_hc -i <postF_hc output> --input2 <step 2's output>` →
   `CaseControlMode.parse_high_confident_somatic` intersects the two → final HC somatic list
   (`..._high_confident_post_filtering_somatic.txt`). **This, not postF_hc, is what should feed A.2's
   `run_shared_tail`/locus table when case-control mode is on.**
4. gene → gVCF run on that final somatic list, writing to the **same output filename** the germline
   gVCF step already used for this sample+type. In the real 2921/L1 run this means the case-control
   `--gVCF` call **silently overwrites** the earlier germline-chain VCF (same `{sample}_{type}.vcf`
   name) — confirmed by inspecting the file on disk: `3762A.recal_LINE1.vcf` contains only the 2 final
   somatic calls, byte-identical to `..._somatic_with_gene.txt`, even though an earlier step in the
   same script also wrote a (germline) VCF to that exact path. **This is master's own de-facto
   behavior, not a bug I'm introducing** — it independently confirms the user's call: one file, somatic
   only, nothing germline survives to disk. xTea2 should do this cleanly (skip writing the germline VCF
   at all when case-control is on, not write-then-overwrite).

### 2. `CaseControlMode.call_somatic_TE_insertion` (`xtea/x_somatic_calling.py:78-194`, confirmed byte-
identical to what's vendored in this repo) is NOT a naive "did clip/disc also fire in control at this
position" diff. Per site it:
- Splits tumor candidates into non-transduction / transduction / transduction-sibling / orphan buckets
  (`XClipDiscFilter.sprt_TEI_to_td_orphan_non_td`).
- Re-collects + realigns clip/disc reads from the **control BAM**, in a window around each tumor
  candidate (`iextnd=400`), and clusters them to a consensus — same machinery the germline chain uses
  on the tumor BAM, just pointed at control and restricted to the tumor's own candidate loci (NOT a
  whole-genome scan of control).
- Extracts per-site genotype features + local read depth from control
  (`extract_sites_raw_gntp_feature_depth` → `XGenotyper.call_genotype` + `ReadDepth`).
- `filter_by_case_ctrl_disc_depth_ratio`: a ~15-cutoff cascade comparing case-vs-control ratios of
  effective-disc/raw-disc/clip/polyA depth-normalized rates (site is dropped as germline if control
  shows comparable-or-stronger signal by several of these ratios). Numeric cutoffs are all hardcoded
  method-local defaults (`f_disc_ef_ratio=0.1`, `f_case_ctrl_cutoff=0.35`, etc.) — **not** currently
  exposed as CLI knobs (matches master; not adding new knobs unless something breaks).
- `parse_high_confident_somatic` step 3 above is a simple exact chrm/pos set-intersection — genuinely
  trivial, unlike step 2.

### 3. Where this slots into xTea2's current architecture (`xtea/run_xtea.py`)
- `build_type_steps` (run_xtea.py:133) already produces both `pf_out` (`postF_all`,
  `candidate_disc_filtered_cns_post_filtering.txt`) and `hc_pf` (`postF_hc`, via `hc_pf_path()`,
  run_xtea.py:125) per type — **both master case-control inputs already exist as named locals in
  today's code**, nothing new to compute there.
- New stage = 2 more steps per type, appended after `postF_hc`, gated behind `--case-control`:
  `case_control` (raw list, step 2 above) then `case_control_hc` (step 3 above, `--somatic_hc`).
  Natural fit as two more entries in `build_type_steps`'s `steps` list (respecting `stop_after`).
- `run_type_downstream` / `run_shared_tail` (run_xtea.py:266-336): when case-control is on, the path
  handed into `hc_pf_paths` (currently always `hc_pf_path(dirs["type"])`) must instead be the
  case-control step's final somatic-HC output. Cleanest: give `case_control_hc`'s output file the
  SAME conceptual role `hc_pf` plays today (a per-type "final raw call list", just somatic-filtered),
  and have `run_type_downstream` return *that* path instead of `hc_pf` when `--case-control` is set.
  `run_shared_tail`'s locus-table merge + gene/genotype/gVCF logic (move A.2) is otherwise **completely
  unaffected** — it doesn't care whether its input list is germline or somatic-filtered.
- `--tumor`/`--purity`: add as new `run_xtea.py` CLI args, threaded into `build_type_steps`' argv exactly
  where master's `gnrt_calling_command` puts `s_tumor`/`s_purity` (clip/disc/cns_filter/transduction/
  sibling get both when tumor; both postF steps get `--tumor` only). Independent of `--case-control`
  (a tumor-only single-sample run without a control is also a valid master mode, though out of scope
  to explicitly support/test right now — just don't block it accidentally).

### 4. bam_list format — verified NOT a risk (read the actual loaders, not just the README)
Every consumer on the case-control path (`XClipDiscFilter.__init__`, `XGenotyper.call_genotype`,
`X_BasicInfo.collect_basic_info_samples`) does `fields = line.split(); sf_bam = fields[0]` — the control
bam_list can be a plain one-path-per-line file (master's `bam_list2.txt` convention) or reuse xTea2's
existing `path\tillumina` two-column convention; only column 0 is ever read. No format-mismatch risk —
build the control bam_list the same way `run_xtea.py:456` already builds the tumor one.

## Validation — two tiers, found real free ground truth this session

1. **Cheap self-vs-self smoke test (user's chosen first gate):** control BAM = the tumor BAM. Every
   tumor candidate should show ≥ as much clip/disc/depth signal in "control" as in "case" →
   `filter_by_case_ctrl_disc_depth_ratio` should drop ~everything → final somatic VCF ≈ empty. Cheap,
   fast, chr1-scale is enough to wire up the plumbing before touching WGS.
2. **Free real-data gate, found this session, not yet run:** the exact LFS `2921/L1` run above has
   **real tumor (`3762A.recal.cram`) and control (`3763.recal.cram`) CRAMs already on disk** at
   `/home/cos689/LFS_light_et_al/analysis/sarek_revised/results/preprocessing/recalibrated/{3762A,3763}/`,
   **plus master's own precomputed final output** already sitting at
   `.../xtea_somatic/results/2921/L1/candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt`
   (2 L1 calls: chr7:5734825, chr19:21769582) — a real byte-level target to diff xTea2's new
   case-control output against, **no new compute needed to get the reference**, only to reproduce it.
   Much stronger than the self-vs-self smoke test; propose running this as gate 2 once gate 1 is green,
   before ever touching the PAUTSB pair or full WGS scale. (Cohort has other samples too —
   `~/LFS_light_et_al/analysis/xtea_somatic/` — 2921 is just the one I happened to inspect; more are
   available as extra gates if 2921 alone proves ambiguous, e.g. its 2-call count is small.)

## Staged implementation (each harness-gated; STOP for approval before starting)

- **CC.0 — `--tumor`/`--purity` plumbing (prerequisite, no case-control logic yet).** DONE, chr1-gated
  2026-07-20 (implementation session). Added the two flags to `run_xtea.py`, threaded into
  `build_type_steps` exactly per master's pattern (§3 above). **Correction to this plan's original gate
  design (written before `x_TEI_locator.py`/`x_post_filter.py` were read line-by-line): `--tumor
  --purity 1.0` is NOT expected to reproduce a purity-less run byte-for-byte, and the chr1 gate run
  correctly does NOT show that (31 calls baseline vs 35 calls with `--tumor --purity 1.0` — same
  automatic coverage-derived cutoffs confirmed identical in both run logs, `Ave coverage is 33.471:
  ... (3, 5, 1)` in both, so the purity-scaling half of `--tumor` genuinely is a no-op at purity=1.0 as
  designed).** The count difference comes from `b_tumor` switching the disc-site AND postF retention
  *formulas* themselves, independent of purity or cutoff value:
  `x_TEI_locator.py:filter_candidate_sites_by_discordant_pairs_multi_alignmts` keeps a site when
  `(n_left+n_right) >= cutoff` under `--tumor` vs `(n_left>half_cutoff AND n_right>half_cutoff) OR
  n_left>cutoff OR n_right>cutoff` otherwise (lines 457-464); `x_post_filter.py` has several more
  `b_tumor`-gated branches (`get_cutoff_by_type`, one-side handling at lines 258/269/103). This is
  master's own real, by-design behavior — not a bug, and not something to make numerically match
  germline mode. **Actual CC.0 gate, corrected:** (a) argv shape matches master's real executed
  `run_xTEA_pipeline.sh` exactly (confirmed by direct diff, not just source reading); (b) both the
  purity-less and `--tumor --purity 1.0` chr1 runs complete end-to-end without crashing; (c) the
  automatic-cutoff VALUE is confirmed identical between the two runs (isolating that purity's specific
  effect — cutoff scaling — is a correct no-op at purity=1.0, decoupled from `--tumor`'s separate
  formula-switching effect). All three hold; CC.0 is complete.
- **CC.1 — case-control steps, chr1 self-vs-self.** DONE, chr1-gated 2026-07-20 (same session as CC.0).
  Added `--case-control`/`--normal-bam` to `run_xtea.py`; extended `build_type_steps` with the two new
  steps (`case_control`, `case_control_hc` — argv shape confirmed against master's real executed
  `run_xTEA_pipeline.sh`, not just source reading); control bam_list built the same way as the tumor
  one; `run_type_downstream` returns the `case_control_hc` path (new `cc_hc_path()` helper, mirrors
  `hc_pf_path()`) instead of `hc_pf` when `--case-control` is on, and the 3 call sites in `main()` now
  pick `stop_after="case_control_hc"` vs `"postF_hc"` based on `opt.case_control`. Also found (reading
  the real script closely): master's somatic path skips `--gntp_classify` entirely (gene output feeds
  `--gVCF` directly, confirmed by the real output VCF showing `GT=./.` for every somatic call) — made
  this conditional in `build_shared_tail_steps` too. **Gate result: chr1, control=tumor BAM, ALL 3
  types' `postF_all` candidates (Alu 39, L1 4, SVA 6) collapsed to 0 after case-control filtering; final
  merged VCF = 0 calls** — full collapse, not just "near-zero" (allowed slack in the original plan
  wording was unnecessary in practice for this run). Note: the gate script's own SLURM job status showed
  FAILED (exit 1) — cosmetic, not a pipeline problem: its trailing `grep -v '^#' "$mvcf" | head -20`
  report line exits 1 when the VCF has 0 data lines, and `set -euo pipefail` treats that as fatal *after*
  the actual run had already completed and printed everything needed (same class of harmless
  reporting-script bug as the one fixed in the move-A.2 WGS timing run, see the 2026-07-20 handoff).
- **CC.2 — wire into A.2's shared tail + WGS self-vs-self.** The *wiring* half turned out to already be
  complete as a side effect of CC.1's implementation (`run_type_downstream`'s generic
  `hc_pf_paths[t] = v` assignment doesn't care whether `v` is a germline `hc_pf` path or a somatic
  `case_control_hc` path — `locus_table.assign_and_write` and the shared gene/gVCF tail are format-
  agnostic, confirmed by reading `parse_high_confident_somatic`: it writes matched lines from the
  hc_pf-format input verbatim, so `cc_hc_path`'s output is byte-format-identical to `hc_pf_path`'s). No
  separate wiring step was needed. What's left of CC.2 is purely the **WGS-scale self-vs-self gate**
  (same near-zero expectation as CC.1, at full scale — also a timing checkpoint since case-control adds
  a full second BAM's worth of targeted clip/disc realignment per type). Submitted 2026-07-20 as job
  `48251912` (`cc2_validate_wgs_selfvself.sh`) — check that job first if picking this up later.
- **CC.3 — real-data gate (LFS 2921/L1).** Not yet run. Script pre-staged at
  `cc3_lfs_2921_L1.sh` (uses master's OWN reference, `Homo_sapiens_assembly38.fa` via `cgap_matches` —
  NOT `hg38_no_alt.fa` used by the other gates — because that's what the real `3762A`/`3763` CRAMs were
  aligned/decoded against). Run xTea2's case-control end-to-end on the real `3762A`/`3763` CRAM pair;
  diff against master's precomputed 2-call output
  (`.../xtea_somatic/results/2921/L1/candidate_disc_filtered_cns.txt.high_confident.post_filtering_somatic_with_gene.txt`,
  chr7:5734825 + chr19:21769582). This is the actual bar-clearing gate, not CC.1/CC.2 (those just prove
  the plumbing doesn't crash/doesn't trivially over/under-filter). WGS-scale real CRAMs (61GB/78GB) —
  budget several hours, not chr1-scale minutes.

**Sequence rationale:** CC.0 isolates a cross-cutting flag change from the case-control logic itself so
a regression is unambiguous about its source. CC.1→CC.2 is chr1→WGS at synthetic zero-signal before any
real biology is on the line. CC.3 is the first point real recall/precision claims can be made, and it's
free (data + reference output already exist) — no reason to skip straight to PAUTSB or a fresh WGS
tumor/normal run before this cheaper real check clears.

## Open items to confirm at CC.0 start
- Master's `--purity` default is 0.45; the real 2921/L1 run used `0.86` (sample-specific, presumably
  from a purity estimator not part of xTea itself) — xTea2 should expose it as a plain user-supplied
  float, not try to estimate it.
- Confirm `automatic_gnrt_parameters_case_control`'s `CaseControlFilterPars` cutoff table
  (`x_parameter.py:114-133`, coverage-keyed like the germline `Parameters` table) is reachable/correct
  as-is via `x_TEA_main.py`'s existing `--case_control` automatic-param path — no evidence it needs
  changing, just hasn't been exercised through xTea2's wrapper yet.
- `--fflank`/`--flklen` are required on the case-control steps (transduction flank realignment) — same
  values `build_type_steps` already threads through today's transduction/sibling steps, just needs
  passing again into the two new steps.
