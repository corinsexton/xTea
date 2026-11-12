#!/usr/bin/env python
"""xTea2 single-process orchestrator  --  Milestone A0 (wrapper-first).

Goal of this milestone: reproduce master's golden call-set from ONE self-contained
process (no SLURM script generation, no per-step hand submission), by driving the
exact `-f 5907` step chain that produced the frozen golden.

Design (A0):
  * This driver runs the 11 pipeline steps *in one process* by shelling out to the
    existing step CLI (`x_TEA_main.py`) once per step -- the same commands the golden
    `run_xTEA_pipeline.sh` scripts issued, minus the SLURM headers and minus job
    generation. It is therefore golden-reproducing *by construction*.
  * L1 / Alu / SVA are run inside a single invocation; clipped-read collection is
    shared across types by the one xtea-scan pass that opens the BAM.
  * All filter/threshold knobs are exposed on the CLI with defaults pinned to the
    golden run (`--nfclip 4 --nfdisc 5 --flklen 3000`), so the defaults reproduce
    the golden and every change is opt-in.

Trim step 4 committed this driver to the scan-only path: candidate generation is
always the C++ xtea-scan pass (clip) plus the native disc-pos adapter (disc), and
master's pysam `-C`/`-D` steps are deleted, not merely bypassed. `--scan-clip`
and `--scan-disc` survive as accepted no-ops so the existing validation scripts
under /n/scratch/users/c/cos689/xtea2_testdata keep running unchanged.

Next milestones (A1..): replace each `x_TEA_main.py` subprocess with a native
in-process call, validating per-step against the frozen golden intermediates
(`candidate_list_from_clip.txt`, `candidate_list_from_disc.txt`, ...), then
integrate/merge steps where that is faster (results stay the contract).

This file intentionally does NOT reimplement detection logic yet -- it is the
always-green scaffold we refactor inward.
"""

import argparse
import concurrent.futures
import contextlib
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import scan_adapter
import disc_adapter
import downstream_unified
import family_params
import global_values
import locus_table
import region_mask

HERE = Path(__file__).resolve().parent
XTEA_MAIN = str(HERE / "x_TEA_main.py")


def _default_scan_bin():
    """Where to find `xtea_scan` when --scan-bin is not given.

    THE REPO-RELATIVE PATH CANNOT BE THE ONLY ANSWER once this is pip-installable: the C++ engine
    is built by `make -C src/csrc` and is deliberately not tracked in git, so after an install
    `HERE.parent/src/csrc/` is a site-packages directory that will never contain it. Resolution
    order, first hit wins:

      1. $XTEA2_SCAN_BIN   -- the same variable the test suite uses, so one export covers both
      2. the in-repo build -- keeps every existing checkout, harness and gate working unchanged
      3. `xtea_scan` on PATH -- the installed case

    Returns the repo path as the fallback even when nothing exists, so the "not found" error names
    the location a developer most likely expected rather than an empty string.
    """
    env = os.environ.get("XTEA2_SCAN_BIN", "")
    if env:
        return env
    in_repo = HERE.parent / "src" / "csrc" / "xtea_scan"
    if in_repo.exists():
        return str(in_repo)
    on_path = shutil.which("xtea_scan")
    return on_path or str(in_repo)


DEFAULT_SCAN_BIN = _default_scan_bin()

# --- per-type annotation layout under <rep_lib>/<subdir>/<genome>/ ---------------
# rtype is the bitmask x_TEA_main uses (L1=1, Alu=2, SVA=4). vcf_alt is the ALT tag
# x_gvcf emits, used to locate the per-type VCF produced by the gVCF step.
REP_TYPES = {
    "L1": {
        "rtype": 1,
        "vcf_alt": "LINE1",
        "copies_with_flank": "LINE/{g}/{g}_L1HS_copies_larger_5K_with_flank.fa",
        "annotation":        "LINE/{g}/{g}_L1_larger_500_with_all_L1HS.out",
        "annotation1":       "LINE/{g}/{g}_L1.fa.out",
        "flank":             "LINE/{g}/{g}_FL_L1_flanks_3k.fa",
        "cns":               "consensus/LINE1.fa",
    },
    "Alu": {
        "rtype": 2,
        "vcf_alt": "ALU",
        "copies_with_flank": "Alu/{g}/{g}_AluJabc_copies_with_flank.fa",
        "annotation":        "Alu/{g}/{g}_Alu.out",
        "annotation1":       "Alu/{g}/{g}_Alu.out",
        "flank":             None,  # no transduction flank for Alu
        "cns":               "consensus/ALU.fa",
    },
    "SVA": {
        "rtype": 4,
        "vcf_alt": "SVA",
        "copies_with_flank": "SVA/{g}/{g}_SVA_copies_with_flank.fa",
        "annotation":        "SVA/{g}/{g}_SVA.out",
        "annotation1":       "SVA/{g}/{g}_SVA.out",
        "flank":             "SVA/{g}/{g}_FL_SVA_flanks_3k.fa",
        "cns":               "consensus/SVA.fa",
    },
}

# accept a few aliases on the command line
TYPE_ALIASES = {"l1": "L1", "line1": "L1", "line": "L1",
                "alu": "Alu", "sva": "SVA"}


# Seconds since the driver started, stamped on every log line. Cheap, and it is the only
# thing that makes CONCURRENCY readable after the fact: the timing file records durations,
# so two rows that overlap and two that abut look identical in it. Under the staged barrier
# the whole claim is about overlap -- one family entering transduction while another is
# still in cns_filter -- and without a stamp that is unfalsifiable from a finished run.
# Nothing greps the `[xtea2]` prefix (checked across docs/handoff/runs and tests); every
# reader matches on the message.
_LOG_T0 = time.monotonic()


def log(msg):
    print(f"[xtea2 +{time.monotonic() - _LOG_T0:8.1f}s] {msg}", flush=True)


def resolve_annot(rep_lib, genome, rel):
    if rel is None:
        return "null"
    return str(Path(rep_lib) / rel.format(g=genome))


# --mate-in-rep-families keys -> REP_TYPES keys. The CLI vocabulary is lowercase.
MATE_IN_REP_KEY_TO_TYPE = {"l1": "L1", "alu": "Alu", "sva": "SVA"}

# Default contributing families for mate_in_rep_clip_cnt: the biologically tightest set.
# The value is consumed only by the SVA cns path, and SVA is a composite element with an
# Alu-derived domain, so a mate landing in Alu is genuine SVA evidence.
MATE_IN_REP_FAMILIES_DEFAULT = "alu,sva"


def mate_in_rep_bed_plan(mate_in_rep_fams):
    """Which families' annotation BEDs feed mate_in_rep_clip_cnt.

    Returns [(rep_type, family, annotation_extend)], ordered stably by family key.

    mate_in_rep_clip_cnt is the number of clipped reads at a site whose mate lands inside
    any selected repeat family's annotation -- a repeat-anchored support signal. It is
    emitted as the 6th column of xtea-scan's realign-out, which is the only reason it was
    ever called `col5`: a column position, not a concept. Master called it
    `mate_within_rep_cnt` and computed it Alu-only (non-deterministically -- whichever
    type won a shared-`pub_clip` race); xTea2 deliberately uses the union over families,
    each family's boundary extend baked into its BED (Alu/L1 = 0, SVA = 200), emitted
    identically for every family at a site. Spec §4 DC1 has the measured justification.

    NOTE THE SIGNATURE: this deliberately does NOT take the requested `--types`. The
    contributing family set is a property of the signal definition, not of what the caller
    asked to be called. Intersecting the two, which is what the inline loop in main() used
    to do, made the BED set `requested_types & mate_in_rep_families`, so the value -- and
    hence the SVA cns consistency check it gates -- varied with the requested type set:

        --types L1,Alu,SVA  ->  ALU + SVA  ->  mate-in-(Alu | SVA)
        --types SVA         ->  SVA only   ->  mate-in-SVA only
        --types L1          ->  no BEDs    ->  0 at every site

    Since union >= subset always, `--types SVA` could emit FEWER SVA calls than the
    same sample run as `--types L1,Alu,SVA`. Asking for less must not change what you
    get for what you did ask for.

    Keeping `types` out of this signature is what stops that being reintroduced;
    tests/test_mate_in_rep_type_partition.py pins it.
    """
    plan = []
    for fam_key in sorted(mate_in_rep_fams):
        t = MATE_IN_REP_KEY_TO_TYPE[fam_key]
        extnd = family_params.for_name(t).annotation_extend
        plan.append((t, scan_adapter.TYPE_TO_FAMILY[t], extnd))
    return plan


def detect_format_and_index(sf_input, reference):
    """Auto-detect BAM vs CRAM by extension; ensure an index exists.

    CRAM requires --reference (threaded to htslib for decode). Returns nothing;
    raises on a missing prerequisite so failures are loud and early.
    """
    p = Path(sf_input)
    if not p.exists():
        sys.exit(f"ERROR: input not found: {sf_input}")
    is_cram = p.suffix.lower() == ".cram"
    if is_cram and not reference:
        sys.exit("ERROR: CRAM input requires --reference (FASTA) for decode.")
    idx_candidates = [p.with_suffix(p.suffix + (".crai" if is_cram else ".bai")),
                      p.with_suffix(".crai" if is_cram else ".bai")]
    if not any(c.exists() for c in idx_candidates):
        log(f"index missing for {p.name}; building with samtools index")
        subprocess.run(["samtools", "index", str(p)], check=True)


def resolve_inputs(input_arg, reference):
    """Parse a (possibly space-delimited) -i/--normal-bam value into N resolved paths.

    Functionally identical for one path or many: each is checked for existence and
    indexed via detect_format_and_index. Multiple libraries of the same sample are
    pooled downstream via bam_list.txt (master's own pre-existing mechanism).
    """
    paths = input_arg.split()
    if not paths:
        sys.exit("ERROR: empty --input/--normal-bam value")
    resolved = []
    for tok in paths:
        p = Path(tok)
        if not p.exists():
            sys.exit(f"ERROR: input not found: {tok}")
        detect_format_and_index(str(p), reference)
        resolved.append(str(p.resolve()))
    return resolved


def check_rg_consistency(paths, reference, ignore_rg):
    """Validate all input files share one @RG SM tag; fail loud (with specifics) if not.

    Pooling libraries from different sequencing centers under one sample is the
    intended use of multi-file --input, but pooling files that turn out to be
    DIFFERENT samples by accident is a real bug class -- so this is fail-loud by
    default. --ignore-rg explicitly opts out (e.g. when centers legitimately use
    different SM tags for the same donor/tissue, as with GCC-split SMaHT data).
    """
    import pysam
    by_path = {}
    for p in paths:
        bf = pysam.AlignmentFile(p, "rb", reference_filename=reference)
        rgs = bf.header.get("RG", [])
        bf.close()
        if not rgs:
            if ignore_rg:
                log(f"WARNING (--ignore-rg): no @RG header in {p}")
                continue
            sys.exit(f"ERROR: no @RG header in {p} (pass --ignore-rg to proceed anyway)")
        sms = sorted({rg.get("SM") for rg in rgs if rg.get("SM")})
        if not sms:
            if ignore_rg:
                log(f"WARNING (--ignore-rg): @RG present but no SM tag in {p}")
                continue
            sys.exit(f"ERROR: @RG present but no SM tag in {p} "
                     f"(pass --ignore-rg to proceed anyway)")
        if len(sms) > 1:
            if ignore_rg:
                log(f"WARNING (--ignore-rg): file {p} itself has multiple SM tags: {sms}")
            else:
                sys.exit(f"ERROR: file {p} itself has multiple SM tags: {sms} "
                         f"(pass --ignore-rg to proceed anyway)")
        by_path[p] = sms

    if ignore_rg:
        all_sms = sorted({sm for sms in by_path.values() for sm in sms})
        log(f"--ignore-rg: proceeding; SM tags found across input: {all_sms or '(none)'}")
        return

    all_sms = sorted({sm for sms in by_path.values() for sm in sms})
    if len(all_sms) > 1:
        by_tag = {sm: sorted(p for p, sms in by_path.items() if sm in sms) for sm in all_sms}
        detail = ", ".join(f"SM:{sm} ({', '.join(ps)})" for sm, ps in by_tag.items())
        sys.exit(f"ERROR: {len(all_sms)} distinct SM tags across input: {detail} "
                 f"(pass --ignore-rg if pooling different-center libraries of the "
                 f"same sample is intentional)")
    log(f"@RG check: single SM tag across all input files (SM:{all_sms[0]})")


# --- timing / peak-memory instrumentation ----------------------------------------
# Every phase and every step appends one row to <sample_root>/step_timings.tsv, so a run's
# cost profile is a data file rather than something to be reconstructed from log
# timestamps. Two caveats, both real:
#   * peak_rss_gb is getrusage(RUSAGE_CHILDREN).ru_maxrss, i.e. the high-water mark across
#     ALL children so far -- it is monotone, never per-step. The useful reading is *which*
#     row first raised it. Recorded as-is rather than as a bogus per-step delta.
#   * under --parallel-downstream the three types' steps overlap, so their wall times are
#     concurrent and must not be summed. The phase rows (which do not overlap) are what
#     add up to the run.
_TIMING_PATH = None
_TIMING_LOCK = threading.Lock()


def _children_peak_gb():
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / (1024.0 * 1024.0)


def record_timing(kind, name, seconds):
    """Append one timing row. No-op until init_timing() has run."""
    if _TIMING_PATH is None:
        return
    with _TIMING_LOCK:
        with open(_TIMING_PATH, "a") as fh:
            fh.write(f"{kind}\t{name}\t{seconds:.2f}\t{_children_peak_gb():.2f}\n")


def init_timing(sample_root):
    global _TIMING_PATH
    _TIMING_PATH = str(Path(sample_root) / "step_timings.tsv")
    with open(_TIMING_PATH, "w") as fh:
        fh.write("kind\tname\twall_s\tchildren_peak_rss_gb\n")
    return _TIMING_PATH


@contextlib.contextmanager
def phase(name):
    """Time a driver-internal phase (the scan pass, the cns emit, the shared tail).

    THE LINE THIS PRINTS IS NOT A SUCCESS SIGNAL UNLESS IT SAYS SO. The timing used to be
    logged from a bare `finally`, so a phase that DIED still printed the identical
    "phase <name>: <n>s" line as one that succeeded. uwsc (job 50447196) printed
    "phase scan_pass: 25665.4s" at the exact second SLURM OOM-killed it -- a child was
    OOM-killed, the subprocess error propagated out through the finally, and the phase
    announced its duration on the way past. Read as completion that line says a 7-hour scan
    finished; on disk scan_clip_realign.tsv was never written at all. Anything deciding
    whether output is reusable must therefore key on the OUTPUT, not on this log line --
    which is why --adopt-existing verifies every expected file itself.
    """
    t0 = time.monotonic()
    ok = False
    try:
        yield
        ok = True
    finally:
        dt = time.monotonic() - t0
        record_timing("phase" if ok else "phase_FAILED", name, dt)
        peak = _children_peak_gb()
        if ok:
            log(f"phase {name}: {dt:.1f}s (children peak RSS {peak:.2f} GB)")
        else:
            log(f"phase {name}: FAILED after {dt:.1f}s (children peak RSS {peak:.2f} GB)")


# --- resume: skip an expensive phase whose output is already on disk -----------------------
#
# WHY THIS EXISTS. scan_pass and cns_emit are the two phases that cost hours (broad: 1.8h +
# 0.5h; uwsc was still in scan_pass at 6h), and everything that has actually gone wrong in
# this fleet went wrong AFTER them -- the 192G OOMs landed during/just after the scan, the
# `| head` pipefail bug failed every arm at the very end, and a downstream code fix means a
# rerun. Losing a completed 6-hour scan to a fault in a later step is the single most
# expensive failure mode here, and it is pure waste: the scan's inputs did not change.
#
# A marker is written ONLY after its phase returns, so a phase killed mid-write (the OOM
# case) leaves no marker and is correctly redone. Existence of the output files is NOT
# sufficient on its own and is never used alone -- an OOM leaves a plausible-looking
# truncated .tsv behind.
#
# The fingerprint is what makes reuse safe. It covers every input that determines the
# output, INCLUDING the scan binary's size+mtime: resuming across a `make -C src/csrc`
# would otherwise silently reuse output from the previous binary, which is exactly the
# stale-artifact trap conftest.py already guards for the test suite.
#
# INTERACTION WITH --keep-intermediates. cleanup_successful_run wipes everything under
# sample_root except the merged VCF, and .resume/ lives there, so a SUCCESSFUL run without
# --keep-intermediates removes both the markers and the outputs they describe. That is the
# right outcome and needs no special case: there is nothing to resume after a success, and
# a later --resume simply finds no marker and reruns. The case --resume is FOR is a run that
# died -- cleanup never executes on that path, so the scan output and its marker both
# survive. The two flags are complementary, not alternatives; the only combination that
# does nothing useful is --resume after a successful run that was allowed to clean itself.
#
# If a marker somehow outlives its outputs (a partial cleanup, a scratch purge), the
# per-output existence+size check below rejects it. Reuse is never decided by the marker
# alone.
def _resume_dir(sample_root):
    return Path(sample_root) / ".resume"


def _fingerprint(parts):
    """Stable hash of everything that determines a phase's output."""
    h = hashlib.sha256()
    for p in parts:
        h.update(repr(p).encode())
        h.update(b"\0")
    return h.hexdigest()


def _binary_stamp(path):
    try:
        st = os.stat(path)
        return (str(path), st.st_size, int(st.st_mtime))
    except OSError:
        return (str(path), None, None)


def phase_is_done(sample_root, name, fingerprint, outputs):
    """True when `name` completed earlier with this exact fingerprint and its outputs survive.

    Size is compared, not just existence: a resumed run must not build on a file that was
    truncated after the marker was written (a full scratch, a partial copy).
    """
    marker = _resume_dir(sample_root) / f"{name}.json"
    if not marker.is_file():
        return False
    try:
        rec = json.loads(marker.read_text())
    except (OSError, ValueError):
        return False
    if rec.get("fingerprint") != fingerprint:
        log(f"[resume] {name}: marker present but inputs changed -- rerunning")
        return False
    for path, size in rec.get("outputs", []):
        if not os.path.isfile(path) or os.path.getsize(path) != size:
            log(f"[resume] {name}: output missing or resized ({path}) -- rerunning")
            return False
    for path in outputs:
        if not os.path.isfile(path):
            log(f"[resume] {name}: expected output absent ({path}) -- rerunning")
            return False
    return True


def adopt_existing_phase(sample_root, name, fingerprint, outputs, adopt):
    """Write a marker for a phase whose output predates --resume, so it can be reused.

    ONLY for adopting runs made before resume markers existed (or by a different driver).
    It trusts output this process did not produce and cannot verify: a file truncated by an
    OOM looks exactly like a complete one here, which is precisely the check the normal
    marker gives you for free by being written only after the phase returns.

    So it demands every expected output be present and NON-EMPTY, refuses if any is missing,
    and logs loudly what it took on trust.

    DO NOT confirm completion from the old log's "phase <name>: <n>s" line. Before the fix in
    phase(), that line was printed from a bare `finally` and appeared identically whether the
    phase succeeded or died: uwsc (job 50447196) printed "phase scan_pass: 25665.4s" in the
    same second SLURM OOM-killed it, having never written scan_clip_realign.tsv. The
    output-presence check below is what actually caught that, and it is the only check here
    that can.
    """
    if not adopt:
        return False
    missing = [p for p in outputs if not os.path.isfile(p) or os.path.getsize(p) == 0]
    if missing:
        log(f"[resume] --adopt-existing: {name} NOT adopted, {len(missing)} output(s) "
            f"missing or empty (first: {missing[0]})")
        return False
    mark_phase_done(sample_root, name, fingerprint, outputs)
    log(f"[resume] --adopt-existing: ADOPTED pre-existing {name} output on trust "
        f"({len(outputs)} file(s)); this process did not verify it is complete")
    return True


def mark_phase_done(sample_root, name, fingerprint, outputs):
    d = _resume_dir(sample_root)
    d.mkdir(parents=True, exist_ok=True)
    rec = {"phase": name, "fingerprint": fingerprint,
           "outputs": [[p, os.path.getsize(p)] for p in outputs if os.path.isfile(p)]}
    (d / f"{name}.json").write_text(json.dumps(rec, indent=2))


@contextlib.contextmanager
def timed(kind, name):
    """phase(), but the caller says which row kind it is.

    The staged barrier runs work that a global barrier ran between waves -- the
    driver gather in particular -- concurrently with other families' steps. That
    work is the same work, but it is no longer a phase: the timing file's contract
    is that phase rows are disjoint and sum to the run, and step rows overlap.
    Recording an overlapping row as a phase would silently break that sum, so the
    barrier records "step" and the serial path keeps "phase".
    """
    t0 = time.monotonic()
    try:
        yield
    finally:
        dt = time.monotonic() - t0
        record_timing(kind, name, dt)
        log(f"{kind} {name}: {dt:.1f}s")


def run_step(name, argv, logfh):
    """Run one pipeline step as a subprocess to x_TEA_main.py, in this process's env.

    cwd is the package dir so x_TEA_main.py's flat imports resolve. Fails loud.
    """
    cmd = [sys.executable, XTEA_MAIN] + argv
    log(f"step: {name}")
    logfh.write(f"\n=== STEP {name} ===\n{' '.join(cmd)}\n")
    logfh.flush()
    t0 = time.monotonic()
    r = subprocess.run(cmd, cwd=str(HERE), stdout=logfh, stderr=subprocess.STDOUT)
    dt = time.monotonic() - t0
    if r.returncode != 0:
        record_timing("step_FAILED", name, dt)
        sys.exit(f"ERROR: step '{name}' failed (exit {r.returncode}); see log.")
    record_timing("step", name, dt)
    log(f"step {name}: {dt:.1f}s (children peak RSS {_children_peak_gb():.2f} GB)")


def hc_pf_path(type_dir):
    """Path to a type's postF_hc raw-result file (round-2 post-filter output).
    Shared by build_type_steps (which writes it) and run_type_downstream /
    run_shared_tail (which read it back for the Phase C move A.2 locus merge),
    so the two can never disagree on the filename."""
    return os.path.join(str(type_dir), "candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt")


def cc_hc_path(type_dir):
    """Path to a type's case-control high-confident somatic output (step 3 of
    master's case-control chain, x_somatic_calling.CaseControlMode.
    parse_high_confident_somatic). Filename mirrors master's own naming
    (confirmed against a real executed run_xTEA_pipeline.sh). Plays the same
    conceptual role hc_pf_path plays for the germline path -- the final per-type
    raw-result file move A.2's locus merge reads back -- when --case-control is
    on."""
    return os.path.join(str(type_dir),
                         "candidate_disc_filtered_cns_high_confident_post_filtering_somatic.txt")


# The per-type chain's recognised TERMINAL steps -- one per mode, each the step whose output
# the Phase C move A.2 locus merge reads back: postF_hc (germline), case_control_hc
# (--case-control), postF_mosaic (--mosaic). A wave that stops anywhere else is an
# intermediate wave with no artifact to hand back, and must return None rather than warn
# about a VCF nobody asked it to produce.
#
# Named once, here, because it is consulted in two places that must never disagree: the guard
# in run_type_downstream and the mode dispatch in main(). It was previously an inline tuple,
# which meant adding a mode edited a literal that a test was pinning by string.
TERMINAL_TYPE_STEPS = ("case_control_hc", "postF_hc", "postF_mosaic")


def mosaic_path(type_dir):
    """Path to a type's mosaic post-filter output (MOSAIC_PLAN M.2).

    Plays the same role for --mosaic that cc_hc_path plays for --case-control and
    hc_pf_path plays for germline: the final per-type raw-result file move A.2's locus
    merge reads back. Named after master's own MosaicCaller output rather than invented,
    so a mosaic tree is recognisable next to a master mosaic run."""
    return os.path.join(str(type_dir),
                        "candidate_disc_filtered_cns_post_filtering_mosaic.txt")


def build_type_steps(t, cfg, paths, dirs, opt,
                     n_cores=None, stop_after=None, start_from=None,
                     n_cores_transduction=None):
    """Return the ordered list of (step_name, argv) for one MEI type.

    Mirrors the golden `run_xTEA_pipeline.sh` step chain exactly, parameterised,
    with one deliberate departure: the chain **starts at cns_filter (-N)**. Trim
    step 4 committed xTea2 to the scan-only path, so candidate_list_from_clip.txt
    and candidate_list_from_disc.txt are always produced natively by the
    xtea-scan pass + disc adapter in the caller, never by x_TEA_main's -C/-D
    steps (which no longer exist). The former skip_clip/skip_disc parameters are
    gone with them -- they had no False caller left.

    n_cores overrides the per-step worker count (`-n`). When the three types'
    downstreams run concurrently (--parallel-downstream) each type gets a slice
    of the total (cores // ntypes) so the sum stays within the allocation; None
    means use opt.cores (the serial default, unchanged).

    n_cores_transduction overrides `-n` for the transduction step ALONE, and
    exists because the even cores//ntypes split is wrong for that one step: Alu
    has no transduction flank at all (REP_TYPES["Alu"]["flank"] is None -> the
    --transduction branch in x_TEA_main takes the empty path), so Alu's
    transduction is a ~1.5s no-op while L1's and SVA's are ~25 and ~30 minutes.
    Splitting cores three ways there hands a third of the machine to a step that
    does nothing. Which types have a flank is known statically from REP_TYPES, so
    this needs no runtime redistribution -- see the caller for the derivation.
    None means fall back to n_cores (unchanged behaviour).

    stop_after truncates the returned list right after the named step (Phase C
    move A.2 uses stop_after="postF_hc": gene/genotype/gVCF run ONCE over the
    unified locus set instead of once per type -- see run_shared_tail).

    start_from truncates from the front, dropping every step BEFORE the named
    one. The caller is then responsible for having staged that step's inputs.
    Used by tests/test_downstream_from_cns.py to run the post-cns chain
    (transduction onward) against pre-computed cns outputs, which turns an
    ~18-minute chr1 gate into ~2 minutes -- the argv is built by this same
    function either way, so the fast harness cannot drift from the real chain.
    """
    P = dirs["type"] + "/"          # per-type output dir (trailing slash)
    bam_list = dirs["bam_list"]
    ref = cfg["reference"]
    rtype = str(paths["rtype"])
    fflank = paths["flank"]
    nc = str(opt.cores if n_cores is None else n_cores)
    nc_td = nc if n_cores_transduction is None else str(n_cores_transduction)
    steps = []
    # SVA needs --sva on the clip/disc/cns steps (it selects the SVA family bundle ->
    # annotation_extend=200 for the mate-in-repeat test), exactly as the golden
    # run_xTEA_pipeline.sh does. Omitting it uses extend 0, which silently drops
    # ~18% of SVA at whole-genome scale (chr1 happened to tolerate it). The golden
    # puts --sva on -C/-D/-N only (not transduction/sibling/postF), so we match that.
    #
    # That omission has a consequence worth knowing before you "fix" it: transduction
    # and sibling therefore resolve family_params.GENERIC even for SVA, so SVA's
    # leftmost_cns=5000 has never executed. See xtea/family_params.py -- waking it moves
    # calls in bulk and belongs to Stage 3, not here.
    sva = ["--sva"] if t == "SVA" else []
    # CC.0: master threads --purity <p> --tumor through clip/disc/cns_filter/
    # transduction/sibling, and --tumor alone (no --purity) through both postF
    # steps (confirmed against a real executed run_xTEA_pipeline.sh, not just
    # gnrt_pipeline_local_v38.py source). Empty lists when opt.tumor is False,
    # so this is a no-op for every existing (non-tumor) caller.
    tumor_purity_args = (["--purity", str(opt.purity), "--tumor"] if opt.tumor else [])
    tumor_only_args = (["--tumor"] if opt.tumor else [])
    # THE GERMLINE CUTOFF THAT CAPPED MOSAIC SENSITIVITY LONG BEFORE postF_mosaic RAN.
    #
    # postF_all is the germline post-filter, and without --tumor its four ratio cutoffs are all
    # 0.075 (x_post_filter.get_cutoff_by_type); with --tumor they drop to 0.001. The binding one
    # for mosaic is f_ef_disc > 0.075 -- the fraction of a site's discordant reads that support
    # the insertion, which is inherently LOW at mosaic allele fractions. Measured on broad's real
    # WGS run against the HapMap truth set, of the truth-matching records reaching postF_all,
    # 92% of Alu and 51% of L1 fail that one gate; f_clip_full_map > 0.075 (a >7.5% allele
    # fraction, against a truth-set VAF ladder that tops out at 8.25%) takes another 36%/17%.
    # Only 82 of 1033 Alu and 72 of 158 L1 truth records survive at 0.075; ALL of them survive
    # at 0.001. postF_mosaic never got to see the rest, so no amount of tuning THERE could
    # recover them -- the ceiling was set one step earlier.
    #
    # Measured end-to-end on broad (postF_all --tumor, then the unchanged mosaic filter):
    #     L1   12.1% -> 19.4% sensitivity, precision 82.3% -> 88.2% (FP unchanged at 11)
    #     Alu   3.6% -> 33.1% sensitivity, precision 61.4% -> 74.1%
    # Precision RISES in both families: postF_all was deleting true positives indiscriminately,
    # and the mosaic filter downstream is what actually discriminates. Alu's VAF curve becomes a
    # real detection curve (9.7% / 38.6% / 55.5% / 73.3% up the ladder) instead of flat noise.
    #
    # The relaxation also stops postF_all discarding one-side orphan/transduction records
    # (x_post_filter :249, :260, :103); that is the same "do not assume germline-strength
    # evidence" intent, so it is wanted here too. Applied to postF_all ONLY, not postF_hc: in
    # mosaic mode the per-type chain terminates at postF_mosaic, which reads pf_out, and
    # run_type_downstream returns the postF_mosaic path as the shared tail's input -- postF_hc's
    # output is never read, so widening it would change nothing and only enlarge the blast radius.
    #
    # --somatic, not --tumor: the flag means "this sample is not germline, do not assume germline
    # allele fractions", which is equally true of a tumor and of a mosaic run. Spelling it
    # --tumor on a mosaic command line would misdescribe the run in every log and step argv.
    # x_TEA_main accepts --tumor as an alias, so the tumor/case-control path is unaffected.
    postF_all_somatic_args = (["--somatic"] if (opt.tumor or getattr(opt, "mosaic", False)) else [])

    disc_out = f"{P}candidate_list_from_disc.txt"
    raw_disc = f"{disc_out}.clip_sites_raw_disc.txt"
    cns_out  = f"{P}candidate_disc_filtered_cns.txt"
    cns2_out = f"{P}candidate_disc_filtered_cns2.txt"
    sib_out  = f"{P}candidate_sibling_transduction2.txt"
    pf_out   = f"{P}candidate_disc_filtered_cns_post_filtering.txt"
    hc_pf    = hc_pf_path(dirs["type"])
    gene_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"
    gntp_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"

    # Steps 1 (clip, -C) and 2 (disc, -D) are gone: trim step 4 deleted master's
    # pysam clip/disc implementation. candidate_list_from_clip.txt and disc_out
    # are still the same paths, now written by
    # scan_adapter.write_clip_candidate_list and disc_adapter.run_disc_from_scan
    # in the caller, before this chain starts.
    # THE FLAG WITHOUT WHICH --mosaic's CUTOFFS DO NOT REACH THE STEPS THAT MATTER.
    #
    # x_TEA_main sets b_automatic = not options.user_specific (x_TEA_main.py:313-315), and
    # every step that takes --cr/--nd then OVERWRITES them from the coverage table when
    # b_automatic is true (:425-427 cns_filter, :510 transduction, :608/:661 sibling). So
    # passing --cr 2 --nd 2 is silently ignored, and at 457x the steps use 14 and 20.
    #
    # This was not theoretical: run 50427774 logged "Average coverage is: 457.245 / Filter
    # (on cns) cutoff: 14 and 20 are used!!!" after --mosaic had already reported forcing
    # 2/2. The disc-adapter override worked (it is applied in-process, so 21,731 L1 and
    # 60,052 Alu candidates survived), and then cns_filter threw that away at 14/20 --
    # which no 1-8% VAF insertion can clear. The run would have finished clean and empty.
    #
    # --user is master's own opt-out ("User specific parameters, by default automatically
    # calc the parameters"), so this asks master for exactly what it already supports.
    user_specific = ["--user"] if getattr(opt, "mosaic", False) else []

    # 3. cns filter
    steps.append(("cns_filter", user_specific + [
        "-N"] + sva + ["--cr", str(opt.nfclip), "--nd", str(opt.nfdisc), "-b", bam_list,
        "-p", dirs["tmp_cns"], "--fflank", fflank, "--flklen", str(opt.flklen),
        "-n", nc, "-i", disc_out, "-r", paths["cns"], "--ref", ref,
        "-a", paths["annotation"], "-o", cns_out] + tumor_purity_args))
    # 4. transduction
    steps.append(("transduction", user_specific + [
        "--transduction", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", nc_td, "-i", cns_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "-o", cns2_out] + tumor_purity_args))
    # 5. sibling
    steps.append(("sibling", user_specific + [
        "--sibling", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", nc, "-i", cns2_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "--blacklist", cfg["blacklist"],
        "-o", sib_out] + tumor_purity_args))
    # 6. post-filter (round 1: all)
    steps.append(("postF_all", [
        "--postF", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", nc,
        "-i", cns2_out, "-a", paths["annotation1"], "-o", pf_out] + postF_all_somatic_args))
    # 7. post-filter (round 2: high-confident, with blacklist)
    steps.append(("postF_hc", [
        "--postF", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", nc,
        "-i", f"{cns2_out}.high_confident", "-a", paths["annotation1"],
        "--blacklist", cfg["blacklist"], "-o", hc_pf] + tumor_only_args))
    # 8. case-control filtering (CC.1, somatic mode only -- gated on --case-control).
    # Two steps mirroring master's x_somatic_calling.CaseControlMode exactly (arg
    # shape confirmed against a real executed run_xTEA_pipeline.sh, LFS 2921/L1):
    # step (a) filters the RAW postF_all list (pf_out, not hc_pf) against the
    # control BAM; step (b) intersects that with the type's own hc_pf list for
    # the final high-confident somatic set. --cr/--nd here reuse nfclip/nfdisc
    # (same cutoffs cns_filter/postF use); no --tumor/--purity (the tumor's own
    # candidates are already built -- this only re-examines the control BAM at
    # those loci) and no --rtype/--blacklist (case_control mode takes neither).
    cc_raw_out = f"{P}candidate_disc_filtered_cns_post_filtering_somatic.txt"
    cc_hc_out = cc_hc_path(dirs["type"])
    if getattr(opt, "case_control", False):
        control_bam_list = dirs["control_bam_list"]
        steps.append(("case_control", [
            "--case_control", "-p", dirs["tmp_cns"], "-b", control_bam_list, "-n", nc,
            "-r", paths["cns"], "-i", pf_out, "--ref", ref,
            "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
            "--fflank", fflank, "--flklen", str(opt.flklen), "-o", cc_raw_out]))
        steps.append(("case_control_hc", [
            "--case_control", "--somatic_hc", "-p", dirs["tmp_cns"], "-i", hc_pf,
            "--input2", cc_raw_out, "-n", nc,
            "--fflank", fflank, "--flklen", str(opt.flklen), "-o", cc_hc_out]))
    # 8b. mosaic post-filter (MOSAIC_PLAN M.2, --mosaic only). One step, not two: unlike
    # case-control this re-reads no BAM -- MosaicCaller filters the already-loaded candidate
    # records and a precomputed repeat interval tree, so it does not scale with coverage.
    #
    # It filters pf_out (postF_all, the RAW list) and NOT hc_pf, for exactly the reason
    # case_control step (a) does: hc_pf has already been through the blacklist and the
    # high-confidence cut, and a 1-5% VAF insertion is precisely what those drop. Feeding the
    # hc list here would produce a mosaic caller that only ever confirms germline-strength
    # calls -- plausible output, no mosaic sensitivity, and nothing in the run would say so.
    mosaic_out = mosaic_path(dirs["type"])
    if getattr(opt, "mosaic", False):
        steps.append(("postF_mosaic", [
            "--postF", "--postFmosaic", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", nc,
            "-i", pf_out, "-a", paths["annotation1"],
            "--blacklist", cfg["blacklist"], "-o", mosaic_out]))
    # 9. gene annotation
    steps.append(("gene", [
        "--gene", "-a", cfg["gene_gff3"], "-i", hc_pf, "-n", nc,
        "-o", gene_out]))
    # 10. genotype (x_genotype_lik; genotype not part of the concordance gate). On this
    # per-type path the EM is fitted over this one family's calls only; the default
    # shared tail (build_shared_tail_steps) fits once over the sample's merged call set.
    steps.append(("genotype", [
        "--gntp_classify", "-i", gene_out, "-o", gntp_out]))
    # 11. gVCF
    steps.append(("gvcf", [
        "--gVCF", "-i", gntp_out, "-o", P, "-b", bam_list, "--ref", ref,
        "--rtype", rtype]))

    if stop_after is not None:
        names = [name for name, _ in steps]
        steps = steps[:names.index(stop_after) + 1]
    if start_from is not None:
        names = [name for name, _ in steps]
        steps = steps[names.index(start_from):]
    return steps


def build_shared_tail_steps(cfg, merged_dir, opt, sf_merged_hcpf, sf_type_map, sf_merged_vcf,
                            bam_list, rtypes_bitmask):
    """gene -> genotype -> gVCF, run ONCE over the unified locus set (Phase C move
    A.2), instead of once per type. Gene annotation and genotyping are already
    type-agnostic (no --rtype dependency); gVCF gets --type-map (+ --merged-vcf-out)
    instead of a single --rtype so it can fan each record out to its own resolved
    type (x_gvcf.cvt_raw_rslt_to_gvcf_multi) while still writing the same per-type
    VCF files the concordance gate expects. rtypes_bitmask is the OR of every
    requested type's rtype code (e.g. 7 = L1|Alu|SVA); x_TEA_main.py decomposes it
    back into the set of per-type VCFs to open.

    CC.1: when opt.case_control is set, the genotype step is skipped -- gene's
    output feeds --gVCF directly. Confirmed against master's own real
    case-control run (LFS 2921/L1): its somatic chain never calls --gntp_classify
    (gene -> gVCF only) and the resulting VCF shows GT=./. for every somatic
    call; x_gvcf.py defaults to "./." when the genotype column isn't present
    (l_fields[50] absent), so this is a graceful, already-supported skip, not a
    new code path.
    """
    P = merged_dir + "/"
    ref = cfg["reference"]
    nc = str(opt.cores)
    gene_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"
    gntp_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"
    steps = [
        ("gene", ["--gene", "-a", cfg["gene_gff3"], "-i", sf_merged_hcpf, "-n", nc, "-o", gene_out]),
    ]
    gvcf_in = gene_out
    # Mosaic skips the genotyper for the same reason case-control does, and the reason is
    # worth stating: x_genotype_lik is a GERMLINE het/hom genotyper (two genotypes, biases
    # fitted by EM over the sample's calls). A 2% VAF insertion is neither het nor hom, so it
    # would assign every mosaic call a germline genotype, and a call set made of mosaic
    # calls would also corrupt the EM fit itself. x_gvcf.py already
    # defaults to "./." when the genotype column is absent, so this is a supported skip
    # rather than a new code path.
    if not (getattr(opt, "case_control", False) or getattr(opt, "mosaic", False)):
        # Genotyping runs HERE, once over the merged all-family locus set, because the
        # genotyper's EM needs all of the sample's calls to fit its per-sample biases.
        steps.append(("genotype", ["--gntp_classify", "-i", gene_out, "-o", gntp_out]))
        gvcf_in = gntp_out
    steps.append(("gvcf", ["--gVCF", "-i", gvcf_in, "-o", P, "-b", bam_list, "--ref", ref,
                           "--type-map", sf_type_map, "--merged-vcf-out", sf_merged_vcf,
                           "--rtype", str(rtypes_bitmask)]))
    return steps


def run_type_downstream(ctx, cfg, opt, n_cores, logfh, stop_after=None,
                        n_cores_transduction=None, start_from=None):
    """Run one MEI type's downstream step chain (serial *within* the type).

    The chain (cns_filter -> transduction -> sibling -> postF x2 -> gene ->
    genotype -> gVCF) must run in order, but the three types' chains are mutually
    independent (separate output dirs; the shared scan+cns-fan already ran), so
    the caller may run several of these concurrently. `logfh` receives every
    step's stdout/stderr.

    stop_after="postF_hc" (Phase C move A.2): stop before gene/genotype/gVCF --
    those run ONCE over the unified locus set instead (see run_shared_tail) --
    and return the postF_hc path instead of a VCF path. stop_after="case_control_hc"
    (CC.1, when --case-control is on): stop after the two case-control filtering
    steps instead, and return that final somatic-HC path -- it plays hc_pf's role
    in the move-A.2 locus merge, so the merge/gene/gVCF tail is unaffected by
    whether the input was germline- or somatic-filtered. stop_after=None (default)
    runs the full chain and returns the produced VCF path, unchanged from before.

    start_from (Stage 2c.3) drops every step BEFORE the named one, so the caller
    can split this chain into waves with a cross-family barrier between them --
    see the family-assignment barrier in main(). It is handed straight to
    build_type_steps, which is also what tests/test_downstream_from_cns.py uses,
    so the argv cannot drift between the real chain and either harness.

    A wave that stops somewhere other than a recognised terminal step (i.e. an
    intermediate wave) returns None: there is no artifact to hand back yet, and
    warning about a "missing" VCF that was never meant to exist would be noise.
    """
    (t, spec, paths, dirs, type_dir) = ctx
    log(f"=== type {t} (rtype {spec['rtype']}): downstream (n={n_cores}"
        + (f", transduction n={n_cores_transduction}" if n_cores_transduction is not None else "")
        + (f", from {start_from}" if start_from else "")
        + (f", through {stop_after}" if stop_after else "")
        + ") ===")
    for name, argv in build_type_steps(t, cfg, paths, dirs, opt,
                                       n_cores=n_cores, stop_after=stop_after,
                                       n_cores_transduction=n_cores_transduction,
                                       start_from=start_from):
        run_step(f"{t}:{name}", argv, logfh)
    if stop_after is not None and stop_after not in TERMINAL_TYPE_STEPS:
        return None
    if stop_after == "postF_mosaic":
        mos = mosaic_path(dirs["type"])
        if os.path.exists(mos):
            log(f"{t} postF_mosaic: {mos}")
            return mos
        log(f"WARNING: expected postF_mosaic not found for {t}: {mos}")
        return None
    if stop_after == "case_control_hc":
        cc_hc = cc_hc_path(dirs["type"])
        if os.path.exists(cc_hc):
            log(f"{t} case_control_hc: {cc_hc}")
            return cc_hc
        log(f"WARNING: expected case_control_hc not found for {t}: {cc_hc}")
        return None
    if stop_after == "postF_hc":
        hc_pf = hc_pf_path(dirs["type"])
        if os.path.exists(hc_pf):
            log(f"{t} postF_hc: {hc_pf}")
            return hc_pf
        log(f"WARNING: expected postF_hc not found for {t}: {hc_pf}")
        return None
    vcf = type_dir / f"{opt.sample}_{spec['vcf_alt']}.vcf"
    if vcf.exists():
        log(f"{t} VCF: {vcf}")
        return str(vcf)
    log(f"WARNING: expected VCF not found for {t}: {vcf}")
    return None


def downstream_waves(assign_early, type_stop_after, n_types, locus_universe=False):
    """The per-family downstream, as a list of (start_from, stop_after) waves.

    One wave is the pre-2c chain. Two waves open a cross-family barrier between
    cns_filter and transduction, which is where anything needing every family's
    cns output and no family's transduction output has to run. Two callers want
    that point, and they are independent:

      - `--family-assign-early` (Stage 2c.3) resolves family there;
      - `--locus-universe` enumerates the loop's iteration domain there, because
        the site-creation rule reads the cns candidate POSITIONS (see
        downstream_unified's module docstring for why that pins the placement).

    n_types == 1 keeps one wave for the assignment barrier -- a lone family has
    nothing to contend with, and `--types L1` on its own is the arm the Stage 0
    gate diffs against to attribute the 12 Alu-over-L1 conflicts, so it must not
    change. The universe does NOT get that exemption: one family still has a
    locus universe, and enumerating it is answer-preserving either way.
    """
    if locus_universe or (assign_early and n_types > 1):
        return [(None, "cns_filter"), ("transduction", type_stop_after)]
    return [(None, type_stop_after)]


class StagedBarrier:
    """The cns_filter -> transduction barrier, decomposed into what each family waits for.

    WHY THIS EXISTS. `downstream_waves` splits the per-family chain into two waves with a
    join between them, and a join makes every family wait for the SLOWEST family's
    cns_filter. At WGS that is Alu's, which is ~2x the others' (124.1s against L1's 57.3
    and SVA's 62.5), and it costs **+52.9s measured** (job 49550925). None of that is work
    -- in the unsplit chain L1 and SVA enter transduction as soon as THEIR OWN cns_filter
    returns and overlap Alu's remaining ~60s, and a join throws that overlap away. It also
    retro-explains the Stage 2c gate's unattributed +40s (job 49502115): --family-assign-early
    opens the same split and paid the same price.

    THE JOIN WAS NEVER REQUIRED, because the three things that run at this point have three
    different dependency sets, and only one of them needs every family:

      - **per-family enumeration** (`--locus-universe`) reads that family's own cns list and
        its own raw_disc. It needs NOBODY. It runs the moment that family's cns_filter lands.
      - **the driver gather** (`--td-gather-driver`) reads the FLANKED families' `.slct`
        files, and writes into the flanked families' working folders. It needs L1 and SVA;
        it does not need Alu, and Alu does not need it -- a flankless family's transduction
        takes x_TEA_main's copy-through branch (`os.path.isfile(sf_flank)` is False) and
        never calls `call_candidate_transduction_v3` at all, so it never reads
        XTEA_TD_GATHER_PRECOMPUTED and cannot be raced by the env var being set late.
      - **the union table** and **`--family-assign-early`** need all three. The union table
        still blocks nobody: nothing in the run reads it. Family assignment genuinely is a
        join, so with that flag on this degenerates to the old behaviour -- correctly.

    So the waits are: a flanked family waits for the gather; every family waits for the
    assignment if it is on; and a run with only `--locus-universe` waits for NOTHING, which
    is what makes flipping that flag on cost its ~3s of enumeration and not a barrier.

    WHO RUNS THE SHARED WORK. Whichever thread completes the dependency, rather than a
    fourth thread or the main one: the gather goes to the last FLANKED family to enumerate
    (it would otherwise be sitting in `wait()`), the union table and the assignment to the
    last family overall. Election is under the lock and latched, so each runs exactly once.

    ORDERING IS PRESERVED, and it is the one that matters: the assignment REWRITES each
    family's cns list to its own winners, so enumerating after it would silently measure a
    post-assignment universe. Here every family enumerates before it registers, and the
    assignment cannot be elected until every family has registered, so all three
    enumerations precede it -- the same guarantee the two-wave order gave, without the join.

    FAILURE. A join has one thing going for it: nobody waits on a family that died. Here
    they do, so every path that can fail must call `abort`, which latches the failure and
    releases both events -- otherwise a cns_filter that exits non-zero leaves the other
    families blocked forever inside the executor's shutdown and the run hangs instead of
    failing. `sys.exit` in a worker thread is a SystemExit, hence BaseException.
    """

    def __init__(self, families, flanked, b_universe, b_gather, b_assign,
                 enumerate_fn=None, gather_fn=None, union_fn=None, assign_fn=None):
        # Kept as a LIST as well as a set, and the list is not decoration: the union
        # table's `families`/`claims` columns are written in the iteration order of the
        # map handed to `build_locus_universe`, so handing it completion order makes the
        # table's bytes depend on which family's cns_filter finished first. Job 49596361
        # failed exactly there -- same 3,272 loci, same claims, `SVA,L1` where the global
        # arm wrote `L1,SVA`. The canonical order is the driver's family order.
        self._family_order = list(families)
        self._families = set(families)
        self._flanked = set(flanked) & self._families
        self._b_universe = b_universe
        self._b_gather = b_gather and bool(self._flanked)
        self._b_assign = b_assign
        self._enumerate_fn = enumerate_fn
        self._gather_fn = gather_fn
        self._union_fn = union_fn
        self._assign_fn = assign_fn
        self._lock = threading.Lock()
        self._enumerated = set()
        self._slct = {}
        self._gather_claimed = False
        self._all_claimed = False
        self._gather_done = threading.Event()
        self._all_done = threading.Event()
        self._failed = None

    def arrive(self, t):
        """Called by family `t`'s worker between its cns_filter and its transduction."""
        try:
            if self._b_universe and self._enumerate_fn is not None:
                _n_cand, _n_new, sf_out = self._enumerate_fn(t)
                with self._lock:
                    self._slct[t] = sf_out
            with self._lock:
                self._enumerated.add(t)
                run_gather = (self._b_gather and not self._gather_claimed
                              and self._flanked <= self._enumerated)
                self._gather_claimed = self._gather_claimed or run_gather
                run_all = (not self._all_claimed
                           and self._families <= self._enumerated)
                self._all_claimed = self._all_claimed or run_all
                # Canonical order, NOT completion order -- see __init__.
                slct = {f: self._slct[f] for f in self._family_order if f in self._slct}
            # The gather first when one thread draws both: it is the only one of these
            # anybody is blocked on.
            if run_gather:
                self._gather_fn()
                self._gather_done.set()
            if run_all:
                if self._b_universe and self._union_fn is not None:
                    self._union_fn(slct)
                if self._b_assign and self._assign_fn is not None:
                    self._assign_fn()
                self._all_done.set()
        except BaseException:
            self.abort(t)
            raise
        if self._b_gather and t in self._flanked:
            self._gather_done.wait()
            self._check()
        if self._b_assign:
            self._all_done.wait()
            self._check()

    def abort(self, t):
        """Latch a failure and release every waiter, so a dead family fails the run."""
        with self._lock:
            if self._failed is None:
                self._failed = t
        self._gather_done.set()
        self._all_done.set()

    def _check(self):
        if self._failed is not None:
            raise RuntimeError(
                f"downstream barrier: {self._failed} failed before the barrier; "
                "aborting the families waiting on it")


def run_family_assign_barrier(type_ctx, sample_root):
    """Stage 2c.3: resolve family ONCE, between cns_filter and transduction.

    WHY HERE, and not at the candidate list, which is what "move the decision up
    front" sounds like. Two facts measured on BASELINE-0 (plan sec 4 Stage 2c):

    - `scan_clip_realign.tsv` carries no family signal to decide on. Its nL/nR
      are per-family by construction (one bwa per family in realign.hpp) but
      come out identical: of the 155,615 positions with >=5 clip support, only
      504 differ across families at all. The three families' clip candidate
      lists are the same 577K positions.
    - `two_side_tprt_both` -- which decides 13 of the 22 conflicts the merge
      resolves today -- is written BY cns_filter (x_clip_disc_filter.py:1845).
      Anything earlier decides without it.

    So cns_filter's output is the first file where the existing rule can run at
    all, and running it there is a strictly smaller change than inventing a new
    rule for an earlier one. Each family then carries only its own candidates
    through transduction/sibling/postF instead of all three carrying overlaps.

    Do not expect throughput from this. The families barely overlap: 672 of
    12,895 disc candidates, and 16 of 2,796 cns rows. The reason to do it is
    that family stops being decided three chains later by whichever list a call
    happened to survive on.

    WHAT IT ACTUALLY MOVES, measured (job 49502115, HG002 WGS): NOTHING. The
    call set is identical with and without it. That is not a null result about
    the rule, it is a fact about the topology -- this barrier does not REPLACE
    the postF_hc merge in run_shared_tail, both run. So it can only change
    output where it drops a cns row that would have WON at the merge, and it
    never does: of the 16 conflicts it resolves, the 7 the merge also sees give
    the identical winner, loser and rule, and the other 9 the merge never sees.
    The 15 the merge still resolves have no cns row for the losing family at
    all, because transduction and sibling create that record after this runs.
    Collapsing to ONE decision point therefore needs unified site creation,
    which is Stage 3 -- see plan sec 4 Stage 2c.
    """
    cns_paths = {t: str(Path(dirs["type"]) / "candidate_disc_filtered_cns.txt")
                 for (t, spec, paths, dirs, type_dir) in type_ctx}
    sf_conflicts = str(sample_root / "cns_family_conflicts.tsv")
    log(f"=== family assignment barrier (cns_filter -> transduction, "
        f"win={locus_table.MERGE_WINDOW}) ===")
    stats, table = locus_table.assign_and_filter(cns_paths, sf_conflicts=sf_conflicts)
    for t in sorted(stats):
        n_in, n_kept = stats[t]
        log(f"    {t:4s} {n_in} cns candidates -> {n_kept} kept "
            f"({n_in - n_kept} lost to another family)")
    summary = table.conflict_summary()
    log(f"family conflicts resolved at cns: {len(table.conflicts)} -> {sf_conflicts}")
    for (win, lose, rule) in sorted(summary):
        log(f"    {win:4s} kept over {lose:4s} by {rule:14s} {summary[(win, lose, rule)]}")
    return stats, table


def run_locus_universe(type_ctx, cfg, opt, sample_root, disc_cutoff, cns_extnd):
    """Enumerate the per-locus loop's iteration domain, driver-side (plan sec 3).

    Runs at the same point as the family-assignment barrier and for the same
    structural reason: every family's cns_filter is finished and none is past it.
    That is forced, not chosen -- the site-creation rule takes the cns candidate
    POSITIONS for its proximity exclusion, so the universe is not computable
    before cns_filter, only before transduction.

    DIAGNOSTIC ONLY, and deliberately so. Each family's transduction step still
    derives its own `<raw_disc>.slct` through the same one copy of the rule
    (site_universe.write_slct_file), so these files are a replica written to
    their own directory and no call can move. What that buys is a gate that is a
    byte comparison instead of an argument: after a run,

        cmp <type>/candidate_list_from_disc.txt.clip_sites_raw_disc.txt.slct \\
            merged/locus_universe/<type>.slct

    must be identical for every flanked family. Alu has no `.slct` from the step
    at all (no transduction flank), so a flanked-family count is the completeness
    check -- expect n_flanked comparisons, not n_types.

    THE ARGUMENTS ARE THE POINT. Everything the rule needs is derived here from
    the driver's own state rather than from a step log (which is what
    docs/handoff/runs/gate_site_universe.py had to do offline):

      - the disc cutoff is the driver's `disc_cutoff`, already
        `Parameters.get_par_by_cov(purity-scaled cov)[1]`. The step ignores
        `--nd` because `b_automatic` is True unless `--user_specific`, so
        reading nfdisc here would be wrong by one at HG002 (4, not 5);
      - the annotation is `annotation1` (LINE/{g}/{g}_L1.fa.out for L1), which
        is what build_type_steps hands transduction -- NOT `annotation`, which
        is cns_filter's;
      - the rmsk boundary extend is the step's own 450, asserted rather than
        retyped, and the min copy length its own 225.

    Returns (stats, universe) with stats {family: (n_candidates, n_new)}.

    THIS IS THE ALL-AT-ONCE FORM, and it is only what the serial path and the
    tests want. Nothing here needs every family at once except the union table:
    each family's `.slct` is derived from ITS OWN cns list and raw_disc and reads
    no other family's anything. `enumerate_family_universe` is the per-family
    half, and `StagedBarrier` is what calls it as each family's cns_filter lands
    instead of after the slowest one -- see that class for why that is worth ~53s.
    """
    log(f"=== locus universe (cns_filter -> transduction, driver-side) ===")
    stats = {}
    slct_by_family = {}
    for ctx in type_ctx:
        n_cand, n_new, sf_out = enumerate_family_universe(
            ctx, cfg, opt, sample_root, disc_cutoff)
        stats[ctx[0]] = (n_cand, n_new)
        slct_by_family[ctx[0]] = sf_out
    universe = write_locus_universe_table(slct_by_family, sample_root, cns_extnd)
    return stats, universe


def enumerate_family_universe(ctx, cfg, opt, sample_root, disc_cutoff):
    """One family's `.slct` replica. Returns (n_candidates, n_new, path).

    Depends on nothing but this family's own `cns_filter` output, so it may run
    the moment that step returns rather than after every family's has.

    Timed as a STEP, not a phase, because under the staged barrier three of these
    overlap each other and one another's `cns_filter`s -- and the timing file's
    contract is that phase rows do not overlap and step rows may (see the
    instrumentation block's header).
    """
    (t, spec, paths, dirs, type_dir) = ctx
    P = dirs["type"] + "/"
    raw_disc = f"{P}candidate_list_from_disc.txt.clip_sites_raw_disc.txt"
    cns_out = f"{P}candidate_disc_filtered_cns.txt"
    sf_out = downstream_unified.family_slct_path(sample_root, t)
    b_flank = bool(spec.get("flank"))
    t0 = time.monotonic()
    n_cand, n_new = downstream_unified.enumerate_family_sites(
        raw_disc, cns_out, paths["annotation1"], spec["rtype"], disc_cutoff,
        opt.tumor, sf_out, b_flank=b_flank, n_jobs=opt.cores,
        sf_reference=cfg["reference"], working_folder=dirs["tmp_td"],
        log=lambda _m: None)   # the rule's per-site rejection log is the step's to print
    record_timing("step", f"{t}:locus_universe", time.monotonic() - t0)
    log(f"    {t:4s} {n_cand} cns candidates + {n_new} created "
        f"= {n_cand + n_new} sites"
        f"{'' if b_flank else '  [flankless: cns list only, no step .slct to compare]'}")
    return n_cand, n_new, sf_out


def write_locus_universe_table(slct_by_family, sample_root, cns_extnd):
    """Union every family's `.slct` into the loop's iteration domain.

    THE ONE PART THAT NEEDS ALL THREE, and the reason it is still not a barrier:
    nothing in the run reads this table. It is the loop's index and a diff target,
    so it is written by whichever family enumerates last and blocks nobody.
    """
    universe = downstream_unified.build_locus_universe(slct_by_family)
    sf_universe = str(Path(downstream_unified.universe_dir(sample_root)) /
                      "locus_universe.tsv")
    universe.write(sf_universe)
    by_n = universe.by_n_families()
    log(f"locus universe: {universe.n_loci} loci -> {sf_universe}")
    log("    claimed by " + ", ".join(
        f"{n} family: {by_n[n]}" if n == 1 else f"{n} families: {by_n[n]}"
        for n in sorted(by_n)))
    # The barrier's own ceiling, measured on the universe instead of on the final
    # call set: a locus no family has a cns row for cannot be resolved by any
    # decision point placed before transduction. See plan sec 4's 7-of-22.
    log(f"    reachable only by creation (no family has a cns row): "
        f"{universe.n_created_only()}")
    log(f"    transduction gather extnd would be "
        f"{downstream_unified.transduction_extnd(cns_extnd) if cns_extnd else 'unknown'}"
        f" (max(400, cns_extnd))")
    return universe


def run_td_gather_driver(type_ctx, cfg, opt, sample_root, cns_extnd, logfh,
                         timing_kind="phase"):
    """Gather every family's transduction reads out of ONE BAM decode, here (plan sec 4 2b).

    `xtea_scan --td-manifest` has been gated since 2026-08-02 (job 49516920, 6/6
    byte-identical at WGS) with no production caller, because hoisting the gather meant
    re-implementing the step's site selection. It does not any more: --locus-universe
    already wrote each family's site list at this barrier, and `transduction_extnd` is the
    step's `iextnd`. So this needs nothing the driver does not have.

    It writes to the exact two paths the step would have written in its own working folder
    and sets XTEA_TD_GATHER_PRECOMPUTED, which every transduction step in the next wave
    inherits and honours by skipping its own gather.

    WHAT IT COSTS AND SAVES, measured at WGS (job 49550925) and chr1 (job 49548949). The
    hoist itself is a small WIN; the BARRIER it needs is the loss, and they are separate
    facts that a single wall-clock number hides.

        WGS   step native_gather   51.49s (L1) + 57.22s (SVA) = 108.7s CPU, ~57s of wall
              driver gather        42.3s, serialized      -> ~66s of CPU returned
              SVA transduction     348.1 -> 292.0s        -> -56.0s off the long pole
              net for the hoist    -56.0 + 42.3           =  -13.7s
              the barrier          +52.9s  (see below)
              downstream_parallel  462.7 -> 505.3s        =  +42.6s

    **The +52.9s is not work; it is lost overlap.** Alu's cns_filter is 124.1s against L1's
    57.3 and SVA's 62.5, and in the one-wave chain L1 and SVA start transduction as soon as
    THEIR OWN cns_filter finishes, overlapping Alu's remaining ~60s. Any barrier here
    destroys that. The critical path decomposes exactly: OFF = 62.5 + 348.1 + 49.9 + 2.1 =
    462.6 (measured 462.7); GATHER = 115.5 + 2.7 + 42.3 + 292.0 + 49.5 + 3.1 = 505.1
    (measured 505.3). It also retro-explains the Stage 2c gate's unattributed +40s (job
    49502115): --family-assign-early uses the same barrier.

    Do NOT read this off the total run wall, which went 2263 -> 2228s only because scan_pass
    varied 1567.3 -> 1500.8s -- noise this change cannot touch.

    THE +52.9s IS TAKEN BACK by `StagedBarrier`, which is now the default whenever this
    flag opens the wave split: enumeration runs per family as its cns_filter finishes, this
    gather runs as soon as the FLANKED families' enumerations are done, and only the union
    table and the family-assign barrier ever need all three. `--global-barrier` restores the
    serialized form the numbers above were measured on, which is how they stay comparable.
    Under the staged barrier this runs while a flankless family may still be in cns_filter,
    so it briefly oversubscribes the allocation -- the same trade, and for the same reason,
    as the transduction core split documented in main().

    The reason to do the hoist is not speed. It is that this is the batched evidence stage of
    the per-locus loop, and the only place a single BAM decode can serve every family.

    Flanked families only: a flankless family's transduction step takes the copy-through
    branch and never gathers at all, so an entry for it would decode reads nobody reads.
    """
    scan_bin = str(Path(opt.scan_bin or DEFAULT_SCAN_BIN).resolve())
    extnd = downstream_unified.transduction_extnd(cns_extnd)
    entries = []
    for (t, spec, paths, dirs, type_dir) in type_ctx:
        if not spec.get("flank"):
            continue
        sf_sites = downstream_unified.family_slct_path(sample_root, t)
        if not Path(sf_sites).is_file():
            sys.exit(f"ERROR: --td-gather-driver: no site list for {t} at {sf_sites} "
                     "(it is written by the locus-universe pass, which must run first)")
        n_sites = sum(1 for _ in open(sf_sites))
        if n_sites == 0:
            # An empty manifest entry would produce an empty FASTQ and the step would take
            # it as a real gather that found nothing. Let the step do its own instead.
            log(f"    {t}: 0 sites -- omitted from the manifest; its step will gather")
            continue
        Path(dirs["tmp_td"]).mkdir(parents=True, exist_ok=True)
        clip_fq, disc_fa = downstream_unified.td_gather_paths(dirs["tmp_td"])
        entries.append((t, sf_sites, extnd, clip_fq, disc_fa))
        log(f"    {t}: {n_sites} sites -> {clip_fq}")
    if not entries:
        log("=== td gather (driver): no flanked family with sites; nothing to do ===")
        return 0
    sf_manifest = str(Path(downstream_unified.universe_dir(sample_root)) / "td_manifest.txt")
    downstream_unified.write_td_manifest(entries, sf_manifest)
    # cfg["input"], not bam_list.txt: xtea_scan takes the space-joined path string (the same
    # convention native_td_gather rebuilds by reading column 0 out of bam_list.txt).
    cmd = [scan_bin, "--bam", cfg["input"], "--reference", cfg["reference"],
           "--td-manifest", sf_manifest, "--threads", str(opt.cores)]
    log(f"=== td gather (driver): {len(entries)} families, ONE decode, extnd={extnd} ===")
    log("  " + " ".join(cmd))
    with timed(timing_kind, "td_gather_driver"):
        subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh)
    for (t, _s, _e, clip_fq, disc_fa) in entries:
        log(f"    {t}: clip {Path(clip_fq).stat().st_size} B, "
            f"disc {Path(disc_fa).stat().st_size} B")
    # Inherited by every x_TEA_main transduction subprocess in the next wave.
    os.environ[downstream_unified.ENV_TD_PRECOMPUTED] = "1"
    return len(entries)


def run_shared_tail(cfg, opt, sample_root, hc_pf_paths, bam_list, logfh):
    """Phase C move A.2: resolve the unified locus set from every type's postF_hc
    output (locus_table.assign_and_write), then run gene -> genotype -> gVCF ONCE
    over it instead of once per type. gVCF fans out per-record (via --type-map)
    into the same per-type VCF files the concordance gate expects, plus the
    single deliverable merged VCF at sample_root/{sample}.vcf (same location/
    contract as move A.1's post-hoc merge, now produced directly in this one
    pass instead of a second conflict-resolving merge step).
    """
    merged_dir = str(sample_root / "merged")
    Path(merged_dir).mkdir(parents=True, exist_ok=True)
    sf_merged_hcpf = hc_pf_path(merged_dir)
    sf_type_map = f"{merged_dir}/locus_type_map.tsv"
    sf_conflicts = f"{merged_dir}/locus_family_conflicts.tsv"
    sf_all_calls = f"{merged_dir}/all_calls.tsv"
    log(f"=== merge per-type postF_hc -> unified locus set (locus table, "
        f"win={locus_table.MERGE_WINDOW}) ===")
    loci, table = locus_table.assign_and_write(hc_pf_paths, sf_merged_hcpf, sf_type_map,
                                               sf_conflicts=sf_conflicts,
                                               sf_all_calls=sf_all_calls,
                                               sf_all_calls_bed=f"{merged_dir}/all_calls.bed")
    log(f"unified locus set: {len(loci)} calls -> {sf_merged_hcpf}")
    log(f"every family call, won or lost: {len(table.all_loci)} -> {sf_all_calls} (+ .bed)")
    # The family-contention population, reported rather than resolved in silence. This is
    # the number Stage 2c's gate is checked against: that stage moves family assignment
    # from here (after three independent per-family chains) to one up-front decision, and
    # "calls may move, but only as family-assignment flips" is uncheckable unless the flips
    # that ALREADY happen have been counted first.
    if table.conflicts:
        summary = table.conflict_summary()
        log(f"family conflicts resolved at merge: {len(table.conflicts)} "
            f"-> {sf_conflicts}")
        for (win, lose, rule) in sorted(summary):
            log(f"    {win:4s} kept over {lose:4s} by {rule:14s} {summary[(win, lose, rule)]}")
    else:
        log("family conflicts resolved at merge: 0")

    rtypes_bitmask = 0
    for t in hc_pf_paths:
        rtypes_bitmask |= REP_TYPES[t]["rtype"]

    log(f"=== shared tail (gene -> genotype -> gVCF), run once ===")
    for name, argv in build_shared_tail_steps(cfg, merged_dir, opt, sf_merged_hcpf, sf_type_map,
                                              str(sample_root / f"{opt.sample}.vcf"),
                                              bam_list, rtypes_bitmask):
        run_step(f"shared:{name}", argv, logfh)

    merged_vcf = str(sample_root / f"{opt.sample}.vcf")
    if os.path.exists(merged_vcf):
        n = sum(1 for ln in open(merged_vcf) if ln and not ln.startswith("#"))
        log(f"merged VCF: {merged_vcf} ({n} calls)")
    else:
        log(f"WARNING: merged VCF not produced: {merged_vcf}")


def cleanup_successful_run(sample_root, sample):
    """Retain only the deliverable merged VCF after a successful pipeline run.

    Cleanup is intentionally all-or-nothing: if gVCF did not produce the expected
    merged VCF, leave the entire tree untouched and fail loudly.  This preserves
    the diagnostics needed to investigate a failed or incomplete run.  Call only
    after every pipeline phase has completed successfully.
    """
    sample_root = Path(sample_root).resolve()
    final_vcf = sample_root / f"{sample}.vcf"
    if not final_vcf.is_file():
        raise RuntimeError(
            f"refusing to clean {sample_root}: final merged VCF is missing: {final_vcf}")

    removed = 0
    for child in sample_root.iterdir():
        if child == final_vcf:
            continue
        if child.is_dir() and not child.is_symlink():
            shutil.rmtree(child)
        else:
            child.unlink()
        removed += 1
    log(f"cleanup complete: retained {final_vcf}; removed {removed} intermediate entries "
        "(use --keep-intermediates to retain them)")


def resolve_decoy_structural(v):
    """The transduction decoy test's policy, from `--(no-)td-decoy-structural`'s three states.

    A module-level function and not two lines inside `main` because this is a SHIPPING DEFAULT,
    and a default that can only be exercised by running the whole driver is a default with no
    test. `None` is "neither flag given", which since 2026-08-05 resolves to ON.

    Kept separate from the env write so the two are separately checkable: this answers "what did
    the options ask for", and `main` is what puts the answer somewhere a subprocess can read.
    """
    return True if v is None else bool(v)


def resolve_polya_cutoff(v):
    """Gate [A]'s polyA threshold from the CLI, or None meaning DERIVE.

    `None` -- the flag absent -- resolves to the frozen `n_clip_cutoff / 2 - 1`, computed in
    `TransductionParams` where it has always lived. This is the third flag in the family and the
    first that is a VALUE rather than a boolean, so "absent means the frozen behaviour" is spelled
    as None-means-derive rather than as False.

    **Zero is a real answer here and must not collide with "no opinion"** -- `--td-polya-cutoff 0`
    is a caller asking for a threshold nothing can fail, which is a different run from one that
    never passed the flag (they happen to agree at HG002, where the derived 0.5 is also
    unreachable, and would diverge at `--nclip 5`). A separate function for the same reason as its
    two siblings: a shipping default only exercisable by running the whole driver has no test.
    """
    return None if v is None else float(v)


def resolve_decoy_concat(v):
    """Whether the transduction step CONCATENATES the decoy consensus block, from the flag.

    Sibling of `resolve_decoy_structural`, and since 2026-08-07 it agrees with it: `None` --
    neither flag given -- resolves to **ON**. `--no-td-decoy-concat` is the rollback.

    WHY THE DEFAULT MOVED. This docstring used to argue for OFF on the grounds that the flag
    "moves calls with no baseline to compare against yet", and that "its two-scale arm is what
    would earn a flip". That arm now exists and both scales PASS -- chr1 containment `49879284`
    (100 identical / 0 differing; L1's bwa reference byte-identical at 709 records) and WGS
    `49879677` against `49792012` as the frozen OFF arm. What they measured:

      * hg38's SVA flank library ships **0** of the 45 decoy records against both L1 libraries'
        45, so SVA has had no decoy channel on either the clip or the disc side. Appending is BY
        NAME, so this is a structural no-op on L1 and the fix on SVA.
      * SVA `decoy_cross_family` 0 -> 44 at WGS, while L1's whole ledger stays byte-identical at
        21,158 rows.
      * Exactly one call moves across the genome, `chr17:45585847`, and the truth set counts it a
        FALSE POSITIVE: SVA FP 8 -> 7, TP 77 and FN 13 unchanged, position precision 90.6% ->
        91.7%. The verifier calls it UNDECIDABLE (self-window, decision 6), so the truth set and
        the verifier agree about it.
      * Stage 3's no-regression gate PASSES (`49995743`), with both arms reporting an identical
        54.0%-76.2% precision range.
      * Cost ~55s, ~12% of the transduction phase, measured and mechanistically explained
        (docs/handoff/HANDOFF_2026-08-06_timing_resolved.md).

    A separate function for the same reason as its sibling: a shipping default that can only
    be exercised by running the whole driver is a default with no test.
    """
    return True if v is None else bool(v)


def main():
    ap = argparse.ArgumentParser(
        prog="xtea2",
        description="xTea2: BAM/CRAM -> MEI VCF, single self-contained process "
                    "(Milestone A0: golden-reproducing wrapper).")
    # required inputs
    ap.add_argument("-i", "--input", required=True,
                    help="input BAM or CRAM. One path, or multiple space-delimited "
                         "paths in a single argument value (e.g. -i \"a.cram b.cram "
                         "c.cram\") to pool several libraries of the same sample -- "
                         "functionally identical for one file or many.")
    ap.add_argument("-s", "--sample", required=True, help="sample identifier")
    ap.add_argument("-o", "--output-dir", required=True, help="output directory")
    ap.add_argument("--reference", required=True,
                    help="genome FASTA (required; also used to decode CRAM)")
    ap.add_argument("--rep-lib-dir", required=True,
                    help="path to rep_lib_annotation/ directory")
    ap.add_argument("--gene-gff3", required=True, help="gene annotation .gff3")
    ap.add_argument("--blacklist", default="null",
                    help="blacklist BED (centromere/dup); default: null (off)")
    # selection / resources
    ap.add_argument("--types", default="L1,Alu,SVA",
                    help="comma list of MEI types (default: L1,Alu,SVA)")
    ap.add_argument("--genome", default="hg38",
                    help="genome build subdir under rep-lib (default: hg38)")
    ap.add_argument("-n", "--cores", type=int, default=8, help="worker cores")
    ap.add_argument("--parallel-downstream", action="store_true",
                    help="Phase C speed: run the per-type (L1/Alu/SVA) downstream chains "
                         "concurrently instead of serially. The chains are independent after "
                         "the shared scan+cns-fan (master runs them as separate jobs), so this "
                         "is a pure wall-clock win. Each type gets cores//ntypes workers so the "
                         "sum stays within the allocation; results are unchanged. Each type "
                         "writes its own <sample>_<type>.steps.log.")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="retain all intermediate files and logs after a successful run. "
                         "By default, successful runs retain only <sample>.vcf; failed runs "
                         "are never cleaned.")
    # exposed filter/threshold knobs (defaults pinned to the golden run)
    # NOTE: --nclip/--cr/--nd were removed by trim step 4 along with the -C/-D steps
    # they fed. Their scan-path equivalents: --scan-min-clip replaces --nclip/--cr,
    # and the disc cutoff is derived from coverage in-process (Parameters.
    # get_par_by_cov, see the disc params block below) rather than passed as --nd.
    # No caller passed them; removing rather than silently ignoring them means an old
    # invocation fails loud instead of appearing to set a cutoff it no longer sets.
    ap.add_argument("--nfclip", type=int, default=4,
                    help="min clipped reads at cns/final filter (golden 4)")
    ap.add_argument("--nfdisc", type=int, default=5,
                    help="min discordant reads at cns/final filter (golden 5)")
    ap.add_argument("--flklen", type=int, default=3000,
                    help="flank length for realign/transduction (golden 3000)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the per-step commands and exit (no execution)")
    ap.add_argument("--region-mask", default=None,
                    help="regions whose CANDIDATE POSITIONS are excluded, applied once after "
                         "the scan pass so masked loci never reach cns_filter/transduction. "
                         "Comma-separated mix of built-in names (simple_repeats, segdup, "
                         "centromeres), explicit BED/BED.gz paths, 'default' for all three "
                         "built-ins, or 'none' to disable. "
                         "UNSET means mode-dependent: 'default' under --mosaic, 'none' "
                         "otherwise -- the evidence for masking was gathered entirely on "
                         "mosaic runs, and switching it on for germline would move a "
                         "long-standing baseline on no measurement at all. Pass it "
                         "explicitly to mask a germline run. NOTE this masks candidate "
                         "positions only; the discordant-read evidence array is deliberately "
                         "left whole (see xtea/region_mask.py).")
    ap.add_argument("--adopt-existing", action="store_true",
                    help="with --resume: accept scan_pass/cns_emit output produced BEFORE "
                         "resume markers existed, writing a marker for it. Trusts files this "
                         "process did not produce and cannot verify -- a truncated output "
                         "looks complete. Use only after confirming the phase finished in "
                         "that run's own log.")
    ap.add_argument("--resume", action="store_true",
                    help="reuse a previous run's completed scan_pass / cns_emit output when "
                         "its inputs (and the scan binary) are unchanged. Markers live in "
                         "<sample_root>/.resume and are written only after a phase succeeds, "
                         "so a phase killed mid-write is always redone.")
    # CC.0: tumor/purity plumbing (prerequisite for case-control mode, independent
    # of it -- a tumor-only single-sample run without a control is also a valid
    # master mode). Mirrors master's gnrt_pipeline_local_v38.py exactly: purity
    # scales the coverage-derived clip/disc cutoffs (x_TEA_main.py:
    # automatic_gnrt_parameters, f_cov *= f_purity when --tumor is set).
    ap.add_argument("--tumor", action="store_true",
                    help="working on a tumor sample: scale automatic clip/disc "
                         "cutoffs by --purity (master's --tumor/--purity)")
    ap.add_argument("--purity", type=float, default=0.45,
                    help="tumor purity, used only with --tumor (default 0.45, "
                         "matches master's default)")
    # CC.1: case-control (somatic) mode. --case-control implies --tumor should
    # also be set (the per-type chain's automatic cutoffs are still tumor-BAM
    # coverage-derived) but this is NOT enforced here -- matches master, which
    # also lets the two vary independently.
    ap.add_argument("--case-control", action="store_true", dest="case_control",
                    help="somatic mode: filter -i's candidates against --normal-bam "
                         "(master's --case_control); output is somatic calls only, "
                         "no germline VCF is produced")
    ap.add_argument("--normal-bam", default=None,
                    help="matched control/normal BAM or CRAM (required with "
                         "--case-control). Same space-delimited multi-file "
                         "convention as -i.")
    # MOSAIC_PLAN.md M.1/M.2. Mosaic is a parameter profile plus one scan field, not a fourth
    # pipeline (plan sec 6.4) -- single-sample in, same one-BAM-list contract as every other mode.
    # --joint (cross-sample singleton filtering) is deliberately out of scope: it is a post-hoc
    # set operation on finished per-sample lists, not a statistical joint caller, and it does not
    # touch the speed story this mode exists for.
    ap.add_argument("--mosaic", action="store_true",
                    help="mosaic mode: call low-AF insertions from high-coverage (200x+) "
                         "bulk data. Turns on the xtea-scan polyA tally, low clip/disc "
                         "cutoffs, and the mosaic post-filter. Single-sample; pool "
                         "replicate libraries via the usual multi-file -i.")
    # THE cutoff decision for mosaic, and it is arithmetic rather than a tuning knob.
    #
    # Every other mode derives its disc cutoff from coverage (Parameters.get_par_by_cov).
    # That is right for germline, where a real insertion sits at 50% or 100% VAF and deeper
    # data should demand proportionally more evidence. It is exactly WRONG for mosaic, whose
    # entire target is 1-8% VAF: at 466x the table returns clip=14 / disc=20, while a 2% VAF
    # insertion has ~9 alt reads in TOTAL across clip and disc combined. The coverage-derived
    # cutoff would reject the whole benchmark ladder before the mosaic post-filter ever ran --
    # and the run would look like a clean negative result rather than a misconfiguration.
    #
    # 2/2 is master's own stated mosaic design, not a guess: x_mosaic_calling.py's module
    # docstring says "The data is sequenced at high coverage (say >200X), and low cutoff is
    # set (say 2 disc and 2 clip)". master reaches it via --user_specific plus explicit
    # --lclip/--rclip; here --mosaic implies it, because leaving the automatic path on at
    # this depth is a silent trap rather than a defensible default.
    ap.add_argument("--mosaic-min-clip", type=int, default=2,
                    help="mosaic clip-support cutoff (default 2). Replaces the "
                         "coverage-derived value, which at 200x+ exceeds the total read "
                         "support a low-AF insertion can have.")
    ap.add_argument("--mosaic-min-disc", type=int, default=2,
                    help="mosaic discordant-support cutoff (default 2). Same reasoning as "
                         "--mosaic-min-clip.")
    ap.add_argument("--ignore-rg", action="store_true",
                    help="skip the @RG SM-tag consistency check across multi-file "
                         "-i/--normal-bam input. Off by default: differing SM tags "
                         "across input files fail loud (likely means you're pooling "
                         "different samples by accident). Pass this when pooling "
                         "libraries from different sequencing centers for the same "
                         "donor/tissue is intentional (their SM tags legitimately differ).")
    # A1 refactor-inward: replace master's clip step with the C++ xtea-scan pass.
    # Trim step 4: this is now unconditional. The flag is retained as an accepted
    # no-op purely so the ~15 validation scripts that pass it keep working.
    ap.add_argument("--scan-clip", action="store_true",
                    help="(no-op, always on since trim step 4) the C++ xtea-scan "
                         "single-pass clip+realign IS the clip step; master's "
                         "x_TEA_main.py -C no longer exists")
    ap.add_argument("--scan-bin", default=DEFAULT_SCAN_BIN,
                    help="path to the xtea_scan binary (for --scan-clip)")
    ap.add_argument("--bwa", default="bwa",
                    help="bwa binary used by xtea-scan realign (default: bwa on PATH)")
    ap.add_argument("--scan-min-clip", type=int, default=2,
                    help="min in-repeat clipped reads (nL+nR) to emit a clip candidate "
                         "from the scan (default 2; drops count-1 permissive-mapping noise)")
    # A2 refactor-inward: replace master's disc step with a native pass over the
    # single scan's --disc-pos stream (removes the per-type BAM re-scan machinery).
    ap.add_argument("--scan-disc", action="store_true",
                    help="(no-op, always on since trim step 4) the xtea-scan --disc-pos "
                         "stream, collected in the same pass, IS the disc step; master's "
                         "x_TEA_main.py -D no longer exists")
    ap.add_argument("--scan-disc-threshold", type=int, default=2000,
                    help="xtea-scan same-chrom disc emission threshold (bp). Kept well "
                         "below master's 100kb test so the stream is a superset; the "
                         "adapter re-applies the exact threshold (default 2000)")
    ap.add_argument("--scan-realign-ref", choices=["copies_with_flank", "consensus"],
                    default="copies_with_flank",
                    help="reference xtea-scan realigns clipped parts to. "
                         "copies_with_flank (default) matches master's clip step and is more "
                         "sensitive (clips over the flank/transduced region map); consensus is "
                         "the bare cns (less sensitive).")
    # DEFAULT ON since 2026-08-13 (user's call), for every mode rather than mosaic only.
    #
    # The measurement that had kept it off -- HANDOFF_2026-07-26, "dropping --cns-emit-scan
    # recovers 31.4 min" on HG002 WGS -- was taken against the IMPLEMENTATION THAT NO LONGER
    # EXISTS. That version folded collection into the detection pass, gathering
    # candidate-agnostically and spilling ~1.1 GB to be fanned afterwards, so it paid for
    # every read in the genome; it is deleted (see the two-phase type loop below). The
    # current one fetches only candidate windows, merged across types, so its cost scales
    # with candidates instead of with genome size and the 31 min has no mechanism left.
    #
    # What replaces it, measured on the chr1 mosaic A/B (jobs 50434768 ON vs 50434907 OFF,
    # one commit, same 4 BAMs, identical candidate lists):
    #
    #     L1:cns_filter    19.0 s ON   vs   1833.5 s OFF      (~96x)
    #     + cns_emit phase 954.8 s ON   vs        --
    #     net              ~974 s ON   vs   ~1834 s OFF       (ON wins by ~14 min)
    #
    # ON wins outright even after paying for its own collection phase, at 83,775 candidates.
    # At germline's ~1,500 the two are much closer, but the WGS gate (job 49502115) settled
    # that the arms are ANSWER-IDENTICAL -- 1493 merged, 188/1223/82 per family, same truth
    # score -- so defaulting ON cannot change a germline call set, only its runtime.
    #
    # --no-cns-emit-scan keeps the OFF arm reachable; mosaic_hapmap_chr1.sbatch's CNS_EMIT=0
    # branch passes it, and that A/B must keep working.
    ap.add_argument("--cns-emit-scan", dest="cns_emit_scan",
                    action="store_true", default=True,
                    help="emit the cns_filter (-N) step's temp_clip.sam/temp_disc.sam/"
                         "low_mapq_clip_cnt.txt natively with xtea-scan and reuse them via the "
                         "bypass, instead of collecting in Python. Same access pattern as the "
                         "Python step -- an indexed fetch per candidate window -- but the windows "
                         "are merged, the collection is threaded rather than forked, and nothing "
                         "goes through per-candidate temp files. Uses --scan-bin/--bwa. "
                         "DEFAULT ON; use --no-cns-emit-scan to collect in Python instead.")
    ap.add_argument("--no-cns-emit-scan", dest="cns_emit_scan", action="store_false",
                    help="collect the cns_filter inputs in Python, one indexed fetch and one "
                         "temp file per candidate. The pre-2026-08-13 default; kept so the "
                         "collection-path A/B stays runnable. Expect it to be far slower at "
                         "mosaic candidate counts (~96x on chr1's L1:cns_filter).")
    # Stage 2c.3. DEFAULT OFF, and the WGS gate (job 49502115) is why it STAYS off rather
    # than why it is still unproven: HG002 WGS, both arms at one commit, the ON arm's call
    # set is IDENTICAL to OFF -- 1493 merged, 188/1223/82 per family, identical truth score.
    #
    # An earlier version of this comment predicted 4 moved loci, all toward truth, from a
    # simulation over BASELINE-0's kept intermediates. That prediction was WRONG twice over
    # and the reason is the thing to carry: the barrier does not REPLACE the postF_hc merge,
    # both run, so it can only change output where it drops a row that would have WON at the
    # merge -- and of the 16 conflicts it resolves, the 7 the merge also sees give the
    # identical winner, loser and rule. (The other error was the probe's, matching each late
    # locus to the nearest early locus of ANY family; see docs/.../probe_early_assign.py.)
    #
    # So the flag is answer-preserving and currently buys nothing on its own: it is a SECOND
    # decision point, not a relocated one. THE DECISION, 2026-08-02: leave it off, and flip
    # it together with the change that makes it load-bearing -- unified site creation, which
    # is what would let the merge's resolution be deleted. 15 of the 22 conflicts are
    # unreachable from cns today because transduction and sibling create the losing family's
    # record AFTER the barrier runs, so no decision point before transduction can resolve
    # them. Until that lands, turning this on adds a barrier and keeps the merge.
    ap.add_argument("--family-assign-early", dest="family_assign_early",
                    action="store_true", default=False,
                    help="decide each locus's family ONCE, right after cns_filter, and carry "
                         "only its own candidates through transduction/sibling/postF -- instead "
                         "of running all three families to completion and resolving the family "
                         "at the postF_hc merge. Same rule either way (locus_table: priority "
                         "L1>SVA>Alu with a two_side_tprt_both override); cns_filter is simply "
                         "the earliest step whose output carries that flag. No-op for a "
                         "single --types entry.")
    # The per-locus loop's iteration domain, enumerated once in the driver (plan sec 3).
    # OFF by default because it is diagnostic until the loop consumes it: the families' own
    # transduction steps still derive the same `.slct` through the same one copy of the rule,
    # so this writes a replica and cannot move a call. It opens the same cns_filter ->
    # transduction barrier --family-assign-early does (see downstream_waves) but does not
    # imply it -- knowing a site exists is not knowing which family claims it, so the
    # 7-of-22 assignment ceiling is untouched by turning this on.
    ap.add_argument("--locus-universe", dest="locus_universe",
                    action="store_true", default=False,
                    help="enumerate the union locus universe (every family's cns candidates "
                         "plus the sites its transduction step will create) once in the "
                         "driver, after cns_filter and before transduction, to "
                         "merged/locus_universe/. Diagnostic: each family's chain still "
                         "computes its own site list, and the two must be byte-identical.")
    # The batched evidence stage of the per-locus loop: one BAM decode for every family's
    # transduction reads, in the driver, instead of one per family inside each step. Implies
    # --locus-universe (it needs the site lists that pass writes). OFF by default: it is
    # answer-preserving by construction (the gather is the same C++ code over the same sites
    # at the same extnd -- run_td_gather_multi IS what the single-family form calls) but it
    # moves work across a process boundary, so it stays opt-in until gated per sample type.
    # Do NOT flip it expecting wall clock -- see run_td_gather_driver's docstring.
    ap.add_argument("--td-gather-driver", dest="td_gather_driver",
                    action="store_true", default=False,
                    help="gather every family's transduction reads once in the driver "
                         "(xtea_scan --td-manifest) instead of once per family inside each "
                         "transduction step. Implies --locus-universe. Requires the native "
                         "gather (--td-gather-scan, the default).")
    # How the cns_filter -> transduction split is RUN, once one of the flags above has
    # opened it. Staged is the default and changes no output: same work, same order, but
    # each family waits only for what it reads (StagedBarrier). The global form is the
    # join every family paid for until now -- +52.9s at WGS, all of it lost overlap with
    # Alu's ~2x cns_filter. It stays reachable because it is the arm the numbers in
    # run_td_gather_driver were measured on, so an A/B remains a one-flag change rather
    # than a checkout. No effect at all on a run with no barrier, which is the default one.
    ap.add_argument("--global-barrier", dest="global_barrier",
                    action="store_true", default=False,
                    help="run the cns_filter -> transduction split as one global join "
                         "(every family waits for the slowest cns_filter) instead of the "
                         "default staged form where each family waits only for the work it "
                         "reads. Same output; slower. For A/B measurement.")
    # DEFAULT ON since the Stage 1 flip gate (job 49480961) passed: HG002 WGS, two arms at
    # one commit, 13/13 artifacts identical between them and 4/4 + 6/6 identical to
    # BASELINE-0 on both. --no-td-gather-scan is the rollback and is kept for exactly that.
    ap.add_argument("--td-gather-scan", dest="td_gather_scan",
                    action="store_true", default=True,
                    help="(default, kept for symmetry with --no-td-gather-scan) run the "
                         "transduction step's read gather natively with xtea-scan.")
    ap.add_argument("--no-td-gather-scan", dest="td_gather_scan", action="store_false",
                    help="fall back to the Python transduction read gather. The native gather "
                         "(src/csrc/td_gather.hpp) is 13-14x faster at WGS scale and takes "
                         "downstream_parallel from 1144.6s to 457.4s, so this is a rollback "
                         "switch, not a tuning knob. The Python version forks a Pool per "
                         "chromosome with one task and two temp files per candidate site, then "
                         "scans the whole genome in 50 Mb bins re-parsing its disc-pos file once "
                         "per bin. The two are byte-identical on real data at 49K/54K sites, "
                         "which matters because bwa's tie-break among equal-scoring placements "
                         "depends on input order. ONE KNOWN DIFFERENCE, and it is a fix: with "
                         "N>1 input BAMs the Python gather truncates one shared all_disc_pos per "
                         "BAM, so only the last BAM's positions survive and every sample's mates "
                         "are gathered against them. No baseline exists for that path.")
    # OFF by default, and the default is not the correct behaviour -- it is the FROZEN one.
    # The transduction step decides "this disc read hit a decoy consensus, not a source" with a
    # substring match on the bwa reference name, and every SVA flank record is named
    # `chr13~51077120~51078923~SVA_D~0L`, so 100% of hg38's 5,454 SVA flank records match and are
    # skipped before the source is extracted. Consequence, measured on WGS job 49634290: 0 of
    # SVA's 105 picked transduction sources are reference copies (L1: 114 of 159). This flag makes
    # the test structural (transduction_evidence.decoy_tags), which is a FIX and which moves SVA
    # calls -- so it ships off, BASELINE-0 and every gate frozen against it stay valid, and the
    # flip is its own measured decision. `docs/handoff/runs/probe_td_source_classes.py` is the
    # before/after instrument.
    # ON BY DEFAULT SINCE 2026-08-05. `default=None` and not `default=True`, so the three states
    # stay distinguishable -- not asked / asked for / asked against -- which is what lets the
    # banner below name the channel. A gate that cannot tell "the default armed it" from "a flag
    # armed it" cannot prove which one it measured.
    ap.add_argument("--td-decoy-structural", dest="td_decoy_structural",
                    action="store_true", default=None,
                    help="in the transduction step, tell a decoy consensus record from a flank "
                         "record by the NAME's STRUCTURE (a flank id is "
                         "chrom~start~end~subfamily~rc, with integer coordinates; a consensus id "
                         "is bare) instead of by a substring match on LINE1/SVA/ALU/HERV. Fixes "
                         "SVA reference-copy transductions, which are unreachable without it. "
                         "ON BY DEFAULT since 2026-08-05; this flag is kept as an explicit "
                         "opt-in so existing harnesses keep working and keep saying what they "
                         "meant.")
    # `default=None` is REQUIRED here, not redundant with the option above. `store_false`'s
    # implicit default is True, and argparse seeds the namespace from the FIRST action that names
    # a dest -- so with these two declared the other way round, "neither flag given" would arrive
    # as True instead of None. The resolved policy would be unchanged (True and None both mean
    # ON), but the BANNER would then claim `--td-decoy-structural` armed a run nothing armed, and
    # two gate harnesses read that banner as their only evidence of which configuration ran.
    # Pinned by `test_the_driver_defaults_the_decoy_test_ON`'s order check.
    ap.add_argument("--no-td-decoy-structural", dest="td_decoy_structural",
                    action="store_false", default=None,
                    help="revert the transduction decoy test to the pre-2026-08-05 SUBSTRING "
                         "match. This is the arm that now needs a reason to exist: the substring "
                         "form classifies 100%% of hg38's 5,454 SVA flank records as the SVA "
                         "consensus decoy, so NO reference SVA copy can be a transduction "
                         "source. Use it to reproduce a baseline frozen before the flip -- "
                         "BASELINE-0 needs it -- and not otherwise.")
    # The PACKAGING half of the same defect. `default=None` for the same reason as the pair
    # above -- so the banner can say whether a flag or the default chose.
    ap.add_argument("--td-decoy-concat", dest="td_decoy_concat",
                    action="store_true", default=None,
                    help="(default since 2026-08-07, kept for symmetry with "
                         "--no-td-decoy-concat) in the transduction step, concatenate the "
                         "45-record decoy consensus block into the bwa reference when the "
                         "packaged flank library does not already carry it. hg38's SVA flank "
                         "library carries NONE (both L1 libraries carry all 45), so without this "
                         "SVA has no decoy channel on either the clip or the disc side. The "
                         "block is rebuilt from <rep_lib>/consensus/ and is byte-identical to "
                         "the packaged one; appending is BY NAME, so it is a structural no-op on "
                         "both L1 libraries.")
    ap.add_argument("--no-td-decoy-concat", dest="td_decoy_concat",
                    action="store_false", default=None,
                    help="leave the transduction bwa reference as the packaged flank library "
                         "plus the polymorphic flanks -- pre-2026-08-07 behaviour, under which "
                         "no reference SVA copy has a decoy competitor. Use it to reproduce a "
                         "baseline frozen before the flip -- BASELINE-0 needs it -- and not "
                         "otherwise.")
    # Gate [A]'s SECOND threshold, which until now was a function of the first.
    ap.add_argument("--td-polya-cutoff", dest="td_polya_cutoff", default=None, metavar="N",
                    help="gate [A]'s polyA-signature threshold, decoupled from --nclip. "
                         "DEFAULT (flag absent): derive it as n_clip_cutoff/2-1, the frozen "
                         "behaviour. One dial for two thresholds is why a clip-support sweep is "
                         "single-variable only up to 4 -- at --nclip 5 the polyA arm wakes and "
                         "rejects 286 (L1) / 364 (SVA) loci the shipping configuration keeps, "
                         "and the two moves get reported as one. At HG002's operating point the "
                         "derived value is 0.5 and UNREACHABLE (every polyA record has a "
                         "count >= 1 by construction), so passing 0.5 here must be "
                         "byte-identical to the default and the first call-moving value is one "
                         "above 1. Does NOT touch the sibling step's n_clip_cutoff/2, which is "
                         "a different quantity.")
    ap.add_argument("--mate-in-rep-families", default=None,
                    help="which repeat families' annotations contribute to "
                         "mate_in_rep_clip_cnt, comma-separated from {l1,alu,sva}. That value "
                         "is the # clipped reads at a site whose mate lands in ANY of these "
                         "families' repeats (union), emitted identically for all types as the "
                         f"6th realign-out column. Default '{MATE_IN_REP_FAMILIES_DEFAULT}' "
                         "(SVA + Alu) -- the biologically tightest set: the value is consumed "
                         "only by the SVA cns path, and SVA is a composite element with an "
                         "Alu-derived domain, so mate-in-Alu is genuine SVA evidence. "
                         "'l1,alu,sva' = full union; 'sva' = per-family-strict. Independent of "
                         "--types by design. See spec §4 DC1.")
    ap.add_argument("--col5-families", default=None,
                    help="DEPRECATED alias for --mate-in-rep-families. 'col5' named a column "
                         "position rather than the concept; kept working so existing run "
                         "scripts do not break.")
    opt = ap.parse_args()

    # Mosaic's low-cutoff profile, applied HERE rather than inside build_type_steps so the
    # whole chain -- cns_filter's --cr/--nd, the scan's clip emission, the disc adapter --
    # sees one consistent set of numbers, and so build_type_steps stays a pure function of
    # opt that the tests can drive directly. See --mosaic-min-clip for the arithmetic.
    if opt.mosaic:
        log(f"[mosaic] clip/disc support cutoffs forced to "
            f"{opt.mosaic_min_clip}/{opt.mosaic_min_disc} (were --cr {opt.nfclip} / "
            f"--nd {opt.nfdisc}); the coverage-derived path is disabled for this mode")
        opt.nfclip = opt.mosaic_min_clip
        opt.nfdisc = opt.mosaic_min_disc
        opt.scan_min_clip = min(opt.scan_min_clip, opt.mosaic_min_clip)

    # Reconcile the new flag with its deprecated alias before anything reads the value.
    if opt.col5_families is not None:
        log("WARNING: --col5-families is deprecated; use --mate-in-rep-families. "
            "'col5' named the 6th realign-out column, not a concept -- the value is "
            "mate_in_rep_clip_cnt, the number of clipped reads at a site whose mate lands "
            "inside any selected repeat family's annotation.")
        if (opt.mate_in_rep_families is not None
                and opt.mate_in_rep_families != opt.col5_families):
            sys.exit("ERROR: --mate-in-rep-families and its deprecated alias --col5-families "
                     f"were both given, with different values ({opt.mate_in_rep_families!r} "
                     f"vs {opt.col5_families!r}). Pass one.")
        if opt.mate_in_rep_families is None:
            opt.mate_in_rep_families = opt.col5_families
    if opt.mate_in_rep_families is None:
        opt.mate_in_rep_families = MATE_IN_REP_FAMILIES_DEFAULT

    opt._mate_in_rep_fams = {s.strip().lower()
                             for s in opt.mate_in_rep_families.split(",") if s.strip()}
    _bad = opt._mate_in_rep_fams - set(MATE_IN_REP_KEY_TO_TYPE)
    if _bad:
        sys.exit(f"ERROR: --mate-in-rep-families: unknown families {sorted(_bad)} "
                 f"(choose from {','.join(sorted(MATE_IN_REP_KEY_TO_TYPE))})")
    if not opt._mate_in_rep_fams:
        sys.exit("ERROR: --mate-in-rep-families: empty family set; pass at least one of "
                 f"{','.join(sorted(MATE_IN_REP_KEY_TO_TYPE))}")

    # Trim step 4: scan-only is the only path. Both flags are forced on regardless of
    # what the caller passed, so the two attributes below still read True everywhere
    # they are logged or reported; nothing branches on them any more.
    opt.scan_clip = True
    opt.scan_disc = True

    if opt.case_control and not opt.normal_bam:
        sys.exit("ERROR: --case-control requires --normal-bam")
    # Checked HERE, not at the barrier where the flag is consumed: the barrier is on the far
    # side of the scan pass, so failing there costs the whole detection run before saying no.
    if getattr(opt, "td_gather_driver", False) and not opt.td_gather_scan:
        sys.exit("ERROR: --td-gather-driver needs the native gather; it precomputes exactly "
                 "what the native path reads. Drop --no-td-gather-scan.")

    # Exported HERE rather than beside XTEA_SCAN_TD_GATHER below, because that one lives inside
    # `if not opt.dry_run` and names a binary the scan pass resolves; this one is a policy the
    # transduction subprocess reads on every path, dry run included.
    #
    # ALWAYS SET, BOTH WAYS, SINCE THE 2026-08-05 FLIP -- and that is a change of kind, not just
    # of default. Before the flip the var was set only when asked, so that an inherited "0" could
    # not look like an opt-in. With the fix ON by default the risk reverses: an inherited value
    # could now look like an opt-OUT and silently restore the defect. So `run_xtea` writes the
    # value its own options resolve to and is the SINGLE WRITER of this policy; whatever the
    # environment said before this line is overwritten, which is the property that makes the
    # banner below trustworthy evidence.
    #
    # THE ENV VAR'S OWN SEMANTICS ARE UNCHANGED: absent still means OFF
    # (`x_transduction.py:78`). That asymmetry is deliberate, not an oversight. `x_TEA_main` is
    # invoked directly -- without `run_xtea` -- by a dozen frozen gate harnesses under
    # `docs/handoff/runs/`, every one of whose baselines was measured with the substring test. If
    # absent meant ON, this flip would silently move all of them at once. The default is a DRIVER
    # policy with a blast radius of one entry point; a direct `x_TEA_main` call is not the default
    # path and gets the frozen behaviour until it asks otherwise.
    b_decoy_structural = resolve_decoy_structural(getattr(opt, "td_decoy_structural", None))
    os.environ[global_values.ENV_TD_DECOY_STRUCTURAL] = "1" if b_decoy_structural else "0"
    if b_decoy_structural:
        # The `[decoy] ... is STRUCTURAL` prefix is load-bearing: two gate harnesses grep for it
        # as proof the policy reached the driver. Keep it verbatim; vary only the tail.
        print("[decoy] the transduction decoy test is STRUCTURAL "
              "({0}): SVA reference-copy sources are reachable, and this run "
              "does NOT reproduce BASELINE-0".format(
                  "--td-decoy-structural" if opt.td_decoy_structural
                  else "DEFAULT since 2026-08-05"))
    else:
        # A banner for the OFF arm too, because otherwise "no STRUCTURAL banner" cannot be told
        # from "the driver never got this far" -- and an off-arm gate whose evidence is the
        # ABSENCE of a line is a gate that passes when the run dies early.
        print("[decoy] the transduction decoy test is the SUBSTRING match "
              "(--no-td-decoy-structural): pre-2026-08-05 behaviour, NO reference SVA copy can "
              "be a transduction source. This arm exists to reproduce a pre-flip baseline")

    b_decoy_concat = resolve_decoy_concat(getattr(opt, "td_decoy_concat", None))
    os.environ[global_values.ENV_TD_DECOY_CONCAT] = "1" if b_decoy_concat else "0"
    if b_decoy_concat:
        # As with the structural banner, the `[decoy-concat] ...` prefix is a gate's evidence.
        # It says the POLICY reached the driver and nothing more -- how many records each family
        # actually gained is x_transduction's own `[td-decoy-concat]` line, which is the number
        # that separates "armed" from "changed anything".
        print("[decoy-concat] the transduction bwa reference GAINS the decoy consensus block "
              "where the packaged flank library lacks it ({0}): this run does NOT reproduce a "
              "baseline measured without it".format(
                  "--td-decoy-concat" if opt.td_decoy_concat is not None else "DEFAULT"))
    else:
        print("[decoy-concat] the transduction bwa reference is the packaged flank library plus "
              "the polymorphic flanks only (--no-td-decoy-concat): pre-2026-08-07 behaviour, so "
              "hg38's SVA has NO decoy channel on either side. This arm exists to reproduce a "
              "pre-flip baseline -- BASELINE-0 needs it")

    # Gate [A]'s polyA threshold. A VALUE, not a boolean, so the env channel carries "derived"
    # rather than "0" for the frozen behaviour -- 0 is a threshold someone may legitimately want.
    f_polya_cutoff = resolve_polya_cutoff(getattr(opt, "td_polya_cutoff", None))
    os.environ[global_values.ENV_TD_POLYA_CUTOFF] = (
        "derived" if f_polya_cutoff is None else repr(f_polya_cutoff))
    if f_polya_cutoff is None:
        print("[td-polya-cutoff] gate [A]'s polyA threshold is DERIVED from --nclip as "
              "n_clip_cutoff/2-1 (DEFAULT): the frozen behaviour")
    else:
        # As with the two decoy banners, this says the POLICY reached the driver and nothing
        # more. What the gate actually compared against is x_transduction's own
        # `[td-polya-cutoff]` line, which prints the derived value it replaced beside it.
        print("[td-polya-cutoff] gate [A]'s polyA threshold is DECOUPLED from --nclip and set "
              "to {0} (--td-polya-cutoff): this run does NOT reproduce a baseline measured "
              "with the derived threshold unless that threshold was also {0}".format(
                  f_polya_cutoff))

    # normalise requested types
    types = []
    for tok in opt.types.split(","):
        tok = tok.strip()
        if not tok:
            continue
        canon = TYPE_ALIASES.get(tok.lower(), tok)
        if canon not in REP_TYPES:
            sys.exit(f"ERROR: unknown type '{tok}' (known: L1,Alu,SVA)")
        types.append(canon)

    input_paths = resolve_inputs(opt.input, opt.reference)
    check_rg_consistency(input_paths, opt.reference, opt.ignore_rg)

    sample_root = Path(opt.output_dir).resolve() / opt.sample
    sample_root.mkdir(parents=True, exist_ok=True)
    pub_clip = str(sample_root / "pub_clip") + "/"
    Path(pub_clip).mkdir(parents=True, exist_ok=True)

    # one bam_list.txt for the whole run -- one line per input file/library.
    # Master's own legacy pooling (X_BasicInfo/XClipDiscFilter/TE_Multi_Locator)
    # already sums coverage / averages read-length+insert-size / flat-sums clip+disc
    # counts across N>1 lines; writing N lines here is all that's needed to use it.
    bam_list = str(sample_root / "bam_list.txt")
    with open(bam_list, "w") as fh:
        for p in input_paths:
            fh.write(f"{p}\tillumina\n")

    # CC.1: control bam_list, built the same way as the tumor one -- every
    # consumer on the case-control path only ever reads column 0 (verified by
    # reading XClipDiscFilter.__init__/XGenotyper.call_genotype/
    # X_BasicInfo.collect_basic_info_samples), so the two-column convention here
    # is a style choice, not a format requirement.
    control_bam_list = None
    if opt.case_control:
        control_paths = resolve_inputs(opt.normal_bam, opt.reference)
        check_rg_consistency(control_paths, opt.reference, opt.ignore_rg)
        control_bam_list = str(sample_root / "control_bam_list.txt")
        with open(control_bam_list, "w") as fh:
            for p in control_paths:
                fh.write(f"{p}\tillumina\n")

    cfg = {
        # space-joined multi-path string -- handed straight to xtea_scan's --bam,
        # which parses the same convention. Downstream Python-side pooling reads
        # bam_list.txt (above) instead of this string directly.
        "input": " ".join(input_paths),
        "reference": str(Path(opt.reference).resolve()),
        "blacklist": ("null" if opt.blacklist in (None, "null", "NULL")
                      else str(Path(opt.blacklist).resolve())),
        "gene_gff3": str(Path(opt.gene_gff3).resolve()),
    }

    logpath = sample_root / f"{opt.sample}_xtea2.log"
    log(f"run start: sample={opt.sample} types={types} -> {sample_root}")
    log(f"full step log: {logpath}")
    if not opt.dry_run:
        log(f"step timings: {init_timing(sample_root)}")

    # resolve each requested type's consensus + paths once
    type_paths = {}
    for t in types:
        spec = REP_TYPES[t]
        type_paths[t] = {
            "rtype": spec["rtype"],
            "copies_with_flank": resolve_annot(opt.rep_lib_dir, opt.genome, spec["copies_with_flank"]),
            "annotation":  resolve_annot(opt.rep_lib_dir, opt.genome, spec["annotation"]),
            "annotation1": resolve_annot(opt.rep_lib_dir, opt.genome, spec["annotation1"]),
            "flank":       resolve_annot(opt.rep_lib_dir, opt.genome, spec["flank"]),
            "cns":         resolve_annot(opt.rep_lib_dir, opt.genome, spec["cns"]),
        }

    # When --cns-emit-scan is on, ONE xtea-scan invocation collects the cns clip+disc reads for
    # ALL types after their candidate lists exist, and emits each type's
    # temp_clip.sam/temp_disc.sam/low_mapq_clip_cnt.txt; we then set XTEA_CNS_PRECOMPUTED_ALIGN so
    # every cns_filter (-N) step reuses them via the bypass. It has to run after the candidate
    # lists are built, which is why the type loop below is two-phase.
    #
    # This used to be folded into the detection pass instead, collecting candidate-agnostically and
    # spilling ~1.1 GB to disk to be fanned afterwards. That is deleted: collecting before the
    # candidates are known meant paying for every read in the genome, and it measured 31 min slower
    # than the Python step it was meant to replace. See src/csrc/cns_emit.hpp.

    hc_pf_paths = {}  # type -> postF_hc path (Phase C move A.2: gene/genotype/gVCF run once, after the loop)
    with open(logpath, "w") as logfh:
        # A1/A2: one xtea-scan pass produces clip candidates (--realign-out) and the
        # discordant-read stream (--disc-pos) for ALL requested types -- all BAM
        # reading happens here, once. Unconditional since trim step 4.
        realign_tsv = None
        disc_arrays = None
        disc_params = None       # (iextend, std_var, disc_cutoff, b_with_chr)
        cns_extnd = None         # cns collection window = int(mean_is + 3*std) (== x_TEA_main max_is)
        cns_read_length = None   # READ_LENGTH the cns step sets (x_TEA_main --cns: set_read_length(rlth))
        if not opt.dry_run:
            scan_bin = str(Path(opt.scan_bin).resolve())
            if not Path(scan_bin).exists():
                sys.exit(f"ERROR: xtea_scan binary not found at {scan_bin}\n"
                         f"  build it:  make -C src/csrc     "
                         f"(PKG_CONFIG_PATH=<htslib-prefix>/lib/pkgconfig if htslib is not found)\n"
                         f"  or point at one:  --scan-bin /path/to/xtea_scan, "
                         f"or export XTEA2_SCAN_BIN, or put xtea_scan on PATH")
            # The transduction step runs as its own x_TEA_main subprocess, so the native gather is
            # selected by an inherited env var rather than an argv flag (same mechanism as
            # XTEA_CNS_PRECOMPUTED_ALIGN below). Unsetting it is the rollback.
            if opt.td_gather_scan:
                os.environ["XTEA_SCAN_TD_GATHER"] = scan_bin
            ref_key = ("copies_with_flank" if opt.scan_realign_ref == "copies_with_flank"
                       else "cns")
            cns_by_family = {scan_adapter.TYPE_TO_FAMILY[t]: type_paths[t][ref_key] for t in types}
            realign_tsv = str(sample_root / "scan_clip_realign.tsv")
            disc_pos_out = str(sample_root / "scan_disc_pos.tsv")
            # Mosaic thresholds on RAW per-position clip counts, which the realign-out does
            # not carry (it counts only clips that realigned to a consensus). Written only
            # under --mosaic so a germline pass produces no extra file.
            clip_sites_out = str(sample_root / "scan_clip_sites.tsv")
            # mate_in_rep_clip_cnt: convert each family's .out annotation to an extended
            # BED (SVA +200, L1/Alu 0 -- same as the clip step's annotation) and let
            # xtea-scan compute the count natively via htslib regidx. It counts a clipped
            # read's mate landing in ANY of the families selected by
            # --mate-in-rep-families (union), emitted identically for all types as the 6th
            # realign-out column; xtea-scan ORs whichever BEDs it receives (spec §4 DC1).
            #
            # The contributing family set is a property of the signal definition, not of
            # which families the caller asked for -- see mate_in_rep_bed_plan(), which is
            # where that reasoning and the bug it fixes are written down. Annotation paths
            # resolve from rep_lib_dir/genome alone, independently of `types`, so a
            # non-requested family's BED is available here; for requested families this
            # resolves to exactly the same path as type_paths[t]["annotation"], which is
            # why three-type output does not move.
            bed_by_family = {}
            for t, fam, extnd in mate_in_rep_bed_plan(opt._mate_in_rep_fams):
                bed = str(sample_root / f"rep_{fam}.bed")
                n_iv = scan_adapter.out_to_bed(
                    resolve_annot(opt.rep_lib_dir, opt.genome, REP_TYPES[t]["annotation"]),
                    bed, extnd)
                bed_by_family[fam] = bed
                log(f"mate-in-rep BED {fam}: {n_iv} intervals (extnd={extnd}) -> {bed}"
                    f"{'' if t in types else ' [not in --types; the count is type-set-independent]'}")
            log(f"mate-in-rep union families: {sorted(opt._mate_in_rep_fams)}")
            log(f"=== xtea-scan single pass (clip+realign {list(cns_by_family)}"
                f"+disc-pos+mate-in-rep) ===")
            scan_outputs = [realign_tsv, disc_pos_out] + \
                ([clip_sites_out] if getattr(opt, "mosaic", False) else [])
            scan_fp = _fingerprint([
                "scan_pass", cfg["input"], cfg["reference"], sorted(cns_by_family.items()),
                sorted(bed_by_family.items()), opt.scan_disc_threshold,
                getattr(opt, "mosaic", False), _binary_stamp(scan_bin)])
            if getattr(opt, "resume", False) and \
                    (phase_is_done(sample_root, "scan_pass", scan_fp, scan_outputs) or
                     adopt_existing_phase(sample_root, "scan_pass", scan_fp, scan_outputs,
                                          getattr(opt, "adopt_existing", False))):
                log(f"[resume] scan_pass: reusing completed output ({', '.join(scan_outputs)})")
                record_timing("phase", "scan_pass_RESUMED", 0.0)
            else:
                with phase("scan_pass"):
                    scan_adapter.run_scan_realign(
                        scan_bin, cfg["input"], cfg["reference"], cns_by_family, realign_tsv,
                        threads=opt.cores, bwa=opt.bwa, logfh=logfh,
                        disc_pos_out=disc_pos_out,
                        disc_threshold=opt.scan_disc_threshold,
                        bed_by_family=bed_by_family,
                        clip_sites_out=(clip_sites_out if getattr(opt, "mosaic", False) else None),
                        mosaic_polya=getattr(opt, "mosaic", False))
                mark_phase_done(sample_root, "scan_pass", scan_fp, scan_outputs)

            # insert-size params (read length, mean IS, std) + coverage-derived disc
            # cutoff -- computed once, exactly as master's disc step does. Unconditional
            # since trim step 4 (was gated on --scan-disc).
            from x_basic_info import X_BasicInfo
            from x_parameter import Parameters
            xbi = X_BasicInfo(str(sample_root / "tmp_basic") + "/", opt.cores, cfg["reference"])
            Path(str(sample_root / "tmp_basic")).mkdir(parents=True, exist_ok=True)
            rcd = xbi.get_cov_is_rlth(bam_list, cfg["reference"], 500, True)
            f_cov, rlth, mean_is, std_var = rcd[0], rcd[1], rcd[2], rcd[3]
            # Purity scaling MUST be applied before the table lookup, exactly as
            # x_TEA_main.py:automatic_gnrt_parameters does (f_cov *= f_purity when
            # b_tumor, then get_par_by_cov on the SCALED value). Threading
            # --tumor/--purity into the master steps' argv is not enough: this
            # cutoff is derived in-process here and bypasses master's automatic
            # path entirely. Omitting it made xTea2 filter disc at 5 where master
            # filtered at 4 on the LFS 3762A CRAM (raw cov 36.6, purity 0.86 ->
            # 31.48), silently dropping real calls the clip step had found -- the
            # CC.3 Alu 0/2 failure. Only the value fed to the lookup is scaled;
            # rlth/mean_is/std_var (and so iextend/cns_extnd below) stay raw,
            # matching master, which never scales those.
            f_cov_raw = f_cov
            if opt.tumor:
                f_cov = f_cov * opt.purity
            disc_cutoff = Parameters().get_par_by_cov(f_cov)[1]
            # Mosaic overrides the coverage-derived cutoff -- see --mosaic-min-disc for why
            # this is arithmetic and not preference. Logged loudly with BOTH numbers, because
            # the failure this prevents (every low-AF call filtered before the mosaic
            # post-filter runs) is indistinguishable in the output from "the sample is clean".
            if getattr(opt, "mosaic", False):
                log(f"[mosaic] disc cutoff {disc_cutoff} (coverage-derived at "
                    f"{f_cov:.1f}x) OVERRIDDEN to {opt.mosaic_min_disc}: a 1-8% VAF "
                    f"insertion cannot produce {disc_cutoff} discordant reads")
                disc_cutoff = opt.mosaic_min_disc
            # bam "chr" prefix convention for the annotation mate-in-rep test,
            # AND a robust read length sampled directly from the BAM. master's
            # get_cov_is_rlth reports rlth=0 here (a long-standing bug: "Read
            # length is: 0.0"); master's own disc step doesn't need it (its
            # bamfile.fetch returns overlapping reads natively), but our disc-pos
            # emulation widens the window by one read length to recover reads
            # that span a breakpoint -- so we must have a real value.
            # Multi-file: pool up to ~2000 read-length samples ACROSS all input
            # files (not just the first), then take one median -- mirrors
            # X_BasicInfo's own averaging convention for read length rather than
            # inventing a new pooling rule, and keeps this correct whether -i
            # named one file or several.
            import pysam
            _rls = []
            b_with_chr = None
            for _p in input_paths:
                _bf = pysam.AlignmentFile(_p, "rb", reference_filename=cfg["reference"])
                if b_with_chr is None:
                    b_with_chr = any(r.startswith("chr") for r in _bf.references[:5])
                for _a in _bf.fetch(until_eof=True):
                    if _a.is_unmapped or _a.is_secondary or _a.is_supplementary:
                        continue
                    if _a.query_length:
                        _rls.append(_a.query_length)
                    if len(_rls) >= 2000:
                        break
                _bf.close()
                if len(_rls) >= 2000:
                    break
            _rls.sort()
            sampled_rlth = _rls[len(_rls) // 2] if _rls else 150   # median; 150 fallback
            b_with_chr = bool(b_with_chr)
            overlap_rlth = int(sampled_rlth)   # real read length, for fetch-overlap emulation
            # iextend (window half-width) MUST match master exactly, incl. its
            # rlth=0 bug -- master builds it as mean_is + 3*std + rcd[1] and rcd[1]
            # is 0 here. Using the real read length would widen the window past
            # master's and over-count disc (over-counted sites trip a downstream
            # upper-bound disc filter and get dropped). So iextend uses master's
            # rlth; only the lower-bound overlap expansion below uses the real one.
            iextend = int(mean_is + 3 * std_var) + int(rlth)
            disc_params = (iextend, std_var, disc_cutoff, overlap_rlth, b_with_chr)
            # cns collection window = max_is = int(mean_is + 3*std) (NO +rlth): this is exactly
            # what x_TEA_main's cns step computes as iextnd, so the fold-in collect matches the
            # cns math's window. (max(400, .) not needed here -- WGS max_is >> 400.)
            cns_extnd = int(mean_is + 3 * std_var)
            # master's cns step places a low-MAPQ right clip at pos + READ_LENGTH - clip length,
            # with READ_LENGTH = this rlth (set_read_length before collecting); the native
            # cns-emit used a fixed 100 until it took --cns-read-length.
            cns_read_length = int(rlth) if rlth and rlth > 0 else overlap_rlth
            # cov_eff is the value the cutoff table was keyed on -- log it separately
            # from the raw estimate so it can be diffed directly against master's
            # "Ave coverage is <cov_eff>: automatic parameters (clip, disc, clip-disc)"
            # line, which also prints the purity-scaled figure.
            log(f"disc params: cov(raw)={f_cov_raw:.5f} cov(eff)={f_cov:.5f} "
                f"tumor={opt.tumor} purity={opt.purity} "
                f"rlth(master)={rlth} rlth(overlap)={overlap_rlth} "
                f"mean_is={mean_is:.0f} std={std_var:.0f} iextend={iextend} "
                f"disc_cutoff={disc_cutoff}")
            log(f"=== loading disc-pos stream {disc_pos_out} ===")
            with phase("disc_pos_load"):
                disc_arrays = disc_adapter.load_disc_pos(disc_pos_out)
            log(f"disc-pos loaded: {sum(len(v[0]) for v in disc_arrays.values())} records "
                f"across {len(disc_arrays)} contigs")

        # --- Phase 1: build every type's candidate lists (needed before the single cns collect) ---
        # Loaded ONCE for all types: the intervals are identical per type and parsing the
        # default set is ~650k intervals. Applied only to clip candidate POSITIONS -- the
        # disc evidence array stays whole, see xtea/region_mask.py for why that distinction
        # is load-bearing rather than an optimisation detail.
        # MODE-DEPENDENT DEFAULT. Every measurement behind this mask came from mosaic runs
        # (FPs enriched 10-20x over TPs in these regions on the HapMap arms). Germline has a
        # different FP profile and a long-standing concordance baseline, so switching it on
        # there by default would move that baseline on no evidence. Unset => mosaic gets the
        # three built-ins, germline gets nothing; an explicit --region-mask always wins.
        _mask_spec = getattr(opt, "region_mask", None)
        if _mask_spec is None:
            _mask_spec = "default" if getattr(opt, "mosaic", False) else "none"
            _why = f"mode default for {'--mosaic' if getattr(opt, 'mosaic', False) else 'germline'}"
        else:
            _why = "explicit --region-mask"
        the_region_mask = region_mask.load(_mask_spec)
        if the_region_mask:
            log(f"region mask ({_why}): {the_region_mask.describe()}")
            for src in the_region_mask.sources:
                log(f"  masking from {src}")
        else:
            log(f"region mask: DISABLED ({_why})")
        type_ctx = []
        for t in types:
            spec = REP_TYPES[t]
            paths = type_paths[t]
            type_dir = sample_root / t
            dirs = {
                "type": str(type_dir),
                "tmp": str(type_dir / "tmp") + "/",
                "tmp_cns": str(type_dir / "tmp" / "cns") + "/",
                "tmp_td": str(type_dir / "tmp" / "transduction") + "/",
                "pub_clip": pub_clip,
                "bam_list": bam_list,
                "control_bam_list": control_bam_list,
            }
            for d in ("type", "tmp", "tmp_cns", "tmp_td"):
                Path(dirs[d]).mkdir(parents=True, exist_ok=True)

            log(f"=== type {t} (rtype {spec['rtype']}): candidate lists ===")
            # A1: write candidate_list_from_clip.txt natively from the scan pass.
            if opt.scan_clip and not opt.dry_run:
                clip_out = str(type_dir / "candidate_list_from_clip.txt")
                nsites, nmasked = scan_adapter.write_clip_candidate_list(
                    realign_tsv, t, clip_out, min_clip=opt.scan_min_clip,
                    region_mask=the_region_mask)
                log(f"{t}:clip (xtea-scan) wrote {nsites} candidate positions -> {clip_out}"
                    + (f" ({nmasked} dropped by --region-mask, "
                       f"{100.0 * nmasked / max(1, nsites + nmasked):.1f}%)" if nmasked else ""))

            # A2: write candidate_list_from_disc.txt natively from the disc-pos stream.
            if opt.scan_disc and not opt.dry_run:
                clip_out = str(type_dir / "candidate_list_from_clip.txt")
                disc_out = str(type_dir / "candidate_list_from_disc.txt")
                iextend, std_var, disc_cutoff, i_rlth, b_with_chr = disc_params
                # b_tumor MUST be threaded through: without it disc_adapter defaults to
                # False and applies master's germline retention formulas to a tumor run
                # (both the effective-disc rule and the raw-disc half-cutoff), leaving
                # its tumor branches dead. Same class of omission as the purity scaling
                # above -- the scan path bypasses master's automatic parameter setup, so
                # every tumor-gated decision has to be handed over explicitly.
                # PEAK_WINDOW 100 -> 30 under --mosaic. Master narrows the same window for
                # BOTH mosaic and somatic (x_TEA_main.py:433-435, a local in the deleted -D
                # branch -- NOT global_values.PEAK_WINDOW, which nothing reads). A 100 bp
                # window merges neighbouring clip peaks; at mosaic's depth that fuses a real
                # low-AF insertion into an adjacent germline peak and the call is lost as a
                # duplicate rather than reported.
                disc_adapter.run_disc_from_scan(
                    clip_out, disc_out, disc_arrays, paths["annotation"], t,
                    iextend, std_var, disc_cutoff, i_rlth, dirs["tmp"],
                    b_with_chr=b_with_chr, b_tumor=opt.tumor,
                    peak_window=(30 if getattr(opt, "mosaic", False) else 100),
                    log=lambda m: log(m))

            type_ctx.append((t, spec, paths, dirs, type_dir))

        # --- Native cns collection: fetch each type's candidate windows with the BAM index,
        # collect + realign, and emit the temp SAMs cns_filter then reuses via the bypass. The
        # candidate windows of all types are unioned and merged, so overlapping ones are read once
        # rather than once per type. ---
        if opt.cns_emit_scan and not opt.dry_run:
            manifest_path = str(sample_root / "cns_manifest.txt")
            with open(manifest_path, "w") as mf:
                for (t, spec, paths, dirs, type_dir) in type_ctx:
                    disc_out = str(type_dir / "candidate_list_from_disc.txt")
                    tc = dirs["tmp_cns"]
                    mf.write(" ".join([t, disc_out, paths["cns"],
                                       tc + "temp_clip.sam", tc + "temp_disc.sam",
                                       tc + "low_mapq_clip_cnt.txt"]) + "\n")
            cns_cmd = [str(Path(opt.scan_bin).resolve()),
                       "--bam", cfg["input"], "--reference", cfg["reference"],
                       "--cns-emit-manifest", manifest_path,
                       "--cns-extnd", str(cns_extnd), "--cns-sample", "0",
                       "--bwa", opt.bwa, "--bwa-threads", str(opt.cores)]
            if cns_read_length:
                cns_cmd += ["--cns-read-length", str(cns_read_length)]
            log(f"=== cns-emit native collection (per-candidate fetch, extnd={cns_extnd}) ===")
            log("  " + " ".join(cns_cmd))
            # The emitted SAMs are what cns_filter reads instead of realigning per candidate;
            # they are the whole point of --cns-emit-scan, and re-emitting them cost broad 27
            # minutes and bcm two hours.
            emit_outputs = []
            for (t, spec, paths, dirs, type_dir) in type_ctx:
                tc = dirs["tmp_cns"]
                emit_outputs += [tc + "temp_clip.sam", tc + "temp_disc.sam"]
            emit_fp = _fingerprint([
                "cns_emit", cfg["input"], cfg["reference"], cns_extnd,
                [(t, str(type_dir / "candidate_list_from_disc.txt"), paths["cns"])
                 for (t, spec, paths, dirs, type_dir) in type_ctx],
                # the mask changes the candidate lists cns_emit reads, so resuming across a
                # mask change must NOT reuse the emitted SAMs
                sorted(the_region_mask.sources), the_region_mask.n_intervals,
                _binary_stamp(Path(opt.scan_bin).resolve())])
            if getattr(opt, "resume", False) and \
                    (phase_is_done(sample_root, "cns_emit", emit_fp, emit_outputs) or
                     adopt_existing_phase(sample_root, "cns_emit", emit_fp, emit_outputs,
                                          getattr(opt, "adopt_existing", False))):
                log(f"[resume] cns_emit: reusing {len(emit_outputs)} emitted SAM(s)")
                record_timing("phase", "cns_emit_RESUMED", 0.0)
            else:
                with phase("cns_emit"):
                    subprocess.run(cns_cmd, check=True)
                mark_phase_done(sample_root, "cns_emit", emit_fp, emit_outputs)
            os.environ["XTEA_CNS_PRECOMPUTED_ALIGN"] = "1"

        # --- Phase 2: per-type steps through postF_hc (cns_filter reuses the emitted SAMs).
        # Phase C move A.2: gene/genotype/gVCF no longer run here -- they run ONCE over the
        # unified locus set after the loop (run_shared_tail), since the type-resolution rule
        # only needs each type's postF_hc output (two_side_tprt_both is set there already).
        # CC.1: when --case-control is on, the per-type chain also runs the two
        # case-control filtering steps, so it stops one step later (after
        # case_control_hc, not postF_hc) -- run_shared_tail's locus merge then
        # reads back that somatic-filtered list instead of the raw germline one. ---
        # Mosaic terminates the per-type chain at its own post-filter, the same way
        # case-control does -- the locus merge downstream reads whichever file is the
        # mode's final per-type raw result.
        if getattr(opt, "mosaic", False):
            type_stop_after = "postF_mosaic"
        elif opt.case_control:
            type_stop_after = "case_control_hc"
        else:
            type_stop_after = "postF_hc"
        # Stage 2c.3: with --family-assign-early the chain runs as TWO waves with a
        # cross-family barrier between them, so family is decided once (over
        # cns_filter's output) instead of three chains later at the postF_hc merge.
        # Off => one wave, the pre-2c chain, byte for byte.
        # One type means nothing to contend with, so the barrier is skipped there too
        # (and --types L1 must stay identical to what it produces today).
        b_assign_early = getattr(opt, "family_assign_early", False)
        b_td_gather_driver = getattr(opt, "td_gather_driver", False)
        # --td-gather-driver consumes the site lists the universe pass writes, so it implies
        # it rather than failing on a missing file.
        b_locus_universe = getattr(opt, "locus_universe", False) or b_td_gather_driver
        waves = downstream_waves(b_assign_early, type_stop_after, len(type_ctx),
                                 locus_universe=b_locus_universe)
        # Both barrier users log their own arrival, since either can open the split.
        if b_assign_early and len(waves) > 1:
            log(f"=== family assignment: EARLY (barrier after cns_filter) ===")
        # How the split is run, once something has opened it. STAGED is the default and
        # is answer-preserving by construction (same work, same order, fewer waits) --
        # see StagedBarrier for which family waits on what and why the join cost ~53s.
        # --global-barrier is the old serialized form, kept because it is the arm the
        # measured numbers in run_td_gather_driver were taken on.
        b_staged = len(waves) > 1 and not getattr(opt, "global_barrier", False)
        if opt.dry_run:
            merged_dir = str(sample_root / "merged")
            rtypes_bitmask = 0
            for t in types:
                rtypes_bitmask |= REP_TYPES[t]["rtype"]
            for wi, (start_from, stop_after) in enumerate(waves):
                if wi and b_locus_universe:
                    print("# barrier: locus universe over each type's "
                          "candidate_disc_filtered_cns.txt + raw_disc "
                          "(downstream_unified.enumerate_family_sites)")
                if wi and b_assign_early:
                    print("# barrier: family assignment over each type's "
                          "candidate_disc_filtered_cns.txt (locus_table.assign_and_filter)")
                for (t, spec, paths, dirs, type_dir) in type_ctx:
                    for name, argv in build_type_steps(t, cfg, paths, dirs, opt,
                                                       stop_after=stop_after,
                                                       start_from=start_from):
                        print(f"# {t}:{name}\npython {XTEA_MAIN} " + " ".join(argv))
            for name, argv in build_shared_tail_steps(
                    cfg, merged_dir, opt, hc_pf_path(merged_dir), f"{merged_dir}/locus_type_map.tsv",
                    str(sample_root / f"{opt.sample}.vcf"), bam_list, rtypes_bitmask):
                print(f"# shared:{name}\npython {XTEA_MAIN} " + " ".join(argv))
        elif opt.parallel_downstream and len(type_ctx) > 1:
            # Phase C speed: run the independent per-type chains concurrently. Each
            # type gets cores//ntypes workers (sum stays within the allocation) and
            # its own step log (a shared handle would interleave). Results are the
            # same as serial -- the chains share nothing after the scan+cns-fan.
            dcores = max(1, opt.cores // len(type_ctx))
            # The transduction step gets its own, wider split. An even cores//ntypes
            # slice is right for every OTHER step, but for transduction it hands a
            # third of the machine to Alu, which has no flank (REP_TYPES["Alu"]
            # ["flank"] is None -- a non-autonomous SINE with no internal Pol II
            # promoter cannot 3'-read-through, so this is biologically correct, not a
            # gap to fill) and therefore finishes that step in ~1.5s while L1 and SVA
            # spend ~25 and ~30 minutes in it. Transduction is ~93% of this whole
            # phase, so the idle third is the largest misallocation in the pipeline.
            #
            # Static, not dynamic: which types have a flank is known from REP_TYPES
            # before anything runs, so this is decided here rather than by
            # redistributing cores at step boundaries (that design was costed at a
            # ~15-20s ceiling and rejected -- a running subprocess cannot widen its
            # own -n mid-flight, which is exactly why deciding up front works).
            #
            # Flankless types keep 1 worker for their no-op. The sum can briefly
            # exceed opt.cores: a flankless type may still be in an EARLIER step
            # (Alu's cns_filter runs to t~116s) while L1/SVA are already in
            # transduction on the wider split. That window is short and bounded by
            # the flankless types' pre-transduction chain, and it is a much better
            # trade than idling those cores for the ~25 minutes that follow.
            td_types = [ctx for ctx in type_ctx if REP_TYPES[ctx[0]].get("flank")]
            dcores_td = max(1, opt.cores // len(td_types)) if td_types else dcores
            log(f"=== Phase 2: {len(type_ctx)} per-type downstreams CONCURRENT "
                f"(n={dcores} each; transduction n={dcores_td} for "
                f"{len(td_types)} flanked type(s), n=1 for the rest) ===")

            def _worker(ctx, start_from, stop_after, mode):
                t = ctx[0]
                nc_td = dcores_td if REP_TYPES[t].get("flank") else 1
                tlog = sample_root / f"{opt.sample}_{t}.steps.log"
                # Append on the second wave: a fresh "w" would truncate wave 1's log.
                with open(tlog, mode) as tfh:
                    return run_type_downstream(ctx, cfg, opt, dcores, tfh,
                                               stop_after=stop_after,
                                               start_from=start_from,
                                               n_cores_transduction=nc_td)

            def _timed_family_assign():
                # "step", not phase, for the same reason as the gather: it runs on a worker
                # thread inside the downstream_parallel phase. (The global path's
                # phase("family_assign") already nests inside phase("downstream_parallel"),
                # so a step row here is the more honest of the two, not a downgrade.)
                with timed("step", "family_assign"):
                    return run_family_assign_barrier(type_ctx, sample_root)

            ctx_by_tag = {ctx[0]: ctx for ctx in type_ctx}
            barrier = StagedBarrier(
                families=list(ctx_by_tag),
                flanked=[ctx[0] for ctx in td_types],
                b_universe=b_locus_universe, b_gather=b_td_gather_driver,
                b_assign=b_assign_early,
                enumerate_fn=lambda t: enumerate_family_universe(
                    ctx_by_tag[t], cfg, opt, sample_root,
                    disc_params[2] if disc_params else None),
                gather_fn=lambda: run_td_gather_driver(
                    type_ctx, cfg, opt, sample_root, cns_extnd, logfh,
                    timing_kind="step"),
                union_fn=lambda slct: write_locus_universe_table(
                    slct, sample_root, cns_extnd),
                assign_fn=_timed_family_assign,
            )

            def _worker_staged(ctx):
                """Both waves on ONE thread, with only this family's own waits between.

                The two-wave form below runs each wave under its own executor, so the
                pool itself is the join. Here the family owns its chain end to end and
                `barrier.arrive` is where it waits -- for the gather if it is flanked,
                for the assignment if that flag is on, and otherwise not at all.

                One log handle for both waves instead of "w" then "a": same bytes, and
                the reopen only ever existed to keep a second executor from truncating
                the first's output.
                """
                t = ctx[0]
                nc_td = dcores_td if REP_TYPES[t].get("flank") else 1
                tlog = sample_root / f"{opt.sample}_{t}.steps.log"
                with open(tlog, "w") as tfh:
                    try:
                        run_type_downstream(ctx, cfg, opt, dcores, tfh,
                                            stop_after=waves[0][1],
                                            start_from=waves[0][0],
                                            n_cores_transduction=nc_td)
                    except BaseException:
                        # Nobody may be left waiting on a family that died -- see
                        # StagedBarrier's FAILURE note.
                        barrier.abort(t)
                        raise
                    barrier.arrive(t)
                    return run_type_downstream(ctx, cfg, opt, dcores, tfh,
                                               stop_after=waves[1][1],
                                               start_from=waves[1][0],
                                               n_cores_transduction=nc_td)

            with phase("downstream_parallel"):
                if b_staged:
                    log("=== barrier: STAGED (per-family; only the assignment joins) ===")
                    if b_td_gather_driver:
                        # Set BEFORE any worker starts, not on the gather thread as the
                        # global form does. Under a join no subprocess exists when the
                        # driver mutates os.environ; here families are spawning steps
                        # throughout, and mutating environ concurrently with fork/exec is
                        # the one race staging introduces. Setting it early is free
                        # because the step's check is `env AND both files present`
                        # (x_transduction.py's b_td_precomputed) -- a family the gather
                        # skipped, or a gather that produced no entries at all, still
                        # falls back to gathering for itself.
                        os.environ[downstream_unified.ENV_TD_PRECOMPUTED] = "1"
                    if b_locus_universe:
                        # The same header the all-at-once form prints, so a run's log
                        # still has one place the universe's per-family lines hang off.
                        log("=== locus universe (cns_filter -> transduction, "
                            "driver-side, per family as each lands) ===")
                    with concurrent.futures.ThreadPoolExecutor(
                            max_workers=len(type_ctx)) as ex:
                        futs = {ex.submit(_worker_staged, ctx): ctx[0] for ctx in type_ctx}
                        for fut in concurrent.futures.as_completed(futs):
                            t = futs[fut]
                            v = fut.result()   # re-raises any step failure (SystemExit) here
                            if v:
                                hc_pf_paths[t] = v
                else:
                    if len(waves) > 1:
                        log("=== barrier: GLOBAL (--global-barrier; every family waits "
                            "for the slowest cns_filter) ===")
                    for wi, (start_from, stop_after) in enumerate(waves):
                        if wi:
                            # The barrier sits between the waves, so both users run with
                            # every family's cns_filter finished and none of them past it.
                            # Universe first: it reads the cns lists, and the assignment
                            # barrier REWRITES them to each family's winners, so the other
                            # order would enumerate a post-assignment universe and silently
                            # change what the number means.
                            if b_locus_universe:
                                with phase("locus_universe"):
                                    run_locus_universe(type_ctx, cfg, opt, sample_root,
                                                       disc_params[2] if disc_params else None,
                                                       cns_extnd)
                            if b_td_gather_driver:
                                run_td_gather_driver(type_ctx, cfg, opt, sample_root,
                                                     cns_extnd, logfh)
                            if b_assign_early:
                                with phase("family_assign"):
                                    run_family_assign_barrier(type_ctx, sample_root)
                        mode = "w" if wi == 0 else "a"
                        with concurrent.futures.ThreadPoolExecutor(
                                max_workers=len(type_ctx)) as ex:
                            futs = {ex.submit(_worker, ctx, start_from, stop_after, mode): ctx[0]
                                    for ctx in type_ctx}
                            for fut in concurrent.futures.as_completed(futs):
                                t = futs[fut]
                                v = fut.result()   # re-raises any step failure (SystemExit)
                                if v:
                                    hc_pf_paths[t] = v
        else:
            with phase("downstream_serial"):
                for wi, (start_from, stop_after) in enumerate(waves):
                    if wi:
                        # Universe before assignment -- see the concurrent branch.
                        if b_locus_universe:
                            with phase("locus_universe"):
                                run_locus_universe(type_ctx, cfg, opt, sample_root,
                                                   disc_params[2] if disc_params else None,
                                                   cns_extnd)
                        if b_td_gather_driver:
                            run_td_gather_driver(type_ctx, cfg, opt, sample_root,
                                                 cns_extnd, logfh)
                        if b_assign_early:
                            with phase("family_assign"):
                                run_family_assign_barrier(type_ctx, sample_root)
                    for ctx in type_ctx:
                        t = ctx[0]
                        v = run_type_downstream(ctx, cfg, opt, opt.cores, logfh,
                                                stop_after=stop_after,
                                                start_from=start_from)
                        if v:
                            hc_pf_paths[t] = v

        # --- Phase 3: shared tail, run ONCE (Phase C move A.2) ---
        if not opt.dry_run and hc_pf_paths:
            with phase("shared_tail"):
                run_shared_tail(cfg, opt, sample_root, hc_pf_paths, bam_list, logfh)

    if not opt.dry_run and not opt.keep_intermediates:
        cleanup_successful_run(sample_root, opt.sample)
    log("run complete")


if __name__ == "__main__":
    main()
