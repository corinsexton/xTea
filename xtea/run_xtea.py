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
    shared across types via one `pub_clip/` dir (as master does).
  * All filter/threshold knobs are exposed on the CLI with defaults pinned to the
    golden run (`--nclip 4 --cr 2 --nd 5 --nfclip 4 --nfdisc 5 --flklen 3000`), so
    the defaults reproduce the golden and every change is opt-in.

Next milestones (A1..): replace each `x_TEA_main.py` subprocess with a native
in-process call, validating per-step against the frozen golden intermediates
(`candidate_list_from_clip.txt`, `candidate_list_from_disc.txt`, ...), then
integrate/merge steps where that is faster (results stay the contract).

This file intentionally does NOT reimplement detection logic yet -- it is the
always-green scaffold we refactor inward.
"""

import argparse
import concurrent.futures
import os
import subprocess
import sys
from pathlib import Path

import scan_adapter
import disc_adapter
import locus_table

HERE = Path(__file__).resolve().parent
XTEA_MAIN = str(HERE / "x_TEA_main.py")
GNTP_MODEL = str(HERE / "genotyping" / "DF21_model_1_2")
DEFAULT_SCAN_BIN = str(HERE.parent / "src" / "csrc" / "xtea_scan")

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


def log(msg):
    print(f"[xtea2] {msg}", flush=True)


def resolve_annot(rep_lib, genome, rel):
    if rel is None:
        return "null"
    return str(Path(rep_lib) / rel.format(g=genome))


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


def run_step(name, argv, logfh):
    """Run one pipeline step as a subprocess to x_TEA_main.py, in this process's env.

    cwd is the package dir so x_TEA_main.py's flat imports resolve. Fails loud.
    """
    cmd = [sys.executable, XTEA_MAIN] + argv
    log(f"step: {name}")
    logfh.write(f"\n=== STEP {name} ===\n{' '.join(cmd)}\n")
    logfh.flush()
    r = subprocess.run(cmd, cwd=str(HERE), stdout=logfh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.exit(f"ERROR: step '{name}' failed (exit {r.returncode}); see log.")


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


def build_type_steps(t, cfg, paths, dirs, opt, skip_clip=False, skip_disc=False,
                     n_cores=None, stop_after=None):
    """Return the ordered list of (step_name, argv) for one MEI type.

    Mirrors the golden `run_xTEA_pipeline.sh` step chain exactly, parameterised.
    When skip_clip is True the clip step is omitted (the caller produced
    candidate_list_from_clip.txt natively via the xtea-scan adapter). When
    skip_disc is True the disc step is omitted too (the caller produced
    candidate_list_from_disc.txt natively via the disc adapter).

    n_cores overrides the per-step worker count (`-n`). When the three types'
    downstreams run concurrently (--parallel-downstream) each type gets a slice
    of the total (cores // ntypes) so the sum stays within the allocation; None
    means use opt.cores (the serial default, unchanged).

    stop_after truncates the returned list right after the named step (Phase C
    move A.2 uses stop_after="postF_hc": gene/genotype/gVCF run ONCE over the
    unified locus set instead of once per type -- see run_shared_tail).
    """
    P = dirs["type"] + "/"          # per-type output dir (trailing slash)
    bam_list = dirs["bam_list"]
    ref = cfg["reference"]
    rtype = str(paths["rtype"])
    fflank = paths["flank"]
    nc = str(opt.cores if n_cores is None else n_cores)
    steps = []
    # SVA needs --sva on the clip/disc/cns steps (turns on IS_CALL_SVA ->
    # SVA_ANNOTATION_EXTND=200 for the mate-in-repeat test), exactly as the golden
    # run_xTEA_pipeline.sh does. Omitting it uses extend 0, which silently drops
    # ~18% of SVA at whole-genome scale (chr1 happened to tolerate it). The golden
    # puts --sva on -C/-D/-N only (not transduction/sibling/postF), so we match that.
    sva = ["--sva"] if t == "SVA" else []
    # CC.0: master threads --purity <p> --tumor through clip/disc/cns_filter/
    # transduction/sibling, and --tumor alone (no --purity) through both postF
    # steps (confirmed against a real executed run_xTEA_pipeline.sh, not just
    # gnrt_pipeline_local_v38.py source). Empty lists when opt.tumor is False,
    # so this is a no-op for every existing (non-tumor) caller.
    tumor_purity_args = (["--purity", str(opt.purity), "--tumor"] if opt.tumor else [])
    tumor_only_args = (["--tumor"] if opt.tumor else [])

    clip_out = f"{P}candidate_list_from_clip.txt"
    disc_out = f"{P}candidate_list_from_disc.txt"
    raw_disc = f"{disc_out}.clip_sites_raw_disc.txt"
    cns_out  = f"{P}candidate_disc_filtered_cns.txt"
    cns2_out = f"{P}candidate_disc_filtered_cns2.txt"
    sib_out  = f"{P}candidate_sibling_transduction2.txt"
    pf_out   = f"{P}candidate_disc_filtered_cns_post_filtering.txt"
    hc_pf    = hc_pf_path(dirs["type"])
    gene_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"
    gntp_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"

    # 1. clip
    if not skip_clip:
        steps.append(("clip", [
            "-C"] + sva + ["-i", bam_list, "--lc", str(opt.nclip), "--rc", str(opt.nclip),
            "--cr", str(opt.cr), "-r", paths["copies_with_flank"], "-a", paths["annotation"],
            "--cns", paths["cns"], "--ref", ref, "-p", dirs["tmp"],
            "-o", clip_out, "-n", nc, "--cp", dirs["pub_clip"]] + tumor_purity_args))
    # 2. disc
    if not skip_disc:
        steps.append(("disc", [
            "-D"] + sva + ["-i", clip_out, "--nd", str(opt.nd), "--ref", ref,
            "-a", paths["annotation"], "-b", bam_list, "-p", dirs["tmp"],
            "-o", disc_out, "-n", nc] + tumor_purity_args))
    # 3. cns filter
    steps.append(("cns_filter", [
        "-N"] + sva + ["--cr", str(opt.nfclip), "--nd", str(opt.nfdisc), "-b", bam_list,
        "-p", dirs["tmp_cns"], "--fflank", fflank, "--flklen", str(opt.flklen),
        "-n", nc, "-i", disc_out, "-r", paths["cns"], "--ref", ref,
        "-a", paths["annotation"], "-o", cns_out] + tumor_purity_args))
    # 4. transduction
    steps.append(("transduction", [
        "--transduction", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", nc, "-i", cns_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "-o", cns2_out] + tumor_purity_args))
    # 5. sibling
    steps.append(("sibling", [
        "--sibling", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", nc, "-i", cns2_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "--blacklist", cfg["blacklist"],
        "-o", sib_out] + tumor_purity_args))
    # 6. post-filter (round 1: all)
    steps.append(("postF_all", [
        "--postF", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", nc,
        "-i", cns2_out, "-a", paths["annotation1"], "-o", pf_out] + tumor_only_args))
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
    # 9. gene annotation
    steps.append(("gene", [
        "--gene", "-a", cfg["gene_gff3"], "-i", hc_pf, "-n", nc,
        "-o", gene_out]))
    # 10. genotype (DF21 deep-forest; genotype not part of the concordance gate)
    steps.append(("genotype", [
        "--gntp_classify", "-i", gene_out, "-n", "1", "--model", cfg["model"],
        "-o", gntp_out]))
    # 11. gVCF
    steps.append(("gvcf", [
        "--gVCF", "-i", gntp_out, "-o", P, "-b", bam_list, "--ref", ref,
        "--rtype", rtype]))

    if stop_after is not None:
        names = [name for name, _ in steps]
        steps = steps[:names.index(stop_after) + 1]
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
    if not getattr(opt, "case_control", False):
        steps.append(("genotype", ["--gntp_classify", "-i", gene_out, "-n", "1",
                                    "--model", cfg["model"], "-o", gntp_out]))
        gvcf_in = gntp_out
    steps.append(("gvcf", ["--gVCF", "-i", gvcf_in, "-o", P, "-b", bam_list, "--ref", ref,
                           "--type-map", sf_type_map, "--merged-vcf-out", sf_merged_vcf,
                           "--rtype", str(rtypes_bitmask)]))
    return steps


def run_type_downstream(ctx, cfg, opt, n_cores, logfh, stop_after=None):
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
    """
    (t, spec, paths, dirs, type_dir) = ctx
    log(f"=== type {t} (rtype {spec['rtype']}): downstream (n={n_cores}) ===")
    for name, argv in build_type_steps(t, cfg, paths, dirs, opt,
                                       skip_clip=opt.scan_clip,
                                       skip_disc=opt.scan_disc,
                                       n_cores=n_cores, stop_after=stop_after):
        run_step(f"{t}:{name}", argv, logfh)
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
    log(f"=== merge per-type postF_hc -> unified locus set (locus table, "
        f"win={locus_table.MERGE_WINDOW}) ===")
    loci = locus_table.assign_and_write(hc_pf_paths, sf_merged_hcpf, sf_type_map)
    log(f"unified locus set: {len(loci)} calls -> {sf_merged_hcpf}")

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


def main():
    ap = argparse.ArgumentParser(
        prog="xtea2",
        description="xTea2: BAM/CRAM -> MEI VCF, single self-contained process "
                    "(Milestone A0: golden-reproducing wrapper).")
    # required inputs
    ap.add_argument("-i", "--input", required=True,
                    help="input BAM or CRAM (single sample)")
    ap.add_argument("-s", "--sample", required=True, help="sample identifier")
    ap.add_argument("-o", "--output-dir", required=True, help="output directory")
    ap.add_argument("--reference", required=True,
                    help="genome FASTA (required; also used to decode CRAM)")
    ap.add_argument("--rep-lib-dir", required=True,
                    help="path to rep_lib_annotation/ directory")
    ap.add_argument("--gene-gff3", required=True, help="gene annotation .gff3")
    ap.add_argument("--blacklist", default="null",
                    help="blacklist BED (centromere/dup); default: null (off)")
    ap.add_argument("--model", default=GNTP_MODEL,
                    help="genotyping model dir (default: bundled DF21)")
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
    # exposed filter/threshold knobs (defaults pinned to the golden run)
    ap.add_argument("--nclip", type=int, default=4,
                    help="min clipped reads per side to seed a candidate (golden 4)")
    ap.add_argument("--cr", type=int, default=2,
                    help="min clipped parts landing in repeat, clip step (golden 2)")
    ap.add_argument("--nd", type=int, default=5,
                    help="min discordant reads, disc step (golden 5)")
    ap.add_argument("--nfclip", type=int, default=4,
                    help="min clipped reads at cns/final filter (golden 4)")
    ap.add_argument("--nfdisc", type=int, default=5,
                    help="min discordant reads at cns/final filter (golden 5)")
    ap.add_argument("--flklen", type=int, default=3000,
                    help="flank length for realign/transduction (golden 3000)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the per-step commands and exit (no execution)")
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
                         "--case-control)")
    # A1 refactor-inward: replace master's clip step with the C++ xtea-scan pass.
    ap.add_argument("--scan-clip", action="store_true",
                    help="use the C++ xtea-scan (single-pass clip+realign) for the "
                         "clip step instead of x_TEA_main.py -C")
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
                    help="use the C++ xtea-scan --disc-pos stream (collected in the same "
                         "pass) for the disc step instead of x_TEA_main.py -D. Implies "
                         "--scan-clip.")
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
    ap.add_argument("--cns-emit-scan", action="store_true",
                    help="Phase C move B.1: have the cns_filter (-N) step emit its "
                         "temp_clip.sam/temp_disc.sam/low_mapq_clip_cnt.txt from a native xtea-scan "
                         "pass (reusing the step's own iextnd) and reuse them via the bypass, instead "
                         "of re-reading the BAM + re-running bwa per candidate. Uses --scan-bin/--bwa.")
    ap.add_argument("--col5-families", default="alu,sva",
                    help="which repeat families' annotations count toward col5 "
                         "(mate-in-repeat clip count), comma-separated from {l1,alu,sva}. "
                         "col5 is the # clipped reads at a site whose mate lands in ANY of "
                         "these families' repeats (union), emitted identically for all types. "
                         "Default 'alu,sva' (SVA + Alu) -- the biologically tightest set: col5 "
                         "is consumed only by the SVA cns path, and SVA is a composite element "
                         "with an Alu-derived domain, so mate-in-Alu is genuine SVA evidence. "
                         "'l1,alu,sva' = full union; 'sva' = per-family-strict. See spec §4 DC1.")
    opt = ap.parse_args()
    opt._col5_fams = {s.strip().lower() for s in opt.col5_families.split(",") if s.strip()}
    _bad = opt._col5_fams - {"l1", "alu", "sva"}
    if _bad:
        sys.exit(f"ERROR: --col5-families: unknown families {sorted(_bad)} "
                 f"(choose from l1,alu,sva)")

    # --scan-disc requires the scan pass (and its clip candidates) to exist.
    if opt.scan_disc:
        opt.scan_clip = True

    if opt.case_control and not opt.normal_bam:
        sys.exit("ERROR: --case-control requires --normal-bam")

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

    detect_format_and_index(opt.input, opt.reference)

    sample_root = Path(opt.output_dir).resolve() / opt.sample
    sample_root.mkdir(parents=True, exist_ok=True)
    pub_clip = str(sample_root / "pub_clip") + "/"
    Path(pub_clip).mkdir(parents=True, exist_ok=True)

    # one bam_list.txt for the whole run (single sample, illumina)
    bam_list = str(sample_root / "bam_list.txt")
    with open(bam_list, "w") as fh:
        fh.write(f"{Path(opt.input).resolve()}\tillumina\n")

    # CC.1: control bam_list, built the same way as the tumor one -- every
    # consumer on the case-control path only ever reads column 0 (verified by
    # reading XClipDiscFilter.__init__/XGenotyper.call_genotype/
    # X_BasicInfo.collect_basic_info_samples), so the two-column convention here
    # is a style choice, not a format requirement.
    control_bam_list = None
    if opt.case_control:
        detect_format_and_index(opt.normal_bam, opt.reference)
        control_bam_list = str(sample_root / "control_bam_list.txt")
        with open(control_bam_list, "w") as fh:
            fh.write(f"{Path(opt.normal_bam).resolve()}\tillumina\n")

    cfg = {
        "input": str(Path(opt.input).resolve()),
        "reference": str(Path(opt.reference).resolve()),
        "blacklist": ("null" if opt.blacklist in (None, "null", "NULL")
                      else str(Path(opt.blacklist).resolve())),
        "gene_gff3": str(Path(opt.gene_gff3).resolve()),
        "model": str(Path(opt.model).resolve()),
    }

    logpath = sample_root / f"{opt.sample}_xtea2.log"
    log(f"run start: sample={opt.sample} types={types} -> {sample_root}")
    log(f"full step log: {logpath}")

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

    # Phase C move B.1 (fold-in): when --cns-emit-scan is on, a SINGLE xtea-scan pass collects the
    # cns clip+disc reads ONCE for all types (2 BAM passes total, not 6) and emits each type's
    # temp_clip.sam/temp_disc.sam/low_mapq_clip_cnt.txt; we then set XTEA_CNS_PRECOMPUTED_ALIGN so
    # every cns_filter (-N) step reuses them via the bypass. Requires --scan-disc (for the insert-size
    # window). The emit happens after candidate lists are built (see the two-phase type loop below).
    if opt.cns_emit_scan and not opt.scan_disc:
        sys.exit("ERROR: --cns-emit-scan requires --scan-disc (needs the insert-size window)")

    hc_pf_paths = {}  # type -> postF_hc path (Phase C move A.2: gene/genotype/gVCF run once, after the loop)
    with open(logpath, "w") as logfh:
        # A1/A2: one xtea-scan pass produces clip candidates (--realign-out) and,
        # when --scan-disc is set, the discordant-read stream (--disc-pos) for ALL
        # requested types -- all BAM reading happens here, once.
        realign_tsv = None
        disc_arrays = None
        disc_params = None       # (iextend, std_var, disc_cutoff, b_with_chr)
        cns_extnd = None         # cns collection window = int(mean_is + 3*std) (== x_TEA_main max_is)
        if opt.scan_clip and not opt.dry_run:
            scan_bin = str(Path(opt.scan_bin).resolve())
            if not Path(scan_bin).exists():
                sys.exit(f"ERROR: --scan-clip: xtea_scan binary not found at {scan_bin}")
            ref_key = ("copies_with_flank" if opt.scan_realign_ref == "copies_with_flank"
                       else "cns")
            cns_by_family = {scan_adapter.TYPE_TO_FAMILY[t]: type_paths[t][ref_key] for t in types}
            realign_tsv = str(sample_root / "scan_clip_realign.tsv")
            disc_pos_out = str(sample_root / "scan_disc_pos.tsv") if opt.scan_disc else None
            # col5 (mate-in-repeat clip count): convert each family's .out annotation
            # to an extended BED (SVA +200, L1/Alu 0 -- same as the clip step's
            # annotation) and let xtea-scan compute col5 natively via htslib regidx.
            # col5 (mate-in-repeat) counts a clipped read's mate landing in ANY of the
            # families selected by --col5-families (union). Only generate/pass those
            # families' BEDs; xtea-scan ORs whichever BEDs it receives (spec §4 DC1).
            _type_to_col5key = {"L1": "l1", "Alu": "alu", "SVA": "sva"}
            bed_by_family = {}
            for t in types:
                if _type_to_col5key[t] not in opt._col5_fams:
                    continue
                fam = scan_adapter.TYPE_TO_FAMILY[t]
                extnd = scan_adapter.SVA_ANNOTATION_EXTND if t == "SVA" else 0
                bed = str(sample_root / f"rep_{fam}.bed")
                n_iv = scan_adapter.out_to_bed(type_paths[t]["annotation"], bed, extnd)
                bed_by_family[fam] = bed
                log(f"col5 BED {fam}: {n_iv} intervals (extnd={extnd}) -> {bed}")
            log(f"col5 union families: {sorted(opt._col5_fams)}")
            log(f"=== xtea-scan single pass (clip+realign {list(cns_by_family)}"
                f"{'+disc-pos' if opt.scan_disc else ''}+col5) ===")
            # full fold-in: spill the cns collection during this same detection pass.
            cns_collect = None
            if opt.cns_emit_scan:
                cns_collect = {
                    "clip":   str(sample_root / "cns_spill_clip.txt"),
                    "anchor": str(sample_root / "cns_spill_anchor.txt"),
                    "lowq":   str(sample_root / "cns_spill_lowq.txt"),
                }
            scan_adapter.run_scan_realign(
                scan_bin, cfg["input"], cfg["reference"], cns_by_family, realign_tsv,
                threads=opt.cores, bwa=opt.bwa, logfh=logfh,
                disc_pos_out=disc_pos_out,
                disc_threshold=(opt.scan_disc_threshold if opt.scan_disc else None),
                bed_by_family=bed_by_family, cns_collect=cns_collect)

            if opt.scan_disc:
                # insert-size params (read length, mean IS, std) + coverage-derived
                # disc cutoff -- computed once, exactly as master's disc step does.
                from x_basic_info import X_BasicInfo
                from x_parameter import Parameters
                xbi = X_BasicInfo(str(sample_root / "tmp_basic") + "/", opt.cores, cfg["reference"])
                Path(str(sample_root / "tmp_basic")).mkdir(parents=True, exist_ok=True)
                rcd = xbi.get_cov_is_rlth(bam_list, cfg["reference"], 500, True)
                f_cov, rlth, mean_is, std_var = rcd[0], rcd[1], rcd[2], rcd[3]
                disc_cutoff = Parameters().get_par_by_cov(f_cov)[1]
                # bam "chr" prefix convention for the annotation mate-in-rep test,
                # AND a robust read length sampled directly from the BAM. master's
                # get_cov_is_rlth reports rlth=0 here (a long-standing bug: "Read
                # length is: 0.0"); master's own disc step doesn't need it (its
                # bamfile.fetch returns overlapping reads natively), but our disc-pos
                # emulation widens the window by one read length to recover reads
                # that span a breakpoint -- so we must have a real value.
                import pysam
                _bf = pysam.AlignmentFile(cfg["input"], "rb", reference_filename=cfg["reference"])
                b_with_chr = any(r.startswith("chr") for r in _bf.references[:5])
                _rls = []
                for _a in _bf.fetch(until_eof=True):
                    if _a.is_unmapped or _a.is_secondary or _a.is_supplementary:
                        continue
                    if _a.query_length:
                        _rls.append(_a.query_length)
                    if len(_rls) >= 2000:
                        break
                _bf.close()
                _rls.sort()
                sampled_rlth = _rls[len(_rls) // 2] if _rls else 150   # median; 150 fallback
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
                log(f"disc params: cov={f_cov:.1f} rlth(master)={rlth} rlth(overlap)={overlap_rlth} "
                    f"mean_is={mean_is:.0f} std={std_var:.0f} iextend={iextend} "
                    f"disc_cutoff={disc_cutoff}")
                log(f"=== loading disc-pos stream {disc_pos_out} ===")
                disc_arrays = disc_adapter.load_disc_pos(disc_pos_out)
                log(f"disc-pos loaded: {sum(len(v[0]) for v in disc_arrays.values())} records "
                    f"across {len(disc_arrays)} contigs")

        # --- Phase 1: build every type's candidate lists (needed before the single cns collect) ---
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
                nsites = scan_adapter.write_clip_candidate_list(
                    realign_tsv, t, clip_out, min_clip=opt.scan_min_clip)
                log(f"{t}:clip (xtea-scan) wrote {nsites} candidate positions -> {clip_out}")

            # A2: write candidate_list_from_disc.txt natively from the disc-pos stream.
            if opt.scan_disc and not opt.dry_run:
                clip_out = str(type_dir / "candidate_list_from_clip.txt")
                disc_out = str(type_dir / "candidate_list_from_disc.txt")
                iextend, std_var, disc_cutoff, i_rlth, b_with_chr = disc_params
                disc_adapter.run_disc_from_scan(
                    clip_out, disc_out, disc_arrays, paths["annotation"], t,
                    iextend, std_var, disc_cutoff, i_rlth, dirs["tmp"],
                    b_with_chr=b_with_chr, log=lambda m: log(m))

            type_ctx.append((t, spec, paths, dirs, type_dir))

        # --- Fold-in FAN: the detection pass already SPILLED the cns collection (cns_spill_*.txt);
        # now fan it to each type's candidates + realign, in ONE more BAM pass (the mate-gather). So
        # detection+cns cost 2 whole-genome passes total. cns_filter then reuses the SAMs via bypass. ---
        if opt.cns_emit_scan and not opt.dry_run:
            manifest_path = str(sample_root / "cns_fan_manifest.txt")
            with open(manifest_path, "w") as mf:
                for (t, spec, paths, dirs, type_dir) in type_ctx:
                    disc_out = str(type_dir / "candidate_list_from_disc.txt")
                    tc = dirs["tmp_cns"]
                    mf.write(" ".join([t, disc_out, paths["cns"],
                                       tc + "temp_clip.sam", tc + "temp_disc.sam",
                                       tc + "low_mapq_clip_cnt.txt"]) + "\n")
            fan_cmd = [str(Path(opt.scan_bin).resolve()),
                       "--bam", cfg["input"], "--reference", cfg["reference"],
                       "--cns-fan-manifest", manifest_path,
                       "--cns-fan-clip", str(sample_root / "cns_spill_clip.txt"),
                       "--cns-fan-anchor", str(sample_root / "cns_spill_anchor.txt"),
                       "--cns-fan-lowq", str(sample_root / "cns_spill_lowq.txt"),
                       "--cns-extnd", str(cns_extnd), "--cns-sample", "0",
                       "--bwa", opt.bwa, "--bwa-threads", str(opt.cores)]
            log(f"=== cns-emit fold-in FAN (from detection-pass spill, extnd={cns_extnd}) ===")
            log("  " + " ".join(fan_cmd))
            subprocess.run(fan_cmd, check=True)
            os.environ["XTEA_CNS_PRECOMPUTED_ALIGN"] = "1"

        # --- Phase 2: per-type steps through postF_hc (cns_filter reuses the emitted SAMs).
        # Phase C move A.2: gene/genotype/gVCF no longer run here -- they run ONCE over the
        # unified locus set after the loop (run_shared_tail), since the type-resolution rule
        # only needs each type's postF_hc output (two_side_tprt_both is set there already).
        # CC.1: when --case-control is on, the per-type chain also runs the two
        # case-control filtering steps, so it stops one step later (after
        # case_control_hc, not postF_hc) -- run_shared_tail's locus merge then
        # reads back that somatic-filtered list instead of the raw germline one. ---
        type_stop_after = "case_control_hc" if opt.case_control else "postF_hc"
        if opt.dry_run:
            merged_dir = str(sample_root / "merged")
            rtypes_bitmask = 0
            for t in types:
                rtypes_bitmask |= REP_TYPES[t]["rtype"]
            for (t, spec, paths, dirs, type_dir) in type_ctx:
                for name, argv in build_type_steps(t, cfg, paths, dirs, opt,
                                                   skip_clip=opt.scan_clip,
                                                   skip_disc=opt.scan_disc,
                                                   stop_after=type_stop_after):
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
            log(f"=== Phase 2: {len(type_ctx)} per-type downstreams CONCURRENT "
                f"(n={dcores} each) ===")

            def _worker(ctx):
                t = ctx[0]
                tlog = sample_root / f"{opt.sample}_{t}.steps.log"
                with open(tlog, "w") as tfh:
                    return run_type_downstream(ctx, cfg, opt, dcores, tfh, stop_after=type_stop_after)

            with concurrent.futures.ThreadPoolExecutor(max_workers=len(type_ctx)) as ex:
                futs = {ex.submit(_worker, ctx): ctx[0] for ctx in type_ctx}
                for fut in concurrent.futures.as_completed(futs):
                    t = futs[fut]
                    v = fut.result()   # re-raises any step failure (SystemExit) here
                    if v:
                        hc_pf_paths[t] = v
        else:
            for ctx in type_ctx:
                t = ctx[0]
                v = run_type_downstream(ctx, cfg, opt, opt.cores, logfh, stop_after=type_stop_after)
                if v:
                    hc_pf_paths[t] = v

        # --- Phase 3: shared tail, run ONCE (Phase C move A.2) ---
        if not opt.dry_run and hc_pf_paths:
            run_shared_tail(cfg, opt, sample_root, hc_pf_paths, bam_list, logfh)

    log("run complete")


if __name__ == "__main__":
    main()
