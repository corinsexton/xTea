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
import os
import subprocess
import sys
from pathlib import Path

import scan_adapter
import disc_adapter

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


def build_type_steps(t, cfg, paths, dirs, opt, skip_clip=False, skip_disc=False):
    """Return the ordered list of (step_name, argv) for one MEI type.

    Mirrors the golden `run_xTEA_pipeline.sh` step chain exactly, parameterised.
    When skip_clip is True the clip step is omitted (the caller produced
    candidate_list_from_clip.txt natively via the xtea-scan adapter). When
    skip_disc is True the disc step is omitted too (the caller produced
    candidate_list_from_disc.txt natively via the disc adapter).
    """
    P = dirs["type"] + "/"          # per-type output dir (trailing slash)
    bam_list = dirs["bam_list"]
    ref = cfg["reference"]
    rtype = str(paths["rtype"])
    fflank = paths["flank"]
    steps = []
    # SVA needs --sva on the clip/disc/cns steps (turns on IS_CALL_SVA ->
    # SVA_ANNOTATION_EXTND=200 for the mate-in-repeat test), exactly as the golden
    # run_xTEA_pipeline.sh does. Omitting it uses extend 0, which silently drops
    # ~18% of SVA at whole-genome scale (chr1 happened to tolerate it). The golden
    # puts --sva on -C/-D/-N only (not transduction/sibling/postF), so we match that.
    sva = ["--sva"] if t == "SVA" else []

    clip_out = f"{P}candidate_list_from_clip.txt"
    disc_out = f"{P}candidate_list_from_disc.txt"
    raw_disc = f"{disc_out}.clip_sites_raw_disc.txt"
    cns_out  = f"{P}candidate_disc_filtered_cns.txt"
    cns2_out = f"{P}candidate_disc_filtered_cns2.txt"
    sib_out  = f"{P}candidate_sibling_transduction2.txt"
    pf_out   = f"{P}candidate_disc_filtered_cns_post_filtering.txt"
    hc_pf    = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"
    gene_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"
    gntp_out = f"{P}candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"

    # 1. clip
    if not skip_clip:
        steps.append(("clip", [
            "-C"] + sva + ["-i", bam_list, "--lc", str(opt.nclip), "--rc", str(opt.nclip),
            "--cr", str(opt.cr), "-r", paths["copies_with_flank"], "-a", paths["annotation"],
            "--cns", paths["cns"], "--ref", ref, "-p", dirs["tmp"],
            "-o", clip_out, "-n", str(opt.cores), "--cp", dirs["pub_clip"]]))
    # 2. disc
    if not skip_disc:
        steps.append(("disc", [
            "-D"] + sva + ["-i", clip_out, "--nd", str(opt.nd), "--ref", ref,
            "-a", paths["annotation"], "-b", bam_list, "-p", dirs["tmp"],
            "-o", disc_out, "-n", str(opt.cores)]))
    # 3. cns filter
    steps.append(("cns_filter", [
        "-N"] + sva + ["--cr", str(opt.nfclip), "--nd", str(opt.nfdisc), "-b", bam_list,
        "-p", dirs["tmp_cns"], "--fflank", fflank, "--flklen", str(opt.flklen),
        "-n", str(opt.cores), "-i", disc_out, "-r", paths["cns"], "--ref", ref,
        "-a", paths["annotation"], "-o", cns_out]))
    # 4. transduction
    steps.append(("transduction", [
        "--transduction", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", str(opt.cores), "-i", cns_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "-o", cns2_out]))
    # 5. sibling
    steps.append(("sibling", [
        "--sibling", "--cr", str(opt.nfclip), "--nd", str(opt.nfdisc),
        "-b", bam_list, "-p", dirs["tmp_td"], "--fflank", fflank,
        "--flklen", str(opt.flklen), "-n", str(opt.cores), "-i", cns2_out,
        "-r", paths["cns"], "--ref", ref, "--input2", raw_disc,
        "--rtype", rtype, "-a", paths["annotation1"], "--blacklist", cfg["blacklist"],
        "-o", sib_out]))
    # 6. post-filter (round 1: all)
    steps.append(("postF_all", [
        "--postF", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", str(opt.cores),
        "-i", cns2_out, "-a", paths["annotation1"], "-o", pf_out]))
    # 7. post-filter (round 2: high-confident, with blacklist)
    steps.append(("postF_hc", [
        "--postF", "--rtype", rtype, "-p", dirs["tmp_cns"], "-n", str(opt.cores),
        "-i", f"{cns2_out}.high_confident", "-a", paths["annotation1"],
        "--blacklist", cfg["blacklist"], "-o", hc_pf]))
    # 8. gene annotation
    steps.append(("gene", [
        "--gene", "-a", cfg["gene_gff3"], "-i", hc_pf, "-n", str(opt.cores),
        "-o", gene_out]))
    # 9. genotype (DF21 deep-forest; genotype not part of the concordance gate)
    steps.append(("genotype", [
        "--gntp_classify", "-i", gene_out, "-n", "1", "--model", cfg["model"],
        "-o", gntp_out]))
    # 10. gVCF
    steps.append(("gvcf", [
        "--gVCF", "-i", gntp_out, "-o", P, "-b", bam_list, "--ref", ref,
        "--rtype", rtype]))
    return steps


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

    produced_vcfs = []
    with open(logpath, "w") as logfh:
        # A1/A2: one xtea-scan pass produces clip candidates (--realign-out) and,
        # when --scan-disc is set, the discordant-read stream (--disc-pos) for ALL
        # requested types -- all BAM reading happens here, once.
        realign_tsv = None
        disc_arrays = None
        disc_params = None       # (iextend, std_var, disc_cutoff, b_with_chr)
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
            scan_adapter.run_scan_realign(
                scan_bin, cfg["input"], cfg["reference"], cns_by_family, realign_tsv,
                threads=opt.cores, bwa=opt.bwa, logfh=logfh,
                disc_pos_out=disc_pos_out,
                disc_threshold=(opt.scan_disc_threshold if opt.scan_disc else None),
                bed_by_family=bed_by_family)

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
                log(f"disc params: cov={f_cov:.1f} rlth(master)={rlth} rlth(overlap)={overlap_rlth} "
                    f"mean_is={mean_is:.0f} std={std_var:.0f} iextend={iextend} "
                    f"disc_cutoff={disc_cutoff}")
                log(f"=== loading disc-pos stream {disc_pos_out} ===")
                disc_arrays = disc_adapter.load_disc_pos(disc_pos_out)
                log(f"disc-pos loaded: {sum(len(v[0]) for v in disc_arrays.values())} records "
                    f"across {len(disc_arrays)} contigs")

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
            }
            for d in ("type", "tmp", "tmp_cns", "tmp_td"):
                Path(dirs[d]).mkdir(parents=True, exist_ok=True)

            log(f"=== type {t} (rtype {spec['rtype']}) ===")
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

            for name, argv in build_type_steps(t, cfg, paths, dirs, opt,
                                               skip_clip=opt.scan_clip,
                                               skip_disc=opt.scan_disc):
                if opt.dry_run:
                    print(f"# {t}:{name}\npython {XTEA_MAIN} " + " ".join(argv))
                    continue
                run_step(f"{t}:{name}", argv, logfh)
            if opt.dry_run:
                continue

            vcf = type_dir / f"{opt.sample}_{spec['vcf_alt']}.vcf"
            if vcf.exists():
                produced_vcfs.append(str(vcf))
                log(f"{t} VCF: {vcf}")
            else:
                log(f"WARNING: expected VCF not found for {t}: {vcf}")

    log("run complete")
    log("per-type VCFs:")
    for v in produced_vcfs:
        log(f"  {v}")
    # (A0) merged single-VCF output is added next; per-type VCFs feed the
    # concordance gate meanwhile.


if __name__ == "__main__":
    main()
