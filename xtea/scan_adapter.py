#!/usr/bin/env python
"""Adapter: C++ xtea-scan  ->  master's clip candidate-list format.

A1 refactor-inward step. Replaces the master clip step (`x_TEA_main.py -C`, slow
pysam scan + per-type bwa realign) with ONE `xtea-scan` pass that:
  * detects clipped reads across the whole BAM in a single streaming pass, and
  * realigns the clipped parts to the L1/Alu/SVA consensi concurrently,
emitting per-site, per-family in-repeat clip counts (`chrom pos family nL nR`).

Key enabling insight (verified against master): the downstream cns-filter step
re-collects clipped/discordant reads *from the BAM by position* and re-counts
them itself (`collect_clipped_disc_reads_by_region` fetches [pos-extnd, pos+extnd]).
So the handoff into the unchanged downstream is only the candidate *positions* +
advisory counts -- xtea-scan does not need to hand off reads or realignments, and
its site list only has to be a superset of master's clip peaks for the pipeline to
converge to the same VCF. Peak-calling itself happens later (in the disc step's
`call_peak_candidate_sites_with_std_derivation`), so this file writes raw
per-position counts with no merging.

Correctness bar (per user): final-VCF concordance vs the golden, not byte-identity
of this intermediate.
"""

import subprocess
import sys
from collections import defaultdict

# MEI type -> the family tag xtea-scan's --realign-out uses.
TYPE_TO_FAMILY = {"L1": "LINE1", "Alu": "ALU", "SVA": "SVA"}

# master's SVA_ANNOTATION_EXTND (global_values.py). SVA repeat intervals are
# widened by this on both sides for the mate-in-repeat test (under --sva); L1/Alu
# use 0. mate_in_rep_clip_cnt uses the SAME extension as the disc step.
SVA_ANNOTATION_EXTND = 200


def out_to_bed(out_path, bed_path, extnd):
    """Convert a RepeatMasker .out to a BED of extended repeat intervals for the
    mate_in_rep_clip_cnt lookup, replicating master's XAnnotation parsing exactly:
    fields[4]=chrom, fields[5]=start, fields[6]=end (min_len=0, no length filter),
    interval [start-extnd, end+extnd). Written as 0-based right-open BED so htslib
    regidx's point stab (begin<=pos<end at the 0-based mate pos) reproduces master's
    half-open interval-tree query -- including master's raw 1-based/0-based mix.
    Returns the number of intervals written.
    """
    rows = []
    with open(out_path) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 7:
                continue
            try:
                start = int(f[5])
                end = int(f[6])
            except ValueError:
                continue  # RepeatMasker header / non-data line
            b = start - extnd
            if b < 0:
                b = 0
            rows.append((f[4], b, end + extnd))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(bed_path, "w") as fout:
        for chrom, b, e in rows:
            fout.write("{0}\t{1}\t{2}\n".format(chrom, b, e))
    return len(rows)


def run_scan_realign(scan_bin, bam, reference, cns_by_family, out_tsv,
                     threads=8, bwa=None, region=None, logfh=None,
                     disc_pos_out=None, disc_threshold=None, bed_by_family=None,
                     bwa_threads=None, clip_sites_out=None, mosaic_polya=False):
    """Run one xtea-scan pass: clip-detect + concurrent realign to all families.

    cns_by_family: {"LINE1": path, "ALU": path, "SVA": path} (any subset).
    Writes `out_tsv` with rows: chrom \t pos \t family \t nL \t nR.

    bwa_threads is the TOTAL bwa-realign core budget that xtea-scan divides across
    the concurrent per-family bwa processes (e.g. 24 -> 8/8/8 for L1/Alu/SVA). It
    defaults to `threads` because the realign runs AFTER the tiling scan finishes,
    so it can use the whole allocation. Leaving it unset makes xtea-scan fall back
    to its small internal default (4 -> 2/1/1), which starves the realign -- the
    bug that made the WGS realign run on ~4 cores instead of the allocated 24.

    If disc_pos_out is given, the SAME pass also streams the discordant-read
    records (chrom read_pos is_rev mate_chrom mate_pos qname is_mate_rev) used by
    the A2 disc adapter -- so all BAM reading happens in this one scan. Set a low
    disc_threshold so the stream is a superset of what master's 100 kb same-chrom
    test counts (the adapter re-applies the exact threshold).
    """
    if bwa_threads is None:
        bwa_threads = threads
    cmd = [scan_bin, "--bam", str(bam), "--reference", str(reference),
           "--realign-out", str(out_tsv), "--threads", str(threads),
           "--bwa-threads", str(bwa_threads)]
    if "LINE1" in cns_by_family:
        cmd += ["--l1-cns", str(cns_by_family["LINE1"])]
    if "ALU" in cns_by_family:
        cmd += ["--alu-cns", str(cns_by_family["ALU"])]
    if "SVA" in cns_by_family:
        cmd += ["--sva-cns", str(cns_by_family["SVA"])]
    if disc_pos_out:
        cmd += ["--disc-pos", str(disc_pos_out)]
    if disc_threshold is not None:
        cmd += ["--disc-threshold", str(disc_threshold)]
    if bed_by_family:
        if "LINE1" in bed_by_family:
            cmd += ["--l1-bed", str(bed_by_family["LINE1"])]
        if "ALU" in bed_by_family:
            cmd += ["--alu-bed", str(bed_by_family["ALU"])]
        if "SVA" in bed_by_family:
            cmd += ["--sva-bed", str(bed_by_family["SVA"])]
    # Mosaic (MOSAIC_PLAN M.1) needs two things this pass does not otherwise produce: the RAW
    # per-position clip counts it thresholds on, and the polyA tally that separates a mosaic
    # insertion from background at 1-5% VAF. Both are opt-in, so a germline pass is unchanged --
    # --mosaic-polya in particular would alter the realign-out column count that
    # write_clip_candidate_list below parses.
    if clip_sites_out:
        cmd += ["--clip-sites", str(clip_sites_out)]
    if mosaic_polya:
        cmd += ["--mosaic-polya"]
    if bwa:
        cmd += ["--bwa", str(bwa)]
    if region:
        cmd += ["--region", region]
    # full fold-in (move B.1): spill the cns collection (clip parts + disc anchors + low-mapq clip
    # positions) DURING this detection pass; a post-candidate fan step turns it into per-type temp
    # SAMs -- so cns costs no extra whole-genome pass beyond the mate-gather.
    where = logfh or sys.stdout
    where.write(f"\n=== xtea-scan (clip+realign) ===\n{' '.join(cmd)}\n")
    where.flush()
    r = subprocess.run(cmd, stdout=where, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.exit(f"ERROR: xtea-scan failed (exit {r.returncode})")


def write_clip_candidate_list(realign_tsv, mei_type, out_path, min_clip=1, region_mask=None):
    """Reformat xtea-scan realign-out into master's candidate_list_from_clip.txt.

    Selects rows for this type's family and writes `chrom pos nL nR 0`
    (third column = master's placeholder; downstream reads it as an int feature).
    Emits any position with nL+nR >= min_clip (superset-safe; downstream filters).
    Returns (n_written, n_masked).

    region_mask (xtea.region_mask.RegionMask or None) drops candidate POSITIONS inside
    excluded regions. This is the single place the mask is applied, and that is deliberate:
    disc_adapter peak-calls from this list and builds its raw-disc rows per candidate peak,
    so a position removed here is removed from every downstream stage -- while the disc
    evidence array itself stays complete, so surviving candidates near a mask boundary keep
    full support on both flanks. See region_mask's module docstring for why filtering the
    disc stream instead would silently drop real calls.
    """
    family = TYPE_TO_FAMILY[mei_type]
    # (chrom,pos) -> [nL, nR, mate_in_rep_clip_cnt]   (realign-out is per-family
    # per-position; mate_in_rep_clip_cnt is the # of clipped reads whose mate lands in any
    # selected family's repeat annotation, xtea-scan's 6th column, computed natively).
    sites = defaultdict(lambda: [0, 0, 0])
    masked = region_mask.contains if region_mask else None
    n_masked_sites = set()
    with open(realign_tsv) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 5 or f[2] != family:
                continue
            chrom, pos, nL, nR = f[0], int(f[1]), int(f[3]), int(f[4])
            # Masked at READ time, not at write time: this dict holds every candidate
            # position in the genome (~57M on a 173x WGS arm) and is the peak-memory
            # structure of this function, so excluding masked positions here keeps them out
            # of it entirely rather than merely omitting them from the output.
            if masked is not None and masked(chrom, pos):
                n_masked_sites.add((chrom, pos))
                continue
            mate_in_rep = int(f[5]) if len(f) >= 6 else 0
            s = sites[(chrom, pos)]
            s[0] += nL
            s[1] += nR
            s[2] += mate_in_rep
    n = 0
    with open(out_path, "w") as fout:
        # deterministic order: chrom then pos
        for (chrom, pos) in sorted(sites, key=lambda k: (k[0], k[1])):
            nL, nR, mate_in_rep = sites[(chrom, pos)]
            if nL + nR < min_clip:
                continue
            fout.write(f"{chrom}\t{pos}\t{nL}\t{nR}\t{mate_in_rep}\n")
            n += 1
    return n, len(n_masked_sites)
