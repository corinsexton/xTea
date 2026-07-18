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
from pathlib import Path

# MEI type -> the family tag xtea-scan's --realign-out uses.
TYPE_TO_FAMILY = {"L1": "LINE1", "Alu": "ALU", "SVA": "SVA"}

# master's SVA_ANNOTATION_EXTND (global_values.py). SVA repeat intervals are
# widened by this on both sides for the mate-in-repeat test (under --sva); L1/Alu
# use 0. col5 (clip mate-in-repeat) uses the SAME extension as the disc step.
SVA_ANNOTATION_EXTND = 200


def out_to_bed(out_path, bed_path, extnd):
    """Convert a RepeatMasker .out to a BED of extended repeat intervals for the
    col5 mate-in-repeat lookup, replicating master's XAnnotation parsing exactly:
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
                     cns_collect=None):
    """Run one xtea-scan pass: clip-detect + concurrent realign to all families.

    cns_by_family: {"LINE1": path, "ALU": path, "SVA": path} (any subset).
    Writes `out_tsv` with rows: chrom \t pos \t family \t nL \t nR.

    If disc_pos_out is given, the SAME pass also streams the discordant-read
    records (chrom read_pos is_rev mate_chrom mate_pos qname is_mate_rev) used by
    the A2 disc adapter -- so all BAM reading happens in this one scan. Set a low
    disc_threshold so the stream is a superset of what master's 100 kb same-chrom
    test counts (the adapter re-applies the exact threshold).
    """
    cmd = [scan_bin, "--bam", str(bam), "--reference", str(reference),
           "--realign-out", str(out_tsv), "--threads", str(threads)]
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
    if bwa:
        cmd += ["--bwa", str(bwa)]
    if region:
        cmd += ["--region", region]
    # full fold-in (move B.1): spill the cns collection (clip parts + disc anchors + low-mapq clip
    # positions) DURING this detection pass; a post-candidate fan step turns it into per-type temp
    # SAMs -- so cns costs no extra whole-genome pass beyond the mate-gather.
    if cns_collect:
        cmd += ["--cns-collect-clip", str(cns_collect["clip"]),
                "--cns-collect-anchor", str(cns_collect["anchor"]),
                "--cns-collect-lowq", str(cns_collect["lowq"])]
    where = logfh or sys.stdout
    where.write(f"\n=== xtea-scan (clip+realign) ===\n{' '.join(cmd)}\n")
    where.flush()
    r = subprocess.run(cmd, stdout=where, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.exit(f"ERROR: xtea-scan failed (exit {r.returncode})")


def write_clip_candidate_list(realign_tsv, mei_type, out_path, min_clip=1):
    """Reformat xtea-scan realign-out into master's candidate_list_from_clip.txt.

    Selects rows for this type's family and writes `chrom pos nL nR 0`
    (third column = master's placeholder; downstream reads it as an int feature).
    Emits any position with nL+nR >= min_clip (superset-safe; downstream filters).
    Returns the number of candidate positions written.
    """
    family = TYPE_TO_FAMILY[mei_type]
    # (chrom,pos) -> [nL, nR, col5]   (realign-out is per-family per-position;
    # col5 = mate-in-repeat clip count, xtea-scan's 6th column, computed natively).
    sites = defaultdict(lambda: [0, 0, 0])
    with open(realign_tsv) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 5 or f[2] != family:
                continue
            chrom, pos, nL, nR = f[0], int(f[1]), int(f[3]), int(f[4])
            col5 = int(f[5]) if len(f) >= 6 else 0
            s = sites[(chrom, pos)]
            s[0] += nL
            s[1] += nR
            s[2] += col5
    n = 0
    with open(out_path, "w") as fout:
        # deterministic order: chrom then pos
        for (chrom, pos) in sorted(sites, key=lambda k: (k[0], k[1])):
            nL, nR, col5 = sites[(chrom, pos)]
            if nL + nR < min_clip:
                continue
            fout.write(f"{chrom}\t{pos}\t{nL}\t{nR}\t{col5}\n")
            n += 1
    return n
