#!/usr/bin/env python
"""Adapter: C++ xtea-scan (with --mosaic-polya)  ->  master's mosaic clip candidate list.

MOSAIC_PLAN.md M.1. Mirrors `scan_adapter.py`'s shape, and replaces master's mosaic clip step
(`x_TEA_main.py --mosaic` -> TE_Multi_Locator.call_TEI_candidate_sites_from_multiple_alignmts_mosaic
-> clip_read.run_cnt_clip_part_aligned_to_rep_by_chrm_mosaic) outright rather than porting it.

WHY REPLACE RATHER THAN PORT, unlike the case/control work
-----------------------------------------------------------
For case/control the whole point was to reproduce master's own step chain command-for-command.
Here the stated priority is the opposite -- "the most important thing with mosaic besides accuracy
is speed! This is prohibitively slow on 300x data" -- and master's mosaic clip path is exactly the
per-read Python loop that does not survive 300x. It was also DELETED from this tree in trim step 4
(x_TEA_main.py:721-726), so there is nothing left to shell out to. The polyA fact that path
existed to compute is now one extra counter on the scan pass we already run (M.0), so this file
only has to join and threshold what the scan already produced.

THE OUTPUT CONTRACT, taken from master rather than invented
------------------------------------------------------------
master's `merge_clip_positions_with_cutoff_polyA` (clip_read.py:1096) writes, per surviving site:

    chrom  pos  nL_raw  nR_raw  nL_realigned  nR_realigned  npL  npR

The last four come from `cnt_clip_part_aligned_to_rep_mosaic`, whose per-chromosome writer
(clip_read.py:993-999) emits `pos` then m_sites_chrm[pos] (the two realigned side counts) then
m_polyA[pos][0], m_polyA[pos][1] -- which is why master's own miss branch writes exactly
"0\t0\t0\t0". A site with raw clip support but no realigned clip is KEPT with zeros, not dropped;
that is master's behaviour and the mosaic post-filter relies on seeing those rows.

TWO INPUTS, because raw and realigned clip counts are different measurements
-----------------------------------------------------------------------------
  * --clip-sites  `chrom pos nL nR`                     -- every clipped read at the position
  * --realign-out `chrom pos family nL nR mate npL npR` -- only those whose clipped part realigned
                                                           to a TE consensus (needs --mosaic-polya)

master thresholds on the RAW counts and reports the realigned ones as evidence, so both are
required and they are not interchangeable. Passing a 6-column realign-out here is a hard error
rather than a silent zero-fill: a mosaic run whose polyA tally is quietly all zeros still produces
a plausible-looking candidate list, and the polyA signal is most of what separates a mosaic
insertion from background at 1-5% VAF.
"""

import sys
from collections import defaultdict

# MEI type -> the family tag xtea-scan's --realign-out uses. Same mapping as scan_adapter.
TYPE_TO_FAMILY = {"L1": "LINE1", "Alu": "ALU", "SVA": "SVA"}


class MosaicPolyAMissing(ValueError):
    """A realign-out was produced without --mosaic-polya (6 columns, not 8)."""


def read_clip_sites(path):
    """--clip-sites -> {(chrom, pos): [nL, nR]}, the RAW per-position clipped-read counts."""
    sites = {}
    with open(path) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 4:
                continue
            try:
                sites[(f[0], int(f[1]))] = [int(f[2]), int(f[3])]
            except ValueError:
                continue          # header or malformed line
    return sites


def read_realign_polya(path, mei_type):
    """--realign-out (8 col) -> {(chrom, pos): [nL, nR, npL, npR]} for one family.

    Raises MosaicPolyAMissing if the file carries the plain 6-column germline format, so a run
    that forgot --mosaic-polya fails loudly instead of reporting zero polyA everywhere.
    """
    family = TYPE_TO_FAMILY[mei_type]
    out = defaultdict(lambda: [0, 0, 0, 0])
    saw_row = False
    for line in open(path):
        f = line.split()
        if len(f) < 3:
            continue
        if len(f) < 8:
            raise MosaicPolyAMissing(
                f"{path}: {len(f)} columns, expected 8. Re-run xtea-scan with --mosaic-polya; "
                "without it the mosaic polyA tally is absent and every site would score 0.")
        if f[2] != family:
            continue
        saw_row = True
        rec = out[(f[0], int(f[1]))]
        # Summed rather than assigned: a (chrom,pos) can legitimately appear more than once
        # per family across the merge, and scan_adapter accumulates for the same reason.
        for i, col in enumerate((3, 4, 6, 7)):
            rec[i] += int(f[col])
    if not saw_row:
        sys.stderr.write(f"WARNING: {path} has no {family} rows -- mosaic list will be raw-clip only\n")
    return out


def write_mosaic_candidate_list(clip_sites_path, realign_tsv, mei_type, out_path,
                                cutoff_left_clip=2, cutoff_right_clip=2, max_cov_cutoff=1000):
    """Join raw + realigned counts, apply master's two cutoffs, write the mosaic candidate list.

    The cutoff rule is master's verbatim (clip_read.py:1121-1124) and its shape matters:

        drop when  nL_raw < cutoff_left AND nR_raw < cutoff_right     -- AND, not OR

    so ONE well-supported side is enough to keep a site. A mosaic insertion at 2% VAF frequently
    has clip support on one side only, and an OR here would discard exactly those.

    max_cov_cutoff drops pileup artefacts (nL+nR above it), which at 300-500x is the coverage
    regime this whole mode exists for -- it is not a formality.

    Returns (n_written, n_dropped_low, n_dropped_high).
    """
    raw = read_clip_sites(clip_sites_path)
    realigned = read_realign_polya(realign_tsv, mei_type)

    n = low = high = 0
    with open(out_path, "w") as fout:
        for (chrom, pos) in sorted(raw, key=lambda k: (k[0], k[1])):
            nl_raw, nr_raw = raw[(chrom, pos)]
            if nl_raw < cutoff_left_clip and nr_raw < cutoff_right_clip:
                low += 1
                continue
            if (nl_raw + nr_raw) > max_cov_cutoff:
                high += 1
                continue
            # Sites with no realigned clip are kept with zeros -- master's own miss branch.
            nl, nr, npl, npr = realigned.get((chrom, pos), (0, 0, 0, 0))
            fout.write(f"{chrom}\t{pos}\t{nl_raw}\t{nr_raw}\t{nl}\t{nr}\t{npl}\t{npr}\n")
            n += 1
    return n, low, high
