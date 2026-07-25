#!/usr/bin/env python
"""Adapter: C++ xtea-scan --disc-pos  ->  master's disc candidate-list format.

A2 refactor-inward step. Replaces the master disc step (`x_TEA_main.py -D`, which
re-opens the BAM per type and re-fetches two windows around every clip peak via
`x_TEI_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts`) with a
Python pass that sources discordant reads from the *single* xtea-scan pass's
`--disc-pos` stream (chrom read_pos is_rev mate_chrom mate_pos qname is_mate_rev).

Redundancy removed (the point of A2 per user): the whole per-type BAM re-scan +
per-chrom candidate-file splitting + `TELocator`/`cnt_discordant_pairs` fetch
machinery collapses to one shared disc-pos read + in-memory window counting. All
BAM reading now happens once, in the C scan.

Fidelity: this reuses master's own helpers verbatim for every output-shaping step
that anything downstream depends on --
  * `XIntemediateSites.call_peak_candidate_sites_with_std_derivation` (peak calling),
  * `XAnnotation.is_within_repeat_region_interval_tree`  (mate-in-repeat test),
  * `DiscCluster.form_one_side_cluster`                  (raw-disc cluster fields),
  * `XIntemediateSites.merge_clip_disc`                  (final candidate_list_from_disc.txt).
Only the *read source* changes (disc-pos stream instead of a per-peak BAM fetch), so
the raw_disc file that feeds transduction keeps master's exact 6-column contract.

Correctness bar (per user): final-VCF concordance vs the golden, not byte-identity
of this intermediate. Positions are the load-bearing output into cns (which
re-collects reads by region); disc counts here are advisory for cns but exact
enough to reproduce master's peak SET and the transduction raw_disc consistency.
"""

import bisect
from collections import defaultdict

import global_values
from x_intermediate_sites import XIntemediateSites
from x_annotation import XAnnotation
from disc_cluster import DiscCluster

# master hard-codes the same-chromosome discordant insert-size threshold to 100 kb
# in the disc step (x_TEA_main.py: `i_is = 100000`). A same-chrom mate counts as
# discordant only if it is farther than this from the *window start*.
I_IS = 100000


def load_disc_pos(disc_pos_file):
    """Load xtea-scan --disc-pos into per-chrom arrays sorted by read position.

    Returns {chrom: (pos_list, end_list, mate_chrom_list, mate_pos_list)} where the
    four lists are parallel and sorted ascending by read start position (for bisect
    windowing). end_list is the read reference END (pysam reference_end), needed to
    reproduce master's fetch-overlap semantics exactly. disc-pos is already MAPQ>=20,
    non-dup/secondary/supplementary, mate-mapped (filtered in the C scan), matching
    master's cnt_discordant_pairs prefilter.
    """
    tmp = defaultdict(list)
    with open(disc_pos_file) as fin:
        for line in fin:
            f = line.split()
            if len(f) < 8:
                continue
            # chrom read_pos is_rev mate_chrom mate_pos qname is_mate_rev read_end
            tmp[f[0]].append((int(f[1]), int(f[7]), f[3], int(f[4])))
    out = {}
    for chrom, recs in tmp.items():
        recs.sort(key=lambda r: r[0])
        out[chrom] = ([r[0] for r in recs],
                      [r[1] for r in recs],
                      [r[2] for r in recs],
                      [r[3] for r in recs])
    return out


def _count_window(disc_arrays, site_chrm, w_start, w_end, i_max_is,
                  xannotation, is_decoy, read_len):
    """Reproduce cnt_discordant_pairs for one [w_start, w_end] window from disc-pos.

    Returns (n_cnt, n_raw_cnt, (b_cluster, c_chrm, c_pos)):
      n_cnt     -- discordant reads whose mate falls in the repeat annotation
                   (the count master filters/merges on).
      n_raw_cnt -- discordant reads farther than i_max_is (or cross-chrom), the
                   raw pool DiscCluster clusters over (feeds transduction).
    Mirrors x_alignments.py:cnt_discordant_pairs exactly, with the read source
    swapped from a bamfile.fetch to the pre-collected disc-pos window.

    disc-pos records a read at its start (reference_start); master's bamfile.fetch
    returns every read *overlapping* [w_start, w_end], including those that start
    just before w_start but extend into the window (reads spanning a breakpoint --
    the dominant disc signal on the far side of a clip). So the lower bound is
    widened by one read length to recover them; the discordance test below still
    uses the true window start (w_start), exactly as master does.
    """
    n_cnt = 0
    n_raw = 0
    m_mate_pos = {}
    arr = disc_arrays.get(site_chrm)
    if arr is not None:
        pos_list, end_list, mc_list, mp_list = arr
        # master's bamfile.fetch(w_start, w_end) returns reads OVERLAPPING [w_start,
        # w_end): ref_start < w_end AND ref_end > w_start. Bisect a generous lower
        # bound on ref_start (read_len covers the longest read span), upper bound at
        # w_end, then filter on the exact ref_end to match fetch read-for-read.
        # ref_start in [w_start-(read_len+buffer), w_end): the buffer covers reads
        # whose ref span exceeds read_len (small indels); the exact ref_end filter
        # below discards any that don't truly overlap, so over-scanning is harmless.
        lo = bisect.bisect_left(pos_list, w_start - read_len - 100)
        hi = bisect.bisect_left(pos_list, w_end)      # ref_start < w_end (strict)
        for i in range(lo, hi):
            if end_list[i] <= w_start:   # ref_end must exceed w_start to overlap [w_start, w_end)
                continue
            mate_chrm = mc_list[i]
            mate_pos = mp_list[i]
            if is_decoy(mate_chrm):
                continue
            # master: discordant if cross-chrom, or same-chrom but far from window start
            if (site_chrm != mate_chrm) or (abs(mate_pos - w_start) > I_IS):
                # raw pool (for clustering): cross-chrom, or same-chrom beyond i_max_is
                if (site_chrm != mate_chrm) or (abs(mate_pos - w_start) > i_max_is):
                    m_mate_pos.setdefault(mate_chrm, []).append(mate_pos)
                    n_raw += 1
                b_in_rep, _ = xannotation.is_within_repeat_region_interval_tree(
                    mate_chrm, mate_pos)
                if b_in_rep:
                    n_cnt += 1
    dc = DiscCluster()
    b_cluster, c_chrm, c_pos = dc.form_one_side_cluster(
        m_mate_pos, I_IS, global_values.MIN_RAW_DISC_CLUSTER_RATIO)
    return n_cnt, n_raw, (b_cluster, c_chrm, c_pos)


def run_disc_from_scan(clip_out, disc_out, disc_arrays, annotation, mei_type,
                       iextend, std_var, disc_cutoff, read_len, wfolder,
                       b_with_chr=True, b_tumor=False, peak_window=100,
                       log=print):
    """Native disc step for one MEI type, sourcing reads from xtea-scan --disc-pos.

    Writes:
      disc_out                     -- candidate_list_from_disc.txt (via merge_clip_disc)
      disc_out + RAW_DISC_TMP_SUFFIX -- raw disc file feeding transduction
    Returns the number of peaks written to disc_out.
    """
    xsites = XIntemediateSites()

    # 1. peak-call the A1 clip candidates (master's exact routine)
    m_original = xsites.load_in_candidate_list(clip_out)
    m_peaks = xsites.call_peak_candidate_sites_with_std_derivation(
        m_original, peak_window)
    m_original.clear()

    # 2. load the type's repeat annotation for the mate-in-repeat test. This must
    #    replicate the `-D --sva` subprocess this adapter replaces: for SVA the
    #    driver now passes --sva (build_type_steps), so IS_CALL_SVA is True and the
    #    boundary extension is SVA_ANNOTATION_EXTND=200, matching the golden. The
    #    cns step (-N --sva) uses the same extend, so disc and cns are consistent
    #    (an inconsistent disc=200/cns=0 is what silently dropped SVA before).
    global_values.IS_CALL_SVA = (mei_type == "SVA")
    boundary_extnd = (global_values.SVA_ANNOTATION_EXTND
                      if global_values.IS_CALL_SVA else 0)
    xannotation = XAnnotation(annotation)
    xannotation.set_with_chr(b_with_chr)
    xannotation.load_rmsk_annotation_with_extnd_with_lenth_cutoff(boundary_extnd, 0)
    xannotation.index_rmsk_annotation_interval_tree()

    i_max_is = 2500
    if int(I_IS + 3 * std_var) > i_max_is:
        i_max_is = int(I_IS + 3 * std_var)

    is_decoy = xsites.is_decoy_contig_chrms

    # 3. per-peak L/R window disc counts + raw/cluster (mirrors
    #    run_filter_by_discordant_pair_by_chrom_non_barcode)
    disc_tmp_rows = []   # (chrm, pos, nL, nR) after the disc cutoff
    raw_rows = []        # (chrm, pos, n_raw_l, n_raw_r, nL, nR, l_ok, r_ok)
    n_half_cutoff = disc_cutoff / 2.0
    for chrm in m_peaks:
        if xsites.is_decoy_contig_chrms(chrm):
            continue
        for pos in m_peaks[chrm]:
            if pos < iextend:
                continue
            nL, nrawL, lclu = _count_window(
                disc_arrays, chrm, pos - iextend, pos, i_max_is, xannotation, is_decoy, read_len)
            nR, nrawR, rclu = _count_window(
                disc_arrays, chrm, pos + 1, pos + iextend, i_max_is, xannotation, is_decoy, read_len)

            # disc-filter cutoff (x_TEI_locator.py:462, germline branch)
            if b_tumor:
                keep = (nL + nR) >= disc_cutoff
            else:
                keep = ((nL > n_half_cutoff and nR > n_half_cutoff)
                        or nL > disc_cutoff or nR > disc_cutoff)
            if keep:
                disc_tmp_rows.append((chrm, pos, nL, nR))
            raw_rows.append((chrm, pos, nrawL, nrawR, nL, nR, lclu[0], rclu[0]))

    # 4. write disc_tmp, then merge with the clip list into candidate_list_from_disc.txt
    #    (reuse master's merge_clip_disc so the column layout is identical)
    sf_disc_tmp = wfolder + "disc_tmp.list"
    with open(sf_disc_tmp, "w") as fout:
        for chrm, pos, nL, nR in disc_tmp_rows:
            fout.write("{0}\t{1}\t{2}\t{3}\n".format(chrm, pos, nL, nR))
    xsites.merge_clip_disc(sf_disc_tmp, clip_out, disc_out)

    # 5. raw_disc file for transduction (x_TEI_locator.py:516-542 + cluster-consistency
    #    gate). Single-sample => no cross-bam raw merge.
    #    b_tumor branch mirrors x_TEI_locator.py:516-520 exactly: tumor uses cutoff/2,
    #    germline raises it to int(cutoff*3/4)+1. This was previously hardcoded to the
    #    germline value regardless of b_tumor, which applied the STRICTER germline
    #    threshold to tumor runs (at cutoff=4: 4 instead of 2) and silently dropped raw
    #    disc sites master would keep. raw_disc feeds the transduction step, so the
    #    dropped sites cost transduction/orphan calls.
    i_half_cutoff = disc_cutoff / 2
    if not b_tumor:
        i_half_cutoff = int(disc_cutoff * 3 / 4) + 1
    if i_half_cutoff < 1:
        i_half_cutoff = 1
    sf_raw_disc = disc_out + global_values.RAW_DISC_TMP_SUFFIX
    n_raw_written = 0
    with open(sf_raw_disc, "w") as fout:
        for chrm, pos, nrawL, nrawR, nL, nR, l_ok, r_ok in raw_rows:
            b_l = (nrawL >= i_half_cutoff) and l_ok
            b_r = (nrawR >= i_half_cutoff) and r_ok
            if b_l or b_r:
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(
                    chrm, pos, nrawL, nrawR, nL, nR))
                n_raw_written += 1

    # count peaks actually landing in disc_out (clip∩disc), for logging
    n_out = 0
    try:
        with open(disc_out) as fin:
            n_out = sum(1 for _ in fin)
    except OSError:
        pass
    log("{0}:disc (xtea-scan) {1} peaks -> {2} pass disc cutoff, "
        "{3} raw-disc rows -> {4}".format(
            mei_type, sum(len(m_peaks[c]) for c in m_peaks),
            len(disc_tmp_rows), n_raw_written, disc_out))
    return n_out
