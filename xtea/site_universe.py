"""The downstream's locus universe: every position the per-family chain can create a
record at, enumerable BEFORE any of that chain runs.

WHY THIS MODULE EXISTS -- the measurement, then the consequence.

`docs/handoff/runs/probe_site_creation.py` measured that **52 of BASELINE-0's 1,515
per-family final records have no `candidate_disc_filtered_cns.txt` row at all** -- 50 L1,
2 SVA, every one created by the transduction step. That is why the Stage 2c family-assign
barrier, which runs on the cns lists, can reach only 7 of the 22 family conflicts: for the
other 15 one side does not exist yet. Plan sec 4 calls the fix "unified site creation" and
names it a Stage 3 dependency.

Two follow-up measurements (this module's reason to exist) say where those sites come from:

1. **They are not moved records, they are new ones.** Every postF_hc position is present at
   distance 0 in `candidate_disc_filtered_cns2.txt`, and the nearest cns row to each of the
   52 is >= 3,293 bp away (median 2.4 Mb). The post-cns chain does not re-anchor
   coordinates, so "created by transduction" is literal.
2. **All 52 are already rows in `candidate_list_from_disc.txt.clip_sites_raw_disc.txt`** --
   the raw-disc file `disc_adapter.run_disc_from_scan` writes BEFORE `cns_filter` runs --
   and all 52 survive into `<raw_disc>.slct`, the selection
   `XTransduction.re_slct_with_clip_raw_disc_sites` derives from it.

So site creation is **a selection over a file that already exists at the barrier**, and it
needs **no BAM, no bwa and no cns evidence** -- only raw_disc, the cns candidate list (for
the proximity exclusion) and the family's rmsk annotation.

**"Up front" means before TRANSDUCTION, not before `cns_filter`.** The cns candidate list
this rule excludes against IS `cns_filter`'s output, so the earliest point the universe is
computable is after every family's `cns_filter` and before any family's transduction -- the
same point the Stage 2c family-assign barrier already occupies. `xtea/downstream_unified.py`
is the driver-side caller and runs exactly there. (raw_disc and the annotation genuinely do
precede `cns_filter`; the cns list is the one input that does not, and it is the one that
fixes the placement.)

    universe(family) = cns candidates(family)  U  new raw-disc sites(family)

`select_new_sites` / `write_slct_file` below are that rule, lifted out of
`x_transduction.py` verbatim so ONE copy serves both the existing step and any driver-side
caller. `XTransduction.re_slct_with_clip_raw_disc_sites` now delegates here and does
nothing else, which is what makes "the lift cannot change the step's output" structural
rather than tested -- the same discipline `locus_table.assign_type` uses.

WHAT THIS DOES **NOT** UNLOCK, so nobody re-derives the disappointment: knowing the site
exists is not knowing which family claims it. `.slct` is a ~34-40K-row superset per family
(L1 34,197, SVA 39,452 on BASELINE-0) whose per-family evidence only comes into being when
the transduction evaluation runs. A barrier placed before that evaluation still cannot
resolve the 15 remaining conflicts, so **the 7-of-22 ceiling is unchanged.** What moves is
the *iteration domain*: the per-locus loop can enumerate its loci before gathering any
evidence, which is what "one pass over loci" in plan sec 3 requires.

The named caller this unblocks: plan sec 4 2b.4 says hoisting the td gather to the driver
would mean "re-implementing step logic the per-locus loop deletes", specifically
`re_slct_with_clip_raw_disc_sites`. It is no longer step logic.
"""
import global_values
from x_intermediate_sites import XIntemediateSites


def is_matched_rep_type(i_rep_type, s_sub_family):
    """Does an rmsk subfamily name belong to the family this rtype bitmask asks for?

    Note there is no Alu arm (bit 2): Alu has no transduction flank, so the transduction
    step never runs for it (x_TEA_main copies its input through instead). Lifted from
    XTransduction.is_matched_rep_type, which now delegates here; `fall_in_low_div_same_type_rep`
    goes through the same delegation, so there is one copy of the name test.
    """
    if (i_rep_type & 1 == 1) and (("L1" in s_sub_family) or ("LINE/L1" in s_sub_family) or "LINE1" in s_sub_family):
        return True
    elif (i_rep_type & 4 == 4) and ("SVA" in s_sub_family):
        return True
    return False


def _sites_close(pos1, pos2, iwin):
    """`XIntemediateSites.are_sites_close(pos1, pos2, iwin)`, without its per-base loop.

    That method walks `range(pos1-iwin, pos1+iwin)` looking for pos2, so the window it
    actually implements is the **half-open** `[pos1-iwin, pos1+iwin)`. Here pos1 is the
    raw-disc site and pos2 the cns candidate, which makes the exclusion window on the
    raw-disc site's offset from the candidate `(-iwin, +iwin]`: a raw-disc site exactly
    TWO_SITES_CLOSE_DIST to the LEFT of a candidate is kept, exactly that far to the RIGHT
    is suppressed. The asymmetry is load-bearing -- it decides whether a site becomes a
    transduction candidate -- so do not "fix" it to a symmetric window. Kept as arithmetic
    rather than a loop because this runs 39,600 raw-disc rows x every candidate on the
    chromosome.
    """
    return (pos1 - iwin) <= pos2 < (pos1 + iwin)


def select_new_sites(sf_raw_disc, m_te_candidates, n_disc_cutoff, xannotation,
                     i_rep_type, b_tumor, log=print):
    """The raw-disc rows transduction re-selects as NEW candidate sites.

    Yields the raw_disc line (rstripped) for each kept row, in raw_disc file order, which
    is the order `write_slct_file` must preserve -- the transduction gather's read
    collection is emitted in site-file order and bwa's tie-break rides on it (see
    src/csrc/td_gather.hpp).

    m_te_candidates is `XIntemediateSites.load_in_candidate_list_str_version(cns_out)`:
    {chrom: {pos: [fields]}}. Only its keys are read here.

    Three rejections, in the original's order (which matters, because the third one logs):
      1. within TWO_SITES_CLOSE_DIST of a cns candidate on the same chromosome -- the cns
         step refines positions, so this is a proximity test, not equality;
      2. neither side reaches the half disc cutoff (int(nd*3/4) germline, int(nd/2) tumor);
      3. the site falls in a low-divergence copy of its OWN family -- a transduction whose
         insertion point is inside a fresh element of the same family is indistinguishable
         from mis-mapped reads out of that element.
    """
    n_half_cutoff = int(n_disc_cutoff / 2)
    if b_tumor == False:  # for non tumor, set a little bit higher cutoff
        n_half_cutoff = int(n_disc_cutoff * 3 / 4)
    iwin = global_values.TWO_SITES_CLOSE_DIST
    with open(sf_raw_disc) as fin_raw:
        for line in fin_raw:
            fields = line.split()
            ins_chrm = fields[0]
            ins_pos = int(fields[1])
            if ins_chrm in m_te_candidates:  # check whether already saved
                b_close = False
                for tmp_pos in m_te_candidates[ins_chrm]:
                    if _sites_close(ins_pos, tmp_pos, iwin):
                        b_close = True
                        break
                if b_close == True:
                    continue
            n_ldisc = int(fields[2])
            n_rdisc = int(fields[3])
            if n_ldisc < n_half_cutoff and n_rdisc < n_half_cutoff:
                continue

            b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
            div_rate, sub_family, family, pos_start, pos_end = xannotation.get_div_subfamily_with_min_div(
                ins_chrm, i_pos)
            b_match_rep = is_matched_rep_type(i_rep_type, sub_family)
            # if fall in same type of repetitive region, then skip
            if (b_in_rep is True) and (b_match_rep is True) and (div_rate < global_values.TD_REP_DIVERGENT_CUTOFF):
                log("[re-select step]:Filtered out: {0}:{1} fall in repetitive region.".format(ins_chrm, ins_pos))
                continue
            yield line.rstrip()


def write_slct_file(sf_raw_disc, sf_candidates, n_disc_cutoff, xannotation,
                    i_rep_type, b_tumor, sf_out, b_with_candidate=True, log=print):
    """Write `<raw_disc>.slct`: the transduction step's full candidate site list.

    Layout, and it is what the file's consumers rely on:
      - when b_with_candidate, every cns candidate first as a bare `chrom<TAB>pos`, in the
        cns file's own order;
      - then each newly selected raw-disc row, full 6 fields, in raw_disc order.

    Returns (n_candidates, n_new).
    """
    xintmdt = XIntemediateSites()
    m_te_candidates = xintmdt.load_in_candidate_list_str_version(sf_candidates)
    n_cand = 0
    n_new = 0
    with open(sf_out, "w") as fout_slct:
        if b_with_candidate == True:
            for ins_chrm in m_te_candidates:
                for ins_pos in m_te_candidates[ins_chrm]:
                    fout_slct.write(ins_chrm + "\t" + str(ins_pos) + "\n")
                    n_cand += 1
        for s_line in select_new_sites(sf_raw_disc, m_te_candidates, n_disc_cutoff,
                                       xannotation, i_rep_type, b_tumor, log=log):
            fout_slct.write(s_line + "\n")
            n_new += 1
    return n_cand, n_new


def read_positions(sf_sites):
    """{chrom: sorted[pos]} from any two-or-more-column site file -- a `.slct` list, a cns
    list, a raw_disc list. Rows with fewer than two fields or a non-integer position are
    skipped rather than raising: raw_disc and cns lists are both written by this pipeline,
    but the cns family also carries `.high_confident` variants whose provenance is looser.
    """
    m = {}
    with open(sf_sites) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            m.setdefault(fields[0], []).append(pos)
    for chrm in m:
        m[chrm].sort()
    return m
