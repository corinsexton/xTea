"""The transduction decision, pinned against the cascade it replaces.

WHY THIS GATE EXISTS
--------------------
`xtea/transduction_evidence.py` moves five groups of decisions out of `x_transduction.py` --
the clip gate, the decoy tally, the source-dominance pick, the candidate conjunction and the
site-sanity filter -- so that every threshold is a named parameter and every drop is a named
reason (DOWNSTREAM_OVERHAUL_PLAN.md §5). The whole safety argument is **moved, not rewritten**,
and this file is what makes that a verified claim: `_legacy_*` below are VERBATIM copies of the
bodies that were deleted, and each is compared against its replacement over randomized inputs.

TWO THINGS THIS FILE ASSERTS THAT AN OUTPUT COMPARISON WOULD MISS.

1. **Key ORDER, not just content.** The step writes its output files by iterating these maps,
   so a replacement that agrees on every locus but builds its map in a different order changes
   the step's BYTES. Job 49596361 failed on exactly that in the locus-universe table. Every
   differential here compares `_ordered()` — nested (key, value) lists — not dicts.

2. **The two preserved defects.** `dominant_side_cluster` reproduces a cluster test that sits
   inside a per-read loop, so a locus can carry a left-side flag from read 3 next to
   boundaries from read 9. `test_the_stale_side_flag_is_reproduced` constructs an input where
   the defect is observable and asserts BOTH implementations show it — if someone "fixes" the
   library, that test fails and says why the fix is a Stage 3 change rather than a cleanup.

The remaining tests pin one gate each, so a failure names the threshold that broke.
"""
import os
import random
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values  # noqa: E402
import transduction_evidence as te  # noqa: E402
from x_cluster_consistency import ClusterChecker  # noqa: E402

L1_RTYPE = 1
SVA_RTYPE = 4

# The step's own arguments at chr1/WGS defaults (`--nclip 4 --nd 5`), so the derived cutoffs
# under test are the production ones: n_polya_cutoff 1.0, n_half_disc 2.5, n_decoy 5.0.
NCLIP = 4
NDISC = 5


def params(i_rep_type=L1_RTYPE, n_clip_cutoff=NCLIP, n_disc_cutoff=NDISC, i_max_cov=100,
           ave_cov=30.0, max_indel_ratio=0.2, f_min_uniq_ratio=0.5):
    return te.TransductionParams(
        n_clip_cutoff=n_clip_cutoff, n_disc_cutoff=n_disc_cutoff, i_max_cov=i_max_cov,
        ave_cov=ave_cov, max_indel_ratio=max_indel_ratio,
        f_min_uniq_ratio=f_min_uniq_ratio, i_rep_type=i_rep_type)


def _ordered(m):
    """Nested (key, value) lists, so a comparison sees insertion order at every level."""
    if isinstance(m, dict):
        return [(k, _ordered(v)) for k, v in m.items()]
    if isinstance(m, (list, tuple)):
        return type(m)(_ordered(v) for v in m)
    return m


# ======================================================================================
# VERBATIM copies of the deleted bodies. Do not tidy: their value is being unchanged.
# ======================================================================================

def _legacy_filter_by_clip_algnmt(m_ltsdct, m_rtsdct, m_polyA, n_clip_cutoff, n_polyA_cutoff):
    m_clip_cnt = {}
    for ins_chrm in m_ltsdct:
        if ins_chrm not in m_clip_cnt:
            m_clip_cnt[ins_chrm] = {}
        for ins_pos in m_ltsdct[ins_chrm]:
            n_left = 0
            for ref_mpos in m_ltsdct[ins_chrm][ins_pos]:
                n_left += len(m_ltsdct[ins_chrm][ins_pos][ref_mpos])
            if ins_pos not in m_clip_cnt[ins_chrm]:
                m_clip_cnt[ins_chrm][ins_pos] = 0
            m_clip_cnt[ins_chrm][ins_pos] += n_left

    for ins_chrm in m_rtsdct:
        if ins_chrm not in m_clip_cnt:
            m_clip_cnt[ins_chrm] = {}
        for ins_pos in m_rtsdct[ins_chrm]:
            n_right = 0
            for ref_mpos in m_rtsdct[ins_chrm][ins_pos]:
                n_right += len(m_rtsdct[ins_chrm][ins_pos][ref_mpos])
            if ins_pos not in m_clip_cnt[ins_chrm]:
                m_clip_cnt[ins_chrm][ins_pos] = 0
            m_clip_cnt[ins_chrm][ins_pos] += n_right
    m_slct = {}
    for ins_chrm in m_clip_cnt:
        for ins_pos in m_clip_cnt[ins_chrm]:
            n_cnt = m_clip_cnt[ins_chrm][ins_pos]
            if n_cnt < n_clip_cutoff:
                continue
            if (ins_chrm not in m_polyA) or (ins_pos not in m_polyA[ins_chrm]):
                continue
            n_polyA = m_polyA[ins_chrm][ins_pos][0] + m_polyA[ins_chrm][ins_pos][1] + \
                      m_polyA[ins_chrm][ins_pos][2] + m_polyA[ins_chrm][ins_pos][3] + \
                      m_polyA[ins_chrm][ins_pos][4] + m_polyA[ins_chrm][ins_pos][5]
            if n_polyA < n_polyA_cutoff:
                continue
            if ins_chrm not in m_slct:
                m_slct[ins_chrm] = {}
            m_slct[ins_chrm][ins_pos] = (n_cnt, n_polyA)
    return m_slct


def _legacy_chk_decoy_cns_cnt(m_trsdct, m_Alu, m_L1, m_SVA, n_cutoff, i_rep_type):
    if i_rep_type < 0:
        return m_trsdct
    m_slct = {}
    for ins_chrm in m_trsdct:
        for ins_pos in m_trsdct[ins_chrm]:
            b_slct = True
            n_Alu = 0
            if (ins_chrm in m_Alu) and (ins_pos in m_Alu[ins_chrm]):
                n_Alu = m_Alu[ins_chrm][ins_pos]
            n_L1 = 0
            if (ins_chrm in m_L1) and (ins_pos in m_L1[ins_chrm]):
                n_L1 = m_L1[ins_chrm][ins_pos]
            n_SVA = 0
            if (ins_chrm in m_SVA) and (ins_pos in m_SVA[ins_chrm]):
                n_SVA = m_SVA[ins_chrm][ins_pos]
            if i_rep_type & 1 != 0:  # for L1
                if n_Alu >= n_cutoff:
                    b_slct = False
                if n_SVA >= n_cutoff:
                    b_slct = False
            elif i_rep_type & 4 != 0:  # for SVA
                if n_Alu >= n_cutoff:
                    b_slct = False
                if n_L1 >= n_cutoff:
                    b_slct = False
            if b_slct == True:
                if ins_chrm not in m_slct:
                    m_slct[ins_chrm] = {}
                m_slct[ins_chrm][ins_pos] = m_trsdct[ins_chrm][ins_pos]
    return m_slct


def _legacy_pick_trsdct_rescued(m_transduction, n_disc_cutoff, m_raw_td_src, f_min_ratio):
    m_trsd_picked = {}
    for ins_chrm in m_transduction:
        for ins_pos in m_transduction[ins_chrm]:
            max_cnt = 0
            max_source = ""
            n_tmp_all_cnt = 0
            for tmp_source in m_transduction[ins_chrm][ins_pos]:
                tmp_cnt = len(m_transduction[ins_chrm][ins_pos][tmp_source])
                n_tmp_all_cnt += tmp_cnt
                if tmp_cnt > max_cnt:
                    max_cnt = tmp_cnt
                    max_source = tmp_source
            if n_tmp_all_cnt < n_disc_cutoff:
                continue

            b_trsdct = False
            f_tmp_ratio = float(max_cnt) / float(n_tmp_all_cnt)
            n_non_unique = 0
            if ((ins_chrm in m_raw_td_src) and (ins_pos in m_raw_td_src[ins_chrm])
                    and (max_source in m_raw_td_src[ins_chrm][ins_pos])):
                n_non_unique = m_raw_td_src[ins_chrm][ins_pos][max_source]
            if (max_cnt >= n_disc_cutoff) and (f_tmp_ratio >= global_values.TRANSDCT_MULTI_SOURCE_MIN_RATIO) \
                    and (float(max_cnt) / float(max_cnt + n_non_unique) > f_min_ratio):
                b_trsdct = True

            if b_trsdct == True:
                if ins_chrm not in m_trsd_picked:
                    m_trsd_picked[ins_chrm] = {}
                l_pos = []
                for tmp_record in m_transduction[ins_chrm][ins_pos][max_source]:
                    l_pos.append(tmp_record)
                m_trsd_picked[ins_chrm][ins_pos] = (max_source, l_pos, max_cnt)
    return m_trsd_picked


def _legacy_slct_orphan_transduction(m_transduction, m_slcted, n_half_cutoff):
    m_orphan_transduct = {}
    for ins_chrm in m_transduction:
        for ins_pos in m_transduction[ins_chrm]:
            if (ins_chrm not in m_slcted) or (ins_pos not in m_slcted[ins_chrm]):
                continue
            n_left = 0
            n_right = 0
            tmp_source = m_slcted[ins_chrm][ins_pos][0]
            for rcd in m_transduction[ins_chrm][ins_pos][tmp_source]:
                anchor_pos = int(rcd[1])
                if (anchor_pos + 5) < ins_pos:
                    n_left += 1
                else:
                    n_right += 1
            if (n_left >= n_half_cutoff) and (n_right >= n_half_cutoff):
                if ins_chrm not in m_orphan_transduct:
                    m_orphan_transduct[ins_chrm] = {}
                m_orphan_transduct[ins_chrm][ins_pos] = (tmp_source, n_left, n_right)
    return m_orphan_transduct


def _legacy_call_out_candidate_transduct(m_repsnt_pos, m_disc_transdct, icorcord, ratio,
                                         m_clip_slct, m_disc_polyA, m_orphan):
    m_transduction = {}
    ck_cluster = ClusterChecker()
    for ins_chrm in m_repsnt_pos:
        for ins_pos in m_repsnt_pos[ins_chrm]:
            b_disc_polyA = False

            if (ins_chrm not in m_clip_slct) or (ins_pos not in m_clip_slct[ins_chrm]):
                continue
            if (ins_chrm in m_disc_polyA) and (ins_pos in m_disc_polyA[ins_chrm]):
                b_disc_polyA = True
            if b_disc_polyA == False:
                continue

            clip_pos = m_repsnt_pos[ins_chrm][ins_pos]
            b_l_trsdct = False
            b_r_trsdct = False
            trsdct_l = -1
            trsdct_r = -1
            n_disc_trsdct = 0
            if (ins_chrm in m_disc_transdct) and (ins_pos in m_disc_transdct[ins_chrm]):
                trsdct_source, l_trsdct_pos, n_unique = m_disc_transdct[ins_chrm][ins_pos]
                l_l_trsdct_pos = []
                l_r_trsdct_pos = []
                for (hit_ref_pos, anchor_pos, dist_from_rep) in l_trsdct_pos:
                    if anchor_pos < clip_pos:
                        l_l_trsdct_pos.append(hit_ref_pos)
                    else:
                        l_r_trsdct_pos.append(hit_ref_pos)
                    if len(l_l_trsdct_pos) > len(l_r_trsdct_pos):
                        b_l_trsdct, trsdct_l, trsdct_r = ck_cluster._is_disc_cluster(l_l_trsdct_pos, icorcord, ratio)
                        n_disc_trsdct = len(l_l_trsdct_pos)
                    else:
                        b_r_trsdct, trsdct_l, trsdct_r = ck_cluster._is_disc_cluster(l_r_trsdct_pos, icorcord, ratio)
                        n_disc_trsdct = len(l_r_trsdct_pos)
                if b_l_trsdct == True or b_r_trsdct == True:
                    if ins_chrm not in m_transduction:
                        m_transduction[ins_chrm] = {}
                    trsdct_info = "{0}~{1}-{2}".format(trsdct_source, trsdct_l, trsdct_r)
                    m_transduction[ins_chrm][ins_pos] = (b_l_trsdct, b_r_trsdct, n_disc_trsdct, trsdct_info)

    for ins_chrm in m_clip_slct:
        for ins_pos in m_clip_slct[ins_chrm]:
            if (ins_chrm in m_orphan) and (ins_pos in m_orphan[ins_chrm]):
                if ins_chrm not in m_transduction:
                    m_transduction[ins_chrm] = {}
                s_source = m_orphan[ins_chrm][ins_pos][0]
                trsdct_info = "{0}~orphan".format(s_source)
                n_disc_trsdct = m_orphan[ins_chrm][ins_pos][1] + m_orphan[ins_chrm][ins_pos][2]
                m_transduction[ins_chrm][ins_pos] = (False, True, n_disc_trsdct, trsdct_info)
    return m_transduction


def _legacy_filter_transduction(m_ori_ins, f_in_rep, m_read_depth, i_max_cov, m_clip_polyA,
                                n_half_clip):
    """`filter_transduction`'s loop, with the annotation call and the batched ReadDepth
    replaced by a predicate and a precomputed map (which is what the library takes)."""
    m_new_ins = {}
    for ins_chrm in m_ori_ins:
        for ins_pos in m_ori_ins[ins_chrm]:
            if f_in_rep(ins_chrm, ins_pos) == True:
                continue
            i_lcov1 = m_read_depth[ins_chrm][ins_pos][0]
            i_rcov1 = m_read_depth[ins_chrm][ins_pos][1]
            i_lcov2 = m_read_depth[ins_chrm][ins_pos][2]
            i_rcov2 = m_read_depth[ins_chrm][ins_pos][3]
            if i_lcov1 > i_max_cov or i_rcov1 > i_max_cov or i_lcov2 > i_max_cov or i_rcov2 > i_max_cov:
                continue
            if (ins_chrm in m_clip_polyA) and (ins_pos in m_clip_polyA[ins_chrm]):
                n_lpolyA = m_clip_polyA[ins_chrm][ins_pos][0]
                n_rpolyA = m_clip_polyA[ins_chrm][ins_pos][1]
                if n_lpolyA > n_half_clip and n_rpolyA > n_half_clip:
                    continue
                n_lpolyT = m_clip_polyA[ins_chrm][ins_pos][2]
                n_rpolyT = m_clip_polyA[ins_chrm][ins_pos][3]
                if n_lpolyT > n_half_clip and n_rpolyT > n_half_clip:
                    continue

            if ins_chrm not in m_new_ins:
                m_new_ins[ins_chrm] = {}
            m_new_ins[ins_chrm][ins_pos] = m_ori_ins[ins_chrm][ins_pos]
    return m_new_ins


# ======================================================================================
# Randomized input generators. Ranges are chosen so every gate both fires and passes --
# a generator that never trips a threshold makes the differential vacuous.
# ======================================================================================

CHRMS = ("chr1", "chr2", "chrX")


def _rand_clip_maps(rnd, n_loci=60):
    """(m_lclip, m_rclip, m_polyA) in the shapes `parse_transduction_clip_algnmt` returns."""
    m_l, m_r, m_polya = {}, {}, {}
    for _ in range(n_loci):
        chrm = rnd.choice(CHRMS)
        pos = rnd.randrange(1000, 9000)
        for m_side in (m_l, m_r):
            if rnd.random() < 0.75:
                per = m_side.setdefault(chrm, {}).setdefault(pos, {})
                for _k in range(rnd.randrange(1, 4)):
                    per[rnd.randrange(100, 999)] = ["r"] * rnd.randrange(0, 4)
        if rnd.random() < 0.8:
            m_polya.setdefault(chrm, {})[pos] = [rnd.randrange(0, 4) for _ in range(6)]
    return m_l, m_r, m_polya


def _rand_source_maps(rnd, n_loci=60):
    """(m_transduction, m_raw_td_src, decoy maps) as `parse_transduction_disc_algnmt` builds."""
    m_td, m_raw, m_alu, m_l1, m_sva = {}, {}, {}, {}, {}
    for _ in range(n_loci):
        chrm = rnd.choice(CHRMS)
        pos = rnd.randrange(1000, 9000)
        srcs = {}
        for i in range(rnd.randrange(1, 4)):
            name = "chr{0}:{1}-{2}~{3}".format(i, 1000 * i, 2000 * i + 500, i % 2)
            srcs[name] = [(rnd.randrange(1, 9000), rnd.randrange(1000, 9000),
                           rnd.randrange(0, 500)) for _ in range(rnd.randrange(0, 6))]
        m_td.setdefault(chrm, {})[pos] = srcs
        if rnd.random() < 0.5:
            m_raw.setdefault(chrm, {})[pos] = {
                s: rnd.randrange(0, 8) for s in srcs if rnd.random() < 0.7}
        for m_decoy in (m_alu, m_l1, m_sva):
            if rnd.random() < 0.4:
                m_decoy.setdefault(chrm, {})[pos] = rnd.randrange(0, 9)
    return m_td, m_raw, m_alu, m_l1, m_sva


def _rand_candidate_maps(rnd, n_loci=60):
    m_repsnt, m_disc_td, m_clip_slct, m_disc_polya, m_orphan = {}, {}, {}, {}, {}
    for _ in range(n_loci):
        chrm = rnd.choice(CHRMS)
        pos = rnd.randrange(2000, 8000)
        m_repsnt.setdefault(chrm, {})[pos] = pos + rnd.randrange(-30, 30)
        if rnd.random() < 0.8:
            l_pos = [(pos + rnd.randrange(-600, 600), pos + rnd.randrange(-200, 200),
                      rnd.randrange(0, 400)) for _ in range(rnd.randrange(0, 9))]
            m_disc_td.setdefault(chrm, {})[pos] = ("src{0}".format(pos), l_pos, len(l_pos))
        if rnd.random() < 0.8:
            m_clip_slct.setdefault(chrm, {})[pos] = (rnd.randrange(4, 20), rnd.randrange(0, 9))
        if rnd.random() < 0.7:
            m_disc_polya.setdefault(chrm, {})[pos] = (rnd.randrange(0, 5), rnd.randrange(0, 5))
        if rnd.random() < 0.2:
            m_orphan.setdefault(chrm, {})[pos] = ("osrc", rnd.randrange(1, 6), rnd.randrange(1, 6))
    return m_repsnt, m_disc_td, m_clip_slct, m_disc_polya, m_orphan


def _rand_sanity_maps(rnd, n_loci=80, i_max_cov=100):
    m_ins, m_depth, m_polya, in_rep = {}, {}, {}, set()
    for _ in range(n_loci):
        chrm = rnd.choice(CHRMS)
        pos = rnd.randrange(1000, 9000)
        m_ins.setdefault(chrm, {})[pos] = (rnd.random() < 0.5, True, 3, "src~1-2")
        m_depth.setdefault(chrm, {})[pos] = [rnd.randrange(0, i_max_cov + 40) for _ in range(4)]
        if rnd.random() < 0.7:
            m_polya.setdefault(chrm, {})[pos] = [rnd.randrange(0, 6) for _ in range(6)]
        if rnd.random() < 0.2:
            in_rep.add((chrm, pos))
    return m_ins, m_depth, m_polya, in_rep


# ======================================================================================
# The differentials
# ======================================================================================

def test_clip_gate_reproduces_the_replaced_implementation():
    p = params()
    for seed in range(40):
        rnd = random.Random(seed)
        m_l, m_r, m_polya = _rand_clip_maps(rnd)
        got = te.select_clip_sites(m_l, m_r, m_polya, p)
        want = _legacy_filter_by_clip_algnmt(m_l, m_r, m_polya, p.n_clip_cutoff,
                                            p.n_polya_cutoff)
        assert _ordered(got) == _ordered(want), "seed {0}".format(seed)


def test_decoy_gate_reproduces_the_replaced_implementation():
    for i_rep_type in (L1_RTYPE, SVA_RTYPE, 2, -1):
        p = params(i_rep_type=i_rep_type)
        for seed in range(30):
            rnd = random.Random(seed)
            m_td, _m_raw, m_alu, m_l1, m_sva = _rand_source_maps(rnd)
            got = te.select_by_decoy(m_td, m_alu, m_l1, m_sva, p)
            want = _legacy_chk_decoy_cns_cnt(m_td, m_alu, m_l1, m_sva, p.n_decoy_cutoff,
                                            i_rep_type)
            assert _ordered(got) == _ordered(want), \
                "rtype {0} seed {1}".format(i_rep_type, seed)


def test_source_gate_reproduces_the_replaced_implementation():
    p = params()
    for seed in range(40):
        rnd = random.Random(seed)
        m_td, m_raw, _a, _l, _s = _rand_source_maps(rnd)
        got = te.pick_dominant_sources(m_td, m_raw, p)
        want = _legacy_pick_trsdct_rescued(m_td, p.n_half_disc, m_raw, p.f_min_uniq_ratio)
        assert _ordered(got) == _ordered(want), "seed {0}".format(seed)


def test_orphan_selection_reproduces_the_replaced_implementation():
    p = params()
    for seed in range(40):
        rnd = random.Random(seed)
        m_td, m_raw, _a, _l, _s = _rand_source_maps(rnd)
        picked = te.pick_dominant_sources(m_td, m_raw, p)
        got = te.select_orphans(m_td, picked, p)
        want = _legacy_slct_orphan_transduction(m_td, picked, p.n_half_disc)
        assert _ordered(got) == _ordered(want), "seed {0}".format(seed)


def test_candidate_gate_reproduces_the_replaced_implementation():
    """Includes the in-loop cluster defect: both sides run the same `_is_disc_cluster`, on
    lists it sorts in place, in the same interleaving."""
    p = params()
    i_concord_dist, f_ratio = 450, 0.4
    for seed in range(40):
        rnd = random.Random(seed)
        m_repsnt, m_disc_td, m_clip_slct, m_disc_polya, m_orphan = _rand_candidate_maps(rnd)
        ck = ClusterChecker()

        def f_is_cluster(l_pos, ck=ck):
            return ck._is_disc_cluster(l_pos, i_concord_dist, f_ratio)

        got = te.call_candidates(m_repsnt, m_disc_td, m_clip_slct, m_disc_polya, m_orphan,
                                 f_is_cluster, p)
        want = _legacy_call_out_candidate_transduct(m_repsnt, m_disc_td, i_concord_dist,
                                                   f_ratio, m_clip_slct, m_disc_polya,
                                                   m_orphan)
        assert _ordered(got) == _ordered(want), "seed {0}".format(seed)


def test_sanity_gate_reproduces_the_replaced_implementation():
    i_max_cov = 100
    p = params(i_max_cov=i_max_cov)
    for seed in range(40):
        rnd = random.Random(seed)
        m_ins, m_depth, m_polya, in_rep = _rand_sanity_maps(rnd, i_max_cov=i_max_cov)

        def f_in_rep(chrm, pos, in_rep=in_rep):
            return (chrm, pos) in in_rep

        got = te.apply_site_filters(m_ins, m_depth, m_polya, f_in_rep, p)
        want = _legacy_filter_transduction(m_ins, f_in_rep, m_depth, i_max_cov, m_polya,
                                          p.n_half_clip)
        assert _ordered(got) == _ordered(want), "seed {0}".format(seed)


def test_the_whole_chain_reproduces_the_cascade():
    """A-through-D composed, on one randomized input, the way the step composes them.

    Each differential above holds its stage's inputs fixed; this one feeds each stage the
    PREVIOUS stage's real output, which is where an off-by-one in what a stage consumes
    (rather than in what it decides) would show.
    """
    p = params()
    i_concord_dist, f_ratio = 450, 0.4
    for seed in range(25):
        rnd = random.Random(seed)
        m_l, m_r, m_clip_polya = _rand_clip_maps(rnd, n_loci=80)
        m_td, m_raw, m_alu, m_l1, m_sva = _rand_source_maps(rnd, n_loci=80)
        ck = ClusterChecker()

        def f_is_cluster(l_pos, ck=ck):
            return ck._is_disc_cluster(l_pos, i_concord_dist, f_ratio)

        # library
        clip_slct = te.select_clip_sites(m_l, m_r, m_clip_polya, p)
        decoyed = te.select_by_decoy(m_td, m_alu, m_l1, m_sva, p)
        picked = te.pick_dominant_sources(decoyed, m_raw, p)
        orphan = te.select_orphans(decoyed, picked, p)
        repsnt = {c: {po: po for po in picked[c]} for c in picked}
        disc_polya = {c: {po: (1, 1) for po in picked[c]} for c in picked}
        cands = te.call_candidates(repsnt, picked, clip_slct, disc_polya, orphan,
                                  f_is_cluster, p)
        m_depth = {c: {po: [10, 10, 10, 10] for po in cands[c]} for c in cands}
        got = te.apply_site_filters(cands, m_depth, m_clip_polya, lambda c, po: False, p)

        # legacy
        clip_slct_l = _legacy_filter_by_clip_algnmt(m_l, m_r, m_clip_polya, p.n_clip_cutoff,
                                                   p.n_polya_cutoff)
        decoyed_l = _legacy_chk_decoy_cns_cnt(m_td, m_alu, m_l1, m_sva, p.n_decoy_cutoff,
                                             p.i_rep_type)
        picked_l = _legacy_pick_trsdct_rescued(decoyed_l, p.n_half_disc, m_raw,
                                              p.f_min_uniq_ratio)
        orphan_l = _legacy_slct_orphan_transduction(decoyed_l, picked_l, p.n_half_disc)
        repsnt_l = {c: {po: po for po in picked_l[c]} for c in picked_l}
        disc_polya_l = {c: {po: (1, 1) for po in picked_l[c]} for c in picked_l}
        cands_l = _legacy_call_out_candidate_transduct(repsnt_l, picked_l, i_concord_dist,
                                                      f_ratio, clip_slct_l, disc_polya_l,
                                                      orphan_l)
        m_depth_l = {c: {po: [10, 10, 10, 10] for po in cands_l[c]} for c in cands_l}
        got_l = _legacy_filter_transduction(cands_l, lambda c, po: False, m_depth_l,
                                           p.i_max_cov, m_clip_polya, p.n_half_clip)
        assert _ordered(got) == _ordered(got_l), "seed {0}".format(seed)


# ======================================================================================
# The preserved defect, made explicit
# ======================================================================================

def test_the_stale_side_flag_is_reproduced():
    """A locus can come out with BOTH sides flagged, from different reads.

    Construction: three tight left-side hits make the left list cluster while it is the
    larger one; then four right-side hits take over the majority, so the right side is
    evaluated last. `b_l_trsdct` keeps the True it got at read 3 and the boundaries belong to
    the right side. `collect_ins_info` reads exactly this pair to decide which side's cns
    interval to blank, so the defect is load-bearing and must not be tidied away without a
    truth-set measurement (see `dominant_side_cluster`'s docstring).
    """
    clip_pos = 5000
    l_trsdct_pos = ([(1000, 4000, 0)] * 3          # left of clip_pos, tightly clustered
                    + [(9000 + 10 * i, 6000, 0) for i in range(4)])  # right, also clustered
    ck = ClusterChecker()

    def f_is_cluster(l_pos):
        return ck._is_disc_cluster(l_pos, 450, 0.4)

    b_l, b_r, _tl, _tr, n_disc = te.dominant_side_cluster(
        list(l_trsdct_pos), clip_pos, f_is_cluster)
    assert b_l is True and b_r is True, "the stale-flag defect is what the original does"
    assert n_disc == 4, "n_disc comes from the LAST side evaluated, not the larger one"

    # ...and the legacy path agrees, which is the point.
    m_repsnt = {"chr1": {clip_pos: clip_pos}}
    m_disc_td = {"chr1": {clip_pos: ("src", list(l_trsdct_pos), 7)}}
    m_clip = {"chr1": {clip_pos: (9, 3)}}
    m_polya = {"chr1": {clip_pos: (1, 1)}}
    want = _legacy_call_out_candidate_transduct(m_repsnt, m_disc_td, 450, 0.4, m_clip,
                                               m_polya, {})
    rcd = want["chr1"][clip_pos]
    assert rcd[0] is True and rcd[1] is True and rcd[2] == 4


# ======================================================================================
# One test per threshold: a failure names what broke
# ======================================================================================

def test_derived_cutoffs_keep_the_step_s_float_arithmetic():
    p = params(n_clip_cutoff=4, n_disc_cutoff=5)
    assert p.n_polya_cutoff == 1.0
    assert p.n_half_disc == 2.5, "halving must not round -- 2.5 and 3 admit different loci"
    assert p.n_decoy_cutoff == 5.0
    assert p.n_half_clip == 2.0
    # The `<= 1` fallback: halving 2 gives 1, so the step uses the UNHALVED cutoff, and the
    # decoy cutoff is then 4 rather than 2.
    p2 = params(n_clip_cutoff=1, n_disc_cutoff=2)
    assert p2.n_half_disc == 2
    assert p2.n_decoy_cutoff == 4
    assert p2.n_polya_cutoff == 0, "a negative polyA cutoff is clamped to 0"


def test_clip_gate_names_its_three_reasons():
    p = params()
    assert te.gate_clip_evidence(
        te.SiteEvidence(n_clip=3, n_polya_clip=5), p).reason == te.R_LOW_CLIP
    assert te.gate_clip_evidence(
        te.SiteEvidence(n_clip=9, n_polya_clip=None), p).reason == te.R_NO_POLYA_CLIP
    assert te.gate_clip_evidence(
        te.SiteEvidence(n_clip=9, n_polya_clip=0), p).reason == te.R_LOW_POLYA_CLIP
    assert te.gate_clip_evidence(te.SiteEvidence(n_clip=4, n_polya_clip=1), p).keep, \
        "both cutoffs are inclusive at the boundary"


def test_decoy_gate_has_no_alu_arm():
    """Alu (rtype 2) is not rejected on any decoy count, because this step never runs for it
    -- `is_matched_rep_type` has the same asymmetry and `test_site_universe` pins it there."""
    ev = te.SiteEvidence(n_decoy_alu=99, n_decoy_l1=99, n_decoy_sva=99)
    assert te.gate_decoy_cross_family(ev, params(i_rep_type=2)).keep
    assert te.gate_decoy_cross_family(ev, params(i_rep_type=-1)).keep, \
        "an unknown family disables the gate entirely"
    assert not te.gate_decoy_cross_family(ev, params(i_rep_type=L1_RTYPE)).keep
    assert not te.gate_decoy_cross_family(ev, params(i_rep_type=SVA_RTYPE)).keep
    # L1 rejects on Alu and SVA but NOT on its own decoy count.
    only_l1 = te.SiteEvidence(n_decoy_alu=0, n_decoy_l1=99, n_decoy_sva=0)
    assert te.gate_decoy_cross_family(only_l1, params(i_rep_type=L1_RTYPE)).keep
    assert not te.gate_decoy_cross_family(only_l1, params(i_rep_type=SVA_RTYPE)).keep


def test_source_gate_needs_dominance_and_uniqueness_together():
    p = params()
    hits = [(1, 2, 3)] * 4
    # 4 of 4 from one source, no non-unique reads: passes.
    ev = te.SiteEvidence(sources={"a": hits}, n_non_unique_by_source={})
    assert te.gate_disc_source(ev, p).keep
    # Same 4, but 4 non-unique reads on the same source: 4/8 == 0.5 is NOT > 0.5.
    ev = te.SiteEvidence(sources={"a": hits}, n_non_unique_by_source={"a": 4})
    assert te.gate_disc_source(ev, p).reason == te.R_NO_DOMINANT_SOURCE
    # Spread across three sources: total passes, dominance (2/6 = 0.33 < 0.45) does not.
    ev = te.SiteEvidence(sources={"a": [(1, 2, 3)] * 2, "b": [(1, 2, 3)] * 2,
                                  "c": [(1, 2, 3)] * 2}, n_non_unique_by_source={})
    assert te.gate_disc_source(ev, p).reason == te.R_NO_DOMINANT_SOURCE
    # Two reads total is below n_half_disc (2.5).
    ev = te.SiteEvidence(sources={"a": [(1, 2, 3)] * 2}, n_non_unique_by_source={})
    assert te.gate_disc_source(ev, p).reason == te.R_LOW_DISC


def test_dominant_source_breaks_ties_on_insertion_order():
    a_first = {"a": [1, 2], "b": [3, 4]}
    b_first = {"b": [3, 4], "a": [1, 2]}
    assert te.dominant_source(a_first)[0] == "a"
    assert te.dominant_source(b_first)[0] == "b", \
        "strictly-greater comparison means the first source wins a tie"


def test_sanity_gate_names_its_four_reasons():
    p = params(i_max_cov=100)
    assert te.gate_site_sanity(
        te.SiteEvidence(b_in_same_family_repeat=True), p).reason == te.R_SAME_FAMILY_REPEAT
    assert te.gate_site_sanity(
        te.SiteEvidence(i_rcov2=101), p).reason == te.R_COVERAGE_ISLAND
    assert te.gate_site_sanity(te.SiteEvidence(i_rcov2=100), p).keep, \
        "the coverage test is strict: equal to i_max_cov survives"
    assert te.gate_site_sanity(
        te.SiteEvidence(n_lpolya=3, n_rpolya=3), p).reason == te.R_TWO_SIDED_POLYA
    assert te.gate_site_sanity(
        te.SiteEvidence(n_lpolyt=3, n_rpolyt=3), p).reason == te.R_TWO_SIDED_POLYT
    assert te.gate_site_sanity(te.SiteEvidence(n_lpolya=3, n_rpolya=2), p).keep, \
        "n_half_clip is 2.0 and the test is strict on BOTH sides"
    assert te.gate_site_sanity(te.SiteEvidence(n_lpolya=None), p).keep, \
        "a locus with no polyA record at all is not subject to the veto"


def test_sibling_params_disable_the_polya_veto():
    """`x_orphan_transduction` calls the sanity filter with an empty polyA map and 0 -- the
    empty map is what disables the veto, and `for_site_sanity` carries the 0 explicitly."""
    p = te.TransductionParams.for_site_sanity(i_max_cov=100, n_half_clip=0,
                                             i_rep_type=L1_RTYPE)
    assert p.n_half_clip == 0
    assert te.gate_site_sanity(te.SiteEvidence(n_lpolya=None), p).keep


def test_a_partial_params_object_fails_loud_on_a_gate_it_does_not_own():
    """A None threshold must raise, not admit everything: the sibling reaches gate D with no
    clip/disc cutoffs at all, and a 0 default there would silently pass every locus."""
    import pytest
    p = te.TransductionParams.for_site_sanity(i_max_cov=100, n_half_clip=0, i_rep_type=1)
    with pytest.raises(TypeError):
        te.gate_clip_evidence(te.SiteEvidence(n_clip=1, n_polya_clip=1), p)
    with pytest.raises(TypeError):
        te.gate_genotype(te.SiteEvidence(n_disc_large_indel=1), p)


def test_genotype_gate_uses_sample_coverage_not_local():
    p = params(ave_cov=30.0, max_indel_ratio=0.2)   # threshold 6.0
    assert te.gate_genotype(te.SiteEvidence(n_disc_large_indel=6), p).keep
    assert te.gate_genotype(
        te.SiteEvidence(n_disc_large_indel=7), p).reason == te.R_LARGE_INDEL


def test_ledger_records_drops_only_and_carries_the_phase():
    p = params()
    ledger = te.DecisionLedger()
    m_l = {"chr1": {100: {5: ["r"] * 9}, 200: {5: ["r"]}}}
    m_polya = {"chr1": {100: [1, 0, 0, 0, 0, 0], 200: [1, 0, 0, 0, 0, 0]}}
    kept = te.select_clip_sites(m_l, {}, m_polya, p, ledger=ledger)
    assert list(kept["chr1"]) == [100], "the 1-read locus is below n_clip_cutoff"
    assert [r[:4] for r in ledger.rows] == [("chr1", 200, "A_clip", te.R_LOW_CLIP)]
    assert ledger.counts == {te.R_LOW_CLIP: 1}
    assert ledger.summary() == "low_clip_support=1"


def test_ledger_writes_a_tsv_with_a_header(tmp_path):
    ledger = te.DecisionLedger()
    ledger.record(te.SiteEvidence(chrm="chr1", pos=42),
                  te.Verdict(False, te.R_COVERAGE_ISLAND, "cov=9"))
    out = ledger.write(str(tmp_path / "led.tsv"))
    lines = open(out).read().splitlines()
    assert lines[0] == "#chrom\tpos\tphase\treason\tdetail"
    assert lines[1] == "chr1\t42\tD_sanity\tcoverage_island\tcov=9"


def test_every_reason_has_a_phase():
    """A reason with no phase writes "?" into the ledger, which is a silent hole in the
    attribution the ledger exists to provide."""
    reasons = [v for k, v in vars(te).items() if k.startswith("R_") and isinstance(v, str)]
    assert reasons, "the reason constants moved -- this test cannot see them"
    missing = [r for r in reasons if r not in te.PHASE_OF_REASON]
    assert not missing, "reasons with no phase: {0}".format(missing)


def test_a_repeat_dropped_locus_is_never_looked_up_in_the_depth_map():
    """The original tests the repeat context BEFORE reading the coverage map, so a locus in a
    same-family repeat is dropped without a depth lookup. Preserving that short-circuit is
    what keeps `apply_site_filters` from raising where the original did not."""
    p = params(i_max_cov=100)
    m_ins = {"chr1": {10: "payload", 20: "payload"}}
    m_depth = {"chr1": {20: [1, 1, 1, 1]}}      # 10 is deliberately ABSENT
    got = te.apply_site_filters(m_ins, m_depth, {}, lambda c, po: po == 10, p)
    assert _ordered(got) == [("chr1", [(20, "payload")])]


# ======================================================================================
# WIRING. The differentials above test the LIBRARY; these test the two methods of
# x_transduction.py that the sibling step shares, in process, with the BAM work stubbed.
#
# They exist because of a real bug this file did not catch: `filter_transduction` was left
# returning THREE values, which every differential above passes through and which only the
# end-to-end gate saw (`ValueError: too many values to unpack`, job 49633492). Nothing in the
# fast suite called these methods at all, so their contract with their callers was untested.
# ======================================================================================

class _FakeDepth:
    """`ReadDepth` with the BAM removed -- returns a fixed coverage per requested locus."""

    def __init__(self, *_a, **_k):
        pass

    def calc_coverage_of_two_regions(self, m_sites, _bams, _search, _focal, _focal2):
        return {c: {p: [7, 7, 7, 7] for p in m_sites[c]} for c in m_sites}


class _FakeAnnotationNoRepeat:
    def is_within_repeat_region_interval_tree2(self, _chrm, _pos):
        return []


def _xtransduction(monkeypatch, tmp_path):
    import x_transduction
    monkeypatch.setattr(x_transduction, "ReadDepth", _FakeDepth)
    return x_transduction.XTransduction(str(tmp_path) + "/", 1, None)


def test_filter_transduction_returns_exactly_the_two_values_its_callers_unpack(monkeypatch,
                                                                              tmp_path):
    """Both callers do `m_new_ins, m_depth = self.filter_transduction(...)`, and
    collect_ins_info reads cov2 out of that second value."""
    xtd = _xtransduction(monkeypatch, tmp_path)
    m_ori = {"chr1": {100: ("payload", 1), 200: ("payload", 2)}}
    rslt = xtd.filter_transduction(m_ori, _FakeAnnotationNoRepeat(), 1, "bam_list.txt",
                                  100, {}, 2.0)
    assert isinstance(rslt, tuple) and len(rslt) == 2, \
        "the callers unpack exactly two values"
    m_new_ins, m_depth = rslt
    assert _ordered(m_new_ins) == _ordered(m_ori), "nothing here should have been dropped"
    assert m_depth["chr1"][100][2] == 7, "the depth map must carry cov2 for collect_ins_info"


def test_filter_transduction_drops_on_coverage_and_records_the_reason(monkeypatch, tmp_path):
    xtd = _xtransduction(monkeypatch, tmp_path)
    ledger = te.DecisionLedger()
    m_ori = {"chr1": {100: "payload"}}
    m_new_ins, _m_depth = xtd.filter_transduction(
        m_ori, _FakeAnnotationNoRepeat(), 1, "bam_list.txt", 6, {}, 2.0, ledger=ledger)
    assert m_new_ins == {}, "coverage 7 is above i_max_cov 6"
    assert [r[3] for r in ledger.rows] == [te.R_COVERAGE_ISLAND]


def test_add_new_rcd_with_cutoff_still_writes_records_and_gates_on_large_indel(monkeypatch,
                                                                              tmp_path):
    """The genotype gate's caller contract: one output line per surviving site, and the
    large-indel rejection still fires. `ave_cov` 30 * 0.2 -> a threshold of 6."""
    xtd = _xtransduction(monkeypatch, tmp_path)
    gntp = [0, 0, 0, 0, 0, 0, 0, "", 0, 0]            # index 6 is n_disc_large_indel
    keep = list(gntp); keep[6] = 6                     # 6 is NOT > 6.0
    drop = list(gntp); drop[6] = 7
    m_gntp = {"chr1": {100: keep, 200: drop}}
    site = (1, -1, -1, 1, -1, -1, 2, 1, 7, 7, "chr9:1-2~0")
    m_site_info = {"chr1": {100: site, 200: site}}
    out = str(tmp_path / "new_rcds.txt")
    ledger = te.DecisionLedger()
    xtd.add_new_rcd_with_cutoff(m_site_info, m_gntp, 30.0, out, ledger=ledger)
    lines = [l for l in open(out).read().splitlines() if l.strip()]
    assert len(lines) == 1 and lines[0].split()[1] == "100"
    assert [r[3] for r in ledger.rows] == [te.R_LARGE_INDEL]


def test_a_locus_absent_from_the_source_map_is_skipped_not_recorded():
    """`call_candidates`' one silent skip, pinned.

    A locus that passes clip + disc-polyA but has no entry in the picked-source map produces
    no record and no ledger row. That branch is UNREACHABLE from the transduction step --
    `m_repsnt_pos` is built as `prep_representative_pos(m_disc_filtered)` and `m_disc_filtered`
    is a filtered subset of `m_disc_transduct`'s keys, measured 0 exceptions for both families
    on the WGS gate (job 49634290) -- so it is a guard for a future caller with its own site
    list, not a population. This test exists so it stays a deliberate skip: without it, a
    refactor that turned the membership test into a plain lookup would raise KeyError, ~40
    minutes into a WGS run, on a path nothing exercises.
    """
    p = params()
    m_repsnt = {"chr1": {500: 500}}
    m_clip = {"chr1": {500: (9, 3)}}
    m_polya = {"chr1": {500: (1, 1)}}
    got = te.call_candidates(m_repsnt, {}, m_clip, m_polya, {},
                            lambda l_pos: (True, 1, 2), p, ledger=(led := te.DecisionLedger()))
    assert got == {}, "no record for a locus with no picked source"
    assert led.rows == [], "and no ledger row -- there is no gate to attribute it to"
