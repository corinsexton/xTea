"""Split-read AF counts supplementary junction reads once per read (B16).

`collect_features_one_site` counts a read as split-read (alt) evidence for AF_SR when it is
clipped within CLIP_EXACT_CLIP_SLACK of the breakpoint (n_af_clip, the AF_CLIP INFO key). It
used to see only soft clips (CIGAR op 4) on primary records. In repeats a junction read's
primary often maps to a TE copy elsewhere, and the piece aligned at the site is a SUPPLEMENTARY
record: hard-clipped (op 5) with bwa defaults, soft-clipped with bwa -Y and similar. Those
reads were missed, so AF_SR was undercounted on the alt side (HG002 NovaSeq het median 0.433;
0.542 with them).

Now a supplementary clipped at the breakpoint counts toward its side's AF clip count, once per
read (query name + R1/R2): not again when the same read is already counted through a
soft-clipped primary at the site. Secondary records never count. The raw clip counts (within
TSD_CUTOFF; LRAWCLIP/RRAWCLIP) and the pair counts are unchanged.
"""
import os
import sys

import pytest

pytest.importorskip("pysam")
sys.path.insert(0, os.path.dirname(__file__))
from test_gntp_pair_counts import (CHRM, INS, MREV, P, R1, R2, REV, SEC, SUPP,  # noqa: E402,F401
                                   _rec, run_site)
import global_values  # noqa: E402  (on sys.path via test_gntp_pair_counts)


@pytest.fixture
def features(run_site, tmp_path):
    def _run(records):
        run_site(records)
        sf = str(tmp_path) + "/{0}_{1}".format(CHRM, INS) + global_values.GNTP_FEATURE_SUFFIX
        f = open(sf).read().split()
        return dict(af_clip=int(f[2]), full=int(f[3]), lraw=int(f[4]), rraw=int(f[5]),
                    disc=int(f[6]), concd=int(f[7]), lens=f[9])
    return _run


def _te_primary(name, flag=P | R1):
    """The read's primary record on chr2 (a TE copy elsewhere); its mate is anchored at 9700."""
    return _rec(name, flag, 5000, "100M", 9700, chrm="chr2")


def test_hard_clipped_supplementary_left_clip_counts(features):
    recs = [_te_primary("HL"), _rec("HL", P | R1 | SUPP, 10000, "60H40M", 9700),
            _rec("HL", P | R2 | MREV, 9700, "100M", 5000, mate_chrm="chr2")]
    out = features(recs)
    assert (out["af_clip"], out["lens"]) == (1, "60")        # old: 0, none
    assert (out["lraw"], out["rraw"]) == (0, 0)              # raw clip counts unchanged
    assert out["disc"] == 1                                  # pair side unchanged


def test_hard_clipped_supplementary_right_clip_counts(features):
    # 40M60H at 9960: aligned part ends at 10000
    recs = [_te_primary("HR"), _rec("HR", P | R1 | SUPP, 9960, "40M60H", 9700)]
    assert features(recs)["af_clip"] == 1                     # old: 0


def test_soft_clipped_supplementary_counts(features):
    # bwa -Y style: the supplementary keeps its clipped bases as a soft clip
    recs = [_te_primary("SL"), _rec("SL", P | R1 | SUPP, 10000, "60S40M", 9700)]
    assert features(recs)["af_clip"] == 1                     # old: 0


def test_supplementary_away_from_breakpoint_does_not_count(features):
    # clipped 3 bp from INS: not < CLIP_EXACT_CLIP_SLACK (3)
    recs = [_te_primary("FAR"), _rec("FAR", P | R1 | SUPP, 10003, "60H40M", 9700)]
    assert global_values.CLIP_EXACT_CLIP_SLACK == 3
    assert features(recs)["af_clip"] == 0


def test_read_counted_once_through_primary_and_supplementary(features):
    # A soft-clipped primary and a supplementary record of the SAME read (R1) are both clipped
    # at INS on the left: one read, counted once (counting records would give 2).
    recs = [_rec("ONE", P | R1, 10000, "40S60M", 5000, mate_chrm="chr2"),
            _rec("ONE", P | R1 | SUPP, 10000, "30H70M", 5000, mate_chrm="chr2"),
            _rec("ONE", P | R2, 5000, "100M", 10000, chrm="chr2")]
    assert features(recs)["af_clip"] == 1


def test_mates_of_one_template_are_two_reads(features):
    # Overlapping mates both cross the junction: R1 as a soft-clipped primary, R2 as a
    # hard-clipped supplementary. Two reads, like two soft-clipped primaries would be.
    recs = [_rec("TPL", P | R1, 10000, "40S60M", 5000, mate_chrm="chr2"),
            _rec("TPL", P | R2 | SUPP, 10000, "30H70M", 10000),
            _rec("TPL", P | R2, 5000, "100M", 10000, chrm="chr2")]
    assert features(recs)["af_clip"] == 2                     # old: 1


def test_secondary_never_counts(features):
    recs = [_te_primary("SECC"), _rec("SECC", P | R1 | SEC, 10000, "60S40M", 9700),
            _rec("SECC", P | R1 | SEC, 9960, "40M60H", 9700)]
    assert features(recs)["af_clip"] == 0


def test_supplementaries_can_make_their_side_dominant(features):
    # left: one soft-clipped primary. right: two hard-clipped supplementaries. The right side
    # is dominant, so n_full_map drops the reads that START near INS (n_r_full_map).
    recs = [_rec("L1", P | R1, 10000, "30S70M", 5000, mate_chrm="chr2"),
            _te_primary("RA"), _rec("RA", P | R1 | SUPP, 9960, "40M60H", 9700),
            _te_primary("RB"), _rec("RB", P | R1 | SUPP, 9950, "50M50H", 9700),
            # ref reads: one starting at INS - 10 (n_r_full_map), one ending at INS + 10
            _rec("FR", P | R1, 9990, "100M", 5000, mate_chrm="chr2"),
            _rec("FL", P | R1, 9910, "100M", 5000, mate_chrm="chr2")]
    out = features(recs)
    assert (out["af_clip"], out["full"]) == (2, 1)            # old: (1, 1) via the left side
    assert sorted(out["lens"].split(":")) == ["50", "60"]
