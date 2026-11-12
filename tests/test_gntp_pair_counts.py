"""`collect_features_one_site` must count read PAIRS, once each, on both sides of AF_PE.

n_concd_pairs (AF_concod, the ref side of AF_PE) and n_disc_pairs (AF_disc, its alt side)
feed the beta-binomial genotyper in x_genotype_lik through AF_PE = disc / (disc + concd).

THE BUG. The concordant tally walked `l_check_concord`, which holds one entry per
ALIGNMENT, and added 1 for every entry that passed `is_concrdant`. `is_concrdant` requires
the pair's span (leftmost start .. rightmost start + READ_LENGTH) to contain ins_pos with
both ends within DFT_IS of it, and the fetch window is ins_pos +- DFT_IS, so a qualifying
pair almost always has BOTH mates in the window and was counted twice. On HG002 30x
NovaSeq (mean insert 419, max_is 954) the counted/unique ratio was ~2.0 at every insert
size, which halved AF_PE's alt fraction for a het. The discordant tally was per alignment
too; its mate is far away, but a hard-clipped supplementary record of the anchor read
(CIGAR ...H, so not skipped as a soft clip) counted the same template a second time.

Both tallies are now sets of query names. Each test below fails on the old per-alignment
code. The BAM is synthetic (two 20 kb contigs, arbitrary sequence), so the suite needs no
external data.
"""
import os
import sys

import pytest

pysam = pytest.importorskip("pysam")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))
import global_values  # noqa: E402
from x_genotype_feature import XGenotyper  # noqa: E402

CHRM = "chr1"
INS = 10000
IS = 500          # DFT_IS: fetch window is INS +- IS
RLEN = 100        # READ_LENGTH used by is_concrdant
SEQ = "ACGTTGCA" * 13  # 104 bp, sliced to each record's query length

# flags
P, PROPER, REV, MREV, R1, R2, SEC, DUP, SUPP = 0x1, 0x2, 0x10, 0x20, 0x40, 0x80, 0x100, 0x400, 0x800


def _rec(name, flag, pos, cigar, mate_pos, mate_chrm=CHRM, chrm=CHRM):
    return dict(name=name, flag=flag, chrm=chrm, pos=pos, cigar=cigar,
                mate_chrm=mate_chrm, mate_pos=mate_pos)


def _pair(name, pos1, pos2, cigar1="100M", cigar2="100M", dup2=False):
    """A forward/reverse pair on CHRM. dup2 marks mate 2 a duplicate (skipped by the code)."""
    return [_rec(name, P | PROPER | MREV | R1, pos1, cigar1, pos2),
            _rec(name, P | PROPER | REV | R2 | (DUP if dup2 else 0), pos2, cigar2, pos1)]


def _write_bam(path, records):
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": "chr1", "LN": 20000}, {"SN": "chr2", "LN": 20000}]}
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        tid = {"chr1": 0, "chr2": 1}
        for r in sorted(records, key=lambda x: (tid[x["chrm"]], x["pos"])):
            a = pysam.AlignedSegment(out.header)
            a.query_name = r["name"]
            a.flag = r["flag"]
            a.reference_id = tid[r["chrm"]]
            a.reference_start = r["pos"]
            a.mapping_quality = 60
            a.cigarstring = r["cigar"]
            qlen = sum(l for op, l in a.cigartuples if op in (0, 1, 4))
            a.query_sequence = SEQ[:qlen]
            a.query_qualities = pysam.qualitystring_to_array("I" * qlen)
            a.next_reference_id = tid[r["mate_chrm"]]
            a.next_reference_start = r["mate_pos"]
            out.write(a)
    pysam.index(path)


@pytest.fixture
def run_site(tmp_path, monkeypatch):
    monkeypatch.setattr(global_values, "DFT_IS", IS)
    monkeypatch.setattr(global_values, "READ_LENGTH", RLEN)
    ref = str(tmp_path / "ref.fa")
    with open(ref, "w") as fout:
        for c in ("chr1", "chr2"):
            fout.write(">{0}\n".format(c))
            for i in range(0, 20000, 80):
                fout.write("ACGT" * 20 + "\n")
    pysam.faidx(ref)

    def _run(records, full=False, extnd=IS):
        """extnd is the fetch half-window; a test that changes DFT_IS passes the same value."""
        bam = str(tmp_path / "reads.bam")
        _write_bam(bam, records)
        wd = str(tmp_path) + "/"
        XGenotyper(ref, wd, 1).collect_features_one_site(((CHRM, INS, extnd), bam, wd))
        sf = wd + "{0}_{1}".format(CHRM, INS) + global_values.GNTP_FEATURE_SUFFIX
        with open(sf) as fin:
            f = fin.read().split()
        out = {"disc": int(f[6]), "concd": int(f[7])}
        if full:
            out["full"] = int(f[3])  # n_full_map
        return out
    return _run


# --- scenarios ------------------------------------------------------------------------------
def both_mates_in_window():
    # span [9800, 10200] covers INS; both records lie inside INS +- IS
    return _pair("BOTH", 9800, 10100)


def single_mate_in_window():
    # a spanning pair with only one usable record (mate 2 is a duplicate), and a pair with
    # one mate in the window that does not span INS (its span ends at 9550)
    return _pair("ONE_SPAN", 9850, 10050, dup2=True) + _pair("ONE_NOSPAN", 9200, 9450)


def clipped_with_concordant_mate():
    # CLIP_AT_BRK: mate 1 is right-clipped at 9950, within TSD_CUTOFF of INS -> the pair is
    #   split-read evidence and must not also be a concordant (ref) pair.
    # CLIP_AWAY: mate 1 is right-clipped at 9890 (110 bp from INS, outside TSD_CUTOFF), so the
    #   pair is an ordinary concordant pair and counts once.
    return (_pair("CLIP_AT_BRK", 9900, 10200, cigar1="50M50S")
            + _pair("CLIP_AWAY", 9800, 10100, cigar1="90M10S"))


def discordant_pair():
    # anchor read near INS with its mate on chr2, plus a hard-clipped supplementary record of
    # the same read (bwa's default for supplementaries) also inside the window
    return [_rec("DISC", P | MREV | R1, 9900, "100M", 5000, mate_chrm="chr2"),
            _rec("DISC", P | MREV | R1 | SUPP, 10050, "40M60H", 5000, mate_chrm="chr2"),
            _rec("DISC", P | REV | R2, 5000, "100M", 9900, chrm="chr2")]


# --- tests (old per-alignment counts in the comments) ---------------------------------------
def test_concordant_pair_with_both_mates_in_window_counts_once(run_site):
    assert run_site(both_mates_in_window()) == {"disc": 0, "concd": 1}      # old: concd 2


def test_single_mate_in_window_counts_once_and_nonspanning_not_at_all(run_site):
    # with a two-record pair alongside, so one-per-pair is checked against a pair that the
    # old code would double: old concd = 2 + 1 = 3
    assert run_site(single_mate_in_window() + both_mates_in_window()) == {"disc": 0, "concd": 2}


def test_clipped_read_with_concordant_mate(run_site):
    assert run_site(clipped_with_concordant_mate()) == {"disc": 0, "concd": 1}  # old: concd 2


def test_discordant_pair_counts_once(run_site):
    assert run_site(discordant_pair()) == {"disc": 1, "concd": 0}           # old: disc 2


def test_all_scenarios_together(run_site):
    recs = (both_mates_in_window() + single_mate_in_window()
            + clipped_with_concordant_mate() + discordant_pair())
    # concordant pairs: BOTH, ONE_SPAN, CLIP_AWAY; discordant: DISC  (old: disc 2, concd 5)
    assert run_site(recs) == {"disc": 1, "concd": 3}


# --- one template is alt OR ref, from primary records only (old/first-fix counts in comments) -
def test_alt_template_is_not_also_a_ref_pair(run_site):
    # R2 is anchored left of the site with R1's primary in a TE on chr2 (alt evidence). R1 also
    # has a hard-clipped supplementary piece at the breakpoint whose mate fields point at R2,
    # so that record alone looks like a pair spanning INS ([9700, 10050]).
    recs = [_rec("ALT", P | R2, 9700, "100M", 5000, mate_chrm="chr2"),
            _rec("ALT", P | R1 | REV | SUPP, 9950, "50M50H", 9700),
            _rec("ALT", P | R1 | REV, 5000, "100M", 9700, chrm="chr2")]
    assert run_site(recs) == {"disc": 1, "concd": 0}       # HEAD and first fix: concd 1


def test_pair_counts_if_any_mate_passes(run_site):
    # R2 is right-clipped at 9900 (100 bp from INS: not < TSD_CUTOFF, so not clip evidence).
    # Its record uses the moved clip position: span [9600, 10000] passes. R1's record uses
    # R2's start (9850): span [9600, 9950] does not. One passing record is enough.
    recs = _pair("ANY", 9600, 9850, cigar2="50M50S")
    assert run_site(recs) == {"disc": 0, "concd": 1}


def test_clip_exclusion_does_not_depend_on_record_order(run_site):
    # the clipped record (R2, clip at 10000) sorts AFTER its concordant mate (R1 at 9700)
    recs = _pair("LATE_CLIP", 9700, 9950, cigar2="50M50S")
    assert run_site(recs) == {"disc": 0, "concd": 0}


def test_secondary_alignment_is_not_ref_evidence(run_site):
    # PRI is an ordinary concordant pair; SECX has only a 100M secondary record over INS
    # (9950-10050, more than BWA_HALF_READ_MIN_SCORE from both ends, so it lands in n_full_map)
    recs = both_mates_in_window() + [
        _rec("SECX", P | PROPER | MREV | R1 | SEC, 9950, "100M", 10100)]
    assert run_site(recs, full=True) == {"disc": 0, "concd": 1, "full": 0}
    # HEAD: concd 3, full 1; first fix: concd 2, full 1


def test_same_chromosome_discordant_pair_inside_window_counts_once(run_site, monkeypatch):
    # mates 2050 bp apart (> DISC_THRESHOLD 2000), both within DFT_IS of INS
    monkeypatch.setattr(global_values, "DFT_IS", 1100)
    assert global_values.DISC_THRESHOLD == 2000
    recs = [_rec("FAR", P | MREV | R1, 8950, "100M", 11000),
            _rec("FAR", P | REV | R2, 11000, "100M", 8950)]
    assert run_site(recs, extnd=1100) == {"disc": 1, "concd": 0}  # HEAD: disc 2
