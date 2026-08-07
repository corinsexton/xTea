"""The transduction decoy block: the recipe, the append rule, and the flag that gates them.

WHY THIS FILE'S CENTRAL TEST IS A BYTE COMPARISON AGAINST SHIPPED DATA.

`xtea/td_decoy_library.py` claims the 45-record decoy consensus block can be rebuilt from
`<rep_lib>/consensus/` on any build. That claim is the whole reason the packaging half of the
SVA decoy defect needs no new data and no library repackaging -- so it is not a claim to
assert with a hand-written fixture. `test_the_recipe_rebuilds_the_packaged_L1_block_byte_for_byte`
compares the rebuilt block against the TAIL OF THE SHIPPED hg38 AND hg19 L1 FLANK LIBRARIES,
character for character, wrapping and header comment fields included. It skips when the
rep_lib directory is absent; it does not weaken.

The hermetic tests below pin the parts that must hold on a build nobody here has: that
appending is by NAME (a duplicate name in a bwa reference is the failure being avoided), that
a missing source record raises rather than silently emitting a short block, and that the OFF
arm builds the same `cat` command it always did.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

import global_values                                                    # noqa: E402
import td_decoy_library as tdl                                          # noqa: E402
import transduction_evidence as te                                      # noqa: E402


# The libraries, and how many decoy records each one PACKAGES. Measured 2026-08-05 against
# ~/xTea/rep_lib_annotation; the census is what makes "hg38 SVA has no decoy channel" a fact
# about the data rather than an inference from the SVA arm's zero.
PACKAGED_LIBRARIES = (
    ("LINE/hg38/hg38_FL_L1_flanks_3k.fa",    707, 45),
    ("LINE/hg19/hg19_FL_L1_flanks_3k.fa",    705, 45),
    ("LINE/chm13/chm13_FL_L1_flanks_3k.fa",  818,  0),
    ("SVA/hg38/hg38_FL_SVA_flanks_3k.fa",   5454,  0),
    ("SVA/hg19/hg19_FL_SVA_flanks_3k.fa",   5158, 44),   # all but HERV-K, in another order
    ("SVA/chm13/chm13_FL_SVA_flanks_3k.fa", 6072,  0),
)


def _consensus(rep_lib_dir):
    return os.path.join(rep_lib_dir, "consensus")


def _require(path):
    if not os.path.isfile(path):
        pytest.skip("not in this rep_lib: {0}".format(path))
    return path


# ======================================================================================
# The recipe, against the shipped libraries
# ======================================================================================

def test_the_recipe_rebuilds_the_packaged_L1_block_byte_for_byte(rep_lib_dir):
    """The claim plan sec 5 rests the packaging fix on, checked as bytes on two builds.

    Both L1 libraries carry the block as their final 45 records, so the rebuilt text must be
    the file's tail from `>LINE1` onward EXACTLY -- 60-column wrapping from
    `transduction_decoy.fa`, 70 from `SVA.fa`, 71 from `ALU.fa`, 120 from `HERV.fa`, and the
    `SINE    Homo sapiens` comment fields the headers carry. Anything that normalises the
    sequence (pysam.FastaFile does) or the header (taking the first token does) fails here.
    """
    text = tdl.decoy_block_text(tdl.build_decoy_block(_consensus(rep_lib_dir)))
    for s_rel, _, n_decoy in PACKAGED_LIBRARIES:
        if n_decoy != tdl.N_DECOY_RECORDS:
            continue
        s_lib = open(_require(os.path.join(rep_lib_dir, s_rel))).read()
        i = s_lib.index(">LINE1\n")
        assert s_lib[i:] == text, s_rel


def test_the_hg19_SVA_library_packages_the_same_block_minus_HERV_K(rep_lib_dir):
    """The 44 are not a different set, they are this set less one record -- so the recipe is
    the libraries' recipe and not the L1 libraries'. Compared as an unordered (name, sequence)
    mapping because hg19's SVA library orders the block differently (the 35 Alu before the 6
    SVA), which is also why the appended layout cannot claim to be THE packaged one."""
    m_recipe = {tdl.record_name(h): "".join(l)
                for h, l in tdl.build_decoy_block(_consensus(rep_lib_dir))}
    s_lib = _require(os.path.join(rep_lib_dir, "SVA/hg19/hg19_FL_SVA_flanks_3k.fa"))
    m_packaged = {tdl.record_name(h): "".join(l) for h, l in tdl.read_fasta_raw(s_lib)
                  if not te.is_flank_record(tdl.record_name(h))}
    assert m_packaged == {k: v for k, v in m_recipe.items() if k != "HERV-K"}


def test_the_library_census_is_what_the_packaging_argument_says_it_is(rep_lib_dir):
    """The table in `td_decoy_library`'s docstring, as an assertion. If a rep_lib update ever
    gives hg38's SVA library a decoy block of its own, this fails and the whole item is moot --
    which is a thing worth being told rather than discovering from a null result."""
    for s_rel, n_total, n_decoy in PACKAGED_LIBRARIES:
        recs = tdl.read_fasta_raw(_require(os.path.join(rep_lib_dir, s_rel)))
        names = [tdl.record_name(h) for h, _ in recs]
        assert len(recs) == n_total, s_rel
        assert sum(1 for n in names if not te.is_flank_record(n)) == n_decoy, s_rel


def test_appending_to_a_packaged_L1_library_appends_nothing(rep_lib_dir, tmp_path):
    """The no-op leg of the gate, stated at the level that makes it structural: not "the output
    compared equal" but "zero records were selected for appending"."""
    sf_out = str(tmp_path / "decoy_only.fa")
    added = tdl.write_missing_decoys(
        _require(os.path.join(rep_lib_dir, "LINE/hg38/hg38_FL_L1_flanks_3k.fa")),
        _consensus(rep_lib_dir), sf_out)
    assert added == []
    assert open(sf_out).read() == ""


def test_appending_to_the_hg38_SVA_library_appends_all_45(rep_lib_dir, tmp_path):
    """And the call-moving leg: on the one library this project measures, every record is new."""
    sf_out = str(tmp_path / "decoy_only.fa")
    added = tdl.write_missing_decoys(
        _require(os.path.join(rep_lib_dir, "SVA/hg38/hg38_FL_SVA_flanks_3k.fa")),
        _consensus(rep_lib_dir), sf_out)
    assert len(added) == 45 == tdl.N_DECOY_RECORDS
    assert added[:4] == ["LINE1", "SVA", "ALU", "HERV-K"]
    # And every appended name is one the step's own classifier calls a decoy, under BOTH
    # policies -- otherwise the block would be added and then read as 45 new SOURCES.
    for s_name in added:
        assert te.decoy_tags(s_name, True), s_name
        assert te.decoy_tags(s_name, False), s_name


def test_the_block_partitions_cleanly_into_the_four_decoy_TAGS(rep_lib_dir):
    """What the block MEANS to the disc gate, and it is the whole point of adding it to SVA.

    Every record must carry exactly ONE tag -- an ambiguous name would be counted into two
    cross-family tallies from one hit -- and the partition must be the one
    `gate_decoy_cross_family`'s SVA arm can act on: it rejects on ALU and L1 counts and NOT on
    its own, so the 36 ALU records and the 1 LINE1 record are what make plan sec 5's
    SVA-is-SINE-VNTR-Alu hazard testable, while the 7 SVA records only stop a hit from being
    read as a specific reference copy's flank.
    """
    from collections import Counter
    counts = Counter()
    for s_header, _ in tdl.build_decoy_block(_consensus(rep_lib_dir)):
        tags = te.decoy_tags(tdl.record_name(s_header), True)
        assert len(tags) == 1, (s_header, tags)
        counts[tags[0]] += 1
    assert counts == Counter({global_values.TD_DECOY_ALU: 36,   # ALU + 35 ALU_Alu*
                              global_values.TD_DECOY_SVA: 7,    # SVA + SVA_A..F
                              global_values.TD_DECOY_LINE: 1,
                              global_values.TD_DECOY_HERV: 1})


def test_appending_to_the_hg19_SVA_library_appends_only_HERV_K(rep_lib_dir, tmp_path):
    """The partial case, which is the one a by-name append exists for: 44 present, 1 missing."""
    sf_out = str(tmp_path / "decoy_only.fa")
    added = tdl.write_missing_decoys(
        _require(os.path.join(rep_lib_dir, "SVA/hg19/hg19_FL_SVA_flanks_3k.fa")),
        _consensus(rep_lib_dir), sf_out)
    assert added == ["HERV-K"]


# ======================================================================================
# Hermetic: the parts that must hold on a build this machine does not have
# ======================================================================================

def _fake_consensus(tmp_path):
    """A consensus directory with the four source files, tiny but structurally real: the
    header comment fields and the irregular wrapping are the properties under test."""
    d = tmp_path / "consensus"
    d.mkdir()
    (d / "transduction_decoy.fa").write_text(
        ">LINE1\nAAAA\nAA\n>SVA\nCCCC\n>ALU\nGGGG\n")
    (d / "HERV.fa").write_text(">HERV-K\nTTTTTT\n")
    (d / "SVA.fa").write_text(
        "".join(">SVA_{0}  SINE    Homo sapiens\nAC{0}\n".format(c) for c in "ABCDEF")
        + ">SVA\nNNNN\n>SVA_S1\nNNNN\n")
    (d / "ALU.fa").write_text(
        "".join(">{0}\tSINE1/7SL\nGT{1}\n".format(n, i)
                for i, n in enumerate(tdl.DECOY_SOURCES[3][1]))
        + ">ALU\nNNNN\n>ALU_v1\nNNNN\n")
    return str(d)


def test_the_block_is_45_records_in_the_packaged_order(tmp_path):
    block = tdl.build_decoy_block(_fake_consensus(tmp_path))
    assert len(block) == tdl.N_DECOY_RECORDS == 45
    names = [tdl.record_name(h) for h, _ in block]
    assert names[:10] == ["LINE1", "SVA", "ALU", "HERV-K",
                          "SVA_A", "SVA_B", "SVA_C", "SVA_D", "SVA_E", "SVA_F"]
    assert names[10] == "ALU_AluJo" and names[-1] == "ALU_AluYb8a1"


def test_records_are_taken_by_name_so_a_reordered_consensus_file_cannot_move_the_block(tmp_path):
    """The reason `DECOY_SOURCES` lists 45 names instead of slicing the first N: the block ends
    up in a bwa reference nothing downstream re-checks, so a positional cut would go wrong
    silently. Here `ALU.fa` is reversed and the block must be unchanged."""
    d = _fake_consensus(tmp_path)
    want = tdl.decoy_block_text(tdl.build_decoy_block(d))
    recs = tdl.read_fasta_raw(os.path.join(d, "ALU.fa"))
    with open(os.path.join(d, "ALU.fa"), "w") as fout:
        fout.write(tdl.decoy_block_text(list(reversed(recs))))
    assert tdl.decoy_block_text(tdl.build_decoy_block(d)) == want


def test_a_missing_record_raises_instead_of_emitting_a_short_block(tmp_path):
    """A partial block would index cleanly, move some calls, and leave nothing in the run that
    says which records were there. Both failure modes name the thing that is missing."""
    d = _fake_consensus(tmp_path)
    recs = [(h, l) for h, l in tdl.read_fasta_raw(os.path.join(d, "SVA.fa"))
            if tdl.record_name(h) != "SVA_D"]
    with open(os.path.join(d, "SVA.fa"), "w") as fout:
        fout.write(tdl.decoy_block_text(recs))
    with pytest.raises(ValueError, match="SVA_D"):
        tdl.build_decoy_block(d)
    os.remove(os.path.join(d, "HERV.fa"))
    with pytest.raises(ValueError, match="HERV.fa"):
        tdl.build_decoy_block(d)


def test_the_append_is_by_name_and_not_by_sequence(tmp_path):
    """A library that carries `>SVA` with a DIFFERENT sequence must still not gain a second
    record called `SVA`: `reference_name` is the string the entire decoy classification tests,
    so a duplicate name is worse than a stale sequence."""
    d = _fake_consensus(tmp_path)
    sf_lib = str(tmp_path / "lib.fa")
    open(sf_lib, "w").write(">chr1~100~200~SVA_D~0L\nACGT\n>SVA\nTTTTTTTT\n")
    added = tdl.write_missing_decoys(sf_lib, d, str(tmp_path / "out.fa"))
    assert "SVA" not in added
    assert len(added) == 44
    # ...and the flank record's own subfamily name never collides, whatever it is called.
    assert "chr1~100~200~SVA_D~0L" not in added


def test_the_name_is_the_first_header_token(tmp_path):
    """bwa reports `reference_name`, not the header. The packaged block keeps `>SVA_A  SINE
    Homo sapiens` verbatim, so these are different strings and only one of them is comparable
    against a SAM record."""
    assert tdl.record_name("SVA_A  SINE    Homo sapiens") == "SVA_A"
    assert tdl.record_name("AluYk13\tSINE1/7SL\tHomo sapiens") == "AluYk13"
    assert tdl.record_name("LINE1") == "LINE1"
    assert tdl.record_name("") == ""


# ======================================================================================
# The seam: construct_novel_flanks, and the flag
# ======================================================================================

class _FakeXRef(object):
    """`gnrt_flank_regions_of_polymerphic_insertions` needs a real genome; the seam under test
    does not care what it writes, only that its output is the last file in the `cat`."""

    def gnrt_flank_regions_of_polymerphic_insertions(self, l_fl, i_len, sf_ref, sf_out):
        open(sf_out, "w").write(">5~77881127~77881127~polymerphic~0L\nACGTACGT\n")


def _transduction(tmp_path, monkeypatch, b_concat):
    import x_reference
    import x_transduction
    monkeypatch.setenv(global_values.ENV_TD_DECOY_CONCAT, "1" if b_concat else "0")
    monkeypatch.setattr(x_reference, "XReference", _FakeXRef)
    monkeypatch.setattr(x_transduction, "XReference", _FakeXRef, raising=False)
    return x_transduction.XTransduction(str(tmp_path) + "/", 1, "ref.fa")


def _flank_library(tmp_path):
    sf_lib = str(tmp_path / "flanks.fa")
    open(sf_lib, "w").write(">chr13~51077120~51078923~SVA_D~0L\nACGTACGTAC\n"
                            ">chr13~51077120~51078923~SVA_D~0R\nTGCATGCATG\n")
    return sf_lib


def test_the_env_var_reads_absent_as_OFF(tmp_path, monkeypatch):
    """Same rule as `ENV_TD_DECOY_STRUCTURAL` and for the same reason -- a direct `x_TEA_main`
    invocation, which every frozen harness under docs/handoff/runs/ uses, must not be moved by
    a driver-level policy. Unlike its sibling this one is also OFF on the driver path."""
    import x_transduction
    for s_val, b_want in (("", False), ("0", False), ("1", True), ("yes", True)):
        monkeypatch.setenv(global_values.ENV_TD_DECOY_CONCAT, s_val)
        assert x_transduction.XTransduction(
            str(tmp_path) + "/", 1, "r.fa").b_decoy_concat is b_want, repr(s_val)
    monkeypatch.delenv(global_values.ENV_TD_DECOY_CONCAT, raising=False)
    assert x_transduction.XTransduction(
        str(tmp_path) + "/", 1, "r.fa").b_decoy_concat is False


def test_the_OFF_arm_builds_the_reference_it_always_did(tmp_path, monkeypatch, rep_lib_dir):
    """The rollback, as bytes: with the flag off the new reference is the library followed by
    the polymorphic flanks and nothing else, whether or not `sf_cns` was passed."""
    xt = _transduction(tmp_path, monkeypatch, False)
    sf_lib = _flank_library(tmp_path)
    sf_cns = os.path.join(_consensus(rep_lib_dir), "SVA.fa")
    sf_out = str(tmp_path / "novel.fa")
    xt.construct_novel_flanks([], sf_lib, 3000, sf_out, sf_cns=sf_cns)
    assert open(sf_out).read() == (open(sf_lib).read()
                                   + open(sf_out + ".polymerphic_only.fa").read())


def test_the_ON_arm_puts_the_block_between_the_library_and_the_polymorphic_flanks(
        tmp_path, monkeypatch, rep_lib_dir):
    """The layout claim: library, then the decoy block, then the polymorphic flanks -- which is
    where hg19's SVA library carries its own 44, so the appended arrangement reproduces a
    packaged one rather than inventing a third."""
    xt = _transduction(tmp_path, monkeypatch, True)
    sf_lib = _flank_library(tmp_path)
    sf_out = str(tmp_path / "novel.fa")
    xt.construct_novel_flanks([], sf_lib, 3000, sf_out,
                              sf_cns=os.path.join(_consensus(rep_lib_dir), "SVA.fa"))
    block = tdl.decoy_block_text(tdl.build_decoy_block(_consensus(rep_lib_dir)))
    assert open(sf_out).read() == (open(sf_lib).read() + block
                                   + open(sf_out + ".polymerphic_only.fa").read())
    # every record name is unique -- the property a by-name append exists to protect
    names = [tdl.record_name(h) for h, _ in tdl.read_fasta_raw(sf_out)]
    assert len(names) == len(set(names))


def test_the_ON_arm_is_a_no_op_on_a_library_that_already_carries_the_block(
        tmp_path, monkeypatch, rep_lib_dir):
    """L1's leg. The ON and OFF arms must produce the same bytes here, which is what lets a
    two-arm gate read an L1 difference as a bug rather than as the change working."""
    sf_lib = _require(os.path.join(rep_lib_dir, "LINE/hg38/hg38_FL_L1_flanks_3k.fa"))
    sf_cns = os.path.join(_consensus(rep_lib_dir), "LINE1.fa")
    l_out = []
    for b_concat in (False, True):
        xt = _transduction(tmp_path, monkeypatch, b_concat)
        sf_out = str(tmp_path / "novel_{0}.fa".format(int(b_concat)))
        xt.construct_novel_flanks([], sf_lib, 3000, sf_out, sf_cns=sf_cns)
        l_out.append(open(sf_out).read())
    assert l_out[0] == l_out[1]
    assert l_out[0].startswith(open(sf_lib).read()[:200])


def test_an_armed_flag_with_no_consensus_dir_is_loud_and_still_builds_the_OFF_reference(
        tmp_path, monkeypatch, capsys):
    """Both shipping call sites pass `sf_cns`, so this is a caller bug -- and the failure mode
    it guards is the expensive one: a run that prints the ON banner and builds the OFF arm's
    reference is a gate measuring the wrong configuration with no way to know."""
    xt = _transduction(tmp_path, monkeypatch, True)
    sf_lib = _flank_library(tmp_path)
    sf_out = str(tmp_path / "novel.fa")
    xt.construct_novel_flanks([], sf_lib, 3000, sf_out)
    assert "REQUESTED BUT NOT APPLIED" in capsys.readouterr().err
    assert open(sf_out).read() == (open(sf_lib).read()
                                   + open(sf_out + ".polymerphic_only.fa").read())


def test_the_driver_defaults_the_concatenation_ON():
    """The sibling of `test_the_driver_defaults_the_decoy_test_ON`, and since 2026-08-07 the
    SAME answer. Both flags address the same defect -- the substring test from the code side,
    this from the packaging side. This one changes the bwa reference the step aligns to, so it
    was held OFF until it had its own two-scale arm; that arm exists and passes (chr1 49879284,
    WGS 49879677, Stage 3 no-regression 49995743), and it moves exactly one call genome-wide,
    which the truth set counts a false positive.

    The assertion that matters is the FIRST one: `None` is "neither flag given", and it is the
    only one of the three a user hits without typing anything."""
    import run_xtea
    assert run_xtea.resolve_decoy_concat(None) is True, "the DEFAULT is ON as of 2026-08-07"
    assert run_xtea.resolve_decoy_concat(True) is True
    assert run_xtea.resolve_decoy_concat(False) is False, (
        "--no-td-decoy-concat must still reach OFF -- it is how BASELINE-0 reproduces")


def test_both_concat_flags_carry_an_explicit_default_None():
    """Same declaration contract as the structural pair, for the same banner reason: with
    `store_true` alone, "neither flag given" and `--no-td-decoy-concat` both arrive as False and
    the banner cannot say which one a harness ran."""
    import inspect
    import run_xtea
    src = inspect.getsource(run_xtea.main)
    for opt in ("--td-decoy-concat", "--no-td-decoy-concat"):
        i = src.index('"{0}"'.format(opt))
        assert "default=None" in src[i:src.index(")\n", i)], opt


def test_the_driver_writes_the_policy_where_the_subprocess_can_read_it(tmp_path):
    """The step is its own `x_TEA_main` subprocess, so the flag has to cross a process boundary.
    Checked at the source level -- the driver must write ENV_TD_DECOY_CONCAT on EVERY path, not
    only when the flag is on, or an inherited value from an earlier arm survives into this run."""
    import inspect
    import run_xtea
    src = inspect.getsource(run_xtea.main)
    assert "os.environ[global_values.ENV_TD_DECOY_CONCAT] = " in src
    assert src.count("ENV_TD_DECOY_CONCAT") == 1, "one writer, unconditional"
