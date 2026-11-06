"""Unit tests for xtea/locus_table.py (Phase C move A.1/A.2).

Pin the type-resolution rule (+-window bp cluster; higher evidence score wins,
tie -> non-orphan_or_sibling_transduction, tie -> L1 > SVA > Alu) in isolation from master's
SampleVCFMerger, whose real implementation has a position-keyed-dict bug (B5,
docs/clipdisc_contract_spec.md sec 4) this module deliberately does not
reproduce: a losing conflict at a different position than the winner would be
written anyway (duplicate), and one at the exact same position silently
clobbers the winner regardless of priority/TPRT. These tests assert the
*correct* (documented, intended) outcome in both situations.

Move A.2 retargeted the resolution point from the three types' final VCFs to
their raw postF_hc candidate-list files (native xTea txt format), so the rows
here are synthetic postF_hc lines, not VCF lines.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "xtea"))

from locus_table import (ALL_CALLS_HEADER, LocusTable, TYPE_TO_RTYPE,  # noqa: E402
                          assign_and_filter, assign_and_write, evidence_components,
                          evidence_score, is_orphan, load_other_families, load_rows,
                          load_type_map)

ORPHAN = "orphan_or_sibling_transduction"


def _row_fields(chrm, lclip_pos, rclip_pos, tprt, ev=None):
    """A 33-column postF_hc-format row. Only the columns the merge reads are set:
    [0] chrm, [2]/[3] clip positions, [4] TSD length, [9]/[10] poly-A reads,
    [19..22] clusters on consensus, [32] subtype. ev: optional dict with
    cls (0-4), polya (reads), tsd (length), subtype. No ev = score 0."""
    ev = ev or {}
    f = ["0"] * 33
    f[0], f[2], f[3] = chrm, str(lclip_pos), str(rclip_pos)
    f[4] = str(ev.get("tsd", -1))
    f[9] = str(ev.get("polya", 0))
    for i in range(19, 23):
        f[i] = "10:20" if i - 19 < ev.get("cls", 0) else "-1:-1"
    f[32] = ev.get("subtype", "two_side_tprt_both" if tprt else "one_end_consistent")
    return f


def _hc_pf(tmp_path, name, rows):
    """rows: list of (chrm, lclip_pos, rclip_pos, tprt[, ev]) -> writes a
    whitespace-separated postF_hc-format file (see _row_fields), returns its path."""
    sf = tmp_path / name
    with open(sf, "w") as fout:
        for r in rows:
            fout.write("\t".join(_row_fields(*r)) + "\n")
    return str(sf)


def test_no_conflict_keeps_everything(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 100, 110, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 5000, 5010, False)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert {(l.chrm, l.pos, l.type) for l in loci} == {("chr1", 100, "L1"), ("chr1", 5000, "Alu")}


def test_load_order_priority_wins_on_equal_evidence(tmp_path):
    # Same score, neither orphan -> L1 (higher priority) wins, Alu is dropped.
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False)])  # 10bp away, within window
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1
    assert (loci[0].chrm, loci[0].pos, loci[0].type) == ("chr1", 1000, "L1")


def test_tprt_flag_alone_no_longer_overrides_priority(tmp_path):
    # The old rule let a two_side_tprt_both call displace a higher-priority family.
    # Rule B reads evidence only: equal scores -> priority -> L1.
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, True)])
    loci = LocusTable().assign_type({"L1": l1, "Alu": alu})
    assert [(l.pos, l.type) for l in loci] == [(1000, "L1")]


def test_higher_score_beats_higher_priority(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False, {"cls": 2})])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False, {"cls": 2, "polya": 3})])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "Alu": alu})
    assert [(l.pos, l.type) for l in loci] == [(1010, "Alu")]
    c = table.conflicts[0]
    assert (c.winner_type, c.loser_type, c.rule, c.winner_score, c.loser_score) == \
        ("Alu", "L1", "score", 3, 2)


def test_lower_score_does_not_displace_placed_call(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False, {"cls": 4, "tsd": 12})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 1000, 1010, False, {"cls": 3, "polya": 5})])
    loci = LocusTable().assign_type({"L1": l1, "SVA": sva})
    assert [l.type for l in loci] == ["L1"]


def test_tie_goes_to_the_non_orphan_call(tmp_path):
    """HG002 chr1:88507142: a transduction-only L1 call and a true SVA call, 2 vs 2.
    Priority alone would keep the L1; the tie-break gives it to SVA."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 88507142, 88507150, False,
                                      {"cls": 2, "subtype": ORPHAN})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 88507142, 88507150, False,
                                        {"cls": 1, "polya": 4})])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "SVA": sva})
    assert [l.type for l in loci] == ["SVA"]
    assert table.conflicts[0].rule == "tie_non_orphan"


def test_orphan_with_more_evidence_still_wins(tmp_path):
    """The orphan tie-break only applies to exact ties; it never outranks score."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False,
                                      {"cls": 3, "subtype": ORPHAN})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 1000, 1010, False, {"cls": 1, "polya": 2})])
    loci = LocusTable().assign_type({"L1": l1, "SVA": sva})
    assert [l.type for l in loci] == ["L1"]


def test_two_orphans_tie_on_priority(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False, {"cls": 1, "subtype": ORPHAN})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 1000, 1010, False, {"cls": 1, "subtype": ORPHAN})])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "SVA": sva})
    assert [l.type for l in loci] == ["L1"]
    assert table.conflicts[0].rule == "tie_priority"


def test_evidence_score_columns():
    f = _row_fields("chr1", 1, 2, False, {"cls": 3, "polya": 7, "tsd": 15})
    assert evidence_components(f) == (3, 1, 1)
    assert evidence_score(f) == 5
    f[9], f[10] = "0", "2"          # right-side poly-A counts too
    assert evidence_components(f)[1] == 1
    f[4] = "0"                       # TSD length 0 is no TSD
    assert evidence_components(f)[2] == 0
    assert not is_orphan(f)
    f[32] = ORPHAN
    assert is_orphan(f)


def test_short_rows_score_zero_and_are_not_orphan():
    """cns_filter output and old fixtures may lack columns: no evidence, no crash."""
    f = "chr1\t0\t1000\t1010\tfiller\tone_end_consistent".split()
    assert evidence_score(f) == 0
    assert not is_orphan(f)


def test_beaten_calls_carry_over_to_the_final_winner(tmp_path):
    """L1 (1) is placed, SVA (2) beats it, Alu (3) beats SVA: Alu must list both."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False, {"cls": 1})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 1005, 1015, False, {"cls": 2})])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False, {"cls": 3})])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "SVA": sva, "Alu": alu})
    assert [l.type for l in loci] == ["Alu"]
    assert sorted(b.type for b in loci[0].beaten) == ["L1", "SVA"]
    assert {table.winner_of()[id(l)].type for l in table.all_loci if l is not loci[0]} == {"Alu"}


def test_exact_position_tie_keeps_priority_winner(tmp_path):
    # Regression for B5: master's real merger silently lets the LAST-loaded type
    # (Alu) win an exact-position tie regardless of priority/TPRT. The fixed rule
    # must keep the higher-priority (or TPRT-both) call instead.
    sva = _hc_pf(tmp_path, "sva.txt", [("chr5", 142000088, 142000098, True)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr5", 142000088, 142000098, True)])
    loci = LocusTable().assign_type({"SVA": sva, "Alu": alu})
    assert len(loci) == 1
    assert loci[0].type == "SVA"


def test_out_of_window_is_not_a_conflict(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1051, 1061, False)])  # just outside 50bp window
    loci = LocusTable(window=50).assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 2


def test_assign_and_write_sorts_and_writes_type_map(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr2", 500, 510, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 100, 110, False)])
    sf_merged = str(tmp_path / "merged.txt")
    sf_type_map = str(tmp_path / "type_map.tsv")
    loci, table = assign_and_write({"L1": l1, "Alu": alu}, sf_merged, sf_type_map)
    assert len(loci) == 2
    assert table.conflicts == []   # two different contigs cannot contend

    with open(sf_merged) as fin:
        lines = [ln.rstrip("\n") for ln in fin]
    assert len(lines) == 2
    assert lines[0].startswith("chr1\t")  # chr1 before chr2 -- sort -k1,1V -k2,2n
    assert lines[1].startswith("chr2\t")

    type_map = load_type_map(sf_type_map)
    assert type_map == {("chr1", 100): TYPE_TO_RTYPE["Alu"], ("chr2", 500): TYPE_TO_RTYPE["L1"]}


# --- the family-contention population (Stage 2c's gate baseline) --------------------------

def test_conflicts_record_the_priority_win(tmp_path):
    """L1 loads before Alu, so L1 keeps the locus -- and the displaced Alu call is
    recorded rather than vanishing. Stage 2c moves this decision up front; "only
    family-assignment flips" is uncheckable unless the flips that already happen have
    been counted."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1020, 1030, False)])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1 and loci[0].type == "L1"
    assert len(table.conflicts) == 1
    c = table.conflicts[0]
    assert (c.chrm, c.winner_type, c.loser_type, c.rule) == ("chr1", "L1", "Alu", "tie_priority")
    assert c.winner_pos == 1000 and c.loser_pos == 1020
    assert table.conflict_summary() == {("L1", "Alu", "tie_priority"): 1}


def test_conflicts_record_the_score_win(tmp_path):
    """A better-evidenced call displaces a higher-priority one: same collision, a
    different rule, and the two must not be conflated -- a Stage 2c flip that changed
    which rule fired is a behaviour change even when the winner is unchanged."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1020, 1030, True, {"cls": 2})])
    table = LocusTable()
    loci = table.assign_type({"L1": l1, "Alu": alu})
    assert len(loci) == 1 and loci[0].type == "Alu"
    assert table.conflict_summary() == {("Alu", "L1", "score"): 1}
    c = table.conflicts[0]
    assert c.winner_tprt is True and c.loser_tprt is False


def test_no_conflict_outside_the_window(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 2000, 2010, False)])
    table = LocusTable()
    assert len(table.assign_type({"L1": l1, "Alu": alu})) == 2
    assert table.conflicts == []


def test_conflicts_are_written_and_do_not_touch_the_two_real_outputs(tmp_path):
    """The sidecar is diagnostic only. If writing it ever changed the merged list or the
    type map, the whole point -- measuring without perturbing -- would be lost."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1020, 1030, False)])
    a_merged, a_map = str(tmp_path / "a.txt"), str(tmp_path / "a.map")
    b_merged, b_map = str(tmp_path / "b.txt"), str(tmp_path / "b.map")
    sf_conf = str(tmp_path / "conf.tsv")
    assign_and_write({"L1": l1, "Alu": alu}, a_merged, a_map)
    assign_and_write({"L1": l1, "Alu": alu}, b_merged, b_map, sf_conflicts=sf_conf)
    assert open(a_merged).read() == open(b_merged).read()
    assert open(a_map).read() == open(b_map).read()

    rows = [ln.rstrip("\n") for ln in open(sf_conf) if not ln.startswith("#")]
    assert len(rows) == 1
    assert rows[0].split("\t") == ["chr1", "L1", "1000", "0", "Alu", "1020", "0",
                                   "tie_priority", "0", "0"]


# --- the record of every call, won or lost --------------------------------------------

def _conflicted_pair(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 88507142, 88507150, False,
                                      {"cls": 2, "subtype": ORPHAN}),
                                     ("chr2", 500, 510, False, {"cls": 4})])
    sva = _hc_pf(tmp_path, "sva.txt", [("chr1", 88507142, 88507150, False,
                                        {"cls": 1, "polya": 4})])
    return {"L1": l1, "SVA": sva}


def test_type_map_lists_the_beaten_calls(tmp_path):
    sf_map = str(tmp_path / "map.tsv")
    assign_and_write(_conflicted_pair(tmp_path), str(tmp_path / "m.txt"), sf_map)
    lines = [ln.split("\t") for ln in open(sf_map).read().splitlines()]
    assert sorted(lines) == [["chr1", "88507142", "4", "L1:88507142"],
                             ["chr2", "500", "1", "."]]
    assert load_type_map(sf_map) == {("chr1", 88507142): 4, ("chr2", 500): 1}
    assert load_other_families(sf_map) == {("chr1", 88507142): "L1:88507142"}


def test_old_three_column_type_map_still_loads(tmp_path):
    sf_map = tmp_path / "old.tsv"
    sf_map.write_text("chr1\t100\t1\n")
    assert load_type_map(str(sf_map)) == {("chr1", 100): 1}
    assert load_other_families(str(sf_map)) == {}


def test_all_calls_tsv_and_bed(tmp_path):
    sf_tsv, sf_bed = str(tmp_path / "all.tsv"), str(tmp_path / "all.bed")
    paths = _conflicted_pair(tmp_path)
    assign_and_write(paths, str(tmp_path / "m.txt"), str(tmp_path / "map.tsv"),
                     sf_all_calls=sf_tsv, sf_all_calls_bed=sf_bed)
    lines = open(sf_tsv).read().splitlines()
    assert lines[0] == ALL_CALLS_HEADER
    rows = [ln.split("\t") for ln in lines[1:]]
    # chr1 (L1 then SVA by family name order), then chr2
    assert [r[:11] for r in rows] == [
        ["chr1", "88507142", "L1", "2", "2", "0", "0", "1", "LOST", "SVA", "88507142"],
        ["chr1", "88507142", "SVA", "2", "1", "1", "0", "0", "WON", "SVA", "88507142"],
        ["chr2", "500", "L1", "4", "4", "0", "0", "0", "WON", "L1", "500"],
    ]
    # the full evidence line rides along, byte-identical to the input row
    l1_line = open(paths["L1"]).read().splitlines()[0]
    assert "\t".join(rows[0][11:]) == l1_line

    bed = open(sf_bed).read().splitlines()
    assert bed[0].startswith("track ")
    assert bed[1:] == ["chr1\t88507141\t88507142\tL1 score=2 LOST",
                       "chr1\t88507141\t88507142\tSVA score=2 WON",
                       "chr2\t499\t500\tL1 score=4 WON"]


def test_all_calls_files_do_not_change_the_merge(tmp_path):
    a = (str(tmp_path / "a.txt"), str(tmp_path / "a.map"))
    b = (str(tmp_path / "b.txt"), str(tmp_path / "b.map"))
    paths = _conflicted_pair(tmp_path)
    assign_and_write(paths, *a)
    assign_and_write(paths, *b, sf_all_calls=str(tmp_path / "x.tsv"),
                     sf_all_calls_bed=str(tmp_path / "x.bed"))
    assert open(a[0]).read() == open(b[0]).read()
    assert open(a[1]).read() == open(b[1]).read()


# --- Stage 2c: the same rule, applied mid-chain (assign_and_filter) -------------------
#
# 2c moves family assignment from the postF_hc merge up to cns_filter's output. The
# thing that must be true for that to be a *relocation* rather than a rewrite is that
# exactly one copy of the rule exists and both callers reach it. These tests pin that,
# plus the file surgery assign_and_filter does on top of it.


def test_assign_rows_and_assign_type_are_the_same_rule(tmp_path):
    """assign_type is defined as load_rows + assign_rows. If they ever diverge, the
    'moved, not rewritten' claim in the plan stops holding."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False), ("chr2", 77, 88, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, True)])
    via_paths = LocusTable().assign_type({"L1": l1, "Alu": alu})
    via_rows = LocusTable().assign_rows({"L1": load_rows(l1, "L1"), "Alu": load_rows(alu, "Alu")})
    assert ([(l.chrm, l.pos, l.type) for l in via_paths]
            == [(l.chrm, l.pos, l.type) for l in via_rows])


def test_assign_and_filter_keeps_only_the_winning_family(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False), ("chr3", 500, 510, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False), ("chr9", 40, 50, False)])
    stats, table = assign_and_filter({"L1": l1, "Alu": alu})
    # L1 wins chr1 by priority; Alu keeps only its uncontended chr9 row.
    assert stats == {"L1": (2, 2), "Alu": (2, 1)}
    assert len(table.conflicts) == 1
    kept_alu = [ln.split()[0] for ln in open(alu)]
    assert kept_alu == ["chr9"]
    assert [ln.split()[0] for ln in open(l1)] == ["chr1", "chr3"]


def test_assign_and_filter_backs_up_the_original(tmp_path):
    """The pre-assignment file is what a gate diffs against to attribute a moved call,
    so losing it would make the movement unexplainable after the fact."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False)])
    before = open(alu).read()
    assign_and_filter({"L1": l1, "Alu": alu})
    assert open(alu + ".pre_family_assign").read() == before
    assert open(alu).read() == ""


def test_assign_and_filter_preserves_row_order_and_content(tmp_path):
    """Surviving rows must come through byte-identical and in the original order --
    the next step parses them positionally."""
    rows = [("chr1", 100 * i, 100 * i + 10, False) for i in range(1, 6)]
    l1 = _hc_pf(tmp_path, "l1.txt", rows)
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 300, 310, True, {"cls": 1})])  # takes row 3 from L1
    original = open(l1).read().splitlines()
    assign_and_filter({"L1": l1, "Alu": alu})
    kept = open(l1).read().splitlines()
    assert kept == [ln for ln in original if ln.split()[2] != "300"]


def test_assign_and_filter_passes_through_unparsed_lines(tmp_path):
    """A line the parser does not recognise as a claim (fewer than 4 fields) is not a
    candidate and must not be silently deleted by a step that only decides family."""
    l1 = tmp_path / "l1.txt"
    l1.write_text("# a comment\nchr1\t0\t1000\t1010\tfiller\tone_end_consistent\n")
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, True, {"cls": 1})])
    stats, _ = assign_and_filter({"L1": str(l1), "Alu": alu})
    assert stats["L1"] == (1, 0)                    # the one real claim lost to Alu
    assert l1.read_text() == "# a comment\n"        # the comment survived


def test_assign_and_filter_is_a_no_op_without_contention(tmp_path):
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr5", 1000, 1010, False)])
    before = open(l1).read(), open(alu).read()
    stats, table = assign_and_filter({"L1": l1, "Alu": alu})
    assert table.conflicts == []
    assert stats == {"L1": (1, 1), "Alu": (1, 1)}
    assert (open(l1).read(), open(alu).read()) == before


# --- the companion files, and the bug that only a real run found -----------------------


def test_assign_and_filter_filters_the_high_confident_companion(tmp_path):
    """cns_filter writes `<candidates>.high_confident` next to its candidate list, and
    the transduction step opens it by DERIVING the path from its own -i rather than
    being handed it (x_TEA_main.py:551). update_high_confident_callset then looks every
    one of its rows up in the transduction output with no membership test, so a row left
    in the companion after being dropped from the candidate list kills the step with a
    bare KeyError. That is exactly how job 49498881 died on SVA:transduction."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    # equal (zero) evidence, so priority gives chr1 to L1 and Alu loses it
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, False),
                                       ("chr9", 40, 50, False)])
    # the companion: a subset of alu.txt's rows, byte-identical lines (as on real data)
    alu_hc = alu + ".high_confident"
    with open(alu_hc, "w") as fout:
        fout.writelines(open(alu).readlines())

    assign_and_filter({"L1": l1, "Alu": alu})

    kept = [ln.split()[0] for ln in open(alu_hc)]
    assert kept == ["chr9"], "the companion still names a locus Alu no longer owns"
    # and the original is preserved for attribution, like the candidate list's
    assert len(open(alu_hc + ".pre_family_assign").read().splitlines()) == 2


def test_companion_is_filtered_on_the_key_downstream_actually_joins_on(tmp_path):
    """update_high_confident_callset keys on int(fields[1]), NOT on the lclip/rclip
    call_position used to cluster loci. Filtering the companion on the wrong one leaves
    exactly the mismatch this is meant to prevent."""
    # fields[1] = 0 in the _hc_pf fixture, so give the companion a row whose fields[1]
    # differs from its call_position and check the join still lands.
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = tmp_path / "alu.txt"
    alu.write_text("chr1\t777\t1010\t1020\tfiller\tone_end_consistent\n"
                   "chr9\t888\t40\t50\tfiller\tone_end_consistent\n")
    (tmp_path / "alu.txt.high_confident").write_text(alu.read_text())

    assign_and_filter({"L1": l1, "Alu": str(alu)})

    kept = [ln.split()[1] for ln in open(str(alu) + ".high_confident")]
    assert kept == ["888"], "companion filtered on the wrong key"


def test_missing_companion_is_not_an_error(tmp_path):
    """--types with a single family, or a step that never wrote one, must not trip this."""
    l1 = _hc_pf(tmp_path, "l1.txt", [("chr1", 1000, 1010, False)])
    alu = _hc_pf(tmp_path, "alu.txt", [("chr1", 1010, 1020, True)])
    assign_and_filter({"L1": l1, "Alu": alu})     # no .high_confident anywhere
    assert not os.path.exists(alu + ".high_confident")
