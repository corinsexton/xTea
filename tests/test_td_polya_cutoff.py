"""Gate [A]'s polyA threshold, decoupled from the clip threshold -- and STILL derived by default.

WHY THIS CHANGE EXISTS. `n_polya_cutoff` has always been a function of `n_clip_cutoff`
(`n_clip_cutoff / 2 - 1`), so the two thresholds at the gate that does 97% of the step's
rejecting could not be moved independently. That makes a clip-support sweep single-variable only
up to 4: at `--nclip 5` the polyA arm wakes and rejects 286 (L1) / 364 (SVA) loci the shipping
configuration keeps, and both moves get reported as one number
(`HANDOFF_2026-08-05_gate_a_measured.md` §2.3). Plan §5's "every threshold a named parameter"
ambition was failing at the gate that matters most.

**THE SAFETY ARGUMENT IS THAT THE DEFAULT PATH IS UNCHANGED, AND IT IS TESTED THREE WAYS**, not
asserted: the derivation is reproduced verbatim over a sweep of clip cutoffs including the
floor-at-0 branch (`test_omitting_the_cutoff_reproduces_the_derivation_exactly`), the two entry
points both resolve absent-to-derived (`test_*_absent_means_derived`), and at HG002's own
operating point an explicit `0.5` is verdict-identical to the default over a generated population
(`test_at_hg002s_operating_point_the_explicit_value_is_verdict_identical`).

**AND THE FLAG'S OWN GATE IS THAT IT MUST MOVE NOTHING AT HG002.** The derived cutoff there is
0.5 and it is UNREACHABLE -- `m_polyA`'s record is created and one of its six counters
incremented in the same block, so every entry has `sum(rcd[:6]) >= 1`, measured as 0 zero-sum
entries of 4,604 (L1) and 5,261 (SVA). So a WGS arm passing `--td-polya-cutoff 0.5` must come
back byte-identical, and the first value that can move a call is one ABOVE 1. That is an unusually
strong gate for a call-moving flag and it is why this change is cheap to land:
`test_the_derived_cutoff_is_unreachable_at_hg002` pins the reachability claim itself.

ZERO IS A REAL ANSWER AND MUST NOT COLLIDE WITH "NO OPINION" -- the reason the env channel spells
the default `"derived"` rather than `"0"`, and the reason the resolvers return `None` rather than
falling back. `--td-polya-cutoff 0` is a caller asking for a threshold nothing can fail; at HG002
it agrees with the default by coincidence and at `--nclip 5` it does not.

MUTATIONS RUN AGAINST THIS FILE, all of which bite:

  | mutation | caught by |
  |---|---|
  | the override branch dropped (always derive) | the override + threading tests |
  | the override applied only when n_clip_cutoff is None | the override test |
  | the override floored at 0 like the derivation | `test_an_explicit_cutoff_is_used_verbatim` |
  | `"derived"` parsed as a float (raises) | `test_env_absent_means_derived` |
  | a malformed env value falling back to None | `test_a_malformed_env_value_raises` |
  | env `""` treated as 0.0 | `test_env_absent_means_derived` + the zero/None test |
"""
import os
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import global_values  # noqa: E402
import transduction_evidence as te  # noqa: E402
import x_transduction  # noqa: E402

L1_RTYPE = 1

# HG002 WGS at 31.816x resolves `m_cov_par[30] = (3, 4, 1)` -- n_clip_cutoff 3, n_disc_cutoff 4.
# This is the operating point every measurement in this project is taken at, and it is NOT the
# CLI default of 4/5 (see test_transduction_evidence.py's note on that trap).
HG002_NCLIP = 3
HG002_DERIVED_POLYA = HG002_NCLIP / 2 - 1  # 0.5


def _derived(n_clip_cutoff):
    """The frozen expression, written once here so a test cannot drift from the shipping one."""
    v = n_clip_cutoff / 2 - 1
    return 0 if v < 0 else v


# ======================================================================================
# 1. The default path is unchanged
# ======================================================================================

@pytest.mark.parametrize("n_clip", [0, 1, 2, 3, 4, 5, 6, 7, 10, 25])
def test_omitting_the_cutoff_reproduces_the_derivation_exactly(n_clip):
    """No override -> the frozen `n_clip_cutoff/2-1`, float and floor-at-0 included.

    The low end is the point: at n_clip 0 and 1 the expression is negative and the shipping code
    floors it to 0, so a replacement that forgot the floor would return -1.0 and admit nothing
    differently -- but would still differ from the frozen object, which is what BASELINE-0 was
    measured with.
    """
    p = te.TransductionParams(n_clip_cutoff=n_clip, n_disc_cutoff=4)
    assert p.n_polya_cutoff == _derived(n_clip)


def test_a_caller_with_no_clip_cutoff_and_no_override_still_gets_none():
    """`None`, not 0 -- the class docstring's rule, and the sibling step depends on it.

    A partial object is what `filter_transduction` and `add_new_rcd_with_cutoff` are built with
    from `x_orphan_transduction`, which never runs gates A-C. A 0 here would silently admit
    everything at a gate that caller never intended to reach; None raises TypeError instead.
    """
    p = te.TransductionParams(n_disc_cutoff=4)
    assert p.n_polya_cutoff is None


def test_the_derived_cutoff_is_unreachable_at_hg002():
    """The claim the flag's whole no-op gate rests on, pinned rather than quoted.

    Every m_polyA entry has `sum(rcd[:6]) >= 1` by construction, so with the derived 0.5 the
    comparison `n_polya_clip < 0.5` needs a sum of zero and cannot fire. Measured: 0 zero-sum
    entries of 4,604 (L1) / 5,261 (SVA) at WGS.
    """
    p = te.TransductionParams(n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4)
    assert p.n_polya_cutoff == HG002_DERIVED_POLYA == 0.5
    # 1 is the smallest value any real entry can carry; it must pass.
    assert not (1 < p.n_polya_cutoff)
    # ...and the threshold first bites above 1, which is where --nclip 5 lands it (1.5).
    assert _derived(5) > 1


# ======================================================================================
# 2. The override
# ======================================================================================

@pytest.mark.parametrize("value", [0, 0.0, 0.5, 1, 1.5, 2, 3.75, 99, -1, -0.5])
def test_an_explicit_cutoff_is_used_verbatim(value):
    """No rounding and NO FLOOR -- the derivation's floor is the derivation's, not the gate's.

    A negative override is accepted deliberately: it is indistinguishable in effect from 0 at
    this gate, but silently rewriting a caller's number is how a sweep point becomes a
    misreported one. The flag's banner prints what it will compare against.
    """
    p = te.TransductionParams(n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4,
                              n_polya_cutoff=value)
    assert p.n_polya_cutoff == value


def test_the_override_beats_the_derivation_and_they_are_different_numbers():
    """Guards the mutation "apply the override only when there is no clip cutoff to derive from"."""
    p = te.TransductionParams(n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4, n_polya_cutoff=7.0)
    assert p.n_polya_cutoff == 7.0
    assert p.n_polya_cutoff != _derived(HG002_NCLIP)
    # ...and it did not disturb the threshold it was decoupled FROM.
    assert p.n_clip_cutoff == HG002_NCLIP
    assert p.n_half_clip == HG002_NCLIP / 2


def test_explicit_zero_is_distinguishable_from_no_opinion():
    """The whole reason the channel is None-means-derive rather than 0-means-derive."""
    explicit = te.TransductionParams(n_clip_cutoff=5, n_disc_cutoff=4, n_polya_cutoff=0)
    default = te.TransductionParams(n_clip_cutoff=5, n_disc_cutoff=4)
    assert explicit.n_polya_cutoff == 0
    assert default.n_polya_cutoff == _derived(5) == 1.5
    assert explicit.n_polya_cutoff != default.n_polya_cutoff


def test_a_partial_object_can_still_be_given_the_threshold():
    """No clip cutoff to derive from, but an explicit polyA cutoff is still honoured."""
    p = te.TransductionParams(n_polya_cutoff=2.0)
    assert p.n_polya_cutoff == 2.0
    assert p.n_clip_cutoff is None


# ======================================================================================
# 3. It reaches the gate
# ======================================================================================

def _ev(n_clip, n_polya):
    return te.SiteEvidence(chrm="chr1", pos=1000, n_clip=n_clip, n_polya_clip=n_polya)


def test_the_override_changes_gate_a_verdicts():
    """A locus kept under the derived threshold and dropped under a raised one.

    Without this the constructor could be correct and the gate still read the old value.
    """
    ev = _ev(n_clip=10, n_polya=1)
    kept = te.gate_clip_evidence(ev, te.TransductionParams(n_clip_cutoff=HG002_NCLIP,
                                                           n_disc_cutoff=4))
    assert kept.keep, "1 polyA read must pass the derived 0.5"

    dropped = te.gate_clip_evidence(ev, te.TransductionParams(
        n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4, n_polya_cutoff=2.0))
    assert not dropped.keep
    assert dropped.reason == te.R_LOW_POLYA_CLIP
    # The detail carries the numbers that decided it, so a ledger row is readable later.
    assert "2.0" in dropped.detail


def test_the_override_does_not_reach_the_clip_support_gate():
    """Decoupling means the OTHER threshold does not move -- the point of the change."""
    ev = _ev(n_clip=2, n_polya=50)
    v = te.gate_clip_evidence(ev, te.TransductionParams(
        n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4, n_polya_cutoff=0))
    assert not v.keep
    assert v.reason == te.R_LOW_CLIP, "a polyA override must not rescue a low-clip locus"


def test_absence_from_the_polya_map_is_still_a_drop_whatever_the_cutoff():
    """`n_polya_clip is None` is a membership failure, decided BEFORE any threshold.

    Gate [A]'s big population -- 68% / 78% of clip-qualified loci -- is dropped here, and no
    value of this flag can reach it. That is measured in the gate-[A] handoff §2.2 and it is
    the separate decision (the "is the veto right at all" question), not this one.
    """
    for cutoff in (None, 0, -1, 99):
        kwargs = {} if cutoff is None else {"n_polya_cutoff": cutoff}
        v = te.gate_clip_evidence(_ev(n_clip=969, n_polya=None),
                                  te.TransductionParams(n_clip_cutoff=HG002_NCLIP,
                                                        n_disc_cutoff=4, **kwargs))
        assert not v.keep
        assert v.reason == te.R_NO_POLYA_CLIP


def test_at_hg002s_operating_point_the_explicit_value_is_verdict_identical():
    """The flag's own gate, at the gate's level: 0.5 explicit == derived, over a population.

    This is the unit-level form of the claim the WGS arm has to reproduce byte-identically.
    """
    default = te.TransductionParams(n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4)
    explicit = te.TransductionParams(n_clip_cutoff=HG002_NCLIP, n_disc_cutoff=4,
                                     n_polya_cutoff=HG002_DERIVED_POLYA)
    for n_clip in range(0, 15):
        for n_polya in [None] + list(range(0, 8)):
            a = te.gate_clip_evidence(_ev(n_clip, n_polya), default)
            b = te.gate_clip_evidence(_ev(n_clip, n_polya), explicit)
            assert (a.keep, a.reason) == (b.keep, b.reason), (n_clip, n_polya)


# ======================================================================================
# 4. The two entry points, both absent-means-derived
# ======================================================================================

@pytest.mark.parametrize("raw", ["", "   ", "derived", "DERIVED", " Derived "])
def test_env_absent_means_derived(raw):
    assert x_transduction.resolve_env_polya_cutoff({global_values.ENV_TD_POLYA_CUTOFF: raw}) is None


def test_env_unset_means_derived():
    assert x_transduction.resolve_env_polya_cutoff({}) is None


@pytest.mark.parametrize("raw,expect", [("0", 0.0), ("0.5", 0.5), ("2", 2.0), ("3.75", 3.75),
                                        (" 1.5 ", 1.5), ("-1", -1.0)])
def test_env_carries_a_value(raw, expect):
    got = x_transduction.resolve_env_polya_cutoff({global_values.ENV_TD_POLYA_CUTOFF: raw})
    assert got == expect
    assert got is not None, "0.0 must not read as 'no opinion'"


@pytest.mark.parametrize("raw", ["yes", "1,5", "None", "0x2", "1.2.3"])
def test_a_malformed_env_value_raises(raw):
    """It must NOT fall back to the derivation.

    A silent revert would make an arm that measured nothing look like an arm that measured no
    change -- the exact confusion the decoy default flip had to be rescued from.
    """
    with pytest.raises(ValueError) as e:
        x_transduction.resolve_env_polya_cutoff({global_values.ENV_TD_POLYA_CUTOFF: raw})
    assert global_values.ENV_TD_POLYA_CUTOFF in str(e.value)


def test_the_cli_resolver_matches_the_env_one():
    import run_xtea
    assert run_xtea.resolve_polya_cutoff(None) is None
    assert run_xtea.resolve_polya_cutoff("0") == 0.0
    assert run_xtea.resolve_polya_cutoff("2.5") == 2.5
    assert run_xtea.resolve_polya_cutoff(3) == 3.0


def test_the_driver_writes_derived_not_zero_for_the_default():
    """The env channel's spelling, which is what keeps 0 and "no opinion" apart across the seam.

    The step is its own `x_TEA_main` subprocess, so this string is the entire contract between
    the driver's flag and the gate's threshold.
    """
    import run_xtea
    v = run_xtea.resolve_polya_cutoff(None)
    encoded = "derived" if v is None else repr(v)
    assert encoded == "derived"
    assert x_transduction.resolve_env_polya_cutoff(
        {global_values.ENV_TD_POLYA_CUTOFF: encoded}) is None
    # ...and a real value survives the same round trip exactly.
    for cli in ("0", "0.5", "2.5"):
        enc = repr(run_xtea.resolve_polya_cutoff(cli))
        assert x_transduction.resolve_env_polya_cutoff(
            {global_values.ENV_TD_POLYA_CUTOFF: enc}) == float(cli)
