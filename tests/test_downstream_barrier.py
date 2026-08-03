"""The cns_filter -> transduction barrier, staged: who waits for whom.

`run_xtea.StagedBarrier` replaces a global join between the two downstream waves with
per-family waits. The join cost +52.9s at WGS (job 49550925) and none of it was work --
every family was waiting for Alu's cns_filter, which is ~2x L1's and SVA's, and the
unsplit chain overlaps L1/SVA transduction with the tail of it.

What can go wrong here is not the work each piece does (the enumeration is pinned by
tests/test_downstream_unified.py, the assignment by tests/test_locus_table.py, the gather
by tests/test_td_gather.py). It is the coordination, and it has four failure modes that
are all silent or fatal in a way a call-set diff would never localize:

  1. **A family waits for something it does not read.** Costs exactly the wall clock this
     change exists to recover, and no output moves -- so only a timing arm catches it.
     The flankless family and the driver gather are the case that matters: Alu never
     calls `call_candidate_transduction_v3`, so it must not block on the gather.
  2. **A family runs before something it DOES read.** Fatal but not obviously: a
     transduction step that starts before the gather it was told to skip finds an absent
     FASTQ, and job 49498881 is what that looks like -- a death inside a step, at WGS,
     invisible on chr1.
  3. **The shared work runs twice, or not at all.** Two threads electing themselves both
     rewrite the same files; nobody electing hangs the run.
  4. **A family dies before the barrier and the rest wait for it forever.** The global
     join had this for free (nothing waits on a dead wave); here it has to be built, and
     the failure is a hang inside the executor's shutdown rather than an error.

Every test below drives real threads with a deadline, because a coordination bug is a
hang, and a hang with no deadline is a test suite that never finishes.
"""
import os
import sys
import threading

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))
import run_xtea  # noqa: E402

# Long enough that a loaded machine does not flake it, short enough that a genuine
# deadlock fails the suite in seconds rather than hanging CI.
DEADLINE = 20.0

FAMILIES = ["L1", "Alu", "SVA"]
FLANKED = ["L1", "SVA"]


class _Recorder:
    """The four callbacks, recording what ran, in what order, on which thread."""

    def __init__(self):
        self.lock = threading.Lock()
        self.events = []            # (label, family-or-None)
        self.enumerated = []
        self.gathers = 0
        self.unions = []            # the slct map each union call saw
        self.assigns = 0

    def _note(self, label, who=None):
        with self.lock:
            self.events.append((label, who))

    def enumerate_fn(self, t):
        self._note("enum", t)
        with self.lock:
            self.enumerated.append(t)
        return (1, 2, f"/tmp/{t}.slct")

    def gather_fn(self):
        self._note("gather")
        with self.lock:
            self.gathers += 1

    def union_fn(self, slct):
        self._note("union")
        with self.lock:
            self.unions.append(dict(slct))

    def assign_fn(self):
        self._note("assign")
        with self.lock:
            self.assigns += 1

    def labels(self):
        with self.lock:
            return [lab for (lab, _who) in self.events]


def _barrier(rec, b_universe=True, b_gather=False, b_assign=False,
             families=FAMILIES, flanked=FLANKED):
    return run_xtea.StagedBarrier(
        families=families, flanked=flanked,
        b_universe=b_universe, b_gather=b_gather, b_assign=b_assign,
        enumerate_fn=rec.enumerate_fn, gather_fn=rec.gather_fn,
        union_fn=rec.union_fn, assign_fn=rec.assign_fn)


def _arrive_all(barrier, families=FAMILIES, before=None, deadline=DEADLINE):
    """One thread per family calling arrive(). Returns {family: exception or None}.

    `before(t)` runs on the family's thread before it arrives -- it stands in for that
    family's cns_filter, and is where a test puts the delay that orders arrivals.
    Anything still running at the deadline is reported as a TimeoutError for that
    family, which is what a missed release looks like.
    """
    out = {}
    threads = {}

    def _one(t):
        try:
            if before is not None:
                before(t)
            barrier.arrive(t)
            out[t] = None
        except BaseException as e:     # noqa: BLE001 -- the point is to report it
            out[t] = e

    for t in families:
        threads[t] = threading.Thread(target=_one, args=(t,), daemon=True)
    for t in families:
        threads[t].start()
    for t in families:
        threads[t].join(deadline)
        if threads[t].is_alive():
            out[t] = TimeoutError(f"{t} still blocked in arrive() after {deadline}s")
    return out


# --- 1. what nobody has to wait for ---------------------------------------------------

def test_universe_alone_makes_nobody_wait():
    """`--locus-universe` on its own must not be a barrier at all.

    Each family's `.slct` is derived from its own cns list and its own raw_disc and
    reads no other family's anything, so the enumeration is per-family work that happens
    to run at the same point -- not a rendezvous. This is what makes flipping the flag on
    cost its ~3s of enumeration instead of the slowest family's cns_filter.

    Proved by starving the pool: one thread at a time, with the others not yet started.
    A real wait would never return.
    """
    rec = _Recorder()
    barrier = _barrier(rec)
    for t in FAMILIES:
        done = threading.Event()

        def _one(tag=t):
            barrier.arrive(tag)
            done.set()

        threading.Thread(target=_one, daemon=True).start()
        assert done.wait(DEADLINE), f"{t} blocked with --locus-universe alone"
    assert rec.enumerated == FAMILIES
    assert rec.gathers == 0 and rec.assigns == 0


def test_flankless_family_never_waits_for_the_driver_gather():
    """Alu must not block on a gather it never reads.

    x_TEA_main takes the copy-through branch when `os.path.isfile(sf_flank)` is False, so
    a flankless family never calls `call_candidate_transduction_v3`, never looks at the
    gathered FASTQ, and never reads XTEA_TD_GATHER_PRECOMPUTED. If it waited here it
    would pay the gather's ~42s for nothing.

    The gather is held open until Alu has come out the other side, so a wait deadlocks.
    """
    rec = _Recorder()
    alu_through = threading.Event()
    release_gather = threading.Event()

    def _slow_gather():
        assert alu_through.wait(DEADLINE), "the flankless family blocked on the gather"
        rec.gather_fn()

    barrier = run_xtea.StagedBarrier(
        families=FAMILIES, flanked=FLANKED, b_universe=True, b_gather=True, b_assign=False,
        enumerate_fn=rec.enumerate_fn, gather_fn=_slow_gather,
        union_fn=rec.union_fn, assign_fn=rec.assign_fn)

    results = {}

    def _one(t):
        try:
            barrier.arrive(t)
            results[t] = None
        except BaseException as e:      # noqa: BLE001
            results[t] = e
        if t == "Alu":
            alu_through.set()

    threads = [threading.Thread(target=_one, args=(t,), daemon=True) for t in FAMILIES]
    for th in threads:
        th.start()
    for th in threads:
        th.join(DEADLINE)
    release_gather.set()
    assert all(not th.is_alive() for th in threads), "a family never left the barrier"
    assert results == {"L1": None, "Alu": None, "SVA": None}
    assert rec.gathers == 1


# --- 2. what everybody has to wait for ------------------------------------------------

def test_the_gather_runs_once_and_after_every_flanked_family_enumerated():
    """The gather reads BOTH flanked families' `.slct`, so both must be written first.

    And exactly once: the manifest names each family's output paths, so a second run
    rewrites the FASTQs the first family may already be reading.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_gather=True)
    errs = _arrive_all(barrier)
    assert all(e is None for e in errs.values()), errs

    labels = rec.labels()
    assert labels.count("gather") == 1
    i_gather = labels.index("gather")
    enums_before = [w for (lab, w) in rec.events[:i_gather] if lab == "enum"]
    for t in FLANKED:
        assert t in enums_before, f"the gather ran before {t} enumerated"


def test_every_flanked_family_transduces_after_the_gather():
    """The other side of the same edge -- a flanked family must not return early.

    This is failure mode 2, and it is the one that killed job 49498881: a transduction
    step told (by XTEA_TD_GATHER_PRECOMPUTED) to skip its own gather, started before the
    driver's gather wrote the file it was going to skip to.
    """
    rec = _Recorder()
    gather_running = threading.Event()
    gather_finished = threading.Event()
    escaped = []

    def _slow_gather():
        gather_running.set()
        # Hold the gather open; any flanked family that returns now returned too early.
        threading.Event().wait(0.5)
        rec.gather_fn()
        gather_finished.set()

    barrier = run_xtea.StagedBarrier(
        families=FAMILIES, flanked=FLANKED, b_universe=True, b_gather=True, b_assign=False,
        enumerate_fn=rec.enumerate_fn, gather_fn=_slow_gather,
        union_fn=rec.union_fn, assign_fn=rec.assign_fn)

    def _one(t):
        barrier.arrive(t)
        if t in FLANKED and not gather_finished.is_set():
            escaped.append(t)

    threads = [threading.Thread(target=_one, args=(t,), daemon=True) for t in FAMILIES]
    for th in threads:
        th.start()
    for th in threads:
        th.join(DEADLINE)
    assert all(not th.is_alive() for th in threads)
    assert escaped == [], f"flanked families left the barrier before the gather: {escaped}"


def test_family_assignment_still_joins_every_family():
    """`--family-assign-early` genuinely needs all three, so it degenerates to a join.

    That is correct, not a regression: the rule compares every family's cns rows against
    every other's. The staging saves nothing here and must not pretend to.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=False, b_assign=True)
    released = []
    all_arrived = threading.Event()
    n_arrived = [0]
    lock = threading.Lock()

    def _before(t):
        with lock:
            n_arrived[0] += 1
            if n_arrived[0] == len(FAMILIES):
                all_arrived.set()

    def _one(t):
        _before(t)
        barrier.arrive(t)
        released.append(t)

    threads = [threading.Thread(target=_one, args=(t,), daemon=True) for t in FAMILIES]
    for th in threads:
        th.start()
    for th in threads:
        th.join(DEADLINE)
    assert all(not th.is_alive() for th in threads)
    assert rec.assigns == 1
    assert sorted(released) == sorted(FAMILIES)
    # Nothing was released before the last family arrived: assign is the last event, and
    # every family's wait is on it.
    assert rec.labels().count("assign") == 1


def test_every_enumeration_precedes_the_assignment():
    """The ordering the two-wave form gave for free, and the one that can silently lie.

    The assignment REWRITES each family's cns list to its own winners. Enumerating after
    it would produce a post-assignment universe -- a smaller number that still looks like
    a universe, with nothing in the output to say so.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=True, b_gather=True, b_assign=True)
    errs = _arrive_all(barrier)
    assert all(e is None for e in errs.values()), errs

    labels = rec.labels()
    i_assign = labels.index("assign")
    assert labels[:i_assign].count("enum") == len(FAMILIES), (
        f"an enumeration ran at or after the assignment: {labels}")
    # And the gather, which reads `.slct` and not the cns lists, is also on the safe side
    # of it -- not required by the data, but it is what "the gather goes first" means.
    assert labels.index("gather") < i_assign


# --- 3. the shared work runs exactly once ---------------------------------------------

def test_the_union_table_is_written_once_with_every_family():
    """Written by whichever family enumerates last, and it must see all three.

    A union built from two families is a valid-looking table with a wrong locus count,
    and the count is the number the barrier's ceiling (38,815 of 41,604) is quoted from.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=True)
    errs = _arrive_all(barrier)
    assert all(e is None for e in errs.values()), errs
    assert len(rec.unions) == 1
    assert sorted(rec.unions[0]) == sorted(FAMILIES)
    assert rec.labels().index("union") == len(rec.labels()) - 1


def test_the_union_sees_families_in_driver_order_not_completion_order():
    """The bytes of the universe table must not depend on who finished cns_filter first.

    `LocusUniverse.rows()` writes the per-locus `families`/`claims` columns in the
    iteration order of the map it was built from. Under the global join that map was built
    by walking `type_ctx`, so it was always the driver's family order. The staged barrier
    populates it as each family lands, and completion order is a race.

    Job 49596361 is this bug: chr1, three arms, everything else identical -- the same
    3,272 loci with the same claims -- and the staged arm wrote `SVA,L1` where the global
    arm wrote `L1,SVA`. Every call was identical, so only a comparison of the universe
    artifacts themselves could see it.

    Here SVA is forced to enumerate first and the map must still come out L1, Alu, SVA.
    """
    rec = _Recorder()
    order_gate = {t: threading.Event() for t in FAMILIES}

    def _staggered_enumerate(t):
        # SVA goes first, then Alu, then L1 -- the reverse of the driver's order.
        for prior, tag in ((None, "SVA"), ("SVA", "Alu"), ("Alu", "L1")):
            if t == tag:
                if prior is not None:
                    assert order_gate[prior].wait(DEADLINE), f"{prior} never enumerated"
                break
        out = rec.enumerate_fn(t)
        order_gate[t].set()
        return out

    barrier = run_xtea.StagedBarrier(
        families=FAMILIES, flanked=FLANKED, b_universe=True, b_gather=False, b_assign=False,
        enumerate_fn=_staggered_enumerate, gather_fn=rec.gather_fn,
        union_fn=rec.union_fn, assign_fn=rec.assign_fn)

    errs = _arrive_all(barrier)
    assert all(e is None for e in errs.values()), errs
    assert rec.enumerated == ["SVA", "Alu", "L1"], "the stagger did not take"
    assert len(rec.unions) == 1
    assert list(rec.unions[0]) == FAMILIES, (
        f"union built in completion order: {list(rec.unions[0])}")


def test_a_lone_flanked_family_elects_itself_for_everything():
    """`--types L1` on its own: one thread draws the gather, the union and the assignment.

    It must not wait on events it is itself responsible for setting -- the degenerate
    case where an election scheme deadlocks on itself.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=True, b_gather=True, b_assign=True,
                       families=["L1"], flanked=["L1"])
    errs = _arrive_all(barrier, families=["L1"])
    assert errs == {"L1": None}
    assert rec.labels() == ["enum", "gather", "union", "assign"]


def test_no_flanked_family_means_no_gather_and_no_wait():
    """`--types Alu`: nothing to gather, so the flag must not leave anyone waiting.

    `--td-gather-driver` with a flankless-only type list is a legitimate command line,
    and the gather election's condition (`flanked <= enumerated`) is vacuously true on an
    empty set -- which would run an empty gather and, worse, look like it worked.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=True, b_gather=True,
                       families=["Alu"], flanked=[])
    errs = _arrive_all(barrier, families=["Alu"])
    assert errs == {"Alu": None}
    assert rec.gathers == 0


# --- 4. a family that dies must not hang the rest -------------------------------------

def test_a_death_before_the_barrier_releases_the_waiters():
    """The global join's one advantage, rebuilt: nothing waits on a dead family.

    Without this, a cns_filter that exits non-zero leaves the other families blocked in
    `arrive()` forever, and the executor's shutdown waits for them -- so the run hangs
    instead of failing, at WGS scale, after ~40 minutes of work.
    """
    rec = _Recorder()
    barrier = _barrier(rec, b_universe=True, b_gather=True, b_assign=True)

    def _before(t):
        if t == "Alu":
            # Alu's cns_filter dies: the worker aborts the barrier and re-raises.
            barrier.abort(t)
            raise SystemExit("Alu:cns_filter failed")

    errs = _arrive_all(barrier, before=_before)
    assert isinstance(errs["Alu"], SystemExit)
    for t in FLANKED:
        assert not isinstance(errs[t], TimeoutError), f"{t} hung on a dead family"
        assert isinstance(errs[t], RuntimeError), errs[t]


def test_a_failing_enumeration_aborts_from_inside_arrive():
    """The same guarantee when the failure is the barrier's own work, not a step's.

    `arrive()` runs the enumeration itself, so it has to abort on its own exception --
    the worker's except clause never sees a thread that is stuck in a callback.
    """
    rec = _Recorder()

    def _boom(t):
        if t == "SVA":
            raise ValueError("annotation1 missing")
        return rec.enumerate_fn(t)

    barrier = run_xtea.StagedBarrier(
        families=FAMILIES, flanked=FLANKED, b_universe=True, b_gather=True, b_assign=True,
        enumerate_fn=_boom, gather_fn=rec.gather_fn,
        union_fn=rec.union_fn, assign_fn=rec.assign_fn)

    errs = _arrive_all(barrier)
    assert isinstance(errs["SVA"], ValueError)
    for t in ("L1", "Alu"):
        assert not isinstance(errs[t], TimeoutError), f"{t} hung on a failed enumeration"
