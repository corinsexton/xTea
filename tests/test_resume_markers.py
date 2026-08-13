"""Resume markers for the two phases that cost hours (scan_pass, cns_emit).

WHY THIS GATE EXISTS
--------------------
Everything that has actually gone wrong in the HapMap fleet went wrong AFTER the scan:
the 192G OOMs landed during/just after it, the `| head` pipefail bug failed every arm at
the very end with the VCF already written, and a downstream code fix means a rerun.
Losing a finished 6-hour scan to a fault in a later step is the most expensive failure
mode here and it is pure waste -- the scan's inputs did not change.

The three properties that make reuse SAFE, rather than merely fast, are what this file
pins. Each is a way the feature could silently corrupt a run instead of saving time:

  1. A phase killed mid-write must be redone. The marker is written only after the phase
     returns, so a crash leaves no marker -- even though it leaves a plausible-looking
     truncated output file behind. Existence of the output is never sufficient on its own.
  2. Changed inputs must invalidate. Including a rebuilt scan binary: resuming across
     `make -C src/csrc` would otherwise reuse the previous binary's output.
  3. A marker that outlives its outputs must invalidate. cleanup_successful_run wipes
     everything under sample_root except the VCF -- .resume/ included -- so the normal
     path removes both together; but a partial cleanup or a scratch purge can leave a
     marker pointing at files that are gone, and reuse must not be decided by the marker
     alone.
"""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "xtea"))

import run_xtea  # noqa: E402


def _outputs(tmp_path, *names_and_sizes):
    paths = []
    for name, size in names_and_sizes:
        p = tmp_path / name
        p.write_text("x" * size)
        paths.append(str(p))
    return paths


def test_marker_absent_means_not_done(tmp_path):
    """The baseline, and the mid-write-crash case: outputs on disk, no marker."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_roundtrip_marker_is_honoured(tmp_path):
    outs = _outputs(tmp_path, ("scan.tsv", 100), ("disc.tsv", 50))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is True


def test_changed_fingerprint_invalidates(tmp_path):
    """Different inputs -- or a rebuilt scan binary -- must force a rerun."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp2", outs) is False


def test_deleted_output_invalidates(tmp_path):
    """cleanup_successful_run / a scratch purge removed the output but left the marker."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    os.unlink(outs[0])
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_truncated_output_invalidates(tmp_path):
    """Size is compared, not just existence: a file truncated after the marker was
    written (full scratch, interrupted copy) must not be built on."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    with open(outs[0], "w") as fh:
        fh.write("x" * 10)
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_corrupt_marker_invalidates(tmp_path):
    """A half-written marker (killed during the write) must read as 'not done'."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    (tmp_path / ".resume" / "scan_pass.json").write_text("{not json")
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_phases_are_independent(tmp_path):
    """A completed scan_pass must not make cns_emit look done."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    assert run_xtea.phase_is_done(tmp_path, "cns_emit", "fp1", outs) is False


def test_binary_stamp_changes_when_binary_is_rebuilt(tmp_path):
    """_binary_stamp feeds the fingerprint; a rebuild must change it. This is the
    stale-artifact trap: same inputs, different code, silently reused output."""
    binp = tmp_path / "xtea_scan"
    binp.write_text("v1")
    before = run_xtea._binary_stamp(str(binp))
    binp.write_text("v2-longer")
    os.utime(str(binp), (0, 0))
    after = run_xtea._binary_stamp(str(binp))
    assert before != after
    assert run_xtea._fingerprint(["scan_pass", before]) != \
           run_xtea._fingerprint(["scan_pass", after])


def test_missing_binary_does_not_raise(tmp_path):
    """A stamp for a binary that is not there must degrade, not explode -- an installed
    deployment may have no source tree beside it."""
    assert run_xtea._binary_stamp(str(tmp_path / "nope")) == (str(tmp_path / "nope"), None, None)


def test_cleanup_removes_resume_dir_so_success_does_not_resume(tmp_path):
    """THE --keep-intermediates INTERACTION, made executable.

    cleanup_successful_run wipes everything under sample_root but the VCF, and .resume/
    lives there. So a successful run without --keep-intermediates leaves nothing to
    resume -- which is correct, not a bug: there is no partial work to recover, and the
    next --resume run must start clean rather than reuse a deleted tree.
    """
    sample = "s"
    (tmp_path / f"{sample}.vcf").write_text("##fileformat=VCFv4.2\n")
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    run_xtea.mark_phase_done(tmp_path, "scan_pass", "fp1", outs)
    assert (tmp_path / ".resume").is_dir()

    run_xtea.cleanup_successful_run(tmp_path, sample)

    assert not (tmp_path / ".resume").exists(), ".resume must not survive cleanup"
    assert (tmp_path / f"{sample}.vcf").is_file(), "the deliverable VCF must survive"
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_adopt_existing_is_off_by_default(tmp_path):
    """Without the flag, pre-existing output is never trusted -- the default must stay
    'a marker or nothing', since that is the property that makes resume safe."""
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    assert run_xtea.adopt_existing_phase(tmp_path, "scan_pass", "fp1", outs, False) is False
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is False


def test_adopt_existing_writes_a_reusable_marker(tmp_path):
    outs = _outputs(tmp_path, ("scan.tsv", 100), ("disc.tsv", 5))
    assert run_xtea.adopt_existing_phase(tmp_path, "scan_pass", "fp1", outs, True) is True
    assert run_xtea.phase_is_done(tmp_path, "scan_pass", "fp1", outs) is True


def test_adopt_existing_refuses_when_an_output_is_missing(tmp_path):
    outs = _outputs(tmp_path, ("scan.tsv", 100))
    outs.append(str(tmp_path / "never_written.tsv"))
    assert run_xtea.adopt_existing_phase(tmp_path, "scan_pass", "fp1", outs, True) is False
    assert not (tmp_path / ".resume").exists()


def test_adopt_existing_refuses_an_empty_output(tmp_path):
    """A zero-byte file is the one incomplete state adoption CAN detect; it must."""
    outs = _outputs(tmp_path, ("scan.tsv", 100), ("empty.tsv", 0))
    assert run_xtea.adopt_existing_phase(tmp_path, "scan_pass", "fp1", outs, True) is False
    assert not (tmp_path / ".resume").exists()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
