"""Clip/disc read detection: golden parsers + an independent pysam oracle.

Two jobs:

1. Parse master's clip/disc *site* tables (``candidate_list_from_clip.txt``,
   ``candidate_list_from_disc.txt``) into position→count maps, for site-level comparison.

2. Provide an **oracle** that re-derives the set of clipped / discordant reads directly from a
   BAM/CRAM region using pysam, implementing exactly the detection rules in
   ``docs/clipdisc_contract_spec.md`` (§1.2–1.3, §2.2). This is the read-level ground truth the
   semantic tests use: it is independent of master's opaque intermediate files, and it is the
   reference the C ``xtea-scan`` output is checked against.

Everything here is short-read only. Constants mirror ``xtea/global_values.py`` (code values, not
the sometimes-wrong comments).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

# --- constants from global_values.py (code values) -------------------------------------------
MINIMUM_CLIP_MAPQ = 12
MINIMUM_DISC_MAPQ = 20
CLIP_PHRED_SCORE_CUTOFF = 15
DISC_THRESHOLD = 2000          # same-chrom insert-size cutoff for "discordant"
FLAG_LEFT_CLIP = "L"
FLAG_RIGHT_CLIP = "R"

# pysam CIGAR op codes
_CIGAR_M, _CIGAR_I, _CIGAR_D, _CIGAR_N, _CIGAR_S, _CIGAR_H = 0, 1, 2, 3, 4, 5
_CLIP_OPS = (_CIGAR_S, _CIGAR_H)
_REF_CONSUMING_SKIP = (_CIGAR_S, _CIGAR_H, _CIGAR_I)  # ops NOT counted toward right-clip ref offset


# ============================================================================================
# 1. Golden site-table parsers
# ============================================================================================
@dataclass(frozen=True)
class ClipSite:
    chrom: str
    pos: int
    left: int
    right: int
    extra: Tuple[int, ...] = ()   # remaining count columns (mate_in_rep, realign_l, realign_r, ...)


@dataclass(frozen=True)
class DiscSite:
    chrom: str
    pos: int
    n_left: int
    n_right: int


def parse_clip_sites(path: str) -> Dict[Tuple[str, int], ClipSite]:
    """Parse ``candidate_list_from_clip.txt`` → {(chrom,pos): ClipSite}.

    Columns: chrom, pos, left_cnt, right_cnt, [mate_in_rep, realign_l, realign_r, ...].
    Robust to trailing columns.
    """
    out: Dict[Tuple[str, int], ClipSite] = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chrom, pos = f[0], int(f[1])
            left, right = int(f[2]), int(f[3])
            extra = tuple(int(x) for x in f[4:] if x.lstrip("-").isdigit())
            out[(chrom, pos)] = ClipSite(chrom, pos, left, right, extra)
    return out


def parse_disc_sites(path: str) -> Dict[Tuple[str, int], DiscSite]:
    """Parse ``candidate_list_from_disc.txt`` → {(chrom,pos): DiscSite}.

    Columns: chrom, pos, n_left, n_right.
    """
    out: Dict[Tuple[str, int], DiscSite] = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 4:
                continue
            chrom, pos = f[0], int(f[1])
            out[(chrom, pos)] = DiscSite(chrom, pos, int(f[2]), int(f[3]))
    return out


# ============================================================================================
# 2. Oracle: re-derive clip/disc reads from a BAM/CRAM region (the executable spec)
# ============================================================================================
@dataclass(frozen=True)
class ClipRead:
    """One clipped read as detected from the alignment record."""
    read_name: str
    chrom: str
    clip_pos: int          # 0-based breakpoint position (see spec §1.3)
    side: str              # "L" or "R"
    is_reverse: bool
    clip_len: int


def _is_qualified_clip(quals, n: int) -> bool:
    """Spec §1.2: reject if > 50% of the n clipped bases are below CLIP_PHRED_SCORE_CUTOFF."""
    if n <= 0:
        return False
    half = n / 2
    low = 0
    for i in range(n):
        if quals[i] < CLIP_PHRED_SCORE_CUTOFF:
            low += 1
            if low > half:
                return False
    return True


def detect_clip_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    reference: Optional[str] = None,
    min_mapq: int = MINIMUM_CLIP_MAPQ,
    apply_baseq: bool = True,
    keep_supplementary: bool = True,
) -> List[ClipRead]:
    """Oracle for the clip step, over a region [start, end) (0-based half-open).

    Implements spec §1.2–1.3: skip dup/secondary/unmapped; MAPQ ≥ min_mapq;
    a read is left-clipped iff CIGAR[0] is S/H and right-clipped iff CIGAR[-1] is S/H;
    base-quality gate over the clipped bases only; left breakpoint = reference_start, right
    breakpoint = reference_start + Σ(ref-consuming ops except the last). A read owned by the
    region iff its breakpoint falls in [start, end) — mirrors the tile POS-ownership rule.

    Supplementary (split-read) alignments carry real breakpoint clips for long elements
    (SVA/L1) and master's clip collection counts them (verified: golden nL=11 = 4 primary +
    7 non-dup supplementary at chr1:168451690), so they are KEPT by default. Set
    keep_supplementary=False to reproduce the old primary-only behavior. Mirrors xtea-scan's
    --no-clip-supplementary toggle.
    """
    import pysam

    reads: List[ClipRead] = []
    kwargs = {}
    if reference:
        kwargs["reference_filename"] = reference
    with pysam.AlignmentFile(bam_path, "rc" if bam_path.endswith(".cram") else "rb", **kwargs) as af:
        for a in af.fetch(chrom, start, end):
            if a.is_duplicate or a.is_secondary or a.is_unmapped:
                continue
            if a.is_supplementary and not keep_supplementary:
                continue
            if a.mapping_quality < min_mapq:
                continue
            cig = a.cigartuples
            if not cig:
                continue
            quals = a.query_qualities
            # --- left clip ---
            if cig[0][0] in _CLIP_OPS:
                clen = cig[0][1]
                pos = a.reference_start
                ok = True
                if apply_baseq and quals is not None and cig[0][0] == _CIGAR_S:
                    ok = _is_qualified_clip(quals[:clen], clen)
                if ok and start <= pos < end:
                    reads.append(ClipRead(a.query_name, chrom, pos, FLAG_LEFT_CLIP,
                                          a.is_reverse, clen))
            # --- right clip ---
            if cig[-1][0] in _CLIP_OPS:
                clen = cig[-1][1]
                pos = a.reference_start
                for op, ln in cig[:-1]:
                    if op not in _REF_CONSUMING_SKIP:
                        pos += ln
                ok = True
                if apply_baseq and quals is not None and cig[-1][0] == _CIGAR_S:
                    ok = _is_qualified_clip(quals[len(quals) - clen:], clen)
                if ok and start <= pos < end:
                    reads.append(ClipRead(a.query_name, chrom, pos, FLAG_RIGHT_CLIP,
                                          a.is_reverse, clen))
    return reads


@dataclass(frozen=True)
class DiscRead:
    read_name: str
    chrom: str
    pos: int               # 0-based read position
    mate_chrom: str
    mate_pos: int
    is_reverse: bool


def detect_disc_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    reference: Optional[str] = None,
    min_mapq: int = MINIMUM_DISC_MAPQ,
    is_threshold: int = DISC_THRESHOLD,
) -> List[DiscRead]:
    """Oracle for the disc step, over a region [start, end) (0-based half-open).

    Implements spec §2.2: skip dup/secondary/supplementary/unmapped/mate-unmapped;
    MAPQ ≥ min_mapq; discordant iff mate on a different chromosome, or same chrom with
    |mate_pos - pos| > is_threshold. Classified purely from the read's own record (RNEXT/PNEXT),
    so it is tile-safe.
    """
    import pysam

    reads: List[DiscRead] = []
    kwargs = {}
    if reference:
        kwargs["reference_filename"] = reference
    with pysam.AlignmentFile(bam_path, "rc" if bam_path.endswith(".cram") else "rb", **kwargs) as af:
        for a in af.fetch(chrom, start, end):
            if (a.is_duplicate or a.is_secondary or a.is_supplementary
                    or a.is_unmapped or a.mate_is_unmapped):
                continue
            if a.mapping_quality < min_mapq:
                continue
            if not (start <= a.reference_start < end):
                continue
            mate_chrom = a.next_reference_name
            mate_pos = a.next_reference_start
            if mate_chrom is None:
                continue
            discordant = (mate_chrom != chrom) or (abs(mate_pos - a.reference_start) > is_threshold)
            if discordant:
                reads.append(DiscRead(a.query_name, chrom, a.reference_start,
                                      mate_chrom, mate_pos, a.is_reverse))
    return reads


# ============================================================================================
# 3. Aggregation + set comparison (the "same reads found" invariant)
# ============================================================================================
def clip_read_keys(reads: List[ClipRead]) -> Set[Tuple[str, str, int, str]]:
    """Read-identity set: (read_name, chrom, clip_pos, side). Used for set-equality tests."""
    return {(r.read_name, r.chrom, r.clip_pos, r.side) for r in reads}


def disc_read_keys(reads: List[DiscRead]) -> Set[Tuple[str, str, int]]:
    return {(r.read_name, r.chrom, r.pos) for r in reads}


@dataclass
class SetDiff:
    n_expected: int
    n_actual: int
    only_expected: set = field(default_factory=set)   # missing from actual
    only_actual: set = field(default_factory=set)      # spurious in actual

    @property
    def equal(self) -> bool:
        return not self.only_expected and not self.only_actual

    def summary(self) -> str:
        return (f"expected={self.n_expected} actual={self.n_actual} "
                f"missing={len(self.only_expected)} spurious={len(self.only_actual)}")


def compare_read_sets(expected: set, actual: set) -> SetDiff:
    """Compare two read-identity sets; the semantic invariant is `equal`."""
    return SetDiff(
        n_expected=len(expected),
        n_actual=len(actual),
        only_expected=expected - actual,
        only_actual=actual - expected,
    )


# ============================================================================================
# 4. Runner for the C++ xtea-scan prototype (parses its TSV into the same read-identity keys)
# ============================================================================================
def _scan_cmd(scan_bin, aln_path, chrom, start, end, reference, threads, tile_size):
    cmd = [scan_bin, "--bam", aln_path, "--region", f"{chrom}:{start}-{end}", "--threads", str(threads)]
    if reference:
        cmd += ["--reference", reference]
    if tile_size:
        cmd += ["--tile-size", str(tile_size)]
    return cmd


def run_xtea_scan_raw(scan_bin, aln_path, chrom, start, end,
                      reference=None, threads=1, tile_size=None) -> str:
    """Return xtea-scan's raw per-read stdout (for byte-level determinism checks).

    Per-read output is opt-in (memory-lean default), so we request it explicitly with --out -.
    """
    import subprocess
    cmd = _scan_cmd(scan_bin, aln_path, chrom, start, end, reference, threads, tile_size) + ["--out", "-"]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def run_xtea_scan(scan_bin: str, aln_path: str, chrom: str, start: int, end: int,
                  reference: Optional[str] = None, threads: int = 1, tile_size=None):
    """Invoke xtea-scan and return (clip_keys, disc_keys) matching the oracle keys.

    Per-read TSV lines:
      CLIP chrom clip_pos side read_name is_rev clip_len -> (read_name, chrom, clip_pos, side)
      DISC chrom pos read_name mate_chrom mate_pos is_rev -> (read_name, chrom, pos)
    """
    out = run_xtea_scan_raw(scan_bin, aln_path, chrom, start, end, reference, threads, tile_size)
    clip_keys: Set[Tuple[str, str, int, str]] = set()
    disc_keys: Set[Tuple[str, str, int]] = set()
    for line in out.splitlines():
        f = line.split("\t")
        if not f:
            continue
        if f[0] == "CLIP":
            clip_keys.add((f[4], f[1], int(f[2]), f[3]))
        elif f[0] == "DISC":
            disc_keys.add((f[3], f[1], int(f[2])))
    return clip_keys, disc_keys


def run_xtea_scan_clip_sites(scan_bin, aln_path, chrom, start, end,
                             reference=None, threads=1, tile_size=None):
    """Run xtea-scan and return its aggregated clip-site table as {(chrom,pos): (n_left,n_right)}."""
    import subprocess, tempfile, os
    fd, path = tempfile.mkstemp(suffix=".clipsites")
    os.close(fd)
    try:
        cmd = _scan_cmd(scan_bin, aln_path, chrom, start, end, reference, threads, tile_size)
        cmd += ["--out", "", "--clip-sites", path]  # suppress per-read stream; emit sites only
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        sites = {}
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 4:
                    sites[(f[0], int(f[1]))] = (int(f[2]), int(f[3]))
        return sites
    finally:
        os.unlink(path)


def aggregate_clip_sites_from_oracle(clip_reads) -> dict:
    """Aggregate oracle ClipReads into {(chrom,pos): (n_left,n_right)} for comparison to --clip-sites."""
    agg = {}
    for r in clip_reads:
        nl, nr = agg.get((r.chrom, r.clip_pos), (0, 0))
        if r.side == FLAG_LEFT_CLIP:
            nl += 1
        else:
            nr += 1
        agg[(r.chrom, r.clip_pos)] = (nl, nr)
    return agg
