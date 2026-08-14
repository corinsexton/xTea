"""Genomic regions excluded from candidate calling, applied once after the scan pass.

WHY IT IS APPLIED TO CANDIDATE POSITIONS AND NOT TO THE DISC STREAM. The obvious reading of
"mask these regions" is to drop every scan record inside them, and for `scan_disc_pos.tsv`
that would be a bug, not an optimisation. Those records are not candidate sites -- they are
EVIDENCE, windowed +/-iextend around each clip peak by disc_adapter.run_disc_from_scan:

    for pos in m_peaks[chrm]:
        nL, ... = _count_window(disc_arrays, chrm, pos - iextend, pos, ...)
        nR, ... = _count_window(disc_arrays, chrm, pos + 1, pos + iextend, ...)

iextend runs 810-874 bp on these samples. Deleting disc rows by their own position would
therefore strip supporting reads from candidates OUTSIDE the mask: a real insertion 500 bp
from a segdup boundary would lose one flank's discordant support, its nL/nR would fall, and
it could drop below disc_cutoff and vanish. With 655k intervals in the default mask, every
candidate within ~874 bp of any of them is exposed.

So the mask is applied to CLIP CANDIDATE POSITIONS only. That is sufficient to skip the work:
`m_peaks` is peak-called from the clip candidate list, and the raw-disc rows that feed
transduction's site list are built per candidate peak, so a position removed here is removed
from every downstream stage. `disc_arrays` stays complete, so surviving candidates keep full
evidence on both flanks.

MEASURED on the HapMap arms (union of the three defaults, 287 Mb, 9.3% of the genome):

    stage                          broad 173x   bcm 450x
    clip sites inside the mask        9.4%       12.0%
    transduction .slct sites          25.5% (L1) 25.3% (SVA)
    disc-pos rows inside the mask     36.2%      23.9%   <- NOT filtered, see above

and on the call sets, false positives are enriched in these regions by 10-20x over true
positives (simple repeats: 0-3% of TPs vs 13-32% of FPs), which is what makes this a
precision fix that also happens to save time rather than a pure speed knob.
"""
import bisect
import gzip
import os

# The SMaHT SNV pipeline's own masks (smaht_snv_pipeline/main.nf params), so an MEI call set
# and an SNV call set from the same project exclude the same regions.
_SNV_NEXTFLOW = "/home/cos689/smaht/test_benchmarking/smaht_snv_pipeline/nextflow"
_SNV_NEXTFLOW_ALT = ("/n/data1/hms/dbmi/park/corinne/smaht/test_benchmarking/"
                     "smaht_snv_pipeline/nextflow")

DEFAULT_MASKS = {
    "simple_repeats": f"{_SNV_NEXTFLOW_ALT}/simple_repeats.bed",
    "segdup": f"{_SNV_NEXTFLOW}/segdup_GRCh38_official.bed.gz",
    "centromeres": f"{_SNV_NEXTFLOW}/centromeres_GRCh38_official.bed.gz",
}


class RegionMask:
    """Merged BED intervals with an O(log n) point-containment test.

    Empty (`RegionMask([])`) is a valid, cheap no-op: `contains` is always False, so callers
    need no None checks.
    """

    __slots__ = ("_iv", "_starts", "n_intervals", "n_bases", "sources")

    def __init__(self, paths=()):
        raw = {}
        self.sources = []
        for p in paths:
            self.sources.append(p)
            opener = gzip.open if str(p).endswith(".gz") else open
            with opener(p, "rt") as fin:
                for line in fin:
                    if not line or line[0] in "#tb":      # comment / track / browser
                        if not line[:1].isalnum():
                            continue
                    f = line.split()
                    if len(f) < 3:
                        continue
                    try:
                        raw.setdefault(f[0], []).append((int(f[1]), int(f[2])))
                    except ValueError:
                        continue                          # header row
        self._iv = {}
        self.n_intervals = 0
        self.n_bases = 0
        for chrom, spans in raw.items():
            spans.sort()
            merged = []
            for s, e in spans:
                if merged and s <= merged[-1][1]:
                    if e > merged[-1][1]:
                        merged[-1] = (merged[-1][0], e)
                else:
                    merged.append((s, e))
            self._iv[chrom] = merged
            self.n_intervals += len(merged)
            self.n_bases += sum(e - s for s, e in merged)
        self._starts = {c: [s for s, _ in v] for c, v in self._iv.items()}

    def __bool__(self):
        return self.n_intervals > 0

    def contains(self, chrom, pos):
        """True when `pos` falls in a masked interval. BED is half-open [start, end)."""
        spans = self._iv.get(chrom)
        if not spans:
            return False
        i = bisect.bisect_right(self._starts[chrom], pos) - 1
        return i >= 0 and spans[i][1] > pos

    def describe(self):
        if not self:
            return "no region mask"
        return (f"{self.n_intervals} merged intervals, {self.n_bases / 1e6:.0f} Mb "
                f"({100.0 * self.n_bases / 3.1e9:.1f}% of genome) from "
                f"{len(self.sources)} file(s)")


def resolve_mask_spec(spec, allow_missing=False):
    """Turn a --region-mask value into a list of paths.

    Accepts a comma-separated mix of DEFAULT_MASKS keys and explicit file paths, so
    `--region-mask simple_repeats,segdup` and `--region-mask /my/own.bed` both work, as does
    mixing them. "none" disables masking entirely.
    """
    if spec is None:
        spec = "default"
    spec = str(spec).strip()
    if spec.lower() in ("none", "off", ""):
        return []
    out = []
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.lower() == "default":
            out.extend(DEFAULT_MASKS.values())
        elif tok in DEFAULT_MASKS:
            out.append(DEFAULT_MASKS[tok])
        else:
            out.append(tok)
    missing = [p for p in out if not os.path.isfile(p)]
    if missing and not allow_missing:
        raise FileNotFoundError(
            "region mask file(s) not found: {0}. Known names: {1}. Use "
            "--region-mask none to disable.".format(
                ", ".join(missing), ", ".join(sorted(DEFAULT_MASKS))))
    return [p for p in out if os.path.isfile(p)]


def load(spec, allow_missing=False):
    return RegionMask(resolve_mask_spec(spec, allow_missing=allow_missing))
