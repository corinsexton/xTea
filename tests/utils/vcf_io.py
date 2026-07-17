"""Parse xTea MEI VCFs into normalized call records.

xTea emits one VCF per MEI type, e.g. ``NA12878.chr1_ALU.vcf``. Records look like::

    chr1  23654637  .  A  <INS:ME:ALU>  .  PASS  SVTYPE=INS:ME:ALU;SVLEN=284;...  GT  1/1

The MEI family is carried both in the symbolic ALT allele (``<INS:ME:ALU>``) and in the
``SVTYPE`` INFO field. We normalize it to a canonical family token so calls from the new
pipeline can be compared to the master golden regardless of which field is populated.

This module is intentionally dependency-free (pure text parsing) so the concordance gate has
no import cost and runs anywhere.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

# Canonical MEI family tokens used everywhere in the harness.
ALU = "ALU"
LINE1 = "LINE1"
SVA = "SVA"
HERV = "HERV"
MT = "MT"
PSEUDOGENE = "PSEUDOGENE"

# Map the various spellings xTea uses (ALT allele, SVTYPE, internal short names) → canonical.
_MEI_ALIASES = {
    "ALU": ALU, "INS:ME:ALU": ALU,
    "LINE1": LINE1, "L1": LINE1, "LINE": LINE1, "INS:ME:LINE1": LINE1,
    "SVA": SVA, "INS:ME:SVA": SVA,
    "HERV": HERV, "HERV-K": HERV, "HERVK": HERV, "INS:ME:HERV-K": HERV,
    "MT": MT, "MITOCHONDRION": MT, "INS:MT": MT,
    "PSEUDOGENE": PSEUDOGENE, "PSDGN": PSEUDOGENE, "INS:PSDGN": PSEUDOGENE,
}


def normalize_mei_type(token: str) -> str:
    """Normalize any known MEI spelling to a canonical family token.

    Raises KeyError on an unrecognized token so a format drift is caught loudly rather than
    silently producing a non-matching type.
    """
    t = token.strip().upper()
    t = t.replace("<", "").replace(">", "")
    if t in _MEI_ALIASES:
        return _MEI_ALIASES[t]
    raise KeyError(f"unrecognized MEI type token: {token!r}")


@dataclass(frozen=True)
class MEICall:
    """One MEI insertion call, normalized for cross-callset comparison."""
    chrom: str
    pos: int              # 1-based VCF POS
    mei_type: str         # canonical family (ALU/LINE1/SVA/...)
    ref: str = ""
    alt: str = ""
    genotype: str = ""
    info: Dict[str, str] = field(default_factory=dict)

    @property
    def key(self):
        return (self.chrom, self.pos, self.mei_type)


def _parse_info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if info_field in (".", ""):
        return out
    for kv in info_field.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = "True"  # flag
    return out


def _mei_type_of(alt: str, info: Dict[str, str]) -> str:
    # Prefer the symbolic ALT allele; fall back to SVTYPE.
    if alt and alt.startswith("<") and ":ME:" in alt or alt.startswith("<INS:"):
        try:
            return normalize_mei_type(alt)
        except KeyError:
            pass
    svtype = info.get("SVTYPE", "")
    if svtype:
        return normalize_mei_type(svtype)
    raise KeyError(f"cannot determine MEI type from ALT={alt!r} INFO.SVTYPE={info.get('SVTYPE')!r}")


def parse_vcf(path: str) -> List[MEICall]:
    """Parse an xTea VCF file into a list of MEICall.

    The sample genotype is taken from the first sample column's GT sub-field when present.
    """
    calls: List[MEICall] = []
    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filter, info_s = f[:8]
            info = _parse_info(info_s)
            genotype = ""
            if len(f) >= 10:
                fmt_keys = f[8].split(":")
                fmt_vals = f[9].split(":")
                if "GT" in fmt_keys:
                    genotype = fmt_vals[fmt_keys.index("GT")]
            calls.append(
                MEICall(
                    chrom=chrom,
                    pos=int(pos),
                    mei_type=_mei_type_of(alt, info),
                    ref=ref,
                    alt=alt,
                    genotype=genotype,
                    info=info,
                )
            )
    return calls


def parse_vcfs(paths) -> List[MEICall]:
    """Parse and concatenate several VCFs (e.g. one per MEI type) into a single call list."""
    out: List[MEICall] = []
    for p in paths:
        out.extend(parse_vcf(p))
    return out
