"""The unified downstream's front half: the locus universe, enumerated by the DRIVER.

Plan sec 3 draws the downstream as "one pass over loci". This module is the part of that
pass which can be built today: the enumeration of the loci themselves, done once in the
driver for every family, before any family's transduction step runs.

WHERE THIS SITS, AND THE THREE FACTS THAT FIX IT THERE
------------------------------------------------------

1. **The iteration domain is knowable up front.** `xtea/site_universe.py` is the
   site-creation rule lifted out of `XTransduction.re_slct_with_clip_raw_disc_sites`; it
   needs no BAM, no bwa and no cns *evidence*, and it is what creates the 52 of BASELINE-0's
   1,515 final records that have no `candidate_disc_filtered_cns.txt` row. Running it for
   every family and taking the union gives **41,604 loci** where the three per-family chains
   between them visit 75,853 (1.82x).

2. **But not before `cns_filter`.** The rule takes the cns candidate *positions* for its
   proximity exclusion, and that file is `cns_filter`'s output. So the earliest point the
   universe is computable is **after every family's `cns_filter` and before any family's
   transduction** -- which is exactly where the Stage 2c family-assignment barrier already
   sits (`run_xtea.run_family_assign_barrier`). `site_universe`'s docstring says the rule's
   inputs are "available before cns_filter runs"; that is true of raw_disc and the
   annotation and false of the cns list, and the difference is the whole placement question.
   Corrected at both docstrings.

3. **Family assignment does NOT move up here with it.** For those same 52 records the only
   per-family evidence is the transduction evaluation's own output, so the 7-of-22 barrier
   ceiling is unchanged by knowing the site exists. `.slct` is a 34-40K-row superset per
   family. The order the data supports is *enumerate loci -> gather per-family evidence ->
   assign family -> filter*, and this module is only the first box.

WHAT THE LOOP CANNOT BE, read before designing its other half
-------------------------------------------------------------
`call_candidate_transduction_v3` is **set-oriented by construction**: it builds one flank
reference and runs four batched passes over the whole `.slct` list (bwa index, clip realign,
disc realign to flanks, disc realign to consensus), then parses the SAMs into per-locus
records. A literal per-locus loop cannot replace that -- one bwa per locus is not a
refactor, it is a different program. The shape the code supports is

    enumerate (here)  ->  batched per-family evidence production  ->  per-locus decision

where only the third box is a loop. The per-locus evidence record it loops over already
exists in today's code, and its fields are exactly what `collect_ins_info` assembles --
`(n_l_disc, lcns_start, lcns_end, n_r_disc, rcns_start, rcns_end, n_cns_clip, n_polyA_clip,
i_lcov2, i_rcov2, s_src)` -- plus the 10-field genotype tuple `_collect_gntp_features`
returns. That answers "what does a per-family evidence record look like at a locus only
`.slct` knows about": the same 11+10 fields as at any other locus, which is why
`add_new_rcd_with_cutoff` can write a created site in the same format as an updated one.

WHAT THIS MODULE DOES NOT DO: it does not change any call. Each family's transduction step
still computes its own `.slct` (through the same one copy of the rule), so the driver's
files are a *replica*, written to their own directory. That is deliberate -- it makes the
production check a byte comparison rather than an argument, and the replica is what the
`--td-manifest` caller will consume next.

SIBLING IS NOT PART OF THE UNIVERSE. The sibling step derives its own `<raw_disc>.slct2`
(`b_with_candidate=False`, annotation boundary extend = mean insert size instead of 450), but
`probe_site_creation.py` measured that it creates **zero** final records -- every one of the
52 comes from transduction. So the universe is the transduction `.slct` per family, and
adding `.slct2` would widen it without adding a locus any record is ever written at.
"""
import os
from pathlib import Path

import site_universe

# NOT a module-level `from x_transduction import XTransduction`: x_transduction imports THIS
# module (for ENV_TD_PRECOMPUTED, the driver-gather bypass), so a top-level import either way
# is a cycle. The step is the authority for the annotation index this module rebuilds, so the
# dependency has to point this way; deferring it to the one function that needs it is the
# smaller of the two evils.

# x_TEA_main.py:519 (transduction) and :614 (sibling) -- the branches' own literal.
I_MIN_COPY_LEN = 225

# XTransduction.__init__'s `_boundary_extnd`. The transduction branch never calls
# set_boundary_extend, so 450 is what that step's annotation index is built at. (The
# sibling branch DOES -- it passes mean insert size - which is one of the reasons its
# `.slct2` is a different file; see the module docstring.)
TD_BOUNDARY_EXTND = 450

# The two `.slct` layouts, as `site_universe.write_slct_file` writes them: a cns candidate
# is a bare `chrom<TAB>pos`, a newly selected raw-disc row carries its full 6 fields.
CLAIM_CNS = "cns"
CLAIM_NEW = "new"


def transduction_extnd(cns_extnd):
    """The transduction step's `iextnd`, from the numbers the driver already has.

    `x_TEA_main.py:473-499` sets `iextnd = 400` and then raises it to
    `max_is = int(mean_is + 3*std_var)` if that is larger. The driver computes the same
    quantity at `run_xtea.py`'s `cns_extnd` (documented there as "exactly what x_TEA_main's
    cns step computes as iextnd"), so the only thing it was missing is the 400 floor -- and
    at WGS insert sizes the floor never binds (HG002: 954). Both steps read the same cached
    `tmp_basic/basic_cov_is_rlth_info.txt`, so the mean/std cannot differ between them.

    This is the last argument a driver-side `--td-manifest` caller needed (plan sec 4 2b.4).

    ONE DEPENDENCY WORTH KNOWING, because `--td-gather-driver` widens its blast radius. The
    step does NOT read the driver's cached basic-info file: `X_BasicInfo` keys that cache on
    its own working folder, so each step writes its own copy under `<type>/tmp/<step>/`. They
    agree because the coverage/insert-size sampling is SEEDED by default
    (`global_values.BASIC_INFO_RANDOM_SEED`; verified -- all seven copies in the 2c gate run
    carry the identical 31.816 / 151 / 418.8067714547886 / 178.46596377929833). The pipeline
    already relies on this identity for the cns-emit collect window. What is new here is the
    consequence of breaking it: with `XTEA_BASIC_INFO_SEED=none` the driver's gather window
    and the step's evaluation could be computed from different estimates, so do not run that
    combination and expect byte-identical output.
    """
    return max(400, int(cns_extnd))


def universe_dir(sample_root):
    return str(Path(sample_root) / "merged" / "locus_universe")


def family_slct_path(sample_root, tag):
    return os.path.join(universe_dir(sample_root), "{0}.slct".format(tag))


def enumerate_family_sites(sf_raw_disc, sf_cns, sf_annotation, i_rep_type, n_disc_cutoff,
                          b_tumor, sf_out, b_flank=True, n_jobs=1, sf_reference=None,
                          working_folder="/tmp", xannotation=None, log=print):
    """Write one family's `.slct` replica -- the sites its transduction step will evaluate.

    Returns (n_candidates, n_new).

    b_flank=False is the flankless family (Alu: no transduction flank, so x_TEA_main copies
    its input through and the rule never runs). Its universe contribution is its cns
    candidate list and nothing else, which `probe_site_creation.py` confirms from the other
    direction: all 1,230 of BASELINE-0's Alu records first exist at cns. The file written
    here is then driver-only -- the step produces no `.slct` for it, so a byte gate must
    expect one fewer comparison than there are families, not treat the absence as a skip.

    xannotation is injectable so tests need no real rmsk `.out`; production leaves it None
    and gets the index the step itself would build, through the step's own method, at the
    step's own boundary extend.
    """
    Path(os.path.dirname(sf_out)).mkdir(parents=True, exist_ok=True)
    if not b_flank:
        n_cand = 0
        with open(sf_out, "w") as fout:
            for chrm, pos in _iter_cns_positions(sf_cns):
                fout.write("{0}\t{1}\n".format(chrm, pos))
                n_cand += 1
        return n_cand, 0
    if xannotation is None:
        from x_transduction import XTransduction   # deferred: see the module's import block
        xtd = XTransduction(working_folder, n_jobs, sf_reference)
        assert xtd._boundary_extnd == TD_BOUNDARY_EXTND, (
            "the transduction step's annotation boundary extend changed; the universe "
            "replica must be built at the same one or its third rejection differs")
        xannotation = xtd.prep_annotation_interval_tree(sf_annotation, I_MIN_COPY_LEN)
    return site_universe.write_slct_file(sf_raw_disc, sf_cns, n_disc_cutoff, xannotation,
                                         i_rep_type, b_tumor, sf_out, log=log)


def _iter_cns_positions(sf_cns):
    """(chrom, pos) per cns candidate, in the cns file's own order.

    `write_slct_file` emits candidates via `load_in_candidate_list_str_version`, a dict of
    dicts; on CPython 3.7+ that preserves insertion order, so file order and dict order are
    the same sequence. Read straight from the file here so the flankless path does not
    depend on that.
    """
    with open(sf_cns) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            yield fields[0], pos


def read_slct_claims(sf_slct):
    """{chrom: {pos: CLAIM_CNS|CLAIM_NEW}} from one family's `.slct`.

    The two layouts are the provenance: a 2-field line is a cns candidate this family
    already had, a longer line is a site the rule CREATED. Nothing else in the file
    distinguishes them, and the distinction is what the 52-of-1,515 finding is about.

    A position present as both (possible only if the cns list itself repeats a position --
    the proximity exclusion suppresses any raw-disc row within 100 bp of a candidate) is
    recorded as CLAIM_CNS: the cns row is the one the barrier can already see.
    """
    claims = {}
    with open(sf_slct) as fin:
        for line in fin:
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            kind = CLAIM_CNS if len(fields) == 2 else CLAIM_NEW
            per_chrm = claims.setdefault(fields[0], {})
            if per_chrm.get(pos) != CLAIM_CNS:
                per_chrm[pos] = kind
    return claims


class LocusUniverse:
    """The union of every family's `.slct`, at exact position.

    Exact position, not a window: the post-cns chain does not re-anchor coordinates (every
    one of the 52 created records sits at distance 0 in `candidate_disc_filtered_cns2.txt`),
    so two families naming the same locus name the same integer. Merging with a window here
    would fold genuinely distinct loci together and hide it inside the universe count --
    cross-family position reconciliation is `locus_table`'s job, at its own MERGE_WINDOW,
    and it happens after evidence exists.

    **THE INPUT MAP'S ITERATION ORDER IS THE OUTPUT FILE'S COLUMN ORDER.** `rows()` emits
    the per-locus `families`/`claims` pairs in `self.families` order, which is the order
    `claims_by_family` was built in. Two runs that agree on every locus and every claim
    still produce different BYTES if they built that map in different orders -- job 49596361
    failed on exactly that, `SVA,L1` against `L1,SVA` over an identical 3,272 loci, once the
    staged barrier started populating the map in cns_filter-completion order. Callers own
    the order; it is not sorted here, because the driver's family order (L1, Alu, SVA as
    configured) is the meaningful one and sorting would silently rewrite every table already
    gated.
    """

    def __init__(self, claims_by_family):
        self.families = list(claims_by_family)
        self.loci = {}          # {chrom: {pos: {family: claim}}}
        for tag, claims in claims_by_family.items():
            for chrm, per_chrm in claims.items():
                dst = self.loci.setdefault(chrm, {})
                for pos, kind in per_chrm.items():
                    dst.setdefault(pos, {})[tag] = kind

    @property
    def n_loci(self):
        return sum(len(v) for v in self.loci.values())

    def per_family(self):
        """{family: n_loci it claims}."""
        out = dict.fromkeys(self.families, 0)
        for per_chrm in self.loci.values():
            for by_fam in per_chrm.values():
                for tag in by_fam:
                    out[tag] += 1
        return out

    def by_n_families(self):
        """{how many families claim this locus: n_loci}."""
        out = {}
        for per_chrm in self.loci.values():
            for by_fam in per_chrm.values():
                out[len(by_fam)] = out.get(len(by_fam), 0) + 1
        return out

    def n_created_only(self):
        """Loci NO family has a cns row for -- the population the barrier cannot see.

        This is the 52-of-1,515 finding's own quantity, measured on the universe instead of
        on the final call set: a locus every claiming family reached only by creation. It is
        the ceiling on what a decision point placed before transduction could ever resolve,
        so it belongs next to the universe count rather than in a probe.
        """
        n = 0
        for per_chrm in self.loci.values():
            for by_fam in per_chrm.values():
                if all(k == CLAIM_NEW for k in by_fam.values()):
                    n += 1
        return n

    def rows(self):
        """(chrom, pos, n_families, [(family, claim)]) sorted by (chrom, pos).

        Plain string sort on the chromosome: this file is an index for the loop and for
        diffing two runs, and any ordering does for that as long as it is deterministic.
        """
        for chrm in sorted(self.loci):
            for pos in sorted(self.loci[chrm]):
                by_fam = self.loci[chrm][pos]
                yield (chrm, pos, len(by_fam),
                       [(t, by_fam[t]) for t in self.families if t in by_fam])

    def write(self, sf_out):
        Path(os.path.dirname(sf_out) or ".").mkdir(parents=True, exist_ok=True)
        with open(sf_out, "w") as fout:
            fout.write("#chrom\tpos\tn_families\tfamilies\tclaims\n")
            for (chrm, pos, n_fam, pairs) in self.rows():
                fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(
                    chrm, pos, n_fam,
                    ",".join(t for (t, _k) in pairs),
                    ",".join(k for (_t, k) in pairs)))
        return sf_out


def build_locus_universe(slct_by_family):
    """{family: path to its `.slct`} -> LocusUniverse. Missing files are skipped."""
    claims = {}
    for tag, path in slct_by_family.items():
        if path and os.path.isfile(path):
            claims[tag] = read_slct_claims(path)
    return LocusUniverse(claims)


# --------------------------------------------------------------------------------------
# The batched evidence stage's first half: gather every family's reads out of ONE BAM
# decode, in the driver, instead of once per family inside each transduction step.
#
# The C++ side has existed since 2026-08-02 (`run_td_gather_multi`, `xtea_scan
# --td-manifest`, gated 6/6 byte-identical at WGS) with no production caller, because
# hoisting it meant re-implementing the step's site selection. It does not any more --
# `site_universe` is that selection and `enumerate_family_sites` above already runs it
# for every family at the barrier. What was still missing was the step's `iextnd`, and
# `transduction_extnd` is it.
#
# WHAT THIS IS AND IS NOT WORTH. Not wall clock, and on chr1 measurably WORSE than not
# doing it -- the per-family gathers already overlap as concurrent subprocesses, this one
# is serialized at the barrier. Measured numbers are in run_xtea.run_td_gather_driver
# (job 49548949): ~25s of CPU returned, ~12s of wall clock lost. The reason to do it is
# that it is the batched evidence stage of the per-locus loop, built once instead of
# twice, and the only place one BAM decode can serve every family.
# --------------------------------------------------------------------------------------

# The two file names XTransduction writes into its working folder
# (x_transduction.call_candidate_transduction_v3). Hard-coded there, so hard-coded here,
# with the step as the authority: a driver-side gather has to land on the exact paths the
# step would have written or the step will not find them.
TD_CLIP_FQ_NAME = "raw_candidate_sites_all_clip.fq"
TD_DISC_FA_NAME = "raw_candidate_sites_all_disc.fa"

# Set by the driver when it has already gathered; read by the step to skip its own gather.
# Same mechanism as XTEA_CNS_PRECOMPUTED_ALIGN, and for the same reason -- the transduction
# step is a separate x_TEA_main subprocess, so an env var is the only channel that does not
# mean threading a new flag through build_type_steps.
ENV_TD_PRECOMPUTED = "XTEA_TD_GATHER_PRECOMPUTED"


def td_gather_paths(tmp_td_dir):
    """(clip FASTQ, disc FASTA) the transduction step expects in its working folder."""
    return (os.path.join(tmp_td_dir, TD_CLIP_FQ_NAME),
            os.path.join(tmp_td_dir, TD_DISC_FA_NAME))


def write_td_manifest(entries, sf_out):
    """`xtea_scan --td-manifest` input: one line per family.

    entries: [(tag, sf_sites, extnd, sf_clip_fq, sf_disc_fa)]. Layout is
    `tag sites extnd clip_fastq disc_fa [disc_pos]` with `-` for an output the entry does
    not want; the disc-pos stream is the legacy `all_disc_pos.disc_pos` file that nothing
    downstream reads, so it is always `-` here.

    The manifest is deliberately strict on the C++ side -- a short line is a hard error and
    two entries sharing an output path are rejected up front -- because a dropped field
    there means a family silently does not run, which on a WGS gate reads as "that family
    found nothing". Nothing is written for a family with no sites: an empty site list would
    produce an empty FASTQ and the step would take it as a real (empty) gather.
    """
    Path(os.path.dirname(sf_out) or ".").mkdir(parents=True, exist_ok=True)
    n = 0
    with open(sf_out, "w") as fout:
        for (tag, sf_sites, extnd, sf_clip_fq, sf_disc_fa) in entries:
            fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t-\n".format(
                tag, sf_sites, int(extnd), sf_clip_fq, sf_disc_fa))
            n += 1
    return n
