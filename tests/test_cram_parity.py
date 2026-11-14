"""The same reads must give the same output whether they come as BAM or as CRAM (bug B17).

WHAT WENT WRONG
---------------
1. mate_in_rep_clip_cnt (xtea-scan's 6th realign-out column) could lose reads at random.
   The repeat BEDs are loaded into htslib regidx indexes that every tile thread shares. regidx
   builds each chromosome's bin index lazily, on the first regidx_overlap() for that
   chromosome, without a lock: it publishes a zero-filled index and only then fills it. A
   second thread querying that chromosome mid-build read the zeros as "no repeat here", so a
   read whose mate DID land in a repeat was not counted. Which reads were lost depended on
   thread timing. CRAM decodes at a different pace from BAM, so the two diverged; the same
   BAM twice could too. The fix builds every index before the threads start.
2. The VCF ##reference line came from the input's @SQ UR tag. bwa writes no UR into a BAM
   (-> "##reference=NULL"), while samtools adds UR to every @SQ when it writes a CRAM, so the
   same reads gave different headers. It now names the --reference FASTA of the run.
3. The sibling-transduction source check could query an empty interval [pos, pos). A BAM
   returns no reads for it and a CRAM returns the reads covering pos-1, so only the CRAM could
   accept the source and add calls. An empty interval now counts no reads for either format.
4. A CRAM written against a FASTA that lacks some @SQ contigs is silently encoded with
   embed_ref=2 by htslib, and its .crai then records a span of 0 for every slice, so region
   queries return almost nothing. run_xtea now refuses a CRAM whose index has span-0 slices.

HOW THE RACE IS PROVOKED
------------------------
A tiny genome built in the test: one 64-tile contig of clipped reads, and 16 repeat contigs
with a dense BED. Every tile starts with 16 clipped reads at one site, one mate on each repeat
contig, inside a repeat. With 64 threads, every tile thread queries every repeat contig at
once, so all 16 lazy index builds are contended. Every site's correct count is 16. On the old
binary about half the runs lose at least one read (measured on a 2-core login node); several
runs per format make a pass on the old code very unlikely. On a single core the race cannot
happen and the test passes trivially, which is correct.
"""
import os
import random
import shutil
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

TILE = 20_000
N_TILES = 64
N_REP_CONTIGS = 16
REP_LEN = 400_000
RUNS_PER_FORMAT = 8


def _synthetic_seq(n, seed):
    """Deterministic ACGT with no run of 5, so bwa sees no low-complexity sequence."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        c = rng.choice("ACGT")
        while len(out) >= 4 and all(x == c for x in out[-4:]):
            c = rng.choice("ACGT")
        out.append(c)
    return "".join(out)


def _tool(name, env):
    path = os.environ.get(env) or shutil.which(name) or f"/home/cos689/software/bin/{name}"
    if not os.path.exists(path):
        pytest.skip(f"{name} not found; set ${env}")
    return path


@pytest.fixture(scope="module")
def parity_case(tmp_path_factory):
    """A tiny reference, a BAM, and the same reads as a CRAM, plus the repeat BED and a consensus."""
    pysam = pytest.importorskip("pysam")
    bwa = _tool("bwa", "XTEA2_BWA")
    samtools = _tool("samtools", "XTEA2_SAMTOOLS")
    d = tmp_path_factory.mktemp("cram_parity")

    cons = _synthetic_seq(300, seed=11)
    cons_fa = d / "cns.fa"
    cons_fa.write_text(">CNS\n" + cons + "\n")
    r = subprocess.run([bwa, "index", str(cons_fa)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]

    clip_contig = _synthetic_seq(TILE * N_TILES, seed=1)
    rep_seq = _synthetic_seq(REP_LEN, seed=2)
    contigs = [("chrA", clip_contig)] + [(f"chrR{j}", rep_seq) for j in range(N_REP_CONTIGS)]
    ref_fa = d / "ref.fa"
    with open(ref_fa, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    subprocess.run([samtools, "faidx", str(ref_fa)], check=True)

    # Dense repeat annotation: 10 bp every 20 bp on every repeat contig.
    bed = d / "rep.bed"
    with open(bed, "w") as f:
        for j in range(N_REP_CONTIGS):
            for p in range(0, REP_LEN - 20, 20):
                f.write(f"chrR{j}\t{p}\t{p + 10}\n")

    # A bwa-written BAM carries no UR tag on @SQ; samtools adds one when it writes the CRAM.
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": n, "LN": len(s)} for n, s in contigs]}
    clip = cons[100:140]
    bam = d / "reads.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=hdr) as bf:
        for t in range(N_TILES):
            pos = t * TILE + 100
            for k in range(N_REP_CONTIGS):
                # Mate inside a repeat interval, near the END of its contig: the last bins are the
                # last ones the lazy build fills, so the window for a false 0 is widest there.
                mpos = ((REP_LEN - 1000 - 40 * (t * N_REP_CONTIGS + k)) // 20) * 20 + 2
                a = pysam.AlignedSegment()
                a.query_name = f"t{t}_m{k}"
                a.flag = 1 | 32 | 64
                a.reference_id = 0
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigartuples = [(4, len(clip)), (0, 60)]
                a.query_sequence = clip + clip_contig[pos:pos + 60]
                a.query_qualities = pysam.qualitystring_to_array("I" * (len(clip) + 60))
                a.next_reference_id = 1 + k
                a.next_reference_start = mpos
                bf.write(a)
    subprocess.run([samtools, "index", str(bam)], check=True)
    cram = d / "reads.cram"
    subprocess.run([samtools, "view", "-C", "-T", str(ref_fa), "-o", str(cram), str(bam)], check=True)
    subprocess.run([samtools, "index", str(cram)], check=True)
    return {"dir": d, "bam": bam, "cram": cram, "ref": ref_fa, "bed": bed, "cons": cons_fa,
            "bwa": bwa}


def _realign_rows(scan_bin, case, inp, out):
    cmd = [scan_bin, "--bam", str(inp), "--reference", str(case["ref"]),
           "--tile-size", str(TILE), "--threads", str(N_TILES),
           "--realign-out", str(out), "--l1-cns", str(case["cons"]),
           "--alu-bed", str(case["bed"]), "--bwa", case["bwa"], "--bwa-threads", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, f"xtea_scan failed:\n{r.stderr[-3000:]}"
    return [ln.split("\t") for ln in out.read_text().splitlines()]


def test_mate_in_rep_count_is_exact_and_identical_for_bam_and_cram(xtea2_scan_bin, parity_case):
    """Every site's mate_in_rep_clip_cnt is 16 in every run, as BAM and as CRAM."""
    expected = None
    for fmt in ("bam", "cram"):
        for i in range(RUNS_PER_FORMAT):
            rows = _realign_rows(xtea2_scan_bin, parity_case, parity_case[fmt],
                                 parity_case["dir"] / f"realign.{fmt}.{i}.tsv")
            assert len(rows) == N_TILES, f"{fmt} run {i}: {len(rows)} sites, want {N_TILES}"
            short = [(r[0], r[1], r[5]) for r in rows if int(r[5]) != N_REP_CONTIGS]
            assert not short, (
                f"{fmt} run {i}: mate_in_rep_clip_cnt below {N_REP_CONTIGS} at {len(short)} "
                f"site(s), e.g. {short[:3]}. Every mate lands in a repeat, so reads were lost: "
                "the shared regidx was queried while another thread was still building it.")
            if expected is None:
                expected = rows
            assert rows == expected, f"{fmt} run {i} differs from the first BAM run"


def test_vcf_reference_line_is_the_run_reference_not_the_container(parity_case):
    """BAM (no UR) and CRAM (UR added by samtools) give the same ##reference line."""
    pysam = pytest.importorskip("pysam")
    from x_gvcf import gVCF
    with pysam.AlignmentFile(str(parity_case["bam"]), "rb") as f:
        assert not any("UR" in sq for sq in f.header.to_dict()["SQ"]), "fixture BAM must lack UR"
    with pysam.AlignmentFile(str(parity_case["cram"]), "rc",
                             reference_filename=str(parity_case["ref"])) as f:
        assert all("UR" in sq for sq in f.header.to_dict()["SQ"]), "samtools should add UR"

    ref = str(parity_case["ref"])
    heads = {fmt: gVCF().gnrt_vcf_head("S", str(parity_case[fmt]), sf_run_ref=ref)
             for fmt in ("bam", "cram")}
    assert heads["bam"] == heads["cram"]
    ref_lines = [ln for ln in heads["bam"].splitlines() if ln.startswith("##reference=")]
    assert ref_lines == ["##reference=ref.fa"]


def test_cram_with_zero_span_index_is_refused(parity_case, tmp_path):
    """A CRAM encoded against a FASTA missing a header contig gets a .crai with span 0 (htslib
    falls back to embed_ref=2), and region queries then drop reads. run_xtea must refuse it, and
    must accept the normal CRAM."""
    pysam = pytest.importorskip("pysam")
    samtools = _tool("samtools", "XTEA2_SAMTOOLS")
    import run_xtea

    # A FASTA holding only chrA: the chrR contigs of the header are missing.
    part_fa = tmp_path / "chrA_only.fa"
    with pysam.FastaFile(str(parity_case["ref"])) as fa, open(part_fa, "w") as f:
        f.write(">chrA\n" + fa.fetch("chrA") + "\n")
    subprocess.run([samtools, "faidx", str(part_fa)], check=True)
    bad = tmp_path / "bad.cram"
    subprocess.run([samtools, "view", "-C", "-T", str(part_fa), "-o", str(bad),
                    str(parity_case["bam"])], check=True, capture_output=True)
    subprocess.run([samtools, "index", str(bad)], check=True)

    n_rows, n_zero = run_xtea.cram_index_zero_span_slices(str(bad) + ".crai")
    assert n_rows > 0 and n_zero == n_rows
    with pytest.raises(SystemExit) as e:
        run_xtea.detect_format_and_index(str(bad), str(part_fa))
    assert "span" in str(e.value) and "embed_ref" in str(e.value)

    assert run_xtea.cram_index_zero_span_slices(str(parity_case["cram"]) + ".crai")[1] == 0
    run_xtea.detect_format_and_index(str(parity_case["cram"]), str(parity_case["ref"]))


def test_empty_source_region_counts_no_reads_for_bam_and_cram(parity_case, tmp_path):
    """The sibling-transduction source check asks for [pos, pos) when both sides' mates cluster
    at one position. htslib returns nothing for that from a BAM but the reads covering pos-1
    from a CRAM, so the CRAM alone could pass the check and gain calls. Both must give 0."""
    pytest.importorskip("pysam")
    from x_orphan_transduction import XOrphanTransduction
    x = XOrphanTransduction(str(tmp_path) + "/", 1, str(parity_case["ref"]))
    x._process_chrm_name = lambda b_with_chr, chrm: chrm   # the toy contig is "chrA", no chr1
    pos = 130                                   # inside the 16 reads at chrA:100-160
    for fmt in ("bam", "cram"):
        f = str(parity_case[fmt])
        assert x._cnt_low_mapq_reads(f, "chrA", pos, pos + 1)[0] == N_REP_CONTIGS, fmt
        assert x._cnt_low_mapq_reads(f, "chrA", pos, pos) == (0, 0), (
            f"{fmt}: an empty interval returned reads")
