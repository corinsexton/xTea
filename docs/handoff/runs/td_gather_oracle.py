"""Byte-compare the legacy Python transduction gather against the native one, on any site list.

WHY THIS EXISTS SEPARATELY FROM THE PYTEST
------------------------------------------
tests/test_td_gather.py covers the ordering rules and runs in seconds, but both of its real-data
fixtures are chr1-only, and the Stage 1 gate is chr1 too. The two things chr1 cannot show are the
two that only appear at whole-genome scale: MULTI-CONTIG emission order (the clip FASTQ follows the
site file's contig order, the disc FASTA the BAM header's -- opposite orders, and a single-contig
run cannot tell them apart) and the ~54K-site memory/time profile. This script does the same
byte-comparison against a real WGS site list, at production thread counts.

It compares the three artifacts the transduction step actually consumes:
  raw_candidate_sites_all_clip.fq   the clip FASTQ  (order-sensitive: bwa tie-breaks ride on it)
  all_disc_pos.disc_pos             the disc anchor stream
  raw_candidate_sites_all_disc.fa   the disc mate FASTA (order-sensitive for the same reason)

usage:
  td_gather_oracle.py --bam BAM --ref FA --scan-bin BIN --out DIR --extnd N --threads N
                      (--sites FILE | --raw-disc FILE --candidates FILE --annotation FILE
                                      [--rtype N] [--ndisc N] [--sva])

--sites takes a ready-made site list. The second form BUILDS the real one the step would use, by
calling re_slct_with_clip_raw_disc_sites -- which matters: it is what drops raw-disc lines close to
a cns candidate, and without that step a hand-made list can contain an exact duplicate position,
which CRASHES the legacy gather (FileNotFoundError out of its merge loop) and proves nothing.
"""
import argparse
import os
import subprocess
import sys
import time

# docs/handoff/runs/<this> -> repo root: four levels up, then xtea/.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(REPO_ROOT, "xtea"))

import global_values                              # noqa: E402
from x_clip_disc_filter import XClipDiscFilter     # noqa: E402
from x_transduction import XTransduction           # noqa: E402

BIN_SIZE = 50_000_000     # x_TEA_main's bin_size for the disc gather
MIN_COPY_LEN = 225        # x_TEA_main's i_min_copy_len for the annotation tree


def build_site_list(a, work):
    """The real sf_raw_slct: what x_TEA_main hands call_candidate_transduction_v3."""
    if a.sva:
        global_values.IS_CALL_SVA = True
    xtd = XTransduction(work, a.threads, a.ref)
    xann = xtd.prep_annotation_interval_tree(a.annotation, MIN_COPY_LEN)
    out = os.path.join(work, "raw_slct.txt")
    xtd.re_slct_with_clip_raw_disc_sites(a.raw_disc, a.candidates, a.ndisc, xann,
                                         a.rtype, False, out)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--scan-bin", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--extnd", type=int, required=True)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sites")
    p.add_argument("--raw-disc")
    p.add_argument("--candidates")
    p.add_argument("--annotation")
    p.add_argument("--rtype", type=int, default=1)
    p.add_argument("--ndisc", type=int, default=5)
    p.add_argument("--sva", action="store_true")
    p.add_argument("--skip-legacy", action="store_true",
                   help="run only the native gather (for a determinism repeat)")
    a = p.parse_args()

    py_dir = os.path.join(a.out, "py") + "/"
    cpp_dir = os.path.join(a.out, "cpp") + "/"
    for d in (py_dir, cpp_dir):
        os.makedirs(d, exist_ok=True)

    sites = a.sites
    if not sites:
        for need in ("raw_disc", "candidates", "annotation"):
            if not getattr(a, need):
                p.error("--sites, or all of --raw-disc/--candidates/--annotation")
        t0 = time.time()
        sites = build_site_list(a, a.out)
        print(f"[sites] built {sites}: {sum(1 for _ in open(sites))} sites "
              f"in {time.time()-t0:.1f}s", flush=True)
    n_sites = sum(1 for _ in open(sites))
    n_contigs = len({ln.split()[0] for ln in open(sites) if ln.strip()})
    print(f"[sites] {n_sites} sites over {n_contigs} contigs  (extnd={a.extnd}, "
          f"threads={a.threads})", flush=True)

    bam_list = os.path.join(a.out, "bam_list.txt")
    with open(bam_list, "w") as f:
        f.write(a.bam + "\n")

    t_py_clip = t_py_disc = float("nan")
    if not a.skip_legacy:
        xf = XClipDiscFilter(bam_list, py_dir, a.threads, a.ref)
        t0 = time.time()
        xf.collect_clipped_reads(sites, a.extnd, py_dir + "clip.fq")
        t1 = time.time()
        xf.collect_disc_reads(sites, a.extnd, BIN_SIZE, py_dir + "disc.fa")
        t2 = time.time()
        t_py_clip, t_py_disc = t1 - t0, t2 - t1
        print(f"[legacy] clip_gather {t_py_clip:.1f}s  disc_gather {t_py_disc:.1f}s  "
              f"total {t_py_clip + t_py_disc:.1f}s", flush=True)

    cmd = [a.scan_bin, "--bam", a.bam, "--reference", a.ref, "--td-sites", sites,
           "--td-extnd", str(a.extnd), "--td-clip-fastq", cpp_dir + "clip.fq",
           "--td-disc-fa", cpp_dir + "disc.fa",
           "--td-disc-pos", cpp_dir + "all_disc_pos.disc_pos",
           "--threads", str(a.threads)]
    print("[native] " + " ".join(cmd), flush=True)
    t0 = time.time()
    subprocess.run(cmd, check=True)
    t_cpp = time.time() - t0
    print(f"[native] total {t_cpp:.1f}s", flush=True)
    if not a.skip_legacy and t_cpp > 0:
        print(f"[speedup] {(t_py_clip + t_py_disc) / t_cpp:.2f}x on the gather", flush=True)

    if a.skip_legacy:
        return 0
    rc = 0
    for name, x, y in [("clip.fq", py_dir + "clip.fq", cpp_dir + "clip.fq"),
                       ("all_disc_pos", py_dir + "all_disc_pos" + global_values.DISC_POS_SUFFIX,
                        cpp_dir + "all_disc_pos.disc_pos"),
                       ("disc.fa", py_dir + "disc.fa", cpp_dir + "disc.fa")]:
        r = subprocess.run(["cmp", x, y], capture_output=True, text=True)
        sx = os.path.getsize(x) if os.path.exists(x) else -1
        sy = os.path.getsize(y) if os.path.exists(y) else -1
        print(f"{'IDENTICAL' if r.returncode == 0 else 'DIFFER   '}  {name:14s} "
              f"py={sx} cpp={sy} {r.stdout.strip()}{r.stderr.strip()}", flush=True)
        rc |= (r.returncode != 0)
    return rc


if __name__ == "__main__":
    sys.exit(main())
