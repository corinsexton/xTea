"""Run postF_mosaic (counting_v2 code copy, B15 inputs) with disc_upper opened to 1.0, so the
ceiling can be applied offline. usage: postf_open.py NAME OUTDIR"""
import os, shlex, subprocess, sys
C = "/n/scratch/users/c/cos689/xtea2_review/counting_v2/mosaic"
CODE = C + "/code/xtea"
name, out = sys.argv[1], sys.argv[2]
root = f"/n/scratch/users/c/cos689/xtea2_testdata/vafmix/runs/{name}/{name}"
for fam in ("L1", "Alu", "SVA"):
    lines = open(f"{root}/{name}_{fam}.steps.log").read().split("\n")
    argv = shlex.split(lines[lines.index(f"=== STEP {fam}:postF_mosaic ===") + 1])[2:]
    os.makedirs(f"{out}/{name}/{fam}", exist_ok=True)
    o = f"{out}/{name}/{fam}/mosaic_open.txt"
    argv[argv.index("-o") + 1] = o
    argv[argv.index("-i") + 1] = f"{C}/inputs/{name}_{fam}.input_cv2.txt"
    env = dict(os.environ, XTEA2_TUNE_DISC_UPPER="1.0", XTEA2_TUNE_AFPE="0.14", XTEA2_TUNE_TD="0.148")
    subprocess.run([sys.executable, CODE + "/x_TEA_main.py"] + argv, cwd=CODE, env=env, check=True,
                   stdout=open(o + ".stdout", "w"), stderr=subprocess.STDOUT)
    print(name, fam, sum(1 for _ in open(o)))
