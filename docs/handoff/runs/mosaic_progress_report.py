#!/usr/bin/env python3
"""
Regenerate the mosaic-fleet progress report (disk footprint, per-arm step
Gantt, realign sub-step timing, memory) from live SLURM state + each arm's
xtea2 log and disk_usage.tsv sampler. Safe to rerun any time; it only reads.

    python3 docs/handoff/runs/mosaic_progress_report.py
    python3 docs/handoff/runs/mosaic_progress_report.py --out /path/to/report.html

Discovers the current job id per arm via squeue (falls back to reading
whatever files exist if an arm isn't currently queued/running -- e.g. it
finished or hasn't been resubmitted yet). Prints a short text summary and
writes a self-contained HTML report next to the run outputs.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time

ARMS = ["bcm", "broad", "nygc", "uwsc", "washu"]
OUT_BASE = "/n/scratch/users/c/cos689/xtea2_testdata"
DEFAULT_REPORT_PATH = f"{OUT_BASE}/mosaic_progress_report.html"

COLOR_VAR = {
    "bcm": "--series-bcm",
    "broad": "--series-broad",
    "nygc": "--series-nygc",
    "uwsc": "--series-uwsc",
    "washu": "--series-washu",
}

STEP_COLOR_VARS = [f"--step-{i}" for i in range(1, 9)]


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout


def parse_slurm_elapsed(s):
    s = s.strip()
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    parts = [int(p) for p in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    h, m, sec = parts[-3], parts[-2], parts[-1]
    return days * 86400 + h * 3600 + m * 60 + sec


def discover_jobs():
    """arm -> {jobid, elapsed_s} for the most recent squeue row matching that job name."""
    out = sh(f"squeue -u $(whoami) -o '%i|%j|%M' -h 2>/dev/null")
    by_arm = {}
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) != 3:
            continue
        jobid, name, elapsed = parts
        name = name.strip()
        if name not in ARMS:
            continue
        jobid = jobid.strip()
        try:
            jobid_i = int(jobid)
        except ValueError:
            continue
        prev = by_arm.get(name)
        if prev is None or jobid_i > prev["jobid"]:
            by_arm[name] = {"jobid": jobid_i, "elapsed_s": parse_slurm_elapsed(elapsed)}
    return by_arm


def sstat_max_rss_gb(jobid):
    out = sh(f"sstat -j {jobid} --format=MaxRSS -a --noheader -P 2>/dev/null")
    best = None
    for line in out.splitlines():
        line = line.strip()
        m = re.match(r"([\d.]+)K", line)
        if not m:
            continue
        kb = float(m.group(1))
        if best is None or kb > best:
            best = kb
    # sstat's MaxRSS is KiB and --mem=256G is 256 GiB, so this MUST divide by 1024^2, not 1e6.
    # Dividing by 1e6 inflates every figure by 7.4% and silently misreports headroom against
    # the cap: it showed uwsc at "233.8G" against 256G (23G to spare) when the true figure was
    # 222.9 GiB (33 GiB to spare). It also mis-told the story of the 192G round -- uwsc's OOM
    # reads 199.4 in the bad units but is 190.2 GiB, i.e. the kernel killed it just UNDER a
    # 192 GiB cap, which is textbook, not the "~4% overshoot" recorded in a7d0def's message.
    return round(best / 1048576, 1) if best is not None else None


def sample_dir(arm):
    return f"{OUT_BASE}/mosaic_hapmap_wgs_{arm}/hapmap_{arm}_wgs"


def out_dir(arm):
    return f"{OUT_BASE}/mosaic_hapmap_wgs_{arm}"


def slurm_stdout_path(arm, jobid):
    """The [xtea2 + Ts] phase/step timing lines are printed to SLURM's stdout
    file (mosaic_wgs_<arm>_<jobid>.out), NOT the internal xtea2 log -- that
    log only has the raw subprocess passthrough (realign, bwa_sam, htslib)."""
    if jobid is not None:
        p = f"{OUT_BASE}/mosaic_wgs_{arm}_{jobid}.out"
        if os.path.exists(p):
            return p
    candidates = glob.glob(f"{OUT_BASE}/mosaic_wgs_{arm}_*.out")
    return max(candidates, key=os.path.getmtime) if candidates else None


def parse_disk_usage(arm):
    path = f"{out_dir(arm)}/disk_usage.tsv"
    rows = []
    try:
        with open(path) as f:
            next(f, None)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                t, total_bytes, _h, n_files = parts
                rows.append({"t": int(t), "gb": round(int(total_bytes) / 1e9, 3), "files": int(n_files)})
    except FileNotFoundError:
        pass
    return rows


BOUNDARY_RE = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*(phase|step) (\S+): ([\d.]+)s")
RUNNING_RE = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*step: (\S+)$")
SECTION_RE = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*=== (.+?) ===\s*$")
CANDLIST_RE = re.compile(r"^\[xtea2 \+\s*([\d.]+)s\]\s*(\w+):(clip|disc) \(xtea-scan\)")
REALIGN_RE = re.compile(r"^\[realign\]\s+(\S+)\s+t=\d+\s+total=([\d.]+)s")


def parse_realign(arm):
    """[realign] family totals are only in the internal xtea2 log (raw subprocess passthrough)."""
    path = f"{sample_dir(arm)}/hapmap_{arm}_wgs_xtea2.log"
    realign = {}
    try:
        with open(path) as f:
            for line in f:
                m = REALIGN_RE.match(line.rstrip("\n"))
                if m:
                    realign[m.group(1)] = float(m.group(2))
    except FileNotFoundError:
        pass
    return realign


def parse_boundaries(arm, jobid, now_elapsed):
    """Gantt segments for one arm from its SLURM stdout file's [xtea2 + Ts] lines."""
    path = slurm_stdout_path(arm, jobid)
    boundaries = []  # (t, label)
    activity_t, activity_name = None, None  # best guess at "what's running now"
    if path:
        try:
            with open(path, errors="replace") as f:
                for line in f:
                    line = line.rstrip("\n")
                    m = BOUNDARY_RE.match(line)
                    if m:
                        t = float(m.group(1))
                        kind, name = m.group(2), m.group(3)
                        boundaries.append((t, f"{kind}:{name}"))
                        continue
                    m = CANDLIST_RE.match(line)
                    if m:
                        t = float(m.group(1))
                        fam, kind = m.group(2), m.group(3)
                        boundaries.append((t, f"{fam}:{kind}"))
                        continue
                    m = RUNNING_RE.match(line)
                    if m:
                        t = float(m.group(1))
                        if activity_t is None or t > activity_t:
                            activity_t, activity_name = t, m.group(2)
                        continue
                    m = SECTION_RE.match(line)
                    if m:
                        t = float(m.group(1))
                        if activity_t is None or t > activity_t:
                            activity_t, activity_name = t, m.group(2)
                        continue
        except FileNotFoundError:
            pass

    boundaries.sort(key=lambda x: x[0])
    # dedupe consecutive identical timestamps, keep the later label
    dedup = []
    for t, label in boundaries:
        if dedup and abs(dedup[-1][0] - t) < 1e-6:
            dedup[-1] = (t, label)
        else:
            dedup.append((t, label))
    boundaries = dedup

    segments = []
    prev_t = 0.0
    for t, label in boundaries:
        if t <= prev_t:
            continue
        segments.append({"name": label, "start": prev_t, "end": t, "dur": t - prev_t, "running": False})
        prev_t = t

    # the currently-running step: best activity marker at/after the last completed boundary,
    # else fall back to "scan_pass" (true before the very first phase line ever appears)
    if now_elapsed > prev_t:
        if activity_t is not None and activity_t >= prev_t - 1e-6:
            label = activity_name
        else:
            label = "scan_pass"
        segments.append({"name": label, "start": prev_t, "end": now_elapsed,
                          "dur": now_elapsed - prev_t, "running": True})

    return segments


def fmt_dur(s):
    if s is None:
        return "n/a"
    if s < 60:
        return f"{s:.0f}s"
    m = s / 60
    if m < 60:
        return f"{m:.1f}m"
    return f"{m/3600*60:.2f}h" if False else f"{m/60:.2f}h"


def build_report(jobs):
    now = time.time()
    disk, gantt, mem, realign_families = {}, {}, {}, set()
    realign_by_arm = {}
    meta = {"generated": time.strftime("%Y-%m-%d %H:%M", time.localtime(now)), "arms": {}}

    for arm in ARMS:
        job = jobs.get(arm)
        running = job is not None
        d_rows = parse_disk_usage(arm)
        if running:
            elapsed = job["elapsed_s"]
        else:
            elapsed = d_rows[-1]["t"] if d_rows else 0

        jobid = job["jobid"] if running else None
        segs = parse_boundaries(arm, jobid, elapsed)
        realign = parse_realign(arm)
        disk[arm] = d_rows
        gantt[arm] = segs
        realign_by_arm[arm] = realign
        realign_families.update(realign.keys())
        mem[arm] = sstat_max_rss_gb(jobid) if running else None
        meta["arms"][arm] = {
            "running": running,
            "jobid": job["jobid"] if running else None,
            "elapsed_s": elapsed,
            "elapsed_h": round(elapsed / 3600, 2),
            "current_step": segs[-1]["name"] if segs else "n/a",
        }

    families = sorted(realign_families) or ["LINE1", "ALU", "SVA"]
    realign = {fam: {arm: realign_by_arm.get(arm, {}).get(fam) for arm in ARMS} for fam in families}

    return {"meta": meta, "disk": disk, "gantt": gantt, "mem": mem, "realign": realign}


def print_summary(report):
    print(f"Mosaic fleet progress -- generated {report['meta']['generated']}")
    for arm in ARMS:
        a = report["meta"]["arms"][arm]
        m = report["mem"][arm]
        status = f"job {a['jobid']}" if a["running"] else "NOT RUNNING"
        mem_s = f"{m:.1f}G" if m is not None else "n/a"
        print(f"  {arm:6s} {status:14s} elapsed {a['elapsed_h']:.2f}h  mem {mem_s:8s}  step: {a['current_step']}")


HTML_TEMPLATE = r"""<title>Mosaic Fleet Progress</title>
<style>
  .viz-root {
    color-scheme: light;
    --surface-1:      #fcfcfb;
    --page:           #f9f9f7;
    --text-primary:   #0b0b0b;
    --text-secondary: #52514e;
    --text-muted:     #898781;
    --grid:           #e1e0d9;
    --axis:           #c3c2b7;
    --border:         rgba(11,11,11,0.10);
    --series-bcm:     #2a78d6;
    --series-broad:   #eb6834;
    --series-nygc:    #1baf7a;
    --series-uwsc:    #eda100;
    --series-washu:   #e87ba4;
    --status-critical:#d03b3b;
    --status-warning: #fab219;
    --step-1: #86b6ef; --step-2: #6da7ec; --step-3: #5598e7; --step-4: #3987e5;
    --step-5: #2a78d6; --step-6: #256abf; --step-7: #1c5cab; --step-8: #184f95;
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) .viz-root {
      color-scheme: dark;
      --surface-1: #1a1a19; --page: #0d0d0d; --text-primary: #ffffff;
      --text-secondary: #c3c2b7; --text-muted: #898781; --grid: #2c2c2a; --axis: #383835;
      --border: rgba(255,255,255,0.10);
      --series-bcm: #3987e5; --series-broad: #d95926; --series-nygc: #199e70;
      --series-uwsc: #c98500; --series-washu: #d55181;
      --status-critical: #e66767; --status-warning: #fab219;
    }
  }
  :root[data-theme="dark"] .viz-root {
    color-scheme: dark;
    --surface-1: #1a1a19; --page: #0d0d0d; --text-primary: #ffffff;
    --text-secondary: #c3c2b7; --text-muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --border: rgba(255,255,255,0.10);
    --series-bcm: #3987e5; --series-broad: #d95926; --series-nygc: #199e70;
    --series-uwsc: #c98500; --series-washu: #d55181;
    --status-critical: #e66767; --status-warning: #fab219;
  }
  * { box-sizing: border-box; }
  body { margin: 0; }
  .viz-root { background: var(--page); color: var(--text-primary); font-family: system-ui, -apple-system, "Segoe UI", sans-serif; padding: 28px 24px 40px; min-height: 100vh; }
  .wrap { max-width: 980px; margin: 0 auto; }
  h1 { font-size: 20px; margin: 0 0 4px; letter-spacing: -0.01em; }
  .subtitle { color: var(--text-secondary); font-size: 13px; margin: 0 0 22px; }
  .card { background: var(--surface-1); border: 1px solid var(--border); border-radius: 12px; padding: 20px 20px 12px; margin-bottom: 20px; }
  .card h2 { font-size: 14px; margin: 0 0 2px; font-weight: 600; }
  .card .desc { color: var(--text-muted); font-size: 12px; margin: 0 0 12px; }
  .legend { display: flex; gap: 16px; flex-wrap: wrap; margin: 4px 0 14px; }
  .legend-item { display: flex; align-items: center; gap: 6px; font-size: 12px; color: var(--text-secondary); }
  .legend-swatch { width: 10px; height: 10px; border-radius: 2px; flex: none; }
  svg { display: block; width: 100%; height: auto; overflow: visible; }
  .axis-label { fill: var(--text-muted); font-size: 10.5px; }
  .grid-line { stroke: var(--grid); stroke-width: 1; }
  .axis-line { stroke: var(--axis); stroke-width: 1; }
  .series-line { fill: none; stroke-width: 2; }
  .end-label { font-size: 11px; font-weight: 600; }
  .hover-dot { r: 4; stroke: var(--surface-1); stroke-width: 2; }
  .crosshair { stroke: var(--text-muted); stroke-width: 1; stroke-dasharray: 3 3; opacity: 0; }
  .tooltip { position: absolute; pointer-events: none; opacity: 0; background: var(--surface-1); border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px; font-size: 11.5px; box-shadow: 0 4px 16px rgba(0,0,0,0.18); transition: opacity 0.08s; min-width: 130px; z-index: 5; }
  .tooltip-time { color: var(--text-muted); font-size: 10.5px; margin-bottom: 4px; }
  .tooltip-row { display: flex; justify-content: space-between; gap: 10px; padding: 1px 0; }
  .tooltip-row .v { font-weight: 600; font-variant-numeric: tabular-nums; }
  .chart-holder { position: relative; }
  .bar-row { display: flex; align-items: center; gap: 10px; margin: 10px 0; }
  .bar-label { width: 58px; font-size: 12px; color: var(--text-secondary); flex: none; }
  .bar-track { position: relative; flex: 1; height: 22px; background: var(--grid); border-radius: 5px; overflow: visible; }
  .bar-fill { position: absolute; top: 0; left: 0; height: 100%; border-radius: 5px; }
  .bar-value { font-size: 12px; font-variant-numeric: tabular-nums; color: var(--text-primary); width: 70px; text-align: right; flex: none; }
  .ref-line { position: absolute; top: -4px; bottom: -4px; width: 0; border-left: 2px dashed var(--status-warning); }
  .ref-line.crit { border-left: 2px dashed var(--status-critical); }
  .ref-label { position: absolute; top: -18px; font-size: 10px; color: var(--text-muted); white-space: nowrap; transform: translateX(-50%); }
  .foot-table { width: 100%; border-collapse: collapse; font-size: 12.5px; margin-top: 4px; }
  .foot-table th { text-align: left; color: var(--text-muted); font-weight: 500; font-size: 11px; padding: 4px 8px; border-bottom: 1px solid var(--grid); }
  .foot-table td { padding: 5px 8px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
  .foot-table tr:last-child td { border-bottom: none; }
  .dot { display: inline-block; width: 8px; height: 8px; border-radius: 2px; margin-right: 6px; vertical-align: -1px; }
  .gantt-row { display: flex; align-items: center; gap: 10px; margin: 8px 0; }
  .gantt-label { width: 58px; font-size: 12px; color: var(--text-secondary); flex: none; }
  .gantt-track { position: relative; flex: 1; height: 26px; border-radius: 6px; overflow: hidden; background: var(--grid); }
  .gantt-seg { position: absolute; top: 0; bottom: 0; border-right: 1px solid var(--surface-1); }
  .gantt-seg.running { background-image: repeating-linear-gradient(45deg, rgba(255,255,255,0.28) 0 4px, transparent 4px 9px); }
  .gantt-not-running { font-size: 11px; color: var(--text-muted); font-style: italic; }
  .rbar-group { margin-bottom: 16px; }
  .rbar-group-title { font-size: 12px; font-weight: 600; color: var(--text-secondary); margin: 0 0 6px; }
  .rbar-row { display: flex; align-items: center; gap: 10px; margin: 6px 0; }
  .rbar-label { width: 58px; font-size: 12px; color: var(--text-secondary); flex: none; }
  .rbar-track { position: relative; flex: 1; height: 18px; background: var(--grid); border-radius: 5px; }
  .rbar-fill { position: absolute; top: 0; left: 0; height: 100%; border-radius: 5px; }
  .rbar-fill.pending { background-image: repeating-linear-gradient(45deg, var(--axis) 0 4px, transparent 4px 9px); background-color: transparent; opacity: 0.6; }
  .rbar-value { font-size: 12px; font-variant-numeric: tabular-nums; color: var(--text-primary); width: 78px; text-align: right; flex: none; }
</style>

<div class="viz-root">
  <div class="wrap">
    <h1>Mosaic Fleet Progress</h1>
    <p class="subtitle" id="subtitle"></p>

    <div class="card">
      <h2>Scratch disk footprint</h2>
      <p class="desc">Total bytes under each arm's output dir, sampled every 60s. Jumps mark phase transitions.</p>
      <div class="legend" id="legend1"></div>
      <div class="chart-holder">
        <svg id="chart1" viewBox="0 0 900 340" preserveAspectRatio="xMidYMid meet"></svg>
        <div class="tooltip" id="tooltip1"></div>
      </div>
    </div>

    <div class="card">
      <h2>Step durations, by arm</h2>
      <p class="desc">Each row is one arm's own timeline since job start, tiled from its xtea2 log's phase/step boundaries. Color darkens with that arm's own pipeline order (not comparable hue-to-hue across arms); hatched segment is still running. Hover a segment for its name and length.</p>
      <div id="ganttRows"></div>
      <div class="tooltip" id="tooltipG"></div>
    </div>

    <div class="card">
      <h2>Realign sub-step time, by arm</h2>
      <p class="desc">Wall time each arm's scan pass has spent realigning clipped/discordant reads against each repeat consensus. Hatched = not reached yet in that arm's log.</p>
      <div id="realignBars"></div>
    </div>

    <div class="card">
      <h2>Current memory (RSS)</h2>
      <p class="desc">Snapshot via <code>sstat</code> at report time, against the 256G cap and the ~199G point where two arms OOM'd under the old 192G cap.</p>
      <div id="bars"></div>
    </div>

    <div class="card">
      <h2>Status</h2>
      <table class="foot-table">
        <thead><tr><th>arm</th><th>job</th><th>elapsed</th><th>current step</th><th>disk</th><th>files</th><th>mem</th></tr></thead>
        <tbody id="statusRows"></tbody>
      </table>
    </div>
  </div>
</div>

<script>
const ARMS = __ARMS_JSON__;
const COLORS = __COLORS_JSON__;
const META = __META_JSON__;
const DISK = __DISK_JSON__;
const GANTT = __GANTT_JSON__;
const MEM = __MEM_JSON__;
const REALIGN = __REALIGN_JSON__;

document.getElementById("subtitle").textContent =
  `WGS mosaic depth-ladder — ${ARMS.join(", ")} · generated ${META.generated}`;

function fmtDur(s) {
  if (s == null) return "n/a";
  if (s < 60) return s.toFixed(0) + "s";
  const m = s / 60;
  if (m < 60) return m.toFixed(1) + "m";
  return (m / 60).toFixed(2) + "h";
}

function downsample(rows, n) {
  if (rows.length <= n) return rows;
  const step = rows.length / n;
  const out = [];
  for (let i = 0; i < n; i++) out.push(rows[Math.floor(i * step)]);
  out.push(rows[rows.length - 1]);
  return out;
}
const diskData = {};
ARMS.forEach(a => diskData[a] = downsample(DISK[a] || [], 110));

// ---------- Legend ----------
const legend1 = document.getElementById("legend1");
ARMS.forEach(a => {
  const el = document.createElement("div");
  el.className = "legend-item";
  el.innerHTML = `<span class="legend-swatch" style="background:${COLORS[a]}"></span>${a}`;
  legend1.appendChild(el);
});

// ---------- Line chart ----------
const svgNS = "http://www.w3.org/2000/svg";
const svg = document.getElementById("chart1");
const W = 900, H = 340;
const M = { top: 14, right: 54, bottom: 32, left: 42 };
const plotW = W - M.left - M.right;
const plotH = H - M.top - M.bottom;
const allT = ARMS.flatMap(a => diskData[a].map(d => d.t));
const allGB = ARMS.flatMap(a => diskData[a].map(d => d.gb));
const maxT = Math.max(1, ...allT);
const maxGB = Math.max(0.001, ...allGB);
const yMax = Math.ceil((maxGB * 1.12) / 5) * 5 || 5;
const xMax = maxT;
function xPos(t) { return M.left + (t / xMax) * plotW; }
function yPos(gb) { return M.top + plotH - (gb / yMax) * plotH; }
function el(tag, attrs) {
  const e = document.createElementNS(svgNS, tag);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  return e;
}
const yTicks = 5;
for (let i = 0; i <= yTicks; i++) {
  const gb = (yMax / yTicks) * i;
  const y = yPos(gb);
  svg.appendChild(el("line", { class: "grid-line", x1: M.left, x2: W - M.right, y1: y, y2: y }));
  const t = el("text", { class: "axis-label", x: M.left - 8, y: y + 3, "text-anchor": "end" });
  t.textContent = gb.toFixed(0) + "G";
  svg.appendChild(t);
}
const hoursMax = xMax / 3600;
const xTickCount = Math.max(2, Math.round(hoursMax * 2));
for (let i = 0; i <= xTickCount; i++) {
  const hr = (hoursMax / xTickCount) * i;
  const x = xPos(hr * 3600);
  const t = el("text", { class: "axis-label", x: x, y: H - M.bottom + 16, "text-anchor": "middle" });
  t.textContent = hr.toFixed(1) + "h";
  svg.appendChild(t);
}
svg.appendChild(el("line", { class: "axis-line", x1: M.left, x2: W - M.right, y1: H - M.bottom, y2: H - M.bottom }));

ARMS.forEach(a => {
  const rows = diskData[a];
  if (!rows.length) return;
  const d = rows.map((r, i) => (i === 0 ? "M" : "L") + xPos(r.t).toFixed(2) + "," + yPos(r.gb).toFixed(2)).join(" ");
  svg.appendChild(el("path", { class: "series-line", d, stroke: COLORS[a] }));
  const last = rows[rows.length - 1];
  const lx = xPos(last.t), ly = yPos(last.gb);
  svg.appendChild(el("circle", { cx: lx, cy: ly, r: 3, fill: COLORS[a] }));
  const label = el("text", { class: "end-label", x: lx + 7, y: ly + 3.5, fill: COLORS[a] });
  label.textContent = a;
  svg.appendChild(label);
});

const crosshair = el("line", { class: "crosshair", y1: M.top, y2: H - M.bottom });
svg.appendChild(crosshair);
const hoverDots = {};
ARMS.forEach(a => {
  const c = el("circle", { class: "hover-dot", fill: COLORS[a], opacity: 0 });
  svg.appendChild(c);
  hoverDots[a] = c;
});
const tooltip1 = document.getElementById("tooltip1");
const overlay = el("rect", { x: M.left, y: M.top, width: plotW, height: plotH, fill: "transparent" });
svg.appendChild(overlay);
function nearestIdx(rows, t) {
  let lo = 0, hi = rows.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (rows[mid].t < t) lo = mid + 1; else hi = mid;
  }
  if (lo > 0 && Math.abs(rows[lo - 1].t - t) < Math.abs(rows[lo].t - t)) return lo - 1;
  return lo;
}
overlay.addEventListener("mousemove", (ev) => {
  const rect = svg.getBoundingClientRect();
  const scale = W / rect.width;
  const mx = (ev.clientX - rect.left) * scale;
  const t = ((mx - M.left) / plotW) * xMax;
  if (t < 0 || t > xMax) return;
  const x = xPos(t);
  crosshair.setAttribute("x1", x);
  crosshair.setAttribute("x2", x);
  crosshair.style.opacity = 1;
  let rows_ = [];
  ARMS.forEach(a => {
    const rows = diskData[a];
    if (!rows.length) return;
    const i = nearestIdx(rows, t);
    const pt = rows[i];
    hoverDots[a].setAttribute("cx", xPos(pt.t));
    hoverDots[a].setAttribute("cy", yPos(pt.gb));
    hoverDots[a].setAttribute("opacity", 1);
    rows_.push({ label: a, color: COLORS[a], gb: pt.gb, t: pt.t });
  });
  rows_.sort((x, y) => y.gb - x.gb);
  const hrs = rows_[0].t / 3600;
  tooltip1.innerHTML = `<div class="tooltip-time">t + ${hrs.toFixed(2)}h</div>` +
    rows_.map(r => `<div class="tooltip-row"><span><span class="dot" style="background:${r.color}"></span>${r.label}</span><span class="v">${r.gb.toFixed(2)} GB</span></div>`).join("");
  const wrapRect = svg.parentElement.getBoundingClientRect();
  const px = (x / W) * rect.width;
  tooltip1.style.left = Math.min(px + 12, wrapRect.width - 150) + "px";
  tooltip1.style.top = "6px";
  tooltip1.style.opacity = 1;
});
overlay.addEventListener("mouseleave", () => {
  crosshair.style.opacity = 0;
  tooltip1.style.opacity = 0;
  Object.values(hoverDots).forEach(d => d.setAttribute("opacity", 0));
});

// ---------- Gantt (one row per arm) ----------
const STEP_VARS = ["--step-1","--step-2","--step-3","--step-4","--step-5","--step-6","--step-7","--step-8"];
const ganttRows = document.getElementById("ganttRows");
const tooltipG = document.getElementById("tooltipG");
const ganttMax = Math.max(1, ...ARMS.flatMap(a => (GANTT[a] || []).map(s => s.end)));
ARMS.forEach(a => {
  const segs = GANTT[a] || [];
  const row = document.createElement("div");
  row.className = "gantt-row";
  const track = document.createElement("div");
  track.className = "gantt-track";
  if (!segs.length) {
    row.innerHTML = `<div class="gantt-label">${a}</div>`;
    const holder = document.createElement("div");
    holder.style.flex = "1";
    holder.innerHTML = `<span class="gantt-not-running">no data yet</span>`;
    row.appendChild(holder);
    ganttRows.appendChild(row);
    return;
  }
  segs.forEach((seg, i) => {
    const leftPct = (seg.start / ganttMax) * 100;
    const widthPct = (seg.dur / ganttMax) * 100;
    const el2 = document.createElement("div");
    el2.className = "gantt-seg" + (seg.running ? " running" : "");
    el2.style.left = leftPct + "%";
    el2.style.width = widthPct + "%";
    el2.style.background = `var(${STEP_VARS[i % STEP_VARS.length]})`;
    el2.addEventListener("mousemove", (ev) => {
      const wrapRect = ganttRows.getBoundingClientRect();
      tooltipG.innerHTML = `<div class="tooltip-time">${a}</div><div class="tooltip-row"><span>${seg.name}${seg.running ? " (running)" : ""}</span><span class="v">${fmtDur(seg.dur)}</span></div>`;
      tooltipG.style.left = Math.min(ev.clientX - wrapRect.left + 12, wrapRect.width - 160) + "px";
      tooltipG.style.top = (ev.clientY - wrapRect.top - 30) + "px";
      tooltipG.style.opacity = 1;
    });
    el2.addEventListener("mouseleave", () => { tooltipG.style.opacity = 0; });
    track.appendChild(el2);
  });
  const label = document.createElement("div");
  label.className = "gantt-label";
  label.textContent = a;
  row.appendChild(label);
  row.appendChild(track);
  ganttRows.appendChild(row);
});

// ---------- Realign sub-step bars ----------
const realignBars = document.getElementById("realignBars");
const families = Object.keys(REALIGN);
const realignMax = Math.max(1, ...families.flatMap(f => ARMS.map(a => REALIGN[f][a]).filter(v => v != null)));
families.forEach(fam => {
  const group = document.createElement("div");
  group.className = "rbar-group";
  const title = document.createElement("div");
  title.className = "rbar-group-title";
  title.textContent = fam;
  group.appendChild(title);
  ARMS.forEach(a => {
    const val = REALIGN[fam][a];
    const row = document.createElement("div");
    row.className = "rbar-row";
    if (val == null) {
      row.innerHTML = `<div class="rbar-label">${a}</div>
        <div class="rbar-track"><div class="rbar-fill pending" style="width:100%"></div></div>
        <div class="rbar-value">not reached</div>`;
    } else {
      const pct = (val / realignMax) * 100;
      row.innerHTML = `<div class="rbar-label">${a}</div>
        <div class="rbar-track"><div class="rbar-fill" style="width:${pct}%;background:${COLORS[a]}"></div></div>
        <div class="rbar-value">${fmtDur(val)}</div>`;
    }
    group.appendChild(row);
  });
  realignBars.appendChild(group);
});

// ---------- Memory bars ----------
const bars = document.getElementById("bars");
const CAP = 256, PRIOR_OOM = 199.4;
const barMax = CAP * 1.08;
ARMS.slice().sort((a, b) => (MEM[b] || 0) - (MEM[a] || 0)).forEach(a => {
  const val = MEM[a];
  const row = document.createElement("div");
  row.className = "bar-row";
  if (val == null) {
    row.innerHTML = `<div class="bar-label">${a}</div>
      <div class="bar-track"></div>
      <div class="bar-value">not running</div>`;
  } else {
    const pct = (val / barMax) * 100;
    row.innerHTML = `<div class="bar-label">${a}</div>
      <div class="bar-track">
        <div class="bar-fill" style="width:${pct}%;background:${COLORS[a]}"></div>
        <div class="ref-line" style="left:${(PRIOR_OOM / barMax) * 100}%"></div>
        <div class="ref-line crit" style="left:${(CAP / barMax) * 100}%"></div>
      </div>
      <div class="bar-value">${val.toFixed(1)}G</div>`;
  }
  bars.appendChild(row);
});
const refWrap = document.createElement("div");
refWrap.style.position = "relative";
refWrap.style.height = "16px";
refWrap.style.marginLeft = "68px";
refWrap.style.marginRight = "62px";
refWrap.innerHTML = `
  <div class="ref-label" style="left:${(PRIOR_OOM / barMax) * 100}%">prior OOM 199G</div>
  <div class="ref-label" style="left:${(CAP / barMax) * 100}%">cap 256G</div>
`;
bars.appendChild(refWrap);

// ---------- Status table ----------
const statusRows = document.getElementById("statusRows");
ARMS.forEach(a => {
  const m = META.arms[a];
  const rows = diskData[a];
  const last = rows.length ? rows[rows.length - 1] : { gb: 0, files: 0 };
  const mem = MEM[a];
  const tr = document.createElement("tr");
  tr.innerHTML = `<td><span class="dot" style="background:${COLORS[a]}"></span>${a}</td>
    <td>${m.running ? "job " + m.jobid : "not running"}</td>
    <td>${m.elapsed_h.toFixed(2)}h</td>
    <td>${m.current_step}</td>
    <td>${last.gb.toFixed(1)} GB</td>
    <td>${last.files}</td>
    <td>${mem != null ? mem.toFixed(1) + " G" : "n/a"}</td>`;
  statusRows.appendChild(tr);
});
</script>
"""


def render_html(report):
    html = HTML_TEMPLATE
    html = html.replace("__ARMS_JSON__", json.dumps(ARMS))
    html = html.replace("__COLORS_JSON__", json.dumps({a: f"var({v})" for a, v in COLOR_VAR.items()}))
    html = html.replace("__META_JSON__", json.dumps(report["meta"]))
    html = html.replace("__DISK_JSON__", json.dumps(report["disk"]))
    html = html.replace("__GANTT_JSON__", json.dumps(report["gantt"]))
    html = html.replace("__MEM_JSON__", json.dumps(report["mem"]))
    html = html.replace("__REALIGN_JSON__", json.dumps(report["realign"]))
    return html


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=DEFAULT_REPORT_PATH)
    ap.add_argument("--json", action="store_true", help="print the raw report data as JSON instead of writing HTML")
    args = ap.parse_args()

    jobs = discover_jobs()
    report = build_report(jobs)

    if args.json:
        print(json.dumps(report, indent=2))
        return

    print_summary(report)
    html = render_html(report)
    with open(args.out, "w") as f:
        f.write(html)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
