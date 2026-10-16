#!/usr/bin/env python3
"""Run the test suite and every reproduction script at one of three tiers.

  python3 run.py                 quick tier: tests and every script at its small default range (measured: 145 s of
                                 wall time on a 2-core machine, nproc = 2)
  python3 run.py --tier medium   larger ranges for some steps
  python3 run.py --tier full     every script with --full, the ranges stated in the paper (measured without the
                                 tests: 630 s on the same machine)

Options: --workers N (default 2), --out DIR (default runs), --only NAME[,NAME...] (run selected steps),
--skip-tests. Each step writes OUT/<step>/<name>.json; OUT/run_summary.json lists every step with its status.
pytest is required (the test suite runs first). Exit status 0 iff every step passed.
"""
import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))

STEPS = ["tests", "two_point_example", "mww_steps", "window_corollary", "window_vertices", "small_atoms", "tail_bounds",
         "staircase_lp", "repeated_thresholds", "density", "density_sharpness", "lattice_tightness", "single_step_density",
         "two_step_density", "nonincreasing_gaps_density", "equally_spaced_density", "density_coupling",
         "hinge_sum_failure", "density_catalogues", "worked_example"]

SCRIPTS = {
    "two_point_example": "check_two_point_example.py", "mww_steps": "check_mww_steps.py",
    "window_corollary": "check_window_corollary.py", "window_vertices": "count_window_vertices.py",
    "small_atoms": "check_small_atoms.py", "tail_bounds": "check_tail_bounds.py", "staircase_lp": "check_staircase_lp.py",
    "repeated_thresholds": "check_repeated_thresholds.py", "density": "check_density_construction.py",
    "density_sharpness": "check_density_sharpness.py", "lattice_tightness": "check_lattice_tightness.py",
    "single_step_density": "check_single_step_density.py", "two_step_density": "check_two_step_density.py",
    "nonincreasing_gaps_density": "check_nonincreasing_gaps_density.py",
    "equally_spaced_density": "check_equally_spaced_density.py", "density_coupling": "check_density_coupling.py",
    "hinge_sum_failure": "check_hinge_sum_failure.py", "density_catalogues": "check_density_catalogues.py",
    "worked_example": "check_worked_example.py"}

QUICK = {}

MEDIUM = {
    "mww_steps": ["--full"], "tail_bounds": ["--full"], "staircase_lp": ["--full"], "density": ["--full"],
    "density_sharpness": ["--full"], "small_atoms": ["--kmax", "6", "--exact-kmax", "5", "--chain-kmax", "10"],
    "single_step_density": ["--grids", "8,12", "--float-rmax", "6", "--float-grid", "30"],
    "two_step_density": ["--full"], "nonincreasing_gaps_density": ["--full"], "equally_spaced_density": ["--full"],
    "density_catalogues": ["--full"]}

WORKER_STEPS = {"single_step_density", "density_catalogues"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", choices=["quick", "medium", "full"], default="quick")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(ROOT, "runs"))
    ap.add_argument("--only", default="")
    ap.add_argument("--skip-tests", action="store_true")
    a = ap.parse_args()
    if sys.version_info[:2] < (3, 12):
        raise SystemExit("Python 3.12 or newer is required")
    steps = [s for s in STEPS if not a.only or s in a.only.split(",")]
    if a.skip_tests and "tests" in steps:
        steps.remove("tests")
    os.makedirs(a.out, exist_ok=True)
    env = dict(os.environ, PYTHONPATH=os.path.join(ROOT, "src") + os.pathsep + os.environ.get("PYTHONPATH", ""),
               PYTHONDONTWRITEBYTECODE="1")
    summary = []; t_all = time.time()
    for s in steps:
        if s == "tests":
            cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", os.path.join(ROOT, "tests")]
        else:
            cmd = [sys.executable, os.path.join(ROOT, "scripts", SCRIPTS[s]), "--out", os.path.join(a.out, s)]
            if s in WORKER_STEPS:
                cmd += ["--workers", str(a.workers)]
            if a.tier == "quick":
                cmd += QUICK.get(s, [])
            elif a.tier == "medium":
                cmd += MEDIUM.get(s, [])
            elif a.tier == "full":
                cmd += ["--full"]
        print(f"\n=== {s}: {' '.join(os.path.relpath(c, ROOT) if c.startswith(ROOT) else c for c in cmd[1:])}", flush=True)
        t = time.time(); rc = subprocess.call(cmd, cwd=ROOT, env=env)
        summary.append({"step": s, "status": "pass" if rc == 0 else "FAIL", "seconds": float(f"{time.time() - t:.1f}")})
        print(f"=== {s}: {'pass' if rc == 0 else 'FAIL'} ({summary[-1]['seconds']} s)", flush=True)
    ok = all(x["status"] == "pass" for x in summary)
    with open(os.path.join(a.out, "run_summary.json"), "w") as f:
        json.dump({"tier": a.tier, "steps": summary, "all_passed": ok}, f, indent=1)
    print(f"\n{'ALL STEPS PASSED' if ok else 'SOME STEPS FAILED'} ({a.tier} tier, {time.time() - t_all:.0f} s)")
    for x in summary:
        print(f"  {x['status']:4s}  {x['seconds']:8.1f} s  {x['step']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
