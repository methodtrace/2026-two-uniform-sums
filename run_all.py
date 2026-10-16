#!/usr/bin/env python3
"""Run every check at the quick tier and compare its results with MANIFEST.sha256, byte for byte.

    python3 run_all.py            # run all checks, write results/, compare with MANIFEST.sha256
    python3 run_all.py --write    # also rewrite MANIFEST.sha256 from the new results

run.py does the work (the test suite, then every script at its quick range) in a scratch folder, runs/run_all/.
Each script's JSON output is then written to results/<step>.json in a canonical form that depends only on the
code, not on the machine: keys sorted, run times, output paths and worker counts removed, the record of which
optional engines ran removed (the C helper for Example 4.2(b), used when a C compiler works, and pycddlib for the
small-atoms vertices, used when installed: each is a cross-check that must agree with the exact computation),
and every floating-point number below 1e-9 in absolute value written as 0.0 and every other one rounded to six
significant digits (floating point enters only in cross-checks; every statement of the paper is checked exactly).
The SHA-256 digests of these files must equal those in MANIFEST.sha256.

Needs CPython 3.12 or newer with numpy, scipy, pytest and sympy (requirements.txt): without sympy the symbolic
checks are skipped, and their results would differ from the manifest. pycddlib is optional here.
"""
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DROP = {"seconds", "out", "workers", "engines", "sources", "pycddlib"}
TIME = re.compile(r"\d+(?:\.\d+)? s\b")
NEEDED = {"numpy": "numpy", "scipy": "scipy", "pytest": "pytest", "sympy": "sympy"}


def canonical(x):
    if isinstance(x, dict):
        return {k: canonical(v) for k, v in sorted(x.items()) if k not in DROP}
    if isinstance(x, list):
        return [canonical(v) for v in x]
    if isinstance(x, float):
        if not math.isfinite(x) or abs(x) < 1e-9:
            return 0.0 if math.isfinite(x) else str(x)
        return float(f"{x:.6g}")
    if isinstance(x, str):
        return TIME.sub("<time> s", x)
    return x


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    if sys.version_info[:2] < (3, 12):
        sys.exit("run_all.py needs CPython 3.12 or newer")
    missing = [pkg for mod, pkg in NEEDED.items() if importlib.util.find_spec(mod) is None]
    if missing:
        sys.exit("missing packages: " + " ".join(missing) + "\ninstall them with: python3 -m pip install -r requirements.txt")
    scratch = ROOT / "runs" / "run_all"
    shutil.rmtree(scratch, ignore_errors=True)
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONHASHSEED"] = "0"
    proc = subprocess.run([sys.executable, str(ROOT / "run.py"), "--out", str(scratch)], cwd=ROOT, env=env)
    if proc.returncode != 0:
        sys.exit("FAILED: run.py reported a failing step (see the output above)")
    results = ROOT / "results"
    shutil.rmtree(results, ignore_errors=True)
    results.mkdir()
    lines = []
    for step_dir in sorted(p for p in scratch.iterdir() if p.is_dir()):
        files = sorted(step_dir.glob("*.json"))
        if len(files) != 1:
            sys.exit(f"FAILED: expected one result file in {step_dir.name}, found {len(files)}")
        data = canonical(json.loads(files[0].read_text(encoding="utf-8")))
        target = results / f"{step_dir.name}.json"
        target.write_bytes((json.dumps(data, indent=1, sort_keys=True, ensure_ascii=True) + "\n").encode("ascii"))
        lines.append(f"{digest(target)}  results/{target.name}")
    new = "\n".join(lines) + "\n"
    manifest = ROOT / "MANIFEST.sha256"
    if "--write" in sys.argv:
        manifest.write_bytes(new.encode("ascii"))
        print("MANIFEST.sha256 written", file=sys.stderr)
        return
    if not manifest.exists():
        sys.exit("MANIFEST.sha256 is missing")
    if manifest.read_bytes().replace(b"\r\n", b"\n") == new.encode("ascii"):
        print("all results match MANIFEST.sha256 byte for byte", file=sys.stderr)
        return
    old = set(manifest.read_bytes().decode("ascii").splitlines())
    for line in lines:
        if line not in old:
            print(f"MISMATCH: {line.split()[1]}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
