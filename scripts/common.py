"""Shared helpers for the scripts: paths, deposited data, output files."""
import gzip
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "src"))


def out_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def _jsonable(o):
    if hasattr(o, "item"):
        return o.item()
    if isinstance(o, (set, frozenset, tuple)):
        return sorted(o) if isinstance(o, (set, frozenset)) else list(o)
    return str(o)


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=_jsonable)
        f.write("\n")


def facets_path(k):
    p = os.path.join(DATA, "facets", f"facets_k{k:02d}.jsonl")
    return p if os.path.exists(p) else p + ".gz"


def load_facet_rows(k, path=None):
    path = path or facets_path(k)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]


class Timer:
    def __init__(self):
        self.t = time.time()

    def __call__(self):
        return float(f"{time.time() - self.t:.1f}")


def report(name, ok, detail=""):
    ok = bool(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""), flush=True)
    return ok


def finish(results, out, name):
    """Write results to OUT/name.json and exit 0 iff every check passed."""
    ok = all(results.get("checks", {}).values()) if results.get("checks") else True
    results["all_checks_passed"] = ok
    write_json(os.path.join(out_dir(out), name + ".json"), results)
    print(("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED") + f" -> {os.path.join(out, name + '.json')}", flush=True)
    sys.exit(0 if ok else 1)


def pmap(fn, items, workers=1, chunksize=64):
    """Ordered parallel map (multiprocessing) or plain map for one worker."""
    items = list(items)
    if workers > 1 and len(items) > chunksize:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            return pool.map(fn, items, chunksize=chunksize)
    return [fn(x) for x in items]


def read_jsonl(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        return [json.loads(line) for line in f if line.strip()]
