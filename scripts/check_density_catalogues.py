#!/usr/bin/env python3
"""Certificate catalogues for staircases at density two: exact verification of the deposited certificates.

Claim: every staircase u = sum_{j in J} 1{x >= j/K} listed below satisfies the staircase inequality
E N(S) <= 2 int u for every S with density at most 2 and S <=cx 2U (Lemma 5.1 with c = 2); each record carries a
convex piecewise linear phi (node values at i/(Km), i = 0..2Km, exact rationals) with
(1/2) int phi + 2 int (N - phi)^+ <= 2 int u. The catalogues:
  data/density/density_two_certificates_K08.jsonl.gz       all J with >= 2 elements, K = 8,  m = 4:    120 records
  data/density/density_two_certificates_K12.jsonl.gz       all J with >= 2 elements, K = 12, m = 4:  2,036 records
  data/density/density_two_certificates_K16.jsonl.gz       all J with >= 2 elements, K = 16, m = 3: 32,752 records
  data/density/density_two_certificates_two_step_K48.jsonl.gz  all J with 2 elements, K = 48, m = 4: 1,081 records
Minimum exact margins 7/512, 7/768, 1/144, 7/3072, each at the two largest thresholds.
The candidates were proposed by a floating-point linear program (as in check_single_step_density.py); only the exact
verification below is evidence. Repeated thresholds are not covered.

Checks, exact rational arithmetic, for every record:
  1. node count 2Km + 1 and convexity (nondecreasing slopes);
  2. N on the cells of width 1/(Km) from the Minkowski sums of the constant-level intervals of u, and equal to
     N_on_cells of density_certificates (the max-plus formula of Lemma 4.1);
  3. J(phi) by exact integration (each linear piece of N - phi split at its zero); the stored margin equals
     budget - J(phi), and is >= 0;
  4. the catalogue lists every required threshold set exactly once (no repetition, nothing missing).
Default: K = 8, 12 and the two-step catalogue K = 48 (about 5 s on 2 cores). --full: also K = 16 (about 1 minute on 2 cores).
"""
import argparse
import os
from fractions import Fraction as Fr
from itertools import combinations

from common import DATA, Timer, finish, out_dir, pmap, read_jsonl, report

from dsxray import density_certificates as DC
from dsxray import density_two as D2

CATALOGUES = {"K08": ("density_two_certificates_K08.jsonl.gz", 8, 4, None, 120, Fr(7, 512), [6, 7]),
              "K12": ("density_two_certificates_K12.jsonl.gz", 12, 4, None, 2036, Fr(7, 768), [10, 11]),
              "K16": ("density_two_certificates_K16.jsonl.gz", 16, 3, None, 32752, Fr(1, 144), [14, 15]),
              "two_step_K48": ("density_two_certificates_two_step_K48.jsonl.gz", 48, 4, 2, 1081, Fr(7, 3072), [46, 47])}


def _check(rec):
    ok, margin = D2.verify_catalogue_record(rec)
    same_N = D2.min_convolution_cells(rec["J"], rec["K"], rec["m"]) == DC.N_on_cells(tuple(rec["J"]), rec["K"], rec["m"])
    return ok and same_N, margin


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/density_catalogues")
    ap.add_argument("--catalogues", default="K08,K12,two_step_K48")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.catalogues = ",".join(CATALOGUES)
    out = out_dir(a.out); res = {}; checks = {}
    for name in a.catalogues.split(","):
        fn, K, m, r, count, mmin, at = CATALOGUES[name]
        t = Timer()
        recs = read_jsonl(os.path.join(DATA, "density", fn))
        shape = all(rec["K"] == K and rec["m"] == m for rec in recs)
        keys = [tuple(rec["J"]) for rec in recs]
        if r is None:
            want = {J for rr in range(2, K) for J in combinations(range(1, K), rr)}
        else:
            want = set(combinations(range(1, K), r))
        complete = len(keys) == len(set(keys)) == count and set(keys) == want
        outs = pmap(_check, recs, a.workers, chunksize=64)
        bad = [list(k) for k, (good, _) in zip(keys, outs) if not good]
        i = min(range(len(outs)), key=lambda j: outs[j][1])
        res[name] = {"file": fn, "K": K, "m": m, "records": len(recs), "complete": complete, "failures": bad[:20],
                     "min_margin": str(outs[i][1]), "min_margin_at": list(keys[i])}
        ok = shape and complete and not bad and outs[i][1] == mmin and list(keys[i]) == at
        checks[name] = report(f"catalogue {fn}: every record certified at c = 2", ok,
                              f"{len(recs)} records, complete {complete}, min margin {outs[i][1]} at {list(keys[i])} ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "density_catalogues")


if __name__ == "__main__":
    main()
