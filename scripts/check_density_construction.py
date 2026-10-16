#!/usr/bin/env python3
"""Theorem 5.1 (densities at most one): the construction of the proof, checked exactly on examples.

For S uniform on [1/2, 3/2], uniform on [0, 2], density 1 on [1/3, 5/6] u [7/6, 5/3] (the example after Theorem 5.1;
it is not unimodal) and density 1 on [1/4, 3/4] u [5/4, 7/4], each with S <=cx 2U (checked exactly), and for
k in {24, 30, 40} (k >= 3 t^(-3/2) = 24 for the margin t = 1/4): the perturbed law S' of Step 1 has density <= 1, its
randomized rounding S'_k to (1/k)Z (Step 2) has atoms <= 1/k and none at 0 or 2, mass and mean 1, the stop-loss gap
(margin) is >= 1/(4k^2) at every m/k, and k S'_k - 1 <=cx 2U_k (Step 3). Exact rational arithmetic.
Proposition 5.1 is checked by check_density_sharpness.py, Remark S4.1 by check_lattice_tightness.py.
Default and --full are the same (seconds).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import tail_bounds as TB

EXAMPLES = {
    "uniform on [1/2, 3/2]": [(Fr(1, 2), Fr(3, 2), Fr(1))],
    "uniform on [0, 2]": [(Fr(0), Fr(2), Fr(1, 2))],
    "density 1 on [1/3, 5/6] u [7/6, 5/3]": [(Fr(1, 3), Fr(5, 6), Fr(1)), (Fr(7, 6), Fr(5, 3), Fr(1))],
    "density 1 on [1/4, 3/4] u [5/4, 7/4]": [(Fr(1, 4), Fr(3, 4), Fr(1)), (Fr(5, 4), Fr(7, 4), Fr(1))]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/density")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}
    t = Timer(); ex = {}
    for name, p in EXAMPLES.items():
        mass = sum(d * (hi - lo) for lo, hi, d in p); mean = sum(d * (hi * hi - lo * lo) / 2 for lo, hi, d in p)
        ex[name] = {"mass_and_mean_one": mass == 1 and mean == 1, "convex_order": TB.cx_below_2U(p),
                    "construction": {k: all(TB.density_construction_check(p, k)) for k in (24, 30, 40)}}
    res["density_at_most_one"] = ex
    checks["construction"] = report("density <= 1: construction checks", all(
        v["mass_and_mean_one"] and v["convex_order"] and all(v["construction"].values()) for v in ex.values()),
        f"({t()} s)")

    res["checks"] = checks
    finish(res, out_dir(a.out), "density")


if __name__ == "__main__":
    main()
