#!/usr/bin/env python3
"""Corollary S3.1 (four equidistant atoms, the 1,104 comparisons of Section S5) against a generic fibered linear
program (floating point, a cross-check).

Spacings 1/b in {2..9, 2.5, 3.3, e, pi, 5.5, 4.7, 1 + sqrt 5} and windows {a, a+b, a+2b, a+3b} with a = 1 - b s,
s = 0, 1/4, ..., 3 and s = 1/3, 2/3, 5/3, kept when 0 < a and a + 3b < 2 (215 windows). For each window the corollary
predicts either infeasibility or the attainable set as a Minkowski combination of explicit vertex lists; the generic
fibered LP (no block structure) is compared through support functions in 8 random directions per feasible window:
77 infeasible windows (prediction agrees), 1,104 comparisons, largest discrepancy below 1e-14 (Section S5; about
4e-15 with the default seed).
Also checks the vertex rule of the window lemma (exact) for m <= 60. Runs in a few seconds.
"""
import argparse
import math

import numpy as np

from common import Timer, finish, out_dir, report

from dsxray.fibered import support_function_on_window
from dsxray.windows import candidate_vertices, four_atom_prediction, support_of_prediction, vertices_exact

LS = [2, 3, 4, 5, 6, 7, 8, 9, 2.5, 3.3, math.e, math.pi, 5.5, 4.7, 1 + math.sqrt(5)]
SHIFTS = [i / 4 for i in range(13)] + [1 / 3, 2 / 3, 5 / 3]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/window_corollary")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--full", action="store_true", help="same as the default (the check is small)")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed); t = Timer()
    windows = comps = infeas = 0; worst = 0.0; bad = []
    for L in LS:
        b = 1 / L
        for s in SHIFTS:
            y0 = 1 - b * s
            if not (y0 > 1e-12 and y0 + 3 * b < 2 - 1e-12):
                continue
            windows += 1; pred = four_atom_prediction(b, y0)
            for _ in range(8):
                w = rng.normal(size=4)
                lp = support_function_on_window(b, y0, [0, 1, 2, 3], w)
                if (lp is None) != (pred is None):
                    bad.append((L, s)); break
                if lp is None:
                    infeas += 1; break
                comps += 1; worst = max(worst, abs(lp - support_of_prediction(pred, w)))
    res = {"windows": windows, "infeasible_windows": infeas, "support_function_comparisons": comps,
           "max_discrepancy": worst, "disagreements": len(bad)}
    checks = {"corollary": report("four-atom window corollary vs generic LP", not bad and windows == 215 and infeas == 77
                                  and comps == 1104 and worst < 1e-14, f"{res} ({t()} s)")}
    checks["vertex_rule"] = report("vertex rule of the window lemma, m <= 60 (exact)",
                                   all(set(vertices_exact(m)) <= candidate_vertices(m) for m in range(1, 61)))
    res["checks"] = checks
    finish(res, out_dir(a.out), "window_corollary")


if __name__ == "__main__":
    main()
