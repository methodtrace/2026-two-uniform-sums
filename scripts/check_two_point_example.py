#!/usr/bin/env python3
"""Example 3.1: the two-point law F = (delta_{3/5} + delta_{7/5})/2, checked exactly.

  - F has mass and mean one;
  - convex order S <=cx 2U: the gap pi_2U - pi_F, with pi_2U(x) = (2 - x)^2/4 and pi_F(x) = E(S - x)^+, equals x^2/4
    on [0, 3/5] and is at least 1/20 on [3/5, 7/5] (it is quadratic there, with minimum 1/20 at x = 1), and it is
    (2 - x)^2/4 >= 0 on [7/5, 2];
  - folded order: the folded law (delta_0 + delta_{2/5} + delta_{3/5} + delta_1)/4 (the law of S - 1{S > 1} and of
    1{S > 1}, averaged) has a stop-loss function at least that of U, (1 - y)^2/2, at each of its atoms, hence
    everywhere (between atoms it is affine and the stop-loss of U is convex);
  - the single step at eta = 3/10: P(S >= 3/5) + P(S >= 13/10) = 3/2 > 7/5 = 2(1 - eta), so F is not attainable;
  - on the grid k = 10 the single step u = 1{i >= 3} defines a facet of hat-P(10) (certified facet test) and its
    inequality I(u) is violated by the count vector of F.
Default and --full are the same (under a second).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray.facet_test import classify, v_of

F = [(Fr(3, 5), Fr(1, 2)), (Fr(7, 5), Fr(1, 2))]


def pi_F(x):
    return sum(m * max(Fr(0), a - x) for a, m in F)


def gap(x):
    return (2 - x) ** 2 / 4 - pi_F(x)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/two_point_example")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}; t = Timer()

    mass = sum(m for _, m in F); mean = sum(x * m for x, m in F)
    checks["mean_one"] = report("F has mass 1 and mean 1", mass == 1 and mean == 1)

    # convex order: the gap is a quadratic on each of [0, 3/5], [3/5, 7/5], [7/5, 2]; fit it exactly by three points
    # per piece, check the fit at further points, and minimize each quadratic exactly on its interval.
    pieces = [(Fr(0), Fr(3, 5)), (Fr(3, 5), Fr(7, 5)), (Fr(7, 5), Fr(2))]
    ok = True; mins = []
    for lo, hi in pieces:
        xs = [lo, (lo + hi) / 2, hi]; ys = [gap(x) for x in xs]
        # Lagrange coefficients c2 x^2 + c1 x + c0
        (x0, x1, x2), (y0, y1, y2) = xs, ys
        c2 = ((y2 - y0) / (x2 - x0) - (y1 - y0) / (x1 - x0)) / (x2 - x1)
        c1 = (y1 - y0) / (x1 - x0) - c2 * (x0 + x1); c0 = y0 - c2 * x0 ** 2 - c1 * x0
        q = lambda x: c2 * x * x + c1 * x + c0
        ok &= all(q(lo + (hi - lo) * Fr(j, 37)) == gap(lo + (hi - lo) * Fr(j, 37)) for j in range(38))
        cand = [lo, hi] + ([-c1 / (2 * c2)] if c2 != 0 and lo <= -c1 / (2 * c2) <= hi else [])
        mins.append(min(q(x) for x in cand))
    ok &= all(gap(Fr(i, 100)) == Fr(i, 100) ** 2 / 4 for i in range(0, 61))
    ok &= mins[1] == Fr(1, 20) == gap(Fr(1)) and min(mins) >= 0
    res["convex_order"] = {"piece_minima": [str(m) for m in mins]}
    checks["convex_order"] = report("S <=cx 2U: gap x^2/4 on [0, 3/5], minimum 1/20 at x = 1 on [3/5, 7/5]", ok)

    M = [Fr(0), Fr(2, 5), Fr(3, 5), Fr(1)]
    piM = lambda y: sum(max(Fr(0), z - y) for z in M) / 4
    piU = lambda y: (1 - y) ** 2 / 2
    ok = sum(M) / 4 == Fr(1, 2) and all(piM(y) >= piU(y) for y in M)
    res["folded_order"] = {"folded_atoms": [str(y) for y in M], "slack_at_atoms": [str(piM(y) - piU(y)) for y in M]}
    checks["folded_order"] = report("folded order: stop-loss of the folded law >= that of U at each atom", ok)

    eta = Fr(3, 10)
    lhs = sum(m for x, m in F if x >= 2 * eta) + sum(m for x, m in F if x >= 1 + eta); rhs = 2 * (1 - eta)
    res["single_step"] = {"eta": str(eta), "lhs": str(lhs), "rhs": str(rhs)}
    checks["single_step"] = report("single step at eta = 3/10 violated: 3/2 > 7/5", lhs == Fr(3, 2) and rhs == Fr(7, 5) and lhs > rhs)

    k = 10; u = [int(i >= 3) for i in range(k)]; v = v_of(u)
    x = [Fr(0)] * (2 * k - 1)
    for s, m in F:
        x[int(k * s) - 1] += k * m
    c = classify(u); excess = sum(p * q for p, q in zip(v, x)) - 2 * sum(u)
    res["grid_k10"] = {"facet": c["facet"], "certificate": c["cert"], "excess": str(excess)}
    checks["grid_facet"] = report("k = 10: 1{i >= 3} is a facet of hat-P(10) and I(u) is violated",
                                  c["facet"] and excess > 0, f"excess {excess} ({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "two_point_example")


if __name__ == "__main__":
    main()
