#!/usr/bin/env python3
"""Proposition 5.1 (the constant in Theorem 5.1 cannot exceed 2) and the floating-point experiments after Problem 5.1.

  1. Exact, for c in {2 + 1/1000, 2 + 1/100, 2 + 1/10, 7/3, 5/2, 3, 4, 5, 8, 10, 50, 1000}: with eta = (2c-1)/(6c),
     h = 1/(2c) and density c on [2 eta, 2 eta + h] u [1 + eta, 1 + eta + h]: each interval has mass 1/2,
     2 eta + h < 1 < 1 + eta, 1 + eta + h < 2, E S = 1; S <=cx 2U (D = pi_2U - pi_S >= 0, exact minimization of a
     piecewise quadratic) with D(1) = (2c-1)/(24c); the folded condition holds, with
     1/4 <= eta < eta + h <= 1/2 <= 2 eta < 2 eta + h <= 3/4 and H(1/3) = H(2/3) = (4c+1)/(144c); the single step
     at eta is violated by exactly (c-2)/(6c). For c in {1, 3/2, 2} the same construction gives violation <= 0
     (0 at c = 2).
  2. Symbolic (sympy, if installed): the mean, the violation, D'(1) = 0, D(1), F_M(1/3) = 1/3, F_M(2/3) = 2/3 and
     H(1/3) = H(2/3) = (4c+1)/(144c) as identities in c.
  3. FLOAT (numerical evidence only, as stated in the paper): for c in {1, 3/2, 2, 5/2}, the largest violation of a
     staircase with at most three thresholds on the grid 1/G, over piecewise-constant densities <= c on K cells of
     [0, 2] with S <=cx 2U (with and without the folded condition), computed by linear programming (scipy HiGHS).
     Expected: no positive violation beyond 1e-9 for c <= 2; a positive violation for c = 5/2.
Defaults: K = 40, G = 20, c in {3/2, 2, 5/2}, convex order only (about 15 s). --full: K = 80, G = 24, all four c, with
and without the folded condition (a few minutes). --no-float skips part 3.
"""
import argparse
import itertools
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import tail_bounds as TB

CS = [Fr(2) + Fr(1, n) for n in (1000, 100, 10)] + [Fr(7, 3), Fr(5, 2)] + [Fr(n) for n in (3, 4, 5, 8, 10, 50, 1000)]


def exact_member(c):
    P, eta = TB.sharpness_family(c); h = 1 / (2 * c)
    mass = [d * (hi - lo) for lo, hi, d in P]
    mean = sum(d * (hi * hi - lo * lo) / 2 for lo, hi, d in P)
    pos = 2 * eta + h < 1 < 1 + eta and 1 + eta + h < 2
    D1 = Fr(1, 4) - TB.pieces_pi(P, Fr(1))          # pi_2U(1) = 1/4
    order = Fr(1, 4) <= eta < eta + h <= Fr(1, 2) <= 2 * eta < 2 * eta + h <= Fr(3, 4)
    # folded law M: atoms 1/4 at 0 and 1, density c/2 on [eta, eta+h] and [2 eta, 2 eta+h]
    piM = lambda x: (Fr(1, 4) * max(Fr(0), -x) + Fr(1, 4) * max(Fr(0), 1 - x)
                     + TB.pieces_pi([(eta, eta + h, c / 2), (2 * eta, 2 * eta + h, c / 2)], x))
    piU = lambda x: (1 - x) ** 2 / 2
    H = [piM(Fr(1, 3)) - piU(Fr(1, 3)), piM(Fr(2, 3)) - piU(Fr(2, 3))]
    v = TB.one_step_violation(P, eta)
    return {"mass_halves": mass == [Fr(1, 2), Fr(1, 2)], "positions": pos, "mean_one": mean == 1,
            "convex_order": TB.cx_below_2U(P), "D(1)": D1, "folded_order": TB.folded_order_ok(P),
            "interval_order": order, "H": H, "violation": v}


def symbolic():
    import sympy as sp
    c, x, y = sp.symbols("c x y", positive=True)
    eta = (2 * c - 1) / (6 * c); h = 1 / (2 * c)
    ok = sp.simplify((2 * eta + h / 2 + 1 + eta + h / 2) / 2 - 1) == 0
    ok &= sp.simplify(sp.Rational(3, 2) - 2 * (1 - eta) - (c - 2) / (6 * c)) == 0
    D = (2 - x) ** 2 / 4 - (1 + eta + h / 2 - x) / 2
    ok &= sp.simplify(sp.diff(D, x).subs(x, 1)) == 0 and sp.simplify(D.subs(x, 1) - (2 * c - 1) / (24 * c)) == 0
    F1 = sp.Rational(1, 4) + c / 2 * (y - eta); F2 = sp.Rational(1, 2) + c / 2 * (y - 2 * eta)
    ok &= sp.simplify(F1.subs(y, sp.Rational(1, 3)) - sp.Rational(1, 3)) == 0
    ok &= sp.simplify(F2.subs(y, sp.Rational(2, 3)) - sp.Rational(2, 3)) == 0
    H13 = sp.Rational(1, 4) * eta + sp.integrate(F1, (y, eta, sp.Rational(1, 3))) - sp.Rational(1, 18)
    H23 = (H13 + sp.integrate(F1, (y, sp.Rational(1, 3), eta + h)) + sp.Rational(1, 2) * (eta - h)
           + sp.integrate(F2, (y, 2 * eta, sp.Rational(2, 3))) - sp.Rational(1, 6))
    target = (4 * c + 1) / (144 * c)
    return ok and sp.simplify(H13 - target) == 0 and sp.simplify(H23 - target) == 0


def float_experiment(K, G, c, fold):
    import numpy as np
    from scipy.optimize import linprog
    edges = np.linspace(0, 2, K + 1); pts = np.linspace(0, 2, 2 * K + 1); h = 2.0 / K

    def sl(j, x):
        a, b = edges[j], edges[j + 1]
        if x <= a:
            return (b - a) * ((a + b) / 2 - x)
        if x >= b:
            return 0.0
        return (b - x) ** 2 / 2
    A = [[sl(j, x) for j in range(K)] for x in pts]; B = list((2 - pts) ** 2 / 4)
    if fold:
        for yv in np.linspace(0, 1, K + 1):
            row = []
            for j in range(K):
                a, b = edges[j], edges[j + 1]
                if b <= 1:
                    lo, hi, extra = a, b, yv * (b - a)
                else:
                    lo, hi, extra = a - 1, b - 1, 0.0
                tt = max(0.0, min(hi, yv) - lo); val = (yv - lo) * tt - tt * tt / 2 if tt > 0 else 0.0
                row.append(-(val + extra))
            A.append(row); B.append(-yv * yv)
    mass = np.full(K, h); mean = np.array([(edges[j] + edges[j + 1]) / 2 * h for j in range(K)])
    best = -9.0; grid = np.arange(1, G) / G
    for r in (1, 2, 3):
        for e in itertools.combinations(grid, r):
            E = list(e) + [1.0]
            X = [max(E[i] + E[l - 1 - i] for i in range(max(0, l - 1 - r), min(r, l - 1) + 1)) for l in range(1, 2 * r + 1)]
            cc = np.zeros(K)
            for xv in X:
                cc += np.clip(np.minimum(edges[1:], 2) - np.maximum(edges[:-1], xv), 0, None)
            res = linprog(-cc, A_ub=A, b_ub=np.array(B) + 1e-12, A_eq=np.vstack([mass, mean]), b_eq=[1, 1],
                          bounds=[(0, c)] * K, method="highs")
            if res.status == 0:
                best = max(best, -res.fun - 2 * sum(1 - xv for xv in e))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/density_sharpness")
    ap.add_argument("--cells", type=int, default=40)
    ap.add_argument("--grid", type=int, default=20)
    ap.add_argument("--no-float", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.cells = max(a.cells, 80); a.grid = max(a.grid, 24)
    res = {"exact": {}}; checks = {}
    t = Timer(); ok = True
    for c in CS:
        m = exact_member(c)
        good = (m["mass_halves"] and m["positions"] and m["mean_one"] and m["convex_order"] and m["folded_order"]
                and m["interval_order"] and m["D(1)"] == (2 * c - 1) / (24 * c)
                and m["H"] == [(4 * c + 1) / (144 * c)] * 2 and m["violation"] == (c - 2) / (6 * c))
        ok &= good
        res["exact"][str(c)] = {key: (str(v) if isinstance(v, Fr) else [str(z) for z in v] if isinstance(v, list) else v)
                               for key, v in m.items()}
    checks["prop_densitysharp"] = report("Proposition 5.1: mean 1, convex and folded order, D(1), H(1/3) = H(2/3), "
                                       "violation (c-2)/(6c)", ok, f"{len(CS)} values of c ({t()} s)")
    low = {str(c): str(TB.one_step_violation(*TB.sharpness_family(c))) for c in (Fr(1), Fr(3, 2), Fr(2))}
    res["c_at_most_2"] = low
    checks["c_at_most_2"] = report("the construction gives violation <= 0 for c = 1, 3/2 and 0 at c = 2",
                                   all(Fr(v) <= 0 for v in low.values()) and Fr(low["2"]) == 0)
    try:
        sym = symbolic(); checks["symbolic"] = report("symbolic identities in c (sympy)", sym)
    except ImportError:
        print("[SKIP] sympy not installed: symbolic identities not checked")
    if not a.no_float:
        t = Timer(); fl = {}
        for cv in ((1.0, 1.5, 2.0, 2.5) if a.full else (1.5, 2.0, 2.5)):
            fl[str(cv)] = {"convex_order": float(float_experiment(a.cells, a.grid, cv, False))}
            if a.full:
                fl[str(cv)]["convex_and_folded"] = float(float_experiment(a.cells, a.grid, cv, True))
        res["float_experiment"] = {"cells": a.cells, "grid": a.grid, "max_violation": fl}
        okf = all(v <= 1e-9 for cv, d in fl.items() if float(cv) <= 2 for v in d.values())
        okf &= all(v > 1e-3 for v in fl["2.5"].values())
        checks["float_experiment"] = report(f"FLOAT: no violation for c <= 2, a violation at c = 5/2 "
                                            f"(K = {a.cells}, grid 1/{a.grid}, <= 3 thresholds)", okf, f"{fl} ({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "density_sharpness")


if __name__ == "__main__":
    main()
