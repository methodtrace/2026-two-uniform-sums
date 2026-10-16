#!/usr/bin/env python3
"""Two-step staircases at density two: the explicit certificate for every pair of thresholds 0 < a < b < 1.

Claim (proved analytically in the paper): with d = b - a, e = 1 - b, rho = d/e, the convex function
  case A (rho <= 3/4): phi(s) = (q/e)(s - t)^+, q = 2/(rho + 1), t = 2 - e alpha/q, alpha = (8/9)(4 + q);
  case B (3/4 < rho <= 1): phi(s) = (s - t)^+/e, t = 2 - e(64 + 16 rho)/17;
  case C (rho > 1): phi(s) = (s - t)^+/d + (1/e - 1/d)(s - 2b)^+, t = 2 - 2e - 8d/3
satisfies (1/2) int phi + 2 int (N - phi)^+ <= 2(1 - a) + 2(1 - b) - min(b - a, 1 - b)/21, with N = u box u, so the
staircase inequality P(S >= 2a) + P(S >= a + b) + P(S >= max(1 + a, 2b)) + P(S >= 1 + b) <= 2(1 - a) + 2(1 - b) holds
for every S with density at most 2 and S <=cx 2U, with that slack.

Checks, exact rational arithmetic, for every pair a = i/K < b = j/K (1 <= i < j < K):
  1. the jumps of N are 2a, a + b, max(1 + a, 2b), 1 + b, and agree with the independent formula tailsums.xi;
  2. J(phi) computed by exact piecewise integration; the actual margin is at least the closed-form margin of the
     case (equal to it when every hinge lies in [0, 2]) and at least min(b - a, 1 - b)/21;
  3. the closed-form margin of each case at rho = 3/4 equals d/21 (the constant 1/21 is attained);
  4. on a sample of pairs, N from direct minimization (tailsums.box_brute) agrees with the jump formula at every
     cell midpoint;
  5. (sympy, if installed) the case identities of the proof as identities in rho (e = 1 - b): case A, the deficits
     h = (13 rho - 3)/(9(1 + rho)), k = (4 + q)/9 and the margin formula in both branches h > 0, h = 0, and the
     factorization M_A - d/21 = e(3/4 - rho)(3252 rho - 252)/(1701(1 + rho)); case B, the heights
     (18 rho - 13, rho + 4, rho + 4, 21 - 16 rho)/17 and 17 J/e = 38 rho^2 - 36 rho + 98, margin
     e(-38 rho^2 + 70 rho - 30)/17.
With K = 120: 7,021 pairs (cases A, B, C: 3,009, 531, 3,481), smallest margin 1/1080.
Default: K = 40 (about 1 s). --full: K = 120 (about 3 s).
"""
import argparse
import random
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import density_two as D2
from dsxray import tailsums as TS

EXPECTED = {120: (7021, {"A": 3009, "B": 531, "C": 3481}, Fr(1, 1080))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/two_step_density")
    ap.add_argument("--denominator", type=int, default=40)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.denominator = max(a.denominator, 120)
    K = a.denominator
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); cases = {"A": 0, "B": 0, "C": 0}; bad = []; worst = None; worst_at = None; jumps_ok = True
    rows = []
    for i in range(1, K):
        for j in range(i + 1, K):
            x, y = Fr(i, K), Fr(j, K)
            jumps_ok &= D2.jumps([x, y]) == TS.xi([x, y])
            r = D2.two_step_check(x, y)
            cases[r["case"]] += 1
            if not r["ok"]:
                bad.append([str(x), str(y)])
            if worst is None or r["margin"] < worst:
                worst, worst_at = r["margin"], [str(x), str(y)]
            ratio = r["margin"] / min(y - x, 1 - y)
            rows.append(ratio)
    n = sum(cases.values())
    res["grid"] = {"K": K, "pairs": n, "cases": cases, "failures": bad[:20], "min_margin": str(worst),
                   "min_margin_at": worst_at, "min_margin_over_min_gap": str(min(rows))}
    ok = not bad and jumps_ok and min(rows) >= Fr(1, 21)
    if K in EXPECTED:
        ok &= (n, cases, worst) == EXPECTED[K]
    checks["two_step_grid"] = report(f"two-step certificate on the grid 1/{K}: margin >= min(b-a, 1-b)/21", ok,
                                     f"{n} pairs, cases {cases}, min margin {worst} at {worst_at} ({t()} s)")

    t = Timer(); ok = True
    for e in (Fr(1, 3), Fr(1, 7), Fr(2, 9)):
        d = Fr(3, 4) * e; b = 1 - e; x = b - d
        h, case, pred = D2.two_step(x, b)
        ok &= case == "A" and pred == d / 21
        dd = d + Fr(1, 10 ** 6) * e
        h2, case2, pred2 = D2.two_step(b - dd, b)
        ok &= case2 == "B" and pred2 > dd / 21 - Fr(1, 10 ** 4) * e
    checks["constant_attained"] = report("the case A margin equals (b - a)/21 at rho = 3/4", ok, f"({t()} s)")

    t = Timer(); rng = random.Random(7); ok = True; tested = 0
    for _ in range(60):
        i, j = sorted(rng.sample(range(1, K), 2)); x, y = Fr(i, K), Fr(j, K)
        xi = D2.jumps([x, y])
        for s2 in range(1, 4 * K, 2):
            s = Fr(s2, 2 * K)
            ok &= TS.box_brute([x, y], s) == sum(1 for z in xi if z <= s); tested += 1
    checks["brute_force_N"] = report("N from direct minimization equals the jump formula (sample)", ok,
                                     f"{tested} evaluations ({t()} s)")
    try:
        import sympy as sp
    except ImportError:
        sp = None
        print("[SKIP] symbolic case identities: sympy not installed", flush=True)
    if sp is not None:
        t = Timer(); r, e = sp.symbols("rho e", positive=True); d = r * e; z = []
        q = 2 / (1 + r); al = sp.Rational(8, 9) * (4 + q); h = 5 - al; k = 4 + q - al
        z += [h - (13 * r - 3) / (9 * (1 + r)), k - (4 + q) / 9, k / q - (2 * r + 3) / 9]
        for hh, xx in ((h, 13 * r - 3), (0, 0)):
            J = (e / q) * (al ** 2 / 4 + 2 * hh ** 2 + 2 * k ** 2)
            z.append(2 * d + 4 * e - J - e * (18 * r * (r + 3) - xx ** 2) / (81 * (1 + r)))
        MA = e * (18 * r * (r + 3) - (13 * r - 3) ** 2) / (81 * (1 + r))
        z.append(MA - d / 21 - e * (sp.Rational(3, 4) - r) * (3252 * r - 252) / (1701 * (1 + r)))
        al = (64 + 16 * r) / 17
        hs = [3 + 2 * r - al, 4 + r - al, 4 + r - al, 5 - al]
        z += [17 * x - y for x, y in zip(hs, [18 * r - 13, r + 4, r + 4, 21 - 16 * r])]
        J = e * (al ** 2 / 4 + sum(x ** 2 for x in hs))
        z += [J * 17 / e - (38 * r ** 2 - 36 * r + 98), (2 * d + 4 * e - J) * 17 / e - (-38 * r ** 2 + 70 * r - 30)]
        ok = all(sp.simplify(x) == 0 for x in z)
        checks["symbolic_identities"] = report("case identities of the proof, symbolically (sympy)", ok, f"{len(z)} identities ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "two_step_density")


if __name__ == "__main__":
    main()
