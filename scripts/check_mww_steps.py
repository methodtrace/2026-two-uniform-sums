#!/usr/bin/env python3
"""Proposition 6.1 (single-step obstructions for the two- and three-point laws of Mao, Wang and Wang), exact.

  1. Two-point laws (Proposition 6.1(a)): every mean-one law on {a, a+b} with a < 1 < a+b, 1/2 < b < 1, a and b on
     the grid 1/100 (3,626 laws), violates the single-step inequality at eta = min(a/2, a+b-1).
  2. Symmetric three-point laws (Proposition 6.1(b)): for 1/4 <= b < 1 (grid 1/100) and middle mass f2 in
     {0, 1/100, ..., 1} (7,575 laws), the single-step inequality at the stated eta holds iff f2 >= dist(1/b, 2Z) b.
  3. The law 1/2 (delta_{3/5} + delta_{7/5}) of Example 3.1 violates the step at 3/10, and on the grid k = 10 the
     single-step facet at 3 by 15 > 14.
  4. The law (1/160, 1/40, 31/32) on {123/280, 29/40, 283/280} (b = 2/7) of the paragraph after Proposition 6.1
     satisfies every single-step inequality (exact minimum slack 0), and the exact maximum violation over staircases
     with at most r steps is 0 for r = 1, 2 and 3/560 for r = 3, in both signs (Proposition S4.1; see also
     check_staircase_lp.py).
  5. With --search (additional, not used in the paper): an arrangement search (bisection to
     2^-25 on the middle mass for b = 0.17..0.24 with two steps; two-point laws with b in {0.35, 0.40, 0.45} excluded
     by two steps and by no single step). This search evaluates exact rational points but is not a certificate.
Default and --full: items 1-4 (seconds).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import mww_steps as MS
from dsxray.tail_bounds import max_violation


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/mww_steps")
    ap.add_argument("--search", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}
    t = Timer(); n = bad = 0
    for bn in range(51, 100):
        b = Fr(bn, 100)
        for an in range(1, 100):
            aa = Fr(an, 100)
            if not (aa < 1 < aa + b):
                continue
            w = (aa + b - 1) / b; law = {aa: w, aa + b: 1 - w}; eta = min(aa / 2, aa + b - 1)
            n += 1; bad += not (0 < eta < 1 and MS.onestep_lhs(law, eta) > 2 * (1 - eta))
    res["two_point"] = {"laws": n, "failures": bad}
    checks["two_point"] = report("two-point laws, 1/2 < b < 1", bad == 0 and n == 3626, f"{n} laws ({t()} s)")

    t = Timer(); n = bad = 0
    for bd in range(25, 100):
        b = Fr(bd, 100); L = 1 / b
        eta = b if (L <= 2 or 3 <= L <= 4) else (1 - b) / 2
        dist = min(abs(L - 2 * m) for m in range(0, 4))
        for fn in range(0, 101):
            f2 = Fr(fn, 100); f = (1 - f2) / 2; law = {1 - b: f, Fr(1): f2, 1 + b: f}
            n += 1; bad += (MS.onestep_lhs(law, eta) <= 2 * (1 - eta)) != (f2 >= dist / L)
            if L == 2:
                bad += (MS.onestep_lhs(law, (1 - b) / 2) <= 2 * (1 - (1 - b) / 2)) != (f2 >= dist / L)
    res["symmetric_three_point"] = {"laws": n, "mismatches": bad}
    checks["three_point"] = report("symmetric three-point laws, 1/4 <= b < 1", bad == 0 and n == 7575, f"{n} laws ({t()} s)")

    F = {Fr(3, 5): Fr(1, 2), Fr(7, 5): Fr(1, 2)}
    u = [1 if i >= 3 else 0 for i in range(10)]
    v = [min(u[i] + u[s - i] for i in range(10) if 0 <= s - i < 10) for s in range(19)]
    checks["example"] = report("two-point example", MS.onestep_lhs(F, Fr(3, 10)) == Fr(3, 2) > 2 * (1 - Fr(3, 10))
                               and 5 * v[5] + 5 * v[13] == 15 > 2 * sum(u))

    t = Timer()
    A = [Fr(123, 280), Fr(29, 40), Fr(283, 280)]; M = [Fr(1, 160), Fr(1, 40), Fr(31, 32)]
    slack = MS.single_step_min_slack(dict(zip(A, M)))[0]
    viol = {r: (max_violation(r, A, M)[0], max_violation(r, [2 - x for x in A], M)[0]) for r in (1, 2, 3)}
    res["asymmetric"] = {"mean": str(sum(x * m for x, m in zip(A, M))), "single_step_min_slack": str(slack),
                         "max_violation": {r: [str(x) for x in v] for r, v in viol.items()}}
    checks["asymmetric"] = report("asymmetric law: all single steps hold, excluded by 3 steps, not by 2",
                                  sum(x * m for x, m in zip(A, M)) == 1 and slack == 0
                                  and viol[1] == (0, 0) and viol[2] == (0, 0) and viol[3] == (Fr(3, 560), Fr(3, 560)),
                                  f"({t()} s)")
    if a.search:
        t = Timer(); sym = []
        for c in range(17, 25):
            b = Fr(c, 100)
            sym.append({"b": str(b), "bound": str(MS.mww_bound(b)), "two_step": str(MS.smallest_middle_mass(b, 2)),
                        "one_step": str(MS.smallest_middle_mass(b, 1))})
        ok = all(abs(Fr(s["two_step"]) - Fr(s["bound"])) < Fr(1, 2**24) and abs(Fr(s["one_step"]) - Fr(s["bound"])) > Fr(1, 100)
                 for s in sym)
        res["search_symmetric"] = sym
        two = []
        for c in (35, 40, 45):
            b = Fr(c, 100); n = n1 = n2 = 0
            for ai in range(1, 100):
                aa = Fr(ai, 100)
                if not (aa < 1 < aa + b):
                    continue
                p = (1 - aa) / b; n += 1
                n1 += MS.best_viol([aa, aa + b], [1 - p, p], 1) > 0
                n2 += MS.best_viol([aa, aa + b], [1 - p, p], 2) > 0
            two.append({"b": str(b), "laws": n, "excluded_by_one_step": n1, "excluded_by_two_steps": n2})
            ok &= n1 == 0 and n2 == n
        res["search_two_point"] = two
        checks["search"] = report("arrangement search (two steps reproduce the three-point bound; two-point laws)", ok, f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "mww_steps")


if __name__ == "__main__":
    main()
