#!/usr/bin/env python3
"""Tail-sum inequalities: Remark 4.1, Lemma 4.3 (counting) and the step counts of Theorem 4.2 and Theorem 4.3.

All checks are exact (rational arithmetic, N(s) = (u box u)(s) by brute force).
  1. Counting lemma: N(s)(2 - s) <= 2 int u for random rational staircases (1 to 7 steps, ties allowed) and s random
     or at every pair sum (222,702 pairs with --trials 8000, the --full setting).
  2. Two-point laws, upper bound: the explicit floor(1/beta)-step potential excludes every law with 1/beta = n + j/D
     (n <= --lmax) or 1/2 < 1/beta < 1, w = i/D, with violation equal to the stated bound (D = 13: 1,482 laws).
  3. Two-point laws, lower bound: random and adversarial m-step staircases never exclude a law with 1/beta >= m + 1
     (2 x --lower-trials staircases; 40,000 with --full).
  4. Two-point laws against the exact maximum over staircases with r and r - 1 steps (r <= 3, D = 6).
  5. Symmetric three-point laws: the explicit r-step potential has N = (r, r, r+1) at 1-b, 1, 1+b and its inequality
     reads t >= dist(1/b, 2Z) b (r <= 8, 4,992 spacings); the exact smallest middle mass t_r(b) for the listed spacings.
  6. Remark 4.1: the thresholds (1/5, 2/5, 1/2, 4/5) and (1/5, 2/5, 3/5, 4/5) give the same xi, (2/5, ..., 9/5), and
     the same u box u on the grid 1/40, with right sides 21/5 and 4.
Example 4.2 (repeated thresholds) is checked by check_repeated_thresholds.py, Proposition S4.1 and the exact
staircase programs by check_staircase_lp.py.
Defaults: --trials 3000, --lower-trials 8000 (about 40 s). --full: 8000 and 20000 (about 90 s).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import tail_bounds as TB

THRESH = {2: ["3/13", "1/5", "2/11", "1/6", "3/19", "2/13"], 3: ["2/13", "1/7", "2/15", "3/23", "1/8", "2/17", "1/9"]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/tail_bounds")
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--lower-trials", type=int, default=8000)
    ap.add_argument("--denominator", type=int, default=13)
    ap.add_argument("--lmax", type=int, default=10)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.trials = max(a.trials, 8000); a.lower_trials = max(a.lower_trials, 20000)
    res = {}; checks = {}
    t = Timer(); tot, bad = TB.counting_lemma_check(a.trials, a.seed)
    res["counting_lemma"] = {"pairs": tot, "failures": bad}
    checks["counting_lemma"] = report("counting lemma", bad == 0, f"{tot} pairs ({t()} s)")

    t = Timer(); cnt, bad, eq = TB.two_point_upper_check(a.denominator, a.lmax)
    res["two_point_upper"] = {"laws": cnt, "failures": bad, "violation_equal_to_bound": eq}
    checks["two_point_upper"] = report("two-point laws: explicit floor(1/beta)-step potential", bad == 0,
                                       f"{cnt} laws, violation equal to the bound in {eq} ({t()} s)")

    t = Timer(); tot, bad = TB.two_point_lower_check(a.lower_trials, a.seed)
    res["two_point_lower"] = {"staircases": tot, "violations": bad}
    checks["two_point_lower"] = report("two-point laws: fewer steps never exclude", bad == 0, f"{tot} staircases ({t()} s)")

    t = Timer(); cnt, bad = TB.two_point_lp_check(6, 3, 3)
    res["two_point_exact_lp"] = {"laws": cnt, "failures": bad}
    checks["two_point_exact_lp"] = report("two-point laws vs exact staircase maximum (r <= 3)", bad == 0, f"{cnt} laws ({t()} s)")

    t = Timer(); cnt, bad = TB.three_point_family_check(8)
    res["three_point_family"] = {"spacings": cnt, "failures": bad}
    checks["three_point_family"] = report("symmetric three-point laws: explicit r-step potential", bad == 0 and cnt == 4992,
                                          f"{cnt} spacings ({t()} s)")
    t = Timer(); tr = {}; ok = True
    for r, bs in THRESH.items():
        for s in bs:
            b = Fr(s); L = 1 / b; val = TB.smallest_middle_mass_exact(b, r)
            target = TB.dist2(L) / L if L <= 2 * r + 2 else Fr(0)
            tr[f"r={r}, b={s}"] = str(val); ok &= val == target
    res["t_r"] = tr
    checks["t_r"] = report("exact t_r(b) = dist(L,2Z)/L for L <= 2r+2 and 0 beyond", ok, f"({t()} s)")

    res["checks"] = checks
    # Remark 4.1: two threshold tuples with the same xi but different right sides; the bound of the first is never attained.
    from dsxray import tailsums as TS
    t = Timer()
    e1 = [Fr(1, 5), Fr(2, 5), Fr(1, 2), Fr(4, 5)]; e2 = [Fr(1, 5), Fr(2, 5), Fr(3, 5), Fr(4, 5)]
    want = [Fr(n, 5) for n in (2, 3, 4, 5, 6, 7, 8, 9)]
    grid = [Fr(j, 40) for j in range(0, 81)]
    same_box = all(TS.box_brute(e1, x) == TS.box_brute(e2, x) for x in grid)
    rhs1 = 2 * sum(1 - e for e in e1); rhs2 = 2 * sum(1 - e for e in e2)
    checks["remark_notsharp"] = report("Remark 4.1: (1/5,2/5,1/2,4/5) and (1/5,2/5,3/5,4/5) give xi = (2/5,...,9/5), "
                                       "right sides 21/5 and 4", TS.xi(e1) == want and TS.xi(e2) == want and same_box
                                       and rhs1 == Fr(21, 5) and rhs2 == 4, f"({t()} s)")
    finish(res, out_dir(a.out), "tail_bounds")


if __name__ == "__main__":
    main()
