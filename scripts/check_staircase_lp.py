#!/usr/bin/env python3
"""Exact staircase linear programs and tail-sum identities (Sections 4 and 5): Proposition S4.1 and checks of Lemma 4.1, Lemma 4.2 and Lemma 4.3,
Remark 4.2 and Theorem 4.2, Theorem 4.3. All arithmetic is exact (fractions).

  1. Proposition S4.1: F_* = (1/160, 1/40, 31/32) on {123/280, 29/40, 283/280} has mean 1 and is the law of
     Corollary 6.1(iv) with b = 2/7, a = 29/40. The thresholds (3/280, 29/80, 363/560) give
     xi = (3/140, 209/560, 29/40, 283/280, 109/80, 923/560), (u box u) = 2, 3, 4 at the atoms, left side 317/80, right
     side 277/70, violation 3/560. The exact supremum of the violation over staircases with r thresholds is 0 for
     r = 1, 2 and 3/560 for r = 3, for F_* and for its reflection (so, by Lemma 4.2, no potential with at most two unit
     steps of either sign excludes F_*), computed by two independent programs: the count-vector reduction of the
     proof with vertex enumeration (dsxray.tailsums.stair_lp_bound) and dsxray.tail_bounds.max_violation.
  2. Lemma 4.1: the closed formula (u box u)(s) = #{l : xi_l <= s} agrees with direct minimization, for random
     rational thresholds with ties and thresholds equal to 1, at all pair sums and random points; Lemma 4.3 (counting)
     on the same pairs.
  3. Lemma 4.2 (reflection): for random staircases u and laws F, the violation of -u at F is at most the violation of
     the staircase with thresholds {1 - eta_q : eta_q < 1} at the reflected law.
  4. Remark 4.2: on random grids k and integer step multisets, the tail-sum inequality with thresholds j/k counted
     with multiplicity equals I(u)/k for the integer potential u (random lattice laws).
  5. Theorem 4.2 by the exact program: for two-point laws with 1/beta in {3/2, 5/3, 9/4, 7/3, 13/5, 10/3, 17/5} and
     w in {1/5, 1/2, 3/4}, staircases with fewer than max(1, floor(1/beta)) thresholds never exclude the law or its
     reflection, and staircases with that many do.
  6. Theorem 4.3 by the exact program: for symmetric three-point laws with L in {5/2, 7/2, 9/2, 11/3, 13/3, 21/4} and
     middle mass t below dist(L, 2Z)/L, the least r with a positive maximum is max(1, floor(L/2)); at
     t = dist(L, 2Z)/L the maximum with that r is <= 0.
Defaults: --trials 3000 (random checks). --full: --trials 12000.
"""
import argparse
import math
import random
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import tailsums as TS
from dsxray.tail_bounds import max_violation

A_STAR = [Fr(123, 280), Fr(29, 40), Fr(283, 280)]
M_STAR = [Fr(1, 160), Fr(1, 40), Fr(31, 32)]
ETA_STAR = [Fr(3, 280), Fr(29, 80), Fr(363, 560)]
XI_STAR = [Fr(3, 140), Fr(209, 560), Fr(29, 40), Fr(283, 280), Fr(109, 80), Fr(923, 560)]


def dist2(L):
    return min(abs(L - 2 * m) for m in range(0, int(L) + 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/staircase_lp")
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.trials = max(a.trials, 12000)
    res = {}; checks = {}

    # 1. Proposition S4.1
    t = Timer(); law = list(zip(A_STAR, M_STAR)); b = Fr(2, 7); aa = Fr(29, 40)
    X = TS.xi(ETA_STAR); vals = [TS.box_brute(ETA_STAR, s) for s in A_STAR]
    lhs = sum(m * c for m, c in zip(M_STAR, vals)); rhs = 2 * sum(1 - e for e in ETA_STAR)
    ok1 = (sum(M_STAR) == 1 and sum(x * m for x, m in law) == 1 and A_STAR == [aa - b, aa, aa + b]
           and X == XI_STAR and vals == [2, 3, 4] and lhs == Fr(317, 80) and rhs == Fr(277, 70)
           and lhs - rhs == Fr(3, 560) == TS.violation(ETA_STAR, law))
    sup = {}; ok2 = True
    for r in (1, 2, 3):
        p1 = TS.stair_lp_bound(law, r)[0]; p2 = TS.stair_lp_bound(TS.reflect(law), r)[0]
        q1 = max_violation(r, A_STAR, M_STAR)[0]; q2 = max_violation(r, [2 - x for x in A_STAR], M_STAR)[0]
        sup[r] = {"law": str(p1), "reflection": str(p2), "second_program_law": str(q1), "second_program_reflection": str(q2)}
        target = Fr(3, 560) if r == 3 else Fr(0)
        ok2 &= p1 == p2 == q1 == q2 == target
    res["prop_asym"] = {"xi": [str(z) for z in X], "values_at_atoms": vals, "lhs": str(lhs), "rhs": str(rhs),
                              "supremum_by_r": sup}
    checks["prop_asym_example"] = report("Proposition S4.1: thresholds (3/280, 29/80, 363/560), xi, 317/80 - 277/70 = 3/560", ok1)
    checks["prop_asym_two_steps"] = report("Proposition S4.1: supremum 0 for r = 1, 2 and 3/560 for r = 3, law and "
                                          "reflection, two programs", ok2, f"({t()} s)")

    # 2. Lemma 4.1 and Lemma 4.3
    t = Timer(); random.seed(a.seed); n = bad = badL = 0
    for _ in range(a.trials):
        r = random.randint(1, 6); D = random.choice([4, 6, 7, 12, 13])
        eta = sorted(Fr(random.randint(1, D), D) for _ in range(r))
        pts = {Fr(0), Fr(2)} | {x + y for x in eta + [Fr(1)] for y in eta + [Fr(1)]}
        pts |= {Fr(random.randint(0, 6 * D), 3 * D) for _ in range(4)}
        for s in pts:
            n += 1; N = TS.box_brute(eta, s)
            bad += N != TS.box_formula(eta, s)
            badL += N * (2 - s) > 2 * sum(1 - e for e in eta)
    res["lem_tailsgen"] = {"pairs": n, "formula_mismatches": bad, "counting_lemma_failures": badL}
    checks["lem_tailsgen"] = report("Lemma 4.1 formula = direct minimization (ties, thresholds at 1)", bad == 0, f"{n} pairs ({t()} s)")
    checks["lem_counting"] = report("Lemma 4.3: N(s)(2 - s) <= 2 int u", badL == 0, f"{n} pairs")

    # 3. Lemma 4.2
    t = Timer(); random.seed(a.seed + 1); bad = 0; n = 0
    for _ in range(a.trials):
        r = random.randint(1, 5); D = random.choice([5, 6, 8, 9])
        eta = sorted(Fr(random.randint(1, D), D) for _ in range(r))
        law = [(Fr(random.randint(1, 2 * D - 1), D), Fr(random.randint(1, 5))) for _ in range(4)]
        tot = sum(m for _, m in law); law = [(x, m / tot) for x, m in law]
        vneg = sum(m * TS.neg_box_brute(eta, x) for x, m in law) + 2 * sum(1 - e for e in eta)
        th = sorted(1 - e for e in eta if e < 1)
        vw = TS.violation_brute(th, TS.reflect(law)) if th else Fr(0)
        n += 1; bad += vneg > vw
    checks["lem_reflectsteps"] = report("Lemma 4.2 (reflection): violation of -u at F <= that of the reflected staircase at "
                                 "the reflected law", bad == 0, f"{n} cases ({t()} s)")

    # 4. Remark 4.2
    t = Timer(); random.seed(a.seed + 2); bad = 0; n = 0
    for _ in range(a.trials):
        k = random.randint(3, 9); r = random.randint(1, 6)
        J = sorted(random.randint(1, k - 1) for _ in range(r))
        u = [sum(1 for j in J if j <= i) for i in range(k)]
        law = [(Fr(tt + 1, k), Fr(random.randint(0, 4))) for tt in range(2 * k - 1)]
        tot = sum(m for _, m in law)
        if tot == 0:
            continue
        law = [(x, m / tot) for x, m in law]
        v = TS.grid_v(u)
        I = sum(m * v[tt] for tt, (_, m) in enumerate(law)) - Fr(2, k) * sum(u)
        n += 1; bad += I != TS.violation([Fr(j, k) for j in J], law)
    checks["rem_multiplicity"] = report("Remark 4.2: tail sum with multiplicities = I(u)/k on the grid", bad == 0,
                                  f"{n} cases ({t()} s)")

    # 5. Theorem 4.2 by the exact program
    t = Timer(); cnt = bad = 0
    for L in [Fr(3, 2), Fr(5, 3), Fr(9, 4), Fr(7, 3), Fr(13, 5), Fr(10, 3), Fr(17, 5)]:
        beta = 1 / L; rs = max(1, math.floor(L))
        for w in [Fr(1, 5), Fr(1, 2), Fr(3, 4)]:
            lw = [(1 - (1 - w) * beta, w), (1 + w * beta, 1 - w)]
            for F in (lw, TS.reflect(lw)):
                ok = all(TS.stair_lp_bound(F, r)[0] <= 0 for r in range(1, rs)) and TS.stair_lp_bound(F, rs)[0] > 0
                cnt += 1; bad += not ok
    checks["thm_twopoint_lp"] = report("Theorem 4.2: least number of thresholds max(1, floor(1/beta)) (exact program)",
                                      bad == 0, f"{cnt} cases ({t()} s)")

    # 6. Theorem 4.3 by the exact program
    t = Timer(); cnt = bad = 0
    for L in [Fr(5, 2), Fr(7, 2), Fr(9, 2), Fr(11, 3), Fr(13, 3), Fr(21, 4)]:
        b = 1 / L; d = dist2(L) / L; rs = max(1, math.floor(L / 2))
        for tt in [d / 2, d * Fr(9, 10)]:
            f = (1 - tt) / 2; lw = [(1 - b, f), (Fr(1), tt), (1 + b, f)]
            ok = all(TS.stair_lp_bound(lw, r)[0] <= 0 for r in range(1, rs)) and TS.stair_lp_bound(lw, rs)[0] > 0
            cnt += 1; bad += not ok
        f = (1 - d) / 2; lw = [(1 - b, f), (Fr(1), d), (1 + b, f)]
        cnt += 1; bad += TS.stair_lp_bound(lw, rs)[0] > 0
    checks["thm_threepoint_lp"] = report("Theorem 4.3: least number of thresholds max(1, floor(L/2)) (exact program)",
                                      bad == 0, f"{cnt} cases ({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "staircase_lp")


if __name__ == "__main__":
    main()
