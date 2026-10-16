#!/usr/bin/env python3
"""Example 4.1: a law with density 4 that satisfies convex order and every single-step inequality but is excluded by
a two-threshold staircase. Exact rational arithmetic throughout.

The law S has density 4 on the four intervals [7/20, 33/80], [61/80, 33/40], [47/40, 99/80], [127/80, 33/20] (each of
length 1/16, mass 1/4). Checks:
  1. mass and mean one;
  2. S <=cx 2U: the gap D(x) = (2 - x)^2/4 - E(S - x)^+ is >= 0 on [0, 2]. D is quadratic between consecutive
     breakpoints {0, 2, interval ends}; on each piece it is recovered exactly from three values, checked at further
     points, and minimized exactly (endpoints and the vertex of the parabola if it lies inside). Cross-checked with
     tail_bounds.cx_below_2U;
  3. single steps: L(eta) - R(eta) = P(S >= 2 eta) + P(S >= 1 + eta) - 2(1 - eta) is piecewise linear in eta with
     breakpoints at a/2 and a - 1 for the interval ends a; it is evaluated at every breakpoint in (0, 1) and at the
     one-sided limits at 0 and 1. Every single-step inequality holds strictly on (0, 1). The largest value over
     0 < eta <= 33/40 (the range in which the left side is positive) is -3/40, attained at eta = 47/80 only. For
     eta > 33/40 the left side is 0 and L - R = -2(1 - eta), which tends to 0 as eta -> 1: the supremum over the whole
     interval (0, 1) is 0 and is not attained. Both values are reported;
  4. the staircase with thresholds (7/40, 47/80): xi = (7/20, 61/80, 47/40, 127/80) (closed formula of Lemma 4.1,
     tailsums.xi, and by direct minimization of u box u,
     tailsums.box_brute, on the grid 1/320 and at the xi), left side sum_l P(S >= xi_l) = 5/2, right side
     2 sum (1 - eta) = 99/40, violation 1/40;
  5. the folded condition law(U) <=cx (law(S - 1{S > 1}) + law(1{S > 1}))/2 is evaluated exactly and its result is
     printed (reported, not asserted).
Default and --full are the same (under a second).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import tail_bounds as TB
from dsxray import tailsums as TS

PIECES = [(Fr(7, 20), Fr(33, 80), Fr(4)), (Fr(61, 80), Fr(33, 40), Fr(4)), (Fr(47, 40), Fr(99, 80), Fr(4)),
          (Fr(127, 80), Fr(33, 20), Fr(4))]
ETA = [Fr(7, 40), Fr(47, 80)]
XI = [Fr(7, 20), Fr(61, 80), Fr(47, 40), Fr(127, 80)]


def stop_loss(x):
    """E(S - x)^+."""
    t = Fr(0)
    for lo, hi, d in PIECES:
        a = max(lo, x)
        if a < hi:
            t += d * ((hi - x) ** 2 - (a - x) ** 2) / 2
    return t


def prob_ge(s):
    return sum(d * max(Fr(0), hi - max(lo, s)) for lo, hi, d in PIECES)


def quad_through(xs, ys):
    (x0, x1, x2), (y0, y1, y2) = xs, ys
    c2 = ((y2 - y0) / (x2 - x0) - (y1 - y0) / (x1 - x0)) / (x2 - x1)
    c1 = (y1 - y0) / (x1 - x0) - c2 * (x0 + x1); c0 = y0 - c2 * x0 ** 2 - c1 * x0
    return c2, c1, c0


def exact_piecewise_min(f, bps, probes=23):
    """Exact minimum of a function that is quadratic between consecutive breakpoints; returns (min, argmin, fit_ok)."""
    best = None; fit_ok = True
    for lo, hi in zip(bps, bps[1:]):
        c2, c1, c0 = quad_through([lo, (lo + hi) / 2, hi], [f(lo), f((lo + hi) / 2), f(hi)])
        q = lambda x: c2 * x * x + c1 * x + c0
        fit_ok &= all(q(lo + (hi - lo) * Fr(j, probes)) == f(lo + (hi - lo) * Fr(j, probes)) for j in range(probes + 1))
        cand = [lo, hi]
        if c2 != 0 and lo < -c1 / (2 * c2) < hi:
            cand.append(-c1 / (2 * c2))
        for x in cand:
            if best is None or q(x) < best[0]:
                best = (q(x), x)
    return best[0], best[1], fit_ok


def folded_condition():
    """U <=cx (law(S^flat) + law(kappa))/2 with S^flat = S - 1{S > 1}, kappa = 1{S > 1}: the stop-loss
    g(y) = (E(S^flat - y)^+ + E(kappa - y)^+)/2 must be >= E(U - y)^+ = (1 - y)^2/2 on [0, 1] (equal means)."""
    flat = [(lo - (1 if lo >= 1 else 0), hi - (1 if lo >= 1 else 0), d) for lo, hi, d in PIECES]   # no piece straddles 1
    assert all(hi <= 1 or lo >= 1 for lo, hi, d in PIECES)
    p_gt1 = sum(d * (hi - lo) for lo, hi, d in PIECES if lo >= 1)

    def sl_flat(y):
        t = Fr(0)
        for lo, hi, d in flat:
            a = max(lo, y)
            if a < hi:
                t += d * ((hi - y) ** 2 - (a - y) ** 2) / 2
        return t
    mean_ok = (sum(d * (hi * hi - lo * lo) / 2 for lo, hi, d in flat) + p_gt1) / 2 == Fr(1, 2)
    g = lambda y: (sl_flat(y) + p_gt1 * (1 - y)) / 2 - (1 - y) ** 2 / 2
    bps = sorted({Fr(0), Fr(1)} | {x for lo, hi, d in flat for x in (lo, hi)})
    m, at, fit = exact_piecewise_min(g, bps)
    return mean_ok, m, at, fit


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/worked_example")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}; t = Timer()

    mass = sum(d * (hi - lo) for lo, hi, d in PIECES); mean = sum(d * (hi * hi - lo * lo) / 2 for lo, hi, d in PIECES)
    checks["mass_mean"] = report("mass 1 and mean 1; each interval has length 1/16 and mass 1/4",
                                 mass == 1 and mean == 1 and all(hi - lo == Fr(1, 16) and d * (hi - lo) == Fr(1, 4)
                                                                 for lo, hi, d in PIECES))

    D = lambda x: (2 - x) ** 2 / 4 - stop_loss(x)
    bps = sorted({Fr(0), Fr(2)} | {x for lo, hi, d in PIECES for x in (lo, hi)})
    m, at, fit = exact_piecewise_min(D, bps)
    res["convex_order"] = {"min_gap": str(m), "argmin": str(at)}
    checks["convex_order"] = report("S <=cx 2U: exact minimum of (2 - x)^2/4 - E(S - x)^+ on [0, 2] is >= 0",
                                    fit and m >= 0 and TB.cx_below_2U(PIECES), f"min {m} at x = {at}")

    diff = lambda e: prob_ge(2 * e) + prob_ge(1 + e) - 2 * (1 - e)
    ends = sorted({x for lo, hi, d in PIECES for x in (lo, hi)})
    br = sorted({e for x in ends for e in (x / 2, x - 1) if 0 < e < 1})
    # L - R is linear between consecutive points of [0] + br + [1] (one-sided values at the ends)
    lim0 = prob_ge(Fr(0)) + prob_ge(Fr(1)) - 2          # eta -> 0+: no atoms, so P is continuous
    lim1 = Fr(0)                                         # eta -> 1-: both probabilities vanish, 2(1 - eta) -> 0
    vals = {e: diff(e) for e in br}
    lin_ok = all(diff((p + q) / 2) == (diff(p) + diff(q)) / 2 for p, q in zip(br, br[1:]))
    cutoff = Fr(33, 40)
    inner = {e: v for e, v in vals.items() if e <= cutoff}
    best_e = max(inner, key=inner.get)
    holds = all(v < 0 for v in vals.values()) and lim0 < 0
    res["single_step"] = {"breakpoints": {str(e): str(v) for e, v in vals.items()}, "limit_eta_to_0": str(lim0),
                          "max_on_(0,33/40]": str(inner[best_e]), "argmax": str(best_e),
                          "sup_on_(0,1)": str(lim1), "sup_attained": False}
    checks["single_steps_hold"] = report("every single-step inequality holds strictly for 0 < eta < 1", holds and lin_ok)
    checks["single_step_max"] = report("max of left minus right over 0 < eta <= 33/40 is -3/40, only at eta = 47/80",
                                       inner[best_e] == Fr(-3, 40) and best_e == Fr(47, 80)
                                       and sum(v == Fr(-3, 40) for v in inner.values()) == 1)
    print(f"[INFO] for eta > 33/40 the left side is 0 and left minus right = -2(1 - eta) -> 0 as eta -> 1: "
          f"the supremum over (0, 1) is {lim1}, not attained", flush=True)

    xi = TS.xi(ETA)
    grid = [Fr(j, 320) for j in range(0, 641)]
    xi_brute_ok = all(TS.box_brute(ETA, s) == sum(1 for z in XI if z <= s) for s in grid + XI)
    lhs = sum(prob_ge(z) for z in xi); rhs = 2 * sum(1 - e for e in ETA)
    res["staircase"] = {"thresholds": [str(e) for e in ETA], "xi": [str(z) for z in xi], "lhs": str(lhs), "rhs": str(rhs),
                        "violation": str(lhs - rhs)}
    checks["staircase"] = report("thresholds (7/40, 47/80): xi = (7/20, 61/80, 47/40, 127/80), 5/2 > 99/40, violation 1/40",
                                 xi == XI and xi_brute_ok and lhs == Fr(5, 2) and rhs == Fr(99, 40) and lhs - rhs == Fr(1, 40))

    mean_ok, fm, fat, ffit = folded_condition()
    folded = mean_ok and ffit and fm >= 0
    cross = TB.folded_order_ok(PIECES)
    res["folded_condition"] = {"holds": folded, "min_slack": str(fm), "argmin": str(fat), "tail_bounds_agrees": cross == folded}
    print(f"[INFO] folded condition law(U) <=cx (law(S - 1{{S > 1}}) + law(1{{S > 1}}))/2: "
          f"{'HOLDS' if folded else 'FAILS'} (min slack {fm} at y = {fat}; tail_bounds.folded_order_ok = {cross})",
          flush=True)
    checks["folded_evaluated"] = report("folded condition holds (exact), two implementations agree", ffit and folded and cross == folded,
                                        f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "worked_example")


if __name__ == "__main__":
    main()
