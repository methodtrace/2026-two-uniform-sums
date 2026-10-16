#!/usr/bin/env python3
"""Equally spaced staircases: the universal cap r + 1 (certificate and sharpness) and a perturbation bound.

Claims (proved analytically in the paper), for thresholds eta_j = 1 - (r - j + 1)d (j = 1..r), 0 < d < 1/r:
  (a) certificate: for 1/2 <= c <= r + 1, phi(s) = (s - t)^+/d with alpha = 2r(2c - 1)/(4cr + 1),
      t = 2 - (2r + alpha)d gives (1/2) int phi + c int (N - phi)^+ = c r (2r + 1)^2 d/(4cr + 1) - ((-t)^+)^2/(4d)
      and margin r(r + 1 - c)d/(4cr + 1) + ((-t)^+)^2/(4d) >= 0;
  (b) sharpness: for c > r + 1, with d = (4cr + 1)/(2cr(2r + 1)) and delta = 1/(2cr), density c on
      [xi_l, xi_l + delta] (the 2r jumps xi_l of N) is a law with mass and mean one, S <=cx 2U, violating the
      staircase inequality by exactly (c - r - 1)/(2c(2r + 1));
  (c) perturbation: for r >= 2 at cap 2, the same phi certifies every threshold vector within sup-distance
      m/(10r) of the equally spaced one, m the margin in (a), because the margin drops by at most 10 r epsilon.
Checks, exact rational arithmetic:
  1. (a) for r <= --rmax, d = x/r with x in {1/50, 1/5, 1/2, 4/5, 49/50}, c in {1/2, 3/4, 1, 3/2, 2, r + 1,
     (r + 3)/2, r + 2}: integrated value equals the closed form, and the margin is >= 0 when c <= r + 1;
  2. (b) for r <= --rmax and c in {r + 101/100, r + 2, 2r + 3}: interval positions, mass, mean, the convex order
     (exact minimum of the piecewise quadratic (2 - x)^2/4 - E(S - x)^+ >= 0; also by tail_bounds.cx_below_2U) and
     the violation;
  3. (c) for 2 <= r <= 6 and the five spacings: 40 random perturbations with |epsilon_j| <= m/(10r) each keep a
     nonnegative margin, and every margin is at least m - 10 r max|epsilon_j|.
--rmax 24: 960 certificate cases and 72 sharpness laws.
Default: r <= 8 (about 1 s). --full: r <= 24 (about 5 s).
"""
import argparse
import random
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import density_two as D2
from dsxray import tail_bounds as TB


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/equally_spaced_density")
    ap.add_argument("--rmax", type=int, default=8)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 24)
    out = out_dir(a.out); res = {}; checks = {}
    spacings = [Fr(1, 50), Fr(1, 5), Fr(1, 2), Fr(4, 5), Fr(49, 50)]

    t = Timer(); n = 0; bad = []
    for r in range(1, a.rmax + 1):
        for c in (Fr(1, 2), Fr(3, 4), Fr(1), Fr(3, 2), Fr(2), Fr(r + 1), Fr(r + 3, 2), Fr(r + 2)):
            for x in spacings:
                m, mc, J, Jc = D2.equally_spaced(r, x / r, c); n += 1
                if not (m == mc and J == Jc and (c > r + 1 or m >= 0)):
                    bad.append([r, str(c), str(x)])
    res["certificate"] = {"rmax": a.rmax, "cases": n, "failures": bad[:10]}
    checks["certificate"] = report(f"equally spaced certificate, r <= {a.rmax}: closed forms, margin >= 0 for c <= r + 1",
                                   not bad and (a.rmax != 24 or n == 960), f"{n} cases ({t()} s)")

    t = Timer(); laws = []; bad = []
    for r in range(1, a.rmax + 1):
        for c in (Fr(r) + Fr(101, 100), Fr(r + 2), Fr(2 * r + 3)):
            eta, iv, viol, gap = D2.equally_spaced_sharpness(r, c)
            pieces = [(lo, hi, c) for lo, hi in iv]
            good = viol == (c - r - 1) / (2 * c * (2 * r + 1)) and viol > 0 and gap >= 0
            if r <= 8:
                good &= TB.cx_below_2U(pieces)
            if not good:
                bad.append([r, str(c)])
            laws.append({"r": r, "c": str(c), "violation": str(viol), "min_convex_order_gap": str(gap)})
    res["sharpness"] = {"laws": len(laws), "failures": bad, "examples": laws[:6]}
    checks["sharpness"] = report(f"sharpness laws for c > r + 1, r <= {a.rmax}: S <=cx 2U, violation (c-r-1)/(2c(2r+1))",
                                 not bad and (a.rmax != 24 or len(laws) == 72), f"{len(laws)} laws ({t()} s)")

    t = Timer(); rng = random.Random(11); n = 0; ok = True
    for r in range(2, min(6, a.rmax) + 1):
        for x in spacings:
            d = x / r
            m = D2.equally_spaced(r, d, 2)[0]
            bound = m / (10 * r)
            for _ in range(40):
                eps = [bound * Fr(rng.randint(-1000, 1000), 1000) for _ in range(r)]
                mm = D2.perturbed_margin(r, d, eps); n += 1
                ok &= mm >= 0 and mm >= m - 10 * r * max(abs(e) for e in eps)
    res["perturbation"] = {"cases": n}
    checks["perturbation"] = report("perturbed thresholds within m/(10r): margin >= m - 10 r eps >= 0", ok,
                                    f"{n} cases ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "equally_spaced_density")


if __name__ == "__main__":
    main()
