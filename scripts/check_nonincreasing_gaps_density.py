#!/usr/bin/env python3
"""Staircases with nonincreasing threshold gaps: the explicit certificate at every depth and cap.

Claim (proved analytically in the paper): let 0 < eta_0 < ... < eta_{r-1} < 1, eta_r = 1, g_j = eta_j - eta_{j-1}
(1 <= j <= r) with g_1 >= g_2 >= ... >= g_r (the last gap 1 - eta_{r-1} included), D = 1 - eta_0, and
1/2 <= c <= 1 + D/g_1. With alpha = 2D(2c - 1)/(4cD + g_1), t = 2 eta_0 - alpha g_1 and
  phi(s) = (s - t)^+/g_1 + sum_{j=1}^{r-1} (1/g_{j+1} - 1/g_j)(s - 2 eta_j)^+,
the margin 2 int u - [(1/2) int phi + c int (N - phi)^+] equals D(D + (1 - c) g_1)/(4cD + g_1) + ((-t)^+)^2/(4 g_1),
which is >= 0 for c <= 1 + D/g_1, and > 0 at c = 2 when r >= 2. At r = 1 the certificate is that of
Proposition 5.2 (t = (8 eta - 2)/3).

Checks, exact rational arithmetic: for every nonincreasing gap pattern with weights in {1, 2, 3, 4} and depth
r <= --rmax, scaled so that D in {1/10, 1/2, 9/10}, and caps c in {1/2, 1, 2, 1 + D/g_1, 3}:
  1. the jumps of N are the balanced pair sums 2 eta_j, eta_j + eta_{j+1};
  2. the integrated margin equals the closed form;
  3. the margin is >= 0 whenever c <= 1 + D/g_1, and > 0 at c = 2 for r >= 2;
  4. at r = 1 and c = 2, the certificate coincides with Proposition 5.2 (density_certificates.single_step_J).
--rmax 8: 494 gap patterns, 7,410 cases.
Default: r <= 5 (about 1 s). --full: r <= 8 (about 10 s).
"""
import argparse
from fractions import Fraction as Fr
from itertools import combinations_with_replacement

from common import Timer, finish, out_dir, report

from dsxray import density_certificates as DC
from dsxray import density_two as D2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/nonincreasing_gaps_density")
    ap.add_argument("--rmax", type=int, default=5)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.rmax = max(a.rmax, 8)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); n = 0; patterns = 0; bad = []; strict = 0
    for r in range(1, a.rmax + 1):
        for seq in combinations_with_replacement(range(1, 5), r):
            patterns += 1
            w = list(reversed(seq)); tot = sum(w)
            for scale in (Fr(1, 10), Fr(1, 2), Fr(9, 10)):
                gaps = [scale * Fr(x, tot) for x in w]
                cap = 1 + sum(gaps) / gaps[0]
                for c in (Fr(1, 2), Fr(1), Fr(2), cap, Fr(3)):
                    m, closed, _ = D2.nonincreasing_gaps(gaps, c); n += 1
                    good = m == closed and (c > cap or m >= 0) and not (c == 2 and r >= 2 and m <= 0)
                    strict += c == 2 and r >= 2
                    if not good:
                        bad.append([[str(g) for g in gaps], str(c), str(m)])
    res["sweep"] = {"rmax": a.rmax, "gap_patterns": patterns, "cases": n, "cap_two_depth_two_or_more": strict,
                    "failures": bad[:10]}
    ok = not bad and (a.rmax != 8 or (patterns, n) == (494, 7410))
    checks["nonincreasing_gaps"] = report(f"nonincreasing gaps, depth <= {a.rmax}: margin equals the closed form, "
                                          "nonnegative for c <= 1 + D/g_1, positive at c = 2 for r >= 2", ok,
                                          f"{patterns} gap patterns, {n} cases ({t()} s)")

    t = Timer(); ok = True
    for j in range(1, 60):
        eta = Fr(j, 60)
        m, closed, tt = D2.nonincreasing_gaps([1 - eta], 2)
        ok &= tt == (8 * eta - 2) / 3 and m == 2 * (1 - eta) - DC.single_step_J(eta)
    checks["single_step_case"] = report("r = 1: the certificate is that of Proposition 5.2", ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "nonincreasing_gaps_density")


if __name__ == "__main__":
    main()
