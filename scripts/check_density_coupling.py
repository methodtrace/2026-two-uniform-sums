#!/usr/bin/env python3
"""An explicit coupling of two uniforms whose sum has density 2 on [1/2, 3/4] u [5/4, 3/2].

Claim: let X be uniform on [0, 1] and, given X = x,
  Y = 3/4 - 2x                         for x in [0, 1/4),
  Y = 3/8 - x/2 or 9/8 - x/2 (prob. 1/2 each) for x in [1/4, 3/4),
  Y = 9/4 - 2x                         for x in [3/4, 1].
Then Y is uniform on [0, 1] and S = X + Y has density 2 on [1/2, 3/4] u [5/4, 3/2] and 0 elsewhere. S has mass
and mean one, S <=cx 2U, and S attains the single-step inequality at eta = 1/4 with equality:
P(S >= 1/2) + P(S >= 5/4) = 3/2 = 2(1 - 1/4). So this density-two law lies on the boundary and is attainable.

Checks, exact rational arithmetic: the push-forward densities of X, Y and S computed branch by branch (each branch is
affine in x, so its image density is (conditional probability)/|slope|); X and Y have density exactly 1 on [0, 1];
the density of S; mass and mean; the convex order S <=cx 2U by the exact minimum of the piecewise quadratic
(2 - x)^2/4 - E(S - x)^+ (endpoints and stationary points), cross-checked with tail_bounds.cx_below_2U; the
equality in the single-step inequality; and, as a second check of the marginal of Y, P(Y <= y) = y for
y = j/64, j = 0..64, integrated directly over x from the branches.
Seconds.
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import density_two as D2
from dsxray import tail_bounds as TB


def prob_y_le(y):
    """P(Y <= y) directly from the branches (integrate over x)."""
    tot = Fr(0)
    for lo, hi, p, b0, b1 in D2.COUPLING:
        # x in [lo, hi] with b0 + b1 x <= y; b1 < 0 so x >= (b0 - y)/(-b1)
        x0 = (b0 - y) / (-b1)
        a = max(lo, x0)
        if a < hi:
            tot += p * (hi - a)
    return tot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/density_coupling")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer()
    P = D2.coupling_profiles()
    marg = all(D2.profile_moment(P[c], 0) == 1 and D2.profile_moment(P[c], 1) == Fr(1, 2)
               and all(h == 1 for _, _, h in P[c]) and P[c][0][0] == 0 and P[c][-1][1] == 1 for c in ("X", "Y"))
    cond = all(sum(p for lo, hi, p, _, _ in D2.COUPLING if lo <= x < hi) == 1 for x in (Fr(0), Fr(1, 8), Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(7, 8)))
    grid = all(prob_y_le(Fr(j, 64)) == Fr(j, 64) for j in range(65))
    checks["uniform_marginals"] = report("X and Y have density 1 on [0, 1]; conditional probabilities sum to one",
                                         marg and cond and grid, f"({t()} s)")
    t = Timer()
    S = P["S"]
    dens = S == [(Fr(1, 2), Fr(3, 4), Fr(2)), (Fr(3, 4), Fr(5, 4), Fr(0)), (Fr(5, 4), Fr(3, 2), Fr(2))]
    mm = D2.profile_moment(S, 0) == 1 and D2.profile_moment(S, 1) == 1
    gap = D2.profile_min_convex_order_gap(S)
    cx = gap >= 0 and TB.cx_below_2U([(lo, hi, h) for lo, hi, h in S if h])
    eq = D2.profile_tail(S, Fr(1, 2)) + D2.profile_tail(S, Fr(5, 4)) == Fr(3, 2) == 2 * (1 - Fr(1, 4))
    res["sum_density"] = [[str(x) for x in row] for row in S]
    res["min_convex_order_gap"] = str(gap)
    checks["sum_law"] = report("S has density 2 on [1/2, 3/4] u [5/4, 3/2], mass and mean one, S <=cx 2U", dens and mm and cx,
                               f"min of (2-x)^2/4 - E(S-x)^+ = {gap} ({t()} s)")
    checks["single_step_equality"] = report("equality in the single step at eta = 1/4: P(S >= 1/2) + P(S >= 5/4) = 3/2", eq)
    finish({"args": vars(a), "results": res, "checks": checks}, out, "density_coupling")


if __name__ == "__main__":
    main()
