"""Readings of the potential inequalities: the tail-sum form of a step potential equals the grid inequality, and the
single-step inequality for symmetric three-point laws decides membership."""
import random
from fractions import Fraction as Fr

from dsxray.facet_test import v_of
from dsxray.mww_steps import boxconv, mww_bound, onestep_lhs


def test_tail_form_matches_grid_inequality():
    rng = random.Random(7)
    for _ in range(300):
        k = rng.randint(4, 12); r = rng.randint(1, k - 2)
        J = sorted(rng.sample(range(1, k), r)); sig = [Fr(j, k) for j in J]
        ext = sig + [Fr(1)]
        xi = [max(ext[i] + ext[m - 1 - i] for i in range(0, m) if i <= r and m - 1 - i <= r) for m in range(1, 2 * r + 1)]
        u = [sum(1 for j in J if j <= i) for i in range(k)]
        for t in range(2 * k - 1):
            s = Fr(t + 1, k)
            assert boxconv(sig, s) == sum(1 for x in xi if s >= x) == v_of(u)[t]


def test_symmetric_three_point_single_step():
    for bn in range(25, 100):
        b = Fr(bn, 100); L = 1 / b
        sigma = (1 - b) / 2 if 2 <= L <= 3 else b
        for f2 in [mww_bound(b) - Fr(1, 200), mww_bound(b), mww_bound(b) + Fr(1, 200)]:
            if not 0 <= f2 <= 1:
                continue
            f = (1 - f2) / 2
            ok = onestep_lhs({1 - b: f, Fr(1): f2, 1 + b: f}, sigma) <= 2 * (1 - sigma)
            assert ok == (f2 >= mww_bound(b))
