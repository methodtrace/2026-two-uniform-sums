"""The coupling linear programs reproduce published results."""
import numpy as np
import pytest

from dsxray.couplings import (Marginals, attainable, comonotone_sum, countermonotone_sum, makarov_cdf_bounds,
                              reach_bounds, stop_loss)

rng = np.random.default_rng(23)


def rand_pmf(k, denom=12):
    w = rng.integers(0, denom, size=k).astype(float) + (rng.random(k) < 0.8)
    return tuple(w / w.sum())


@pytest.mark.parametrize("trial", range(100))
def test_lu_ding_dasgupta_tau(trial):
    """Lu, Ding and Dasgupta (2018), Prop. 1: sharp bounds on P(Y1 >= Y0) for ordinal outcomes."""
    J = int(rng.integers(2, 7))
    a, b = rand_pmf(J), rand_pmf(J)
    lo, hi = reach_bounds(Marginals(a, tuple(b[::-1])), J - 1)
    ge = lambda v, j: sum(v[j:])
    delta = [ge(a, j) - ge(b, j) for j in range(J)]
    assert abs(lo - max(b[j] + delta[j] for j in range(J))) < 1e-8 and abs(hi - (1 + min(delta))) < 1e-8


@pytest.mark.parametrize("trial", range(100))
def test_makarov_lattice(trial):
    """Pointwise bounds on P(S <= t) equal the lattice Makarov bounds."""
    M = Marginals(rand_pmf(int(rng.integers(2, 6))), rand_pmf(int(rng.integers(2, 6))))
    lo_mk, hi_mk = makarov_cdf_bounds(M)
    for t in range(M.N):
        rmin, rmax = reach_bounds(M, t + 1)
        assert abs((1 - rmax) - lo_mk[t]) < 1e-8 and abs((1 - rmin) - hi_mk[t]) < 1e-8


def test_mao_wang_example():
    """X, Y uniform on {0,1,2}; S uniform on {1,3} lies between the comonotone and countermonotone sums in convex order
    but is not attainable."""
    M = Marginals((1 / 3, 1 / 3, 1 / 3), (1 / 3, 1 / 3, 1 / 3))
    s = np.array([0, .5, 0, .5, 0])
    assert np.all(stop_loss(s) <= stop_loss(comonotone_sum(M)) + 1e-12)
    assert np.all(stop_loss(s) >= stop_loss(countermonotone_sum(M)) - 1e-12)
    assert not attainable(M, s)


def test_bernoulli():
    lo, hi = reach_bounds(Marginals((0.4, 0.6), (0.7, 0.3)), 2)
    assert abs(lo) < 1e-12 and abs(hi - 0.3) < 1e-12


u = lambda k: tuple([1 / k] * k)


def test_two_unattainable_laws_on_four_points():
    assert not attainable(Marginals(u(4), u(4)), np.array([1, 1, 0, 0, 0, 1, 1]) / 4)
    assert not attainable(Marginals(u(4), u(4)), np.array([1, 0, 0, 1.5, 0, 1.5, 0]) / 4)


def test_attainable_law_spacing_three_sevenths():
    assert attainable(Marginals(u(7), u(7)), np.array([1, 0, 0, 2, 0, 0, 0, 0, 0, 4, 0, 0, 0]) / 7)
