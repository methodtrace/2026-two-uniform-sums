"""Window theorem: vertex rule of the window lemma, the four-atom corollary against the generic fibered LP, the
three-point interval (including the a + b = 1 endpoint), and the five-point window vertex counts."""
import math

import numpy as np
import pytest

from dsxray.fibered import optimize_on_support, support_function_on_window
from dsxray.two_uniforms import w3_f2_interval
from dsxray.windows import (candidate_vertices, five_point_window_vertex_count, four_atom_prediction,
                            support_of_prediction, vertices_exact)


@pytest.mark.parametrize("m", range(1, 61))
def test_vertex_rule(m):
    assert set(vertices_exact(m)) <= candidate_vertices(m)


CASES = [(1 / k, 1 - (1 / k) * t / 4) for k in range(2, 8) for t in range(1, 4 * k) if 1 - (1 / k) * t / 4 > 0] + \
        [(1 / L, 1 - j / L) for L in (2.5, math.e, math.pi, 4.7) for j in (1, 2)]


@pytest.mark.parametrize("b,a", CASES)
def test_four_atom_corollary(b, a):
    if a + 3 * b >= 2:
        pytest.skip("window leaves (0, 2)")
    pred = four_atom_prediction(b, a)
    rng = np.random.default_rng(int(1e6 * b + 1e3 * a))
    for _ in range(4):
        w = rng.normal(size=4)
        lp = support_function_on_window(b, a, [0, 1, 2, 3], w)
        assert (lp is None) == (pred is None)
        if lp is None:
            return
        assert abs(lp - support_of_prediction(pred, w)) < 1e-7


@pytest.mark.parametrize("a,b", [(3 / 5, 2 / 5), (0.7, 0.3), (1, 1 / 5), (1, 2 / 7), (13 / 18, 1 / 3), (5 / 6, 1 / 4),
                                 (9 / 10, 1 / 5), (0.9, 0.3), (1, 1 / math.pi)])
def test_three_point_interval(a, b):
    pred = w3_f2_interval(a, b)
    lo = optimize_on_support(b, a - b, [0, 1, 2], 1, 1)
    if pred is None:
        assert lo is None
        return
    hi = abs(optimize_on_support(b, a - b, [0, 1, 2], 1, -1))
    assert abs(pred[0] - lo) < 1e-7 and abs(pred[1] - hi) < 1e-7


@pytest.mark.parametrize("b", [0.2, 0.25, 0.3, 1 / math.pi])
def test_window_starting_at_one(b):
    pred = four_atom_prediction(b, 1.0); w = np.array([0.3, -1.1, 0.7, 0.2])
    assert pred is not None
    assert abs(support_of_prediction(pred, w) - 0.3) < 1e-9
    assert abs(support_function_on_window(b, 1.0, [0, 1, 2, 3], w) - 0.3) < 1e-9


def test_five_point_window_counts():
    assert [five_point_window_vertex_count(m) for m in range(1, 9)] == [1, 2, 5, 5, 9, 6, 9, 9]
