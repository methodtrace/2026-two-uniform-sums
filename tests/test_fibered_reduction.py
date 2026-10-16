"""The fibered reduction recovers the bi-atomic theorem of Mao, Wang and Wang: {a, a+b}, mean 1, attainable iff 1/b is
an integer; and the three-point minimum for non-integer 1/(2b)."""
import math

import numpy as np
import pytest

from dsxray.fibered import feasible_support, fibers, optimize_on_support

rng = np.random.default_rng(3)
CASES = []
for b in [1 / k for k in range(1, 9)] + [2 / 5, 3 / 7, 0.37, 1 / math.pi, 2 / 3, 0.61803, 0.45, 1 / 3.5]:
    for _ in range(8):
        lo = max(0.0, 1 - b); a = lo + (1 - lo) * rng.random()
        if a < 1 < a + b <= 2:
            CASES.append((b, a))


@pytest.mark.parametrize("b,a", CASES)
def test_biatomic(b, a):
    assert feasible_support(b, a, [0, 1]) == (abs(1 / b - math.floor(1 / b + 0.5)) < 1e-9)


@pytest.mark.parametrize("b", [2 / 5, 3 / 7, 0.37, 1 / math.pi, 2 / 3])
def test_support_condition(b):
    L = 1 / b; on = 1 - b * math.floor(1 / b); off = on + 0.3 * b
    ok = lambda a: all(n1 == n2 for _, n1, n2, _ in fibers(L, (a / b) % 1.0))
    assert ok(on) and not ok(off)


def test_nonconstant_example():
    assert feasible_support(0.4, 0.6, [0, 2, 3])


@pytest.mark.parametrize("b", [1 / 2.5, 1 / 3.3, 1 / math.e, 0.3, 0.45, 1 / 4.7, 1 / 5.5, 0.21, 1 / math.pi, 2 / 7, 0.61803])
def test_three_point_minimum(b):
    x = 1 / (2 * b); pred = min(math.ceil(x) / x - 1, 1 - math.floor(x) / x)
    assert abs(optimize_on_support(b, 1 - b, [0, 1, 2], 1, 1) - pred) < 1e-7
    assert abs(optimize_on_support(b, 1 - b, [0, 1, 2], 1, -1) - 1) < 1e-7
