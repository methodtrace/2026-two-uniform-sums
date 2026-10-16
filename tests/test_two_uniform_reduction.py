"""Lattice reduction of U + V reproduces the tri-atomic theorem of Mao, Wang and Wang, with the case (iii) bound
a + b - 1 of their proof (their statement prints a + b - 1/2), and the corrected vertex normalization."""
from fractions import Fraction as Fr

import pytest

from dsxray.exact_constructions import verify
from dsxray.two_uniforms import (mww_theorem4_as_printed, mww_theorem4_corrected, non_integrity,
                                 tri_atomic_f2_range)

CASES = [(k, (1 / k) + (1 - 1 / k) * t / (4 * k)) for k in range(2, 9) for t in range(0, 4 * k + 1)]
CASES = [(k, a) for k, a in CASES if 1 - (1 - a) * k > 1e-9]


@pytest.mark.parametrize("k,a", CASES)
def test_tri_atomic(k, a):
    b = 1 / k; x = 1 / (2 * b)
    if abs(a - 1) < 1e-12:
        lo = non_integrity(x)
    elif abs(x - int(x + 0.5)) < 1e-9:
        lo = 0.0
    else:
        lo = a + b - 1
    r = tri_atomic_f2_range(a, b)
    assert r is not None and abs(r[0] - max(lo, 0)) < 1e-7 and abs(r[1] - (1 - (1 - a) / b)) < 1e-7


def test_printed_statement_fails_on_an_attainable_law():
    lo, hi = tri_atomic_f2_range(13 / 18, 1 / 3)
    assert lo <= 1 / 18 + 1e-9 <= hi + 1e-9
    assert abs(mww_theorem4_as_printed(13 / 18, 1 / 3) - 5 / 9) < 1e-12
    assert abs(mww_theorem4_corrected(13 / 18, 1 / 3) - lo) < 1e-9 and abs(lo - 1 / 18) < 1e-9


def test_vertex_normalization():
    a, b = Fr(13, 18), Fr(1, 3); C = (a + b - 1) / b
    printed = [C * (1 - b) / 2, b * C, (2 - C * (1 - b)) / 2]
    corrected = [C * (1 - b) / 2, b * C, (2 - C * (1 + b)) / 2]
    assert sum(printed) == 1 + b * C and sum(corrected) == 1
    assert corrected == [Fr(1, 18), Fr(1, 18), Fr(8, 9)]


def test_exact_construction_of_the_law():
    law = verify(Fr(3), Fr(13, 18), [(1, 2, 0), (2, 1, 0)])
    assert law == {Fr(7, 18): Fr(1, 18), Fr(13, 18): Fr(1, 18), Fr(19, 18): Fr(8, 9)}


def test_exact_construction_spacing_two_fifths():
    law = verify(Fr(5, 2), Fr(1, 5), [(1, 0, 2), (0, 1)])
    assert law == {Fr(3, 5): Fr(3, 5), Fr(7, 5): Fr(1, 5), Fr(9, 5): Fr(1, 5)}


def test_printed_parameter_leaves_the_simplex():
    """Appendix A, item 2: the printed c = a/b + 1 - 1/(2b) is 2.2 at b = 1/4, a = 0.8; the corrected c = (a+b-1)/b
    lies in (0, 1)."""
    a, b = Fr(4, 5), Fr(1, 4)
    assert a / b + 1 - 1 / (2 * b) == Fr(11, 5)
    assert 0 < (a + b - 1) / b < 1
