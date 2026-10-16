"""Corollary 3.2 and the lemmas of its proof on small instances (exact)."""
import itertools
import random
from fractions import Fraction as Fr

from dsxray import small_atoms as SA


def test_statement_small_k():
    for k, nv in ((2, 2), (3, 7), (4, 44)):
        assert all(SA.satisfies_cx(SA.xray(s), k) for s in itertools.permutations(range(k)))
        V = SA.vertices_exact(k)
        assert len(V) == nv
        for x in V:
            Z = SA.find_Z(list(x), k)
            assert Z is not None and SA.check_Z(Z, list(x), k)


def test_negative_control():
    V = SA.vertices_exact(3, cap=None)
    bad = [x for x in V if SA.find_Z(list(x), 3) is None]
    assert len(V) == 5 and bad == [(0, Fr(3, 2), 0, Fr(3, 2), 0)]
    for x in V:
        assert (SA.separating_potential(list(x), 3) is None) == (list(x) not in [list(b) for b in bad])


def test_find_Z_fast_is_exact():
    rng = random.Random(1)
    for k in (3, 5, 6):
        x = SA.random_Q_point(k, rng)
        Z = SA.find_Z_fast(x, k)
        assert SA.check_Z(Z, x, k) and all(isinstance(v, Fr) for row in Z for v in row)


def test_minorant_lemma_and_chain():
    rng = random.Random(2)
    for k in range(1, 16):
        for _ in range(20):
            SA.minorant_lemma(SA.random_u(k, rng))
    for k in range(1, 9):
        for mask in range(1 << (k - 1)):
            J = [j for j in range(1, k) if mask >> (j - 1) & 1]
            for s in (1, -1):
                SA.minorant_lemma([Fr(s * sum(1 for j in J if i >= j)) for i in range(k)])
    for k in (3, 5):
        x = SA.random_Q_point(k, rng)
        assert SA.in_Q(x, k)
        assert all(SA.chain_check(x, SA.random_u(k, rng)) for _ in range(30))
