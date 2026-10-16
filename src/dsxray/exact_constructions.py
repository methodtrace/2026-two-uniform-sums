"""Exact (rational) verification of explicit constructions from the fibered reduction, U[0,1] marginals.

A construction is: spacing b = 1/L (L rational), offset a, and for each fiber (an interval of the fractional
part r of X/b) a permutation sigma of {0..n-1} (n = number of integer parts available in that fiber).
X' = I + R with I uniform given R; Y' = sigma(I) + phi(R), phi(r) = frac(c - r); S = b (X' + Y').
We verify exactly: (1) fiber sizes match, n(r) = n(phi(r)); (2) Y' has density 1/L on [0, L) (piecewise check);
(3) the law of S.
"""
from fractions import Fraction as Fr
import math

def verify(L, a, perms):
    b = 1 / L
    c = (a / b) % 1
    theta = L - math.floor(L)
    pts = sorted(set([Fr(0), Fr(1), theta % 1, (c - theta) % 1, c % 1]))
    fibers = [(lo, hi) for lo, hi in zip(pts, pts[1:]) if hi > lo]
    assert len(fibers) == len(perms), (fibers, perms)
    n_of = lambda r: math.ceil(L - r)
    law = {}
    ypieces = []                                  # (start, end, density) pieces of Y'
    for (lo, hi), sig in zip(fibers, perms):
        mid = (lo + hi) / 2
        n = n_of(mid); n2 = n_of((c - mid) % 1)
        assert n == n2, "fractional-part laws do not match"
        assert sorted(sig) == list(range(n))
        B = 1 if mid > c else 0
        w = (hi - lo) * Fr(n) / L                 # probability of this fiber
        for i in range(n):
            s_val = b * (i + sig[i] + c + B)
            law[s_val] = law.get(s_val, Fr(0)) + w / n
            # Y' = sig[i] + phi(r), r in (lo, hi): phi maps (lo,hi) onto an interval of length hi-lo
            p1, p2 = (c - hi) % 1, (c - lo) % 1
            if p2 == 0: p2 = Fr(1)
            if p1 < p2:
                ypieces.append((sig[i] + p1, sig[i] + p2))
            else:                                  # wraps around: split
                ypieces.append((sig[i] + p1, sig[i] + 1)); ypieces.append((Fr(sig[i]), sig[i] + p2))
    # each piece carries density (1/L) * (n/n) = 1/L per unit length (R density n/L, I|R uniform 1/n)
    ypieces.sort()
    cover = Fr(0)
    for s, e in ypieces:
        assert s == cover, ("gap or overlap in Y support", s, cover)
        cover = e
    assert cover == L, ("Y support does not reach L", cover)
    total = sum(law.values()); mean = sum(v * p for v, p in law.items())
    assert total == 1 and mean == 1
    return dict(sorted(law.items()))
