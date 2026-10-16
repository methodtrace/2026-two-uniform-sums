"""Lattice reduction for sums of two U[0,1] variables.

Let S = X + Y with X, Y ~ U[0,1], and suppose S takes values in a + (1/k)Z.
Write c = frac(k a). Then kS - c = I + J + B where, with probability c (event R <= c), B = 0 and
(I, J) is a coupling of two uniforms on {0..k-1}; with probability 1 - c, B = 1 and (I, J) is another such
coupling. Conversely every such mixture is attained. Hence the law of kS - c ranges exactly over
  c * P_k + (1 - c) * (P_k shifted by 1),   P_k = aggregation set of two uniforms on {0..k-1}.
"""
import math

import numpy as np
from scipy.optimize import linprog


def _nint(x):
    """Nearest integer (halves upward)."""
    return math.floor(x + 0.5)


def _lp_tri(k, c, target_support, fixed, objective):
    """Optimize a linear objective of the law f of kS - c over the reduction, with support restricted."""
    nv = 2 * k * k
    A, b = [], []
    for part, w in ((0, c), (1, 1 - c)):
        off = part * k * k
        for i in range(k):
            r = np.zeros(nv); r[off + i * k:off + (i + 1) * k] = 1; A.append(r); b.append(w / k)
        for j in range(k):
            r = np.zeros(nv); r[off + j:off + k * k:k] = 1; A.append(r); b.append(w / k)
    L = 2 * k  # values 0..2k-1
    F = np.zeros((L, nv))
    for part in (0, 1):
        off = part * k * k
        for i in range(k):
            for j in range(k):
                F[i + j + part, off + i * k + j] += 1
    for t in range(L):
        if t not in target_support:
            A.append(F[t]); b.append(0.0)
    for t, val in fixed.items():
        A.append(F[t]); b.append(val)
    res = linprog(objective @ F, A_eq=np.array(A), b_eq=np.array(b), bounds=(0, None), method="highs")
    return res


def tri_atomic_f2_range(a, b):
    """Range of the middle mass f2 for S on {a-b, a, a+b} with mean 1, S in D2(U,U).
    Requires 1/b = k integer (lattice spacing 1/k)."""
    k = _nint(1 / b); assert abs(k * b - 1) < 1e-12
    c = (k * a) % 1.0
    if abs(c - 1) < 1e-12: c = 0.0
    base = k * (a - b) - c                      # integer position of atom a-b in kS - c
    assert abs(base - _nint(base)) < 1e-9; base = _nint(base)
    sup = {base, base + 1, base + 2}
    if min(sup) < 0 or max(sup) > 2 * k - 1:
        return None
    # mean constraint: E[S]=1 automatically holds for any coupling (E[X+Y]=1); so only support matters.
    out = []
    for sign in (1, -1):
        obj = np.zeros(2 * k); obj[base + 1] = sign
        res = _lp_tri(k, c, sup, {}, obj)
        if res.status != 0:
            return None
        out.append(sign * res.fun)
    return out[0], out[1]


def non_integrity(x):
    return min(math.ceil(x) / x - 1, 1 - math.floor(x) / x)


def mww_theorem4_as_printed(a, b):
    """Lower bound on f2 exactly AS PRINTED in Mao-Wang-Wang (2019) Theorem 4 (arXiv v2 Theorem 3.5).
    Case (iii) is misprinted (a + b - 1/2); kept only to document the misprint. Do not use for decisions."""
    if abs(a - 1) < 1e-12:
        return non_integrity(1 / (2 * b))
    x = 1 / (2 * b)
    if abs(x - _nint(x)) < 1e-9:
        return 0.0
    if abs((x - 0.5) - _nint(x - 0.5)) < 1e-9:
        return a + b - 0.5
    return None


def mww_theorem4_corrected(a, b):
    """Corrected lower bound on f2 (case (iii): a + b - 1, as in their proof and Lemma 4(iii)(c) with closed endpoint)."""
    if abs(a - 1) < 1e-12:
        return non_integrity(1 / (2 * b))
    x = 1 / (2 * b)
    if abs(x - _nint(x)) < 1e-9:
        return 0.0
    if abs((x - 0.5) - _nint(x - 0.5)) < 1e-9:
        return a + b - 1
    return None


def w3_f2_interval(a, b, tol=1e-9):
    """Three-point corollary of the window theorem, including the a + b = 1 endpoint: exact range of the middle mass f2 for mean-1 laws
    in D2U supported on {a - b, a, a + b}, 0 < b <= a <= 1. Returns (lo, hi), or None if no such law exists."""
    if a + b < 1 - tol:
        return None                                  # support lies strictly below 1: no mean-1 law
    if abs(a - 1) < tol:                             # 1 is the middle atom
        L = 1 / b
        d = abs(L - 2 * _nint(L / 2))                # dist(L, 2Z)
        return (d / L, 1.0)
    if abs(a + b - 1) < tol:                         # 1 is the top atom: the only mean-1 law is delta_1
        return (0.0, 0.0)
    k = 1 / b                                        # 1 strictly inside (a, a + b): need 1/b integer
    if abs(k - _nint(k)) > tol:
        return None
    k = _nint(k); C = (a + b - 1) / b
    return (C * (k % 2) / k, C)
