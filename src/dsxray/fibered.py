"""Fibered reduction, specialized to X, Y ~ U[0,1] with arbitrary lattice spacing b.

Scale by 1/b: X' = X/b ~ U[0, L], L = 1/b. I = floor(X'), R = frac(X'). n(r) = #{i >= 0: i + r < L} = ceil(L - r).
Given R = r, I is uniform on {0..n(r)-1}; R has density n(r)/L. If S' = S/b lies in A + c + Z (c = frac(a/b)),
then R' = phi(R) = frac(c - R), R + R' = c + 1{R > c}, and conditional on R = r, (I, J) couples
Unif{0..n(r)-1} with Unif{0..n(phi(r))-1}; a joint law exists only if n(r) = n(phi(r)) a.e.
Fibers are the intervals between breakpoints {0, theta, frac(c - theta), c, 1}, theta = frac(L).
"""
import math
import numpy as np
from scipy.optimize import linprog


def _nint(x):
    """Nearest integer (halves upward)."""
    return math.floor(x + 0.5)


def n_of(r, L):
    return math.ceil(L - r - 1e-15)


def fibers(L, c):
    theta = L - math.floor(L)
    pts = sorted(set([0.0, 1.0, theta % 1.0, (c - theta) % 1.0, c % 1.0]))
    out = []
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi - lo < 1e-13:
            continue
        mid = 0.5 * (lo + hi)
        n1 = n_of(mid, L); n2 = n_of((c - mid) % 1.0, L)
        out.append((hi - lo, n1, n2, 1 if mid > c else 0))
    return out


def feasible_support(b, a, support_offsets):
    """Is there a coupling of U[0,1], U[0,1] with S in {a + b*o : o in support_offsets}?"""
    L = 1.0 / b
    c = (a / b) % 1.0
    A = _nint(a / b - c)                     # a/b = A + c
    fb = fibers(L, c)
    if any(n1 != n2 for _, n1, n2, _ in fb):
        return False                         # fractional-part laws cannot be matched: no lattice-supported sum
    allowed = {A + o for o in support_offsets}   # allowed values of I + J + B
    var = []; A_eq = []; b_eq = []
    idx = 0; blocks = []
    for (length, n, _, B) in fb:
        w = length * n / L
        blocks.append((idx, n, B, w)); idx += n * n
    nv = idx
    for (off, n, B, w) in blocks:
        for i in range(n):
            r = np.zeros(nv); r[off + i * n: off + (i + 1) * n] = 1; A_eq.append(r); b_eq.append(w / n)
            r = np.zeros(nv); r[off + i: off + n * n: n] = 1; A_eq.append(r); b_eq.append(w / n)
    ub = []
    for (off, n, B, w) in blocks:
        for i in range(n):
            for j in range(n):
                ub.append((0, 0) if (i + j + B) not in allowed else (0, None))
    res = linprog(np.zeros(nv), A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=ub, method="highs")
    return res.status == 0


def optimize_on_support(b, a, support_offsets, target_offset, sense=1):
    """min (sense=1) or max (sense=-1) of P(S = a + b*target_offset) over couplings of U[0,1], U[0,1]
    with S in {a + b*o}. Returns None if infeasible."""
    L = 1.0 / b
    c = (a / b) % 1.0
    if abs(c - 1) < 1e-12: c = 0.0
    A = _nint(a / b - c)
    fb = fibers(L, c)
    if any(n1 != n2 for _, n1, n2, _ in fb):
        return None
    allowed = {A + o for o in support_offsets}; tgt = A + target_offset
    blocks = []; idx = 0
    for (length, n, _, B) in fb:
        blocks.append((idx, n, B, length * n / L)); idx += n * n
    A_eq = []; b_eq = []; cost = np.zeros(idx); bounds = []
    for (off, n, B, w) in blocks:
        for i in range(n):
            r = np.zeros(idx); r[off + i * n: off + (i + 1) * n] = 1; A_eq.append(r); b_eq.append(w / n)
            r = np.zeros(idx); r[off + i: off + n * n: n] = 1; A_eq.append(r); b_eq.append(w / n)
        for i in range(n):
            for j in range(n):
                v = i + j + B
                bounds.append((0, None) if v in allowed else (0, 0))
                if v == tgt: cost[off + i * n + j] = sense
    res = linprog(cost, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    return None if res.status != 0 else sense * res.fun


def support_function_on_window(b, a, window, w):
    """max over couplings of U[0,1], U[0,1] with S in {a + b*o : o in window} of sum_o w[o] P(S = a + b o).
    Returns None if infeasible."""
    L = 1.0 / b
    c = (a / b) % 1.0
    if abs(c - 1) < 1e-12: c = 0.0
    A = _nint(a / b - c)
    fb = fibers(L, c)
    if any(n1 != n2 for _, n1, n2, _ in fb):
        return None
    allowed = {A + o: idx for idx, o in enumerate(window)}
    blocks = []; nv = 0
    for (length, n, _, B) in fb:
        blocks.append((nv, n, B, length * n / L)); nv += n * n
    A_eq = []; b_eq = []; cost = np.zeros(nv); bounds = []
    for (off, n, B, wt) in blocks:
        for i in range(n):
            r = np.zeros(nv); r[off + i * n: off + (i + 1) * n] = 1; A_eq.append(r); b_eq.append(wt / n)
            r = np.zeros(nv); r[off + i: off + n * n: n] = 1; A_eq.append(r); b_eq.append(wt / n)
        for i in range(n):
            for j in range(n):
                v = i + j + B
                if v in allowed:
                    bounds.append((0, None)); cost[off + i * n + j] = -w[allowed[v]]
                else:
                    bounds.append((0, 0))
    res = linprog(cost, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
    return None if res.status != 0 else -res.fun
