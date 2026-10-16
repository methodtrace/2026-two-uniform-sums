"""Tail-sum inequalities with repeated thresholds (Lemma 4.1) in exact rational arithmetic; a second, independent
implementation of the staircase quantities of tail_bounds.py.

A staircase with thresholds 0 < eta_1 <= ... <= eta_r <= 1 is u(x) = #{q : eta_q <= x} on [0, 1].
  xi(eta)              the thresholds xi_1..xi_2r of Lemma 4.1 (closed formula)
  box_formula(eta, s)  (u box u)(s) = #{l : xi_l <= s} (Lemma 4.1)
  box_brute(eta, s)    (u box u)(s) by direct minimization over x (independent of the formula)
  violation(eta, law)  sum_l P(S >= xi_l) - 2 sum_q (1 - eta_q) for a finite law [(atom, mass), ...]
  stair_lp_bound(law, r)  the exact supremum of the violation over staircases with exactly r thresholds in (0, 1]:
                       for every nondecreasing count vector c (values of u box u at the atoms) the thresholds with
                       (u box u)(s_m) >= c_m form a polyhedron R_c, and the supremum is the maximum over c of
                       sum_m f_m c_m - 2r + 2 max{sum eta : eta in closure(R_c)}, each maximum found by enumerating the
                       vertices of closure(R_c) exactly (the reduction in the proof of Proposition S4.1).
"""
from fractions import Fraction as Fr
from itertools import combinations, combinations_with_replacement


def u_of(eta, x):
    return sum(1 for e in eta if e <= x)


def box_brute(eta, s):
    lo, hi = max(Fr(0), s - 1), min(Fr(1), s)
    pts = {lo, hi}
    for e in eta:
        for p in (e, s - e):
            if lo <= p <= hi:
                pts.add(p)
    pts = sorted(pts)
    cand = list(pts) + [(a + b) / 2 for a, b in zip(pts, pts[1:])]
    return min(u_of(eta, x) + u_of(eta, s - x) for x in cand)


def neg_box_brute(eta, s):
    """((-u) box (-u))(s) by direct minimization (used for the reflection lemma)."""
    lo, hi = max(Fr(0), s - 1), min(Fr(1), s)
    pts = {lo, hi}
    for e in eta:
        for p in (e, s - e):
            if lo <= p <= hi:
                pts.add(p)
    pts = sorted(pts)
    cand = pts + [(a + b) / 2 for a, b in zip(pts, pts[1:])]
    return -max(u_of(eta, x) + u_of(eta, s - x) for x in cand)


def xi(eta):
    r = len(eta); E = list(eta) + [Fr(1)]
    return [max(E[i] + E[l - 1 - i] for i in range(0, l) if i <= r and l - 1 - i <= r) for l in range(1, 2 * r + 1)]


def box_formula(eta, s):
    return sum(1 for x in xi(eta) if x <= s)


def violation(eta, law):
    X = xi(eta)
    return sum(m * sum(1 for x in X if x <= a) for a, m in law) - 2 * sum(1 - e for e in eta)


def violation_brute(eta, law):
    return sum(m * box_brute(eta, a) for a, m in law) - 2 * sum(1 - e for e in eta)


def reflect(law):
    return [(2 - a, m) for a, m in law]


def _solve(A, b):
    n = len(A[0]); M = [list(r) + [bb] for r, bb in zip(A, b)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]; M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def stair_lp_bound(law, r):
    """Exact supremum of the violation over staircases with exactly r thresholds in (0, 1] (see module docstring).
    Returns (value, maximizing thresholds)."""
    law = sorted(law); atoms = [a for a, _ in law]; F = [m for _, m in law]
    best = None; bestE = None
    for c in combinations_with_replacement(range(2 * r + 1), len(atoms)):
        cons = []
        for m, cm in enumerate(c):
            for i in range(0, cm):
                j = cm - 1 - i
                if i > r or j > r:
                    continue
                row = [Fr(0)] * r; rhs = atoms[m]
                for idx in (i, j):
                    if idx == r:
                        rhs -= 1
                    else:
                        row[idx] += 1
                cons.append((row, rhs))
        for q in range(r):
            row = [Fr(0)] * r; row[q] = Fr(-1); cons.append((row, Fr(0)))
            row = [Fr(0)] * r; row[q] = Fr(1)
            if q + 1 < r:
                row[q + 1] = Fr(-1); cons.append((row, Fr(0)))
            else:
                cons.append((row, Fr(1)))
        cons = [(list(a), b) for a, b in {(tuple(a), b) for a, b in cons}]
        if any(not any(a) and b < 0 for a, b in cons):
            continue
        opt = None; optE = None
        for S in combinations(range(len(cons)), r):
            E = _solve([cons[s][0] for s in S], [cons[s][1] for s in S])
            if E is None:
                continue
            if all(sum(x * y for x, y in zip(a, E)) <= b for a, b in cons):
                val = sum(E)
                if opt is None or val > opt:
                    opt, optE = val, E
        if opt is None:
            continue
        V = sum(f * cm for f, cm in zip(F, c)) - 2 * r + 2 * opt
        if best is None or V > best:
            best, bestE = V, optE
    return best, bestE


def grid_v(u):
    """v(t) = min_{i+j=t} u(i) + u(j) on {0..k-1}."""
    k = len(u)
    return [min(u[i] + u[t - i] for i in range(max(0, t - k + 1), min(k - 1, t) + 1)) for t in range(2 * k - 1)]
