"""Linear programs over couplings of two lattice marginals.

X on {0..m} with pmf p, Y on {0..n} with pmf q. A coupling pi is an (m+1) x (n+1) matrix with row sums p and
column sums q; the law of S = X + Y is s_t = sum_{i+j=t} pi_ij. These routines use floating-point linear
programming (HiGHS through scipy) and serve as independent cross-checks of the exact results; no claim of the
paper rests on them alone.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, lil_matrix, vstack


@dataclass(frozen=True)
class Marginals:
    p: tuple
    q: tuple

    def __post_init__(self):
        for v in (self.p, self.q):
            if abs(sum(float(x) for x in v) - 1) > 1e-12 or min(float(x) for x in v) < 0:
                raise ValueError("marginals must be probability vectors")

    @property
    def m(self):
        return len(self.p) - 1

    @property
    def n(self):
        return len(self.q) - 1

    @property
    def N(self):
        return self.m + self.n + 1


def _idx(i, j, n):
    return i * (n + 1) + j


def _coupling_constraints(M: Marginals):
    m, n = M.m, M.n
    A = lil_matrix((m + 1 + n + 1, (m + 1) * (n + 1)))
    b = np.zeros(m + 1 + n + 1)
    for i in range(m + 1):
        for j in range(n + 1):
            A[i, _idx(i, j, n)] = 1
            A[m + 1 + j, _idx(i, j, n)] = 1
        b[i] = float(M.p[i])
    for j in range(n + 1):
        b[m + 1 + j] = float(M.q[j])
    return A, b


def _sum_map(M: Marginals):
    m, n = M.m, M.n
    T = lil_matrix((M.N, (m + 1) * (n + 1)))
    for i in range(m + 1):
        for j in range(n + 1):
            T[i + j, _idx(i, j, n)] = 1
    return T


def _quantile_sum(M: Marginals, anti: bool):
    """Law of F^{-1}(U) + G^{-1}(U) (anti=False) or F^{-1}(U) + G^{-1}(1-U) (anti=True)."""
    P = np.cumsum([float(x) for x in M.p])
    Q = np.cumsum([float(x) for x in M.q])
    Qr = np.cumsum([float(x) for x in M.q[::-1]])
    cuts = sorted(set([0.0, 1.0] + list(P) + list(Qr if anti else Q)))
    s = np.zeros(M.N)
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a < 1e-15:
            continue
        u = 0.5 * (a + b)
        x = int(np.searchsorted(P, u, side="left"))
        y = M.n - int(np.searchsorted(Qr, u, side="left")) if anti else int(np.searchsorted(Q, u, side="left"))
        s[x + y] += b - a
    return s


def comonotone_sum(M):
    return _quantile_sum(M, anti=False)


def countermonotone_sum(M):
    return _quantile_sum(M, anti=True)


def stop_loss(s):
    """E(S - d)^+ for d = 0..len(s)-1."""
    N = len(s)
    return np.array([sum((t - d) * s[t] for t in range(d, N)) for d in range(N)])


def reach_bounds(M: Marginals, k: int):
    """Minimum and maximum of P(X + Y >= k) over all couplings."""
    A, b = _coupling_constraints(M)
    c = np.asarray(_sum_map(M)[k:, :].sum(axis=0)).ravel()
    out = []
    for sign in (1, -1):
        res = linprog(sign * c, A_eq=csr_matrix(A), b_eq=b, bounds=(0, None), method="highs")
        if res.status != 0:
            raise RuntimeError(res.message)
        out.append(sign * res.fun)
    return out[0], out[1]


def attainable(M: Marginals, s):
    """Is the pmf s on {0..m+n} the law of X + Y for some coupling?"""
    A, b = _coupling_constraints(M)
    Aeq = vstack([csr_matrix(A), csr_matrix(_sum_map(M))])
    beq = np.concatenate([b, np.asarray(s, dtype=float)])
    res = linprog(np.zeros(Aeq.shape[1]), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    return res.status == 0


def makarov_cdf_bounds(M: Marginals):
    """Lattice Makarov bounds on P(S <= t), t = 0..m+n."""
    F = np.cumsum([float(x) for x in M.p])
    G = np.cumsum([float(x) for x in M.q])
    Fv = lambda i: 0.0 if i < 0 else (1.0 if i > M.m else F[i])
    Gv = lambda j: 0.0 if j < 0 else (1.0 if j > M.n else G[j])
    lo, hi = [], []
    for t in range(M.N):
        lo.append(max(max(Fv(i) + Gv(t - i) - 1, 0) for i in range(-1, t + 2)))
        hi.append(min(min(Fv(i) + Gv(t - i - 1), 1) for i in range(-1, t + 2)))
    return np.array(lo), np.array(hi)
