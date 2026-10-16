"""Step inequalities for laws on finitely many atoms in [0, 2] (exact rational arithmetic).

For thresholds 0 < sigma_1 < ... < sigma_r < 1 let u(x) = #{i : sigma_i <= x} on [0, 1]. The potential inequality
reads sum_s F(s) (u box u)(s) <= 2 int_0^1 u, with (u box u)(s) = min_x u(x) + u(s - x). For a law F with atoms A and
masses M, `viol(A, M, sigma)` is the left side minus the right side (positive: F is excluded).

The violation changes only when a threshold sigma_i + sigma_j (i = j allowed) or 1 + sigma_i crosses an atom.
`best_viol` evaluates it at every vertex of the arrangement of these hyperplanes and at perturbations of each vertex
by 10^-9 in every coordinate (to capture one-sided limits); this is a search over finitely many exact rational
points, not a certificate of the maximum. The exact maximum over r-step staircases is tail_bounds.max_violation.
"""
import itertools
from fractions import Fraction as Fr


def boxconv(sig, s):
    """(u box u)(s) for u(x) = #{q in sig : q <= x}, exactly over the breakpoints and the midpoints between them."""
    lo, hi = max(Fr(0), s - 1), min(Fr(1), s)
    pts = sorted({lo, hi} | {q for q in sig if lo <= q <= hi} | {s - q for q in sig if lo <= s - q <= hi})
    u = lambda x: sum(1 for q in sig if q <= x)
    best = min(u(x) + u(s - x) for x in pts)
    for a, b in zip(pts, pts[1:]):
        m = (a + b) / 2; best = min(best, u(m) + u(s - m))
    return best


def viol(A, M, sig):
    return sum(m * boxconv(sig, a) for a, m in zip(A, M)) - 2 * sum(1 - q for q in sig)


def candidates(A, r):
    """Vertices of the arrangement {2 sigma_i = a, sigma_i + sigma_j = a, 1 + sigma_i = a : a in A} in (0,1)^r."""
    eqs = []
    for a in A:
        for i in range(r):
            eqs.append(({i: 2}, a)); eqs.append(({i: 1}, a - 1))
            for j in range(i + 1, r):
                eqs.append(({i: 1, j: 1}, a))
    out = set()
    for combo in itertools.combinations(eqs, r):
        Mx = [[Fr(e[0].get(i, 0)) for i in range(r)] + [Fr(e[1])] for e in combo]
        ok = True
        for c in range(r):
            p = next((rr for rr in range(c, r) if Mx[rr][c] != 0), None)
            if p is None:
                ok = False; break
            Mx[c], Mx[p] = Mx[p], Mx[c]
            for rr in range(r):
                if rr != c and Mx[rr][c] != 0:
                    f = Mx[rr][c] / Mx[c][c]; Mx[rr] = [x - f * y for x, y in zip(Mx[rr], Mx[c])]
        if not ok:
            continue
        sol = tuple(sorted(Mx[i][r] / Mx[i][i] for i in range(r)))
        if all(0 < s < 1 for s in sol) and len(set(sol)) == r:
            out.add(sol)
    return out


def best_viol(A, M, r):
    """Largest violation found over r-step thresholds (both signs; sign - is the reflected law)."""
    b = Fr(-10)
    for sign in (1, -1):
        AA = A if sign == 1 else [2 - a for a in A]
        for sig in candidates(AA, r):
            for d in itertools.product((Fr(0), Fr(1, 10**9), -Fr(1, 10**9)), repeat=r):
                s2 = [x + y for x, y in zip(sig, d)]
                if all(0 < x < 1 for x in s2):
                    b = max(b, viol(AA, M, s2))
    return b


def mww_bound(b):
    """dist(L, 2Z)/L, L = 1/b: the smallest middle mass of a symmetric three-point law on {1-b, 1, 1+b} in D_U."""
    L = 1 / b; n = L // 2
    return min(L - 2 * n, 2 * (n + 1) - L) / L


def smallest_middle_mass(b, r, steps=25):
    """Bisection (exact rational midpoints) on the middle mass of the symmetric law: the smallest value not excluded
    by the r-step search of best_viol, to within 2^-steps."""
    lo, hi = Fr(0), Fr(1)
    for _ in range(steps):
        mid = (lo + hi) / 2; f = (1 - mid) / 2
        v = max(best_viol([1 - b, Fr(1), 1 + b], [f, mid, f], rr) for rr in range(1, r + 1))
        if v > 0:
            lo = mid
        else:
            hi = mid
    return hi


# ---------------- single-step inequality ----------------
def onestep_lhs(law, eta):
    """P(S >= 2 eta) + P(S >= 1 + eta) for a law given as {atom: mass}."""
    return sum(p for s, p in law.items() if s >= 2 * eta) + sum(p for s, p in law.items() if s >= 1 + eta)


def single_step_min_slack(law):
    """min over eta in (0,1) of 2(1 - eta) - P(S >= 2 eta) - P(S >= 1 + eta), exactly: the left side is piecewise
    constant between the breakpoints s/2 and s - 1, so the binding points are the right ends of the pieces (approached
    from the left) and the breakpoints themselves. Returns (min slack, piece)."""
    bps = sorted({Fr(0), Fr(1)} | {s / 2 for s in law if 0 < s / 2 < 1} | {s - 1 for s in law if 0 < s - 1 < 1})
    worst = None
    for a, b in zip(bps, bps[1:]):
        lhs_piece = onestep_lhs(law, (a + b) / 2)
        cands = [2 * (1 - b) - lhs_piece]
        if b < 1:
            cands.append(2 * (1 - b) - onestep_lhs(law, b))
        m = min(cands)
        if worst is None or m < worst[0]:
            worst = (m, a, b)
    return worst
