"""Joint tail-sum inequalities for U + V (U, V standard uniform, arbitrary dependence); exact rational arithmetic.

A unit staircase with r steps is u(x) = #{i < r : eta_i <= x}, 0 < eta_0 <= ... <= eta_{r-1} < 1 (ties allowed),
N(s) = (u box u)(s) = min_x u(x) + u(s - x), and I = int_0^1 u = sum (1 - eta_i). A law F of S is excluded by u if
sum_s F(s) N(s) > 2 I (the violation is the difference). Contents:

  boxconv / violation          brute-force N(s) for step functions (independent of any closed formula)
  counting lemma               N(s) (2 - s) <= 2 I for every staircase and s in [0, 2]
  two-point laws               explicit floor(1/beta)-step potential and its violation; lower-bound checks
  max_violation                exact maximum of the violation over r-step staircases (vertex enumeration)
  three-point family           explicit r-step potential for symmetric three-point laws, t_r(b)
  densities                    stop-loss transforms of piecewise-constant densities, convex order, folded order,
                               the construction for densities at most 1, the sharpness family for c > 2
  lattice witnesses            laws on (1/k)Z with small atoms that satisfy convex and folded order but not a facet
Repeated thresholds (Example 4.2) and a second implementation of the staircase programs are in tailsums.py.
"""
import itertools
from fractions import Fraction as Fr


# ---------------- brute-force evaluation ----------------
def make_u(pts, vals):
    def u(x):
        return vals[sum(1 for p in pts if p <= x)]
    return u


def boxconv(pts, vals, s):
    """min over x of u(x) + u(s - x) for the right-continuous step function with breakpoints pts and cell values vals."""
    u = make_u(pts, vals)
    lo, hi = max(Fr(0), s - 1), min(Fr(1), s)
    c = sorted({lo, hi} | {p for p in pts if lo <= p <= hi} | {s - p for p in pts if lo <= s - p <= hi})
    best = min(u(x) + u(s - x) for x in c)
    for a, b in zip(c, c[1:]):
        m = (a + b) / 2; best = min(best, u(m) + u(s - m))
    return best


def integral(pts, vals):
    e = [Fr(0)] + list(pts) + [Fr(1)]
    return sum(v * (e[i + 1] - e[i]) for i, v in enumerate(vals))


def staircase(eta):
    """(breakpoints, values) of u(x) = #{i : eta_i <= x}."""
    pts = sorted(set(eta))
    return pts, [0] + [sum(1 for q in eta if q <= p) for p in pts]


def violation(eta, atoms, masses):
    pts, vals = staircase(eta)
    return sum(f * boxconv(pts, vals, s) for s, f in zip(atoms, masses)) - 2 * integral(pts, vals)


# ---------------- counting lemma ----------------
def counting_lemma_check(trials, seed):
    """Random rational staircases (1 to 7 steps, ties allowed) and s random or at all pair sums.
    Returns (pairs checked, failures)."""
    import random
    rng = random.Random(seed); bad = tot = 0
    for _ in range(trials):
        n = rng.randint(1, 7); M = rng.choice([6, 12, 24, 60, 97])
        e = sorted(Fr(rng.randint(1, M), M) for _ in range(n))
        pts, vals = staircase(e); I = sum(1 - x for x in e)
        S = [Fr(rng.randint(0, 2 * M), M) for _ in range(3)] + [Fr(rng.randint(0, 2000), 1000)]
        S += [x + y for x in e for y in e + [Fr(1)]]
        for s in S:
            if s > 2:
                continue
            tot += 1
            if boxconv(pts, vals, s) * (2 - s) > 2 * I:
                bad += 1
    return tot, bad


# ---------------- two-point laws ----------------
def two_point_law(L, w):
    """Mean-one law with mass w at a = 1 - (1-w) beta and 1 - w at a + beta, beta = 1/L."""
    beta = 1 / L; a = 1 - (1 - w) * beta
    return [a, a + beta], [w, 1 - w]


def two_point_potential(L, w):
    """The explicit potential of the two-point theorem: r = floor(L) thresholds centered at a/2 with spacing beta,
    the first capped at w beta; for L < 1 one step at (w + min(w beta, 1))/2."""
    beta = 1 / L; a = 1 - (1 - w) * beta
    if L < 1:
        return [(w + min(w * beta, Fr(1))) / 2]
    r = int(L)
    e = [a / 2 + (i - Fr(r - 1, 2)) * beta for i in range(r)]
    e[0] = min(e[0], w * beta)
    return e


def two_point_bound(L, w):
    """Lower bound of the proof for the violation of two_point_potential."""
    if L < 1:
        return 2 * (two_point_potential(L, w)[0] - w)
    r = int(L)
    return min((1 - w) * (L - r), w * (r + 1 - L)) / L


def two_point_upper_check(D, Lmax):
    """All laws with L in (1/2, 1) or L = n + j/D (n <= Lmax), w = i/D, support in [0, 2]: the potential has valid
    thresholds and violation >= the bound > 0. Returns (laws, failures, equal to the bound)."""
    cnt = bad = eq = 0
    Ls = [Fr(j, D) for j in range(D // 2 + 1, D)] + [n + Fr(j, D) for n in range(1, Lmax + 1) for j in range(1, D)]
    for L in Ls:
        beta = 1 / L
        for wi in range(1, D):
            w = Fr(wi, D); a = 1 - (1 - w) * beta; c = a + beta
            if a < 0 or c > 2:
                continue
            e = two_point_potential(L, w)
            ok = all(0 < x < 1 for x in e) and all(e[i] <= e[i + 1] for i in range(len(e) - 1))
            v = violation(e, [a, c], [w, 1 - w]); bound = two_point_bound(L, w)
            ok = ok and v >= bound > 0
            eq += (v == bound); cnt += 1; bad += (not ok)
    return cnt, bad, eq


def two_point_lower_check(trials, seed):
    """Random staircases with m steps against two-point laws with 1/beta >= m + 1, both w and 1 - w (the reflected
    law), thresholds random or snapped to the critical values (with +-10^-6 shifts): no violation expected.
    Returns (staircases, violations)."""
    import random
    rng = random.Random(seed); tot = bad = 0
    for _ in range(trials):
        m = rng.randint(1, 8); D = rng.choice([5, 7, 11, 17])
        Lb = (m + 1) + Fr(rng.randint(0, 3 * D), D); beta = 1 / Lb; w = Fr(rng.randint(1, D - 1), D)
        for ww in (w, 1 - w):
            a = 1 - (1 - ww) * beta; c = 1 + ww * beta
            M = rng.choice([30, 120, 840]); e = set()
            while len(e) < m:
                if rng.random() < 0.5:
                    x = Fr(rng.randint(1, M - 1), M)
                else:
                    kk = rng.randint(0, 2 * m + 2)
                    x = rng.choice([a / 2 + kk * beta / 2, a / 2 - kk * beta / 2, c / 2 + kk * beta / 2,
                                    c / 2 - kk * beta / 2, ww * beta, kk * beta, c - 1 + kk * beta, a - kk * beta])
                    x += rng.choice([0, 0, Fr(1, 10**6), -Fr(1, 10**6)])
                if 0 < x < 1:
                    e.add(x)
            e = sorted(e); vals = [0] + list(range(1, m + 1))
            lhs = ww * boxconv(e, vals, a) + (1 - ww) * boxconv(e, vals, c)
            tot += 1; bad += (lhs > 2 * sum(1 - x for x in e))
    return tot, bad


# ---------------- exact maximum over r-step staircases ----------------
def _solve_exact(rows, n):
    M = [list(r[0]) + [r[1]] for r in rows]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]; pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]


def staircase_constraints(r, atoms, counts):
    """Constraints coef.eta <= rhs on thresholds eta in [0,1]^r (nondecreasing) under which N(s) >= counts(s) at every
    atom s: N(s) >= l iff eta_i + eta_j <= s for all i + j = l - 1 (eta_r = 1)."""
    cons = []
    vec = lambda d: tuple(Fr(d.get(i, 0)) for i in range(r))
    for i in range(r):
        cons.append((vec({i: -1}), Fr(0))); cons.append((vec({i: 1}), Fr(1)))
    for i in range(r - 1):
        cons.append((vec({i: 1, i + 1: -1}), Fr(0)))
    for s, c in zip(atoms, counts):
        if c == 0:
            continue
        for i in range(0, r + 1):
            j = c - 1 - i
            if j < 0 or j > r:
                continue
            d = {}; rhs = Fr(s)
            for t in (i, j):
                if t == r:
                    rhs -= 1
                else:
                    d[t] = d.get(t, 0) + 1
            if not d:
                if rhs < 0:
                    return None
                continue
            cons.append((vec(d), rhs))
    return list(set(cons))


def lp_max_sum(r, cons):
    best = None; arg = None
    for combo in itertools.combinations(range(len(cons)), r):
        sol = _solve_exact([cons[k] for k in combo], r)
        if sol is None:
            continue
        if all(sum(a * x for a, x in zip(c[0], sol)) <= c[1] for c in cons):
            v = sum(sol)
            if best is None or v > best:
                best = v; arg = sol
    return best, arg


def max_violation(r, atoms, masses):
    """Exact supremum of the violation over staircases with at most r steps (a threshold at 0 removes a step).
    Exponential in r; intended for r <= 4."""
    order = sorted(range(len(atoms)), key=lambda i: atoms[i])
    atoms = [atoms[i] for i in order]; masses = [masses[i] for i in order]
    best = None; arg = None
    for counts in itertools.combinations_with_replacement(range(2 * r + 1), len(atoms)):
        cons = staircase_constraints(r, atoms, counts)
        if cons is None:
            continue
        m, eta = lp_max_sum(r, cons)
        if m is None:
            continue
        v = sum(f * c for f, c in zip(masses, counts)) - 2 * r + 2 * m
        if best is None or v > best:
            best = v; arg = (counts, eta)
    return best, arg


def two_point_lp_check(D, Lmax, rmax):
    """For L = n + j/D (n <= Lmax), w = i/D: with r = floor(L) <= rmax, the exact r-step maximum equals the explicit
    violation, and the (r-1)-step maximum is <= 0 for the law and its reflection. Returns (laws, failures)."""
    cnt = bad = 0
    for n in range(1, Lmax + 1):
        for j in range(1, D):
            L = n + Fr(j, D); b = 1 / L; r = n
            if r > rmax:
                continue
            for wi in range(1, D):
                w = Fr(wi, D); a = 1 - (1 - w) * b
                v = min((1 - w) * (L - r), w * (r + 1 - L)) / L
                ok = max_violation(r, [a, a + b], [w, 1 - w])[0] == v
                if r >= 2:
                    ok = ok and max_violation(r - 1, [a, a + b], [w, 1 - w])[0] <= 0
                    ok = ok and max_violation(r - 1, [2 - a - b, 2 - a], [1 - w, w])[0] <= 0
                cnt += 1; bad += (not ok)
    return cnt, bad


# ---------------- symmetric three-point laws ----------------
def dist2(L):
    n = L // 2
    return min(L - 2 * n, 2 * n + 2 - L)


def three_point_potential(r, b):
    """eta_1 = min(b, (1-(2r-1)b)/2), eta_i = (1-(2r+3)b)/2 + 2ib (2 <= i <= r), for 2r <= 1/b <= 2r + 2."""
    return [min(b, (1 - (2 * r - 1) * b) / 2)] + [(1 - (2 * r + 3) * b) / 2 + 2 * i * b for i in range(2, r + 1)]


def three_point_family_check(R, Dmax=24):
    """For r <= R and L = 2r + j/D (D <= Dmax, 0 <= j <= 2D): valid thresholds, N at 1-b, 1, 1+b equal to
    (r, r, r+1) by brute force, and the inequality reads t >= dist(L, 2Z)/L. Returns (spacings, failures)."""
    bad = cnt = 0
    for r in range(1, R + 1):
        for D in range(1, Dmax + 1):
            for j in range(0, 2 * D + 1):
                L = Fr(2 * r) + Fr(j, D); b = 1 / L; eta = three_point_potential(r, b)
                ok = all(0 < x < 1 for x in eta) and all(eta[i] <= eta[i + 1] for i in range(r - 1))
                pts, vals = staircase(eta)
                ok = ok and [boxconv(pts, vals, s) for s in (1 - b, Fr(1), 1 + b)] == [r, r, r + 1]
                ok = ok and 4 * sum(eta) - (2 * r - 1) == dist2(L) / L
                cnt += 1; bad += (not ok)
    return cnt, bad


def smallest_middle_mass_exact(b, r):
    """t_r(b): exact smallest middle mass t such that ((1-t)/2, t, (1-t)/2) on {1-b, 1, 1+b} satisfies every
    potential inequality with at most r unit steps (0 if none excludes any member)."""
    atoms = [1 - b, Fr(1), 1 + b]; best = Fr(0)
    for rr in range(1, r + 1):
        for c in itertools.combinations_with_replacement(range(2 * rr + 1), 3):
            den = Fr(c[0] + c[2], 2) - c[1]
            if den <= 0:
                continue
            cons = staircase_constraints(rr, atoms, c)
            if cons is None:
                continue
            m, eta = lp_max_sum(rr, cons)
            if m is None:
                continue
            t = (Fr(c[0] + c[2], 2) - 2 * rr + 2 * m) / den
            best = max(best, t)
    return best


# ---------------- laws with a piecewise-constant density ----------------
def pieces_pi(pieces, x):
    """E(S - x)^+ for a density given by pieces (lo, hi, height)."""
    t = Fr(0)
    for lo, hi, d in pieces:
        a = max(lo, x)
        if a < hi:
            t += d * ((hi - x) ** 2 - (a - x) ** 2) / 2
    return t


def _min_quadratic_nonneg(f, bps):
    """Is the piecewise quadratic f (quadratic between consecutive breakpoints) >= 0 on [bps[0], bps[-1]]?"""
    for a, b in zip(bps, bps[1:]):
        fa, fb, fm = f(a), f(b), f((a + b) / 2)
        if min(fa, fb, fm) < 0:
            return False
        A2 = 2 * (fa - 2 * fm + fb) / ((b - a) ** 2)
        if A2 > 0:
            B = (fb - fa) / (b - a) - A2 * (b - a) / 2; xv = a - B / A2
            if a < xv < b and f(xv) < 0:
                return False
    return True


def cx_below_2U(pieces):
    """S <=cx 2U for S with the given density (exact): E(2U - x)^+ - E(S - x)^+ >= 0 on [0, 2]."""
    bps = sorted({Fr(0), Fr(2)} | {p for lo, hi, d in pieces for p in (lo, hi)})
    return _min_quadratic_nonneg(lambda x: (2 - x) ** 2 / 4 - pieces_pi(pieces, x), bps)


def folded_order_ok(pieces):
    """The folded condition: U <=cx (law(S^flat) + law(kappa))/2, S^flat = S - 1{S > 1}, kappa = 1{S > 1}."""
    flat = []
    for lo, hi, d in pieces:
        if hi <= 1:
            flat.append((lo, hi, d))
        elif lo >= 1:
            flat.append((lo - 1, hi - 1, d))
        else:
            flat.append((lo, Fr(1), d)); flat.append((Fr(0), hi - 1, d))
    pS1 = sum(d * max(Fr(0), min(hi, Fr(1)) - lo) for lo, hi, d in pieces)

    def g(y):
        t = Fr(0)
        for lo, hi, d in flat:
            b = min(hi, y)
            if b > lo:
                t += d * ((y - lo) ** 2 - (y - b) ** 2) / 2
        return t + y * pS1 - y * y
    bps = sorted({Fr(0), Fr(1)} | {x for lo, hi, d in flat for x in (lo, hi) if 0 <= x <= 1})
    return _min_quadratic_nonneg(g, bps)


def prob_ge(pieces, s):
    return sum(d * max(Fr(0), hi - max(lo, s)) for lo, hi, d in pieces)


def one_step_violation(pieces, eta):
    return prob_ge(pieces, 2 * eta) + prob_ge(pieces, 1 + eta) - 2 * (1 - eta)


def sharpness_family(c):
    """Two bumps of height c > 2, mass 1/2 each, at [2 eta, 2 eta + h] and [1 + eta, 1 + eta + h], h = 1/(2c),
    eta = (2c - 1)/(6c). Returns (pieces, eta)."""
    eta = (2 * c - 1) / (6 * c); h = 1 / (2 * c)
    return [(2 * eta, 2 * eta + h, c), (1 + eta, 1 + eta + h, c)], eta


def tent_discretization(pieces, k):
    """Randomized rounding of S to (1/k)Z: mass at j/k is int g(s) (1 - |ks - j|)^+ ds."""
    def F2(lo, hi, j):
        tot = Fr(0)
        for A, B, left in ((Fr(j - 1, k), Fr(j, k), True), (Fr(j, k), Fr(j + 1, k), False)):
            x0, x1 = max(lo, A), min(hi, B)
            if x1 <= x0:
                continue
            if left:
                tot += k * (x1 * x1 - x0 * x0) / 2 - (j - 1) * (x1 - x0)
            else:
                tot += (j + 1) * (x1 - x0) - k * (x1 * x1 - x0 * x0) / 2
        return tot
    return {j: sum(d * F2(lo, hi, j) for lo, hi, d in pieces) for j in range(0, 2 * k + 1)}


def density_construction_check(pieces, k, t=Fr(1, 4)):
    """The construction of the density theorem with margin t: S' = (1-t) law(1 + rho(S-1)) + t U[1/4, 7/4],
    rho = 1 - t/3; G_k = tent discretization of S'. Returns the five checks: density of S' <= 1; atoms of G_k <= 1/k
    and none at 0 or 2; mass and mean 1; stop-loss gap >= 1/(4k^2) at every m/k; k G_k - 1 <=cx 2 U_k."""
    rho = 1 - t / 3
    Sp = [(1 + rho * (lo - 1), 1 + rho * (hi - 1), (1 - t) * d / rho) for lo, hi, d in pieces] + [(Fr(1, 4), Fr(7, 4), t * Fr(2, 3))]
    bps = sorted({x for lo, hi, d in Sp for x in (lo, hi)})
    dens_ok = max(sum(d for lo, hi, d in Sp if lo <= (x + y) / 2 < hi) for x, y in zip(bps, bps[1:])) <= 1
    G = tent_discretization(Sp, k)
    atoms_ok = all(v <= Fr(1, k) for v in G.values()) and G[0] == 0 and G[2 * k] == 0
    mass_mean = sum(G.values()) == 1 and sum(v * Fr(j, k) for j, v in G.items()) == 1
    gap_ok = all((2 - Fr(m, k)) ** 2 / 4 - pieces_pi(Sp, Fr(m, k)) >= Fr(1, 4 * k * k) for m in range(1, 2 * k))
    slT = lambda m: sum(v * max(0, (j - 1) - m) for j, v in G.items())
    slU = lambda m: sum(Fr(1, k) * max(0, 2 * i - m) for i in range(k))
    cx = all(slT(m) <= slU(m) for m in range(-1, 2 * k))
    return dens_ok, atoms_ok, mass_mean, gap_ok, cx


# ---------------- lattice witnesses (count coordinates on the grid k) ----------------
def v_of(u):
    k = len(u)
    return [min(u[i] + u[t - i] for i in range(max(0, t - k + 1), min(k - 1, t) + 1)) for t in range(2 * k - 1)]


def lattice_witness_check(k, x, u, c):
    """x in count coordinates (sum x = k, sum t x_t = k(k-1)), atoms <= c: satisfies every hinge inequality of both
    signs (convex and folded order) and violates I(u). Returns (affine hull, atoms bounded, hinges, violation)."""
    T = 2 * k - 1; x = [Fr(v) for v in x]
    hull = sum(x) == k and sum(t * v for t, v in enumerate(x)) == k * (k - 1)
    atoms = min(x) >= 0 and max(x) <= c
    H = [[max(0, i - j) for i in range(k)] for j in range(1, k - 1)]; H += [[-y for y in h] for h in H]
    hinges = all(sum(a * b for a, b in zip(v_of(h), x)) <= 2 * sum(h) for h in H)
    viol = sum(a * b for a, b in zip(v_of(u), x)) - 2 * sum(u)
    return hull, atoms, hinges, viol


def lattice_witness_search(k, c, facet_rays):
    """Search (floating-point LP, scipy) for a law with count coordinates in [0, c] satisfying all hinge
    inequalities and violating a facet I(u), u among facet_rays; the LP vertex is recovered exactly from its active
    constraints and every condition is rechecked in rational arithmetic. Returns (violation, u, x) or None."""
    import numpy as np
    from scipy.optimize import linprog
    T = 2 * k - 1
    H = [[max(0, i - j) for i in range(k)] for j in range(1, k - 1)]; H += [[-y for y in h] for h in H]
    Aub = [v_of(h) for h in H]; bub = [2 * sum(h) for h in H]
    Aeq = [[1] * T, list(range(T))]; beq = [k, k * (k - 1)]

    def exact_vertex(xf):
        rows = [(Aeq[0], Fr(beq[0])), (Aeq[1], Fr(beq[1]))]
        for a, b in zip(Aub, bub):
            if abs(np.dot(a, xf) - b) < 1e-7:
                rows.append((a, Fr(b)))
        for t in range(T):
            e = [0] * T; e[t] = 1
            if abs(xf[t]) < 1e-9:
                rows.append((e, Fr(0)))
            elif abs(xf[t] - float(c)) < 1e-9:
                rows.append((e, c))
        piv = []; R = []
        for a, b in rows:
            r = [Fr(v) for v in a] + [b]
            for (pc, pr) in zip(piv, R):
                if r[pc] != 0:
                    f = r[pc] / pr[pc]; r = [x1 - f * y1 for x1, y1 in zip(r, pr)]
            nz = next((j for j in range(T) if r[j] != 0), None)
            if nz is None:
                continue
            piv.append(nz); R.append(r)
        if len(piv) < T:
            return None
        sol = [Fr(0)] * T
        for pc, pr in reversed(list(zip(piv, R))):
            sol[pc] = (pr[T] - sum(pr[j] * sol[j] for j in range(T) if j != pc)) / pr[pc]
        return sol
    best = None
    for u in facet_rays:
        v = v_of(u)
        res = linprog(-np.array(v, float), A_ub=np.array(Aub, float), b_ub=np.array(bub, float),
                      A_eq=np.array(Aeq, float), b_eq=beq, bounds=[(0, float(c))] * T, method="highs-ds")
        if res.status == 0 and -res.fun - 2 * sum(u) > 1e-7:
            xs = exact_vertex(res.x)
            if xs is None:
                continue
            hull, atoms, hinges, viol = lattice_witness_check(k, xs, u, c)
            if hull and atoms and hinges and viol > 0 and (best is None or viol > best[0]):
                best = (viol, list(u), xs)
    return best
