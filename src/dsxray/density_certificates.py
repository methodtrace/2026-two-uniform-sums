"""Single steps at density two (Proposition 5.2) and dual certificates for staircases (Lemma 5.1), exact.

Lemma 5.1: if phi is convex on [0, 2] and
        J(phi) = (1/2) int_0^2 phi + c int_0^2 (N - phi)^+  <=  2 int_0^1 u,          N = u box u,
then (1) holds for u at every law with a density at most c and S <=cx 2U.

Single step (Proposition 5.2): u = 1{x >= eta}, N = 1{s >= 2 eta} + 1{s >= 1 + eta}, 2 int u = 2(1 - eta), and with
t = (8 eta - 2)/3 the function phi(s) = (s - t)^+ / (1 - eta) gives J(phi) = 2(1 - eta) - (t^-)^2 / (4(1 - eta)) for
c = 2 (t^- = max(0, -t)). `single_step_J` evaluates J(phi) exactly by integrating piecewise linear functions.

Staircases with distinct thresholds J/K (J a subset of {1..K-1}): N(s) = #{l : xi_l <= s} with the thresholds xi_l of
Lemma 4.1 (multiples of 1/K). A convex piecewise linear phi with nodes on the grid 1/(Km) is found by a
floating-point linear program (scipy; it only proposes a candidate), its slopes are rounded to rationals and made
nondecreasing (so the rounded phi is convex by construction), and J(phi) is evaluated in exact rational arithmetic
(`certificate_value`). The certificate is accepted iff J(phi) <= 2 int u exactly.
"""
from fractions import Fraction as Fr
from itertools import combinations


# ---------------- exact piecewise-linear integration ----------------
def int_pos_linear(a, b, y0, y1):
    """int_a^b of (the linear function from y0 at a to y1 at b)^+ , exact."""
    if y0 >= 0 and y1 >= 0:
        return (y0 + y1) / 2 * (b - a)
    if y0 <= 0 and y1 <= 0:
        return Fr(0)
    p, q = (y0, y1) if y0 > 0 else (y1, y0)
    return p * p / (p - q) / 2 * (b - a)


# ---------------- single step ----------------
def single_step_phi(eta):
    t = (8 * eta - 2) / 3; lam = 1 / (1 - eta)
    return t, (lambda s: max(Fr(0), (s - t) * lam))


def single_step_J(eta, c=Fr(2)):
    """Exact J(phi) for the single step at eta with the phi of Proposition 5.2."""
    eta = Fr(eta); t, phi = single_step_phi(eta)
    pts = sorted({Fr(0), Fr(2), 2 * eta, 1 + eta} | ({t} if 0 < t < 2 else set()))
    J = Fr(0)
    for a, b in zip(pts, pts[1:]):
        Nv = (1 if a >= 2 * eta else 0) + (1 if a >= 1 + eta else 0)
        J += Fr(1, 2) * (phi(a) + phi(b)) / 2 * (b - a) + c * int_pos_linear(a, b, Nv - phi(a), Nv - phi(b))
    return J


def single_step_closed_form(eta):
    eta = Fr(eta); t = (8 * eta - 2) / 3; tm = max(Fr(0), -t)
    return 2 * (1 - eta) - tm * tm / (4 * (1 - eta))


def single_step_facts(eta):
    """The intermediate facts of the proof of Proposition 5.2, exact (each entry True)."""
    eta = Fr(eta); t, phi = single_step_phi(eta)
    x1, x2 = (5 * eta + 1) / 3, (2 * eta + 4) / 3
    tri = Fr(0)
    for a, b, lev in ((2 * eta, x1, 1), (1 + eta, x2, 2)):
        tri += int_pos_linear(a, b, lev - phi(a), lev - phi(b))
    P = sorted({Fr(0), Fr(2)} | ({t} if 0 < t < 2 else set()))
    half_int_phi = sum(Fr(1, 2) * (phi(a) + phi(b)) / 2 * (b - a) for a, b in zip(P, P[1:]))
    out = {"t_le_2eta": t <= 2 * eta,
           "phi_values": (phi(2 * eta), phi(x1), phi(1 + eta), phi(x2)) == (Fr(2, 3), 1, Fr(5, 3), 2) if t <= 2 * eta else False,
           "order": 2 * eta <= x1 <= 1 + eta <= x2 <= 2,
           "two_triangles": 2 * tri == 2 * (1 - eta) / 9,
           "t_nonneg_iff_eta_ge_quarter": (t >= 0) == (eta >= Fr(1, 4))}
    if t >= 0:
        out["half_integral"] = half_int_phi == 16 * (1 - eta) / 9 == (2 - t) ** 2 / (4 * (1 - eta))
    else:
        out["half_integral"] = half_int_phi == (1 - t) / (1 - eta) == (2 - t) ** 2 / (4 * (1 - eta)) - t * t / (4 * (1 - eta))
    return out


# ---------------- staircases with distinct grid thresholds ----------------
def xi_int(J, K):
    """Thresholds xi_1..xi_2r of Lemma 4.1, as multiples of 1/K, for thresholds J/K (distinct, sorted)."""
    e = sorted(J) + [K]; r = len(J); out = []
    for l in range(1, 2 * r + 1):
        out.append(max(e[q - 1] + e[l - q] for q in range(1, r + 2) if 1 <= l + 1 - q <= r + 1))
    return out


def N_on_cells(J, K, m):
    """N = u box u on the cells [i/G, (i+1)/G), G = Km, i = 0..2G-1."""
    xi = xi_int(J, K); G = K * m
    return [sum(1 for x in xi if x * m <= i) for i in range(2 * G)]


def lp_candidate(J, K, m, c=2.0):
    """Floating-point LP: node values of a convex phi on the grid 1/(Km) minimizing a discretization of J(phi)."""
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix
    G = K * m; n = 2 * G; nv = n + 1
    Nc = np.array(N_on_cells(J, K, m), float)
    obj = np.zeros(nv + 2 * n)
    for i in range(n):
        obj[i] += 0.25 / G; obj[i + 1] += 0.25 / G
    obj[nv:] = c * 0.5 / G
    A = lil_matrix((2 * n + nv - 2, nv + 2 * n)); b = np.zeros(2 * n + nv - 2); row = 0
    for i in range(n):
        A[row, i] = -1; A[row, nv + i] = -1; b[row] = -Nc[i]; row += 1
        A[row, i + 1] = -1; A[row, nv + n + i] = -1; b[row] = -Nc[i]; row += 1
    for i in range(1, nv - 1):
        A[row, i - 1] = -1; A[row, i] = 2; A[row, i + 1] = -1; row += 1
    res = linprog(obj, A_ub=A.tocsr(), b_ub=b, bounds=[(None, None)] * nv + [(0, None)] * (2 * n), method="highs")
    return list(res.x[:nv])


def rationalize_convex(phi_f, den=10000):
    """Approximate node values by a convex rational piecewise linear function: rational slopes, made nondecreasing."""
    slopes = [Fr(phi_f[i + 1] - phi_f[i]).limit_denominator(den) for i in range(len(phi_f) - 1)]
    for i in range(1, len(slopes)):
        if slopes[i] < slopes[i - 1]:
            slopes[i] = slopes[i - 1]
    phi = [Fr(phi_f[0]).limit_denominator(den)]
    for s in slopes:
        phi.append(phi[-1] + s)
    return phi


def is_convex(phi):
    return all(phi[i + 1] - phi[i] <= phi[i + 2] - phi[i + 1] for i in range(len(phi) - 2))


def certificate_value(phi, J, K, m, c=Fr(2)):
    """Exact J(phi) for node values phi at i/G (G = Km), and the budget 2 int u = 2 sum (1 - j/K)."""
    G = K * m; h = Fr(1, G); Nc = N_on_cells(J, K, m); tot = Fr(0)
    for i in range(2 * G):
        a, b = phi[i], phi[i + 1]
        tot += Fr(1, 2) * (a + b) / 2 * h + c * int_pos_linear(Fr(0), h, Nc[i] - a, Nc[i] - b)
    return tot, 2 * sum(1 - Fr(j, K) for j in J)


def certify_staircase(J, K, m, c=Fr(2)):
    """(certified, exact margin = budget - J(phi)) for the staircase with thresholds J/K."""
    phi = rationalize_convex(lp_candidate(J, K, m, float(c)))
    tot, bud = certificate_value(phi, J, K, m, c)
    return is_convex(phi) and tot <= bud, bud - tot


def staircases(K, rmin=2):
    for r in range(rmin, K):
        yield from combinations(range(1, K), r)


# ---------------- floating-point search (numerical evidence only) ----------------
def _density_constraints(N):
    """Densities g_m <= c on the cells [m/N, (m+1)/N) of [0, 2]: mass and mean equations, and the convex-order
    inequalities pi_S(x) <= pi_2U(x) = (2 - x)^2/4 at the cell endpoints and midpoints (FLOAT)."""
    import numpy as np
    M = 2 * N; a = np.arange(M) / N; b = a + 1.0 / N; ctr = (a + b) / 2
    Aeq = np.vstack([np.ones(M) / N, ctr / N]); beq = np.array([1.0, 1.0])
    rows = []; rhs = []
    for x in np.arange(1, 4 * N) / (2 * N):
        e = np.where(x <= a, ctr - x, np.where(x >= b, 0.0, (b - x) ** 2 / (2 * (b - a))))
        rows.append(e / N); rhs.append((2 - x) ** 2 / 4)
    return Aeq, beq, np.array(rows), np.array(rhs)


def float_max_violation(eta, N, c, cons):
    """FLOAT: max over densities <= c on cells 1/N with S <=cx 2U of sum_l P(S >= xi_l) - 2 sum (1 - eta_q)."""
    import numpy as np
    from scipy.optimize import linprog
    Aeq, beq, Aub, bub = cons
    e = sorted(eta) + [1.0]; r = len(eta)
    xi = [max(e[q - 1] + e[l - q] for q in range(1, r + 2) if 1 <= l + 1 - q <= r + 1) for l in range(1, 2 * r + 1)]
    a = np.arange(2 * N) / N; b = a + 1.0 / N
    w = sum(np.clip(b - np.maximum(a, t), 0, 1.0 / N) for t in xi)
    res = linprog(-w, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=[(0, c)] * (2 * N), method="highs")
    return -res.fun - 2 * sum(1 - x for x in eta)


def float_search(c, N, K, r, seeds, rng):
    """FLOAT: coordinate ascent over staircases with r thresholds in {1/K, ..., K/K} (repetitions allowed), from
    random starts; returns the list of (best violation, thresholds as multiples of 1/K) per start."""
    cons = _density_constraints(N); cache = {}

    def F(et):
        key = tuple(et)
        if key not in cache:
            cache[key] = float_max_violation([x / K for x in et], N, c, cons)
        return cache[key]
    out = []
    for _ in range(seeds):
        et = sorted(int(x) for x in rng.integers(1, K, size=r)); cur = F(et)
        for _sweep in range(8):
            changed = False
            for q in range(r):
                lo = et[q - 1] if q > 0 else 1; hi = et[q + 1] if q < r - 1 else K
                for x in range(lo, hi + 1):
                    if x != et[q]:
                        e2 = et.copy(); e2[q] = x; val = F(e2)
                        if val > cur + 1e-12:
                            cur, et, changed = val, e2, True
            if not changed:
                break
        out.append((cur, et))
    return out
