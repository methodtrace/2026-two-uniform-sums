"""Laws with small atoms on the grid: exact checks of the small-atoms corollary and of its proof.

Count coordinates on the grid k: x(t) = k P(T = t), T = kS - 1 in {0, ..., 2k-2}. The polytope

    Q(k) = {x : 0 <= x(t) <= 1, sum x = k, sum t x(t) = k(k-1), the k-2 hinge inequalities I(h_j)},
    I(h_j): sum_t (t - 2j)^+ x(t) <= sum_i (2i - 2j)^+,  1 <= j <= k-2,

contains hat-P(k) intersected with the box, and Corollary 3.2 says the two are equal. The proof uses
Lemma 3.2 (x >= 0 lies in hat-P(k) iff sum v(t) x(t) <= 2 sum u(i) for every u) and Lemma 3.3: for the
convex minorant phi of u (the lower convex envelope of the points (i, u(i))) and w = u - phi,

    sum_t (v(t) - 2 phi(t/2))^+ <= 2 sum_i w(i),

with the pointwise bounds v(t) <= 2 phi(t/2) + w(t - a) for 2a <= t <= a + b and v(t) <= 2 phi(t/2) + w(t - b) for
a + b <= t <= 2b on every piece [a, b] between consecutive contact points.

Contents (all exact, fractions.Fraction):
  simplex               a small exact simplex method (Bland's rule), feasibility and optimization
  find_Z, check_Z       a doubly stochastic matrix with prescribed anti-diagonal sums, by exact linear programming
  find_Z_fast           the same, with a floating-point proposal whose support is solved and verified exactly
  in_Q, satisfies_cx    membership in Q(k); T <=cx 2U_k at every integer hinge
  vertices_exact        all vertices of Q(k) (or of the polytope without the bound x <= 1) by exhaustive enumeration
                        of bases, in exact arithmetic (complete by construction; practical for k <= 6)
  vertices_cdd          the same with pycddlib (floating point); every vertex is recomputed exactly from its active
                        constraints and checked to lie in Q(k)
  convex_minorant, minorant_lemma, chain_check   the steps of the proof
  random_Q_point        a vertex of Q(k) from an exact LP with a random objective
  duality_necessity, separating_potential   both directions of the grid duality lemma, exactly
"""
import itertools
from fractions import Fraction as Fr


# ---------------- exact simplex ----------------
def simplex(A, b, c=None):
    """max c.x s.t. A x = b, x >= 0 (exact, Bland's rule). Returns (status, x) with status 'optimal', 'infeasible'
    or 'unbounded'. With c = None only feasibility is decided (phase I) and a feasible x is returned."""
    m = len(A); n = len(A[0])
    A = [[Fr(v) for v in row] for row in A]; b = [Fr(v) for v in b]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-v for v in A[i]]; b[i] = -b[i]
    T = [A[i] + [Fr(1) if j == i else Fr(0) for j in range(m)] + [b[i]] for i in range(m)]
    basis = list(range(n, n + m)); N = n + m

    def pivot(r, col):
        pv = T[r][col]; T[r] = [v / pv for v in T[r]]
        for i in range(m):
            if i != r and T[i][col] != 0:
                f = T[i][col]; Ti = T[i]; Tr = T[r]
                T[i] = [Ti[j] - f * Tr[j] for j in range(N + 1)]
        basis[r] = col

    def run(cost, allowed):
        while True:
            cb = [cost[basis[i]] for i in range(m)]
            enter = None
            for j in allowed:
                if j in basis:
                    continue
                rc = cost[j] - sum(cb[i] * T[i][j] for i in range(m) if cb[i] != 0)
                if rc > 0:
                    enter = j; break
            if enter is None:
                return "optimal"
            best = None
            for i in range(m):
                if T[i][enter] > 0:
                    ratio = T[i][N] / T[i][enter]
                    if best is None or ratio < best[0] or (ratio == best[0] and basis[i] < basis[best[1]]):
                        best = (ratio, i)
            if best is None:
                return "unbounded"
            pivot(best[1], enter)
    run([Fr(0)] * n + [Fr(-1)] * m, range(N))
    if sum(T[i][N] for i in range(m) if basis[i] >= n) != 0:
        return "infeasible", None
    for i in range(m):
        if basis[i] >= n:
            for j in range(n):
                if T[i][j] != 0:
                    pivot(i, j); break
    if c is not None:
        if run([Fr(v) for v in c] + [Fr(0)] * m, range(n)) == "unbounded":
            return "unbounded", None
    x = [Fr(0)] * n
    for i in range(m):
        if basis[i] < n:
            x[basis[i]] = T[i][N]
    return "optimal", x


# ---------------- Q(k), convex order, doubly stochastic matrices ----------------
def hinge_rhs(k, j):
    return sum(max(2 * i - 2 * j, 0) for i in range(k))


def on_affine_hull(x, k):
    return sum(x) == k and sum(t * x[t] for t in range(2 * k - 1)) == k * (k - 1)


def hinges_hold(x, k):
    n = 2 * k - 1
    return all(sum(max(t - 2 * j, 0) * x[t] for t in range(n)) <= hinge_rhs(k, j) for j in range(1, k - 1))


def in_Q(x, k, cap=1):
    return (len(x) == 2 * k - 1 and all(0 <= v and (cap is None or v <= cap) for v in x) and on_affine_hull(x, k)
            and hinges_hold(x, k))


def satisfies_cx(x, k):
    """T <=cx 2U_k at every integer m (not only the even hinges), with equal mass and mean."""
    n = 2 * k - 1
    return on_affine_hull(x, k) and all(sum(max(t - m, 0) * x[t] for t in range(n)) <= sum(max(2 * i - m, 0) for i in range(k))
                                        for m in range(0, 2 * k))


def xray(sigma):
    k = len(sigma); x = [0] * (2 * k - 1)
    for i, j in enumerate(sigma):
        x[i + j] += 1
    return x


def find_Z(x, k):
    """A doubly stochastic k x k matrix whose anti-diagonal sums are x (exact LP), or None."""
    idx = lambda i, j: i * k + j
    A = []; b = []
    for i in range(k):
        r = [0] * (k * k)
        for j in range(k):
            r[idx(i, j)] = 1
        A.append(r); b.append(1)
    for j in range(k):
        r = [0] * (k * k)
        for i in range(k):
            r[idx(i, j)] = 1
        A.append(r); b.append(1)
    for t in range(2 * k - 1):
        r = [0] * (k * k)
        for i in range(k):
            if 0 <= t - i < k:
                r[idx(i, t - i)] = 1
        A.append(r); b.append(x[t])
    st, z = simplex(A, b)
    if st != "optimal":
        return None
    return [[z[idx(i, j)] for j in range(k)] for i in range(k)]


def find_Z_fast(x, k):
    """As find_Z, but a floating-point linear program (scipy, dual simplex) proposes a basic solution; its support is
    then solved exactly and the result is verified with check_Z. Falls back to the exact simplex when the proposal
    does not verify. The returned matrix is always exact."""
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError:
        return find_Z(x, k)
    n2 = k * k; A = []; b = []
    for i in range(k):
        A.append([1 if p // k == i else 0 for p in range(n2)]); b.append(1)
    for j in range(k):
        A.append([1 if p % k == j else 0 for p in range(n2)]); b.append(1)
    for t in range(2 * k - 1):
        A.append([1 if p // k + p % k == t else 0 for p in range(n2)]); b.append(x[t])
    r = linprog(np.zeros(n2), A_eq=np.array(A, float), b_eq=np.array([float(v) for v in b]), bounds=(0, None),
                method="highs-ds")
    if r.status == 0:
        supp = [p for p in range(n2) if r.x[p] > 1e-9]
        rows = [[Fr(A[i][p]) for p in supp] for i in range(len(A))]
        z = _least_rows_solve(rows, [Fr(v) for v in b], len(supp))
        if z is not None:
            full = [Fr(0)] * n2
            for p, v in zip(supp, z):
                full[p] = v
            Z = [[full[i * k + j] for j in range(k)] for i in range(k)]
            if check_Z(Z, x, k):
                return Z
    return find_Z(x, k)


def check_Z(Z, x, k):
    return (all(v >= 0 for row in Z for v in row) and all(sum(Z[i]) == 1 for i in range(k))
            and all(sum(Z[i][j] for i in range(k)) == 1 for j in range(k))
            and all(sum(Z[i][t - i] for i in range(k) if 0 <= t - i < k) == x[t] for t in range(2 * k - 1)))


# ---------------- vertices of Q(k) ----------------
def _solve_square(rows, rhs):
    """Unique solution of a square system (Fractions), or None if singular."""
    n = len(rows)
    M = [list(r) + [h] for r, h in zip(rows, rhs)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [v / pv for v in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]


def vertices_exact(k, cap=1):
    """All vertices of Q(k) (cap = 1) or of {x >= 0, affine hull, hinges} (cap = None), exactly.

    A vertex is a feasible point at which the active constraints have rank n = 2k - 1. Every vertex is therefore found
    by choosing n - 2 active inequalities (each coordinate at its lower bound 0, at its upper bound cap, or free, and a
    subset H of the hinges) such that, after substituting the fixed coordinates, the two equations of the affine hull
    and the hinges in H determine the free coordinates uniquely; the solution is kept if it is feasible."""
    n = 2 * k - 1; hinges = list(range(1, k - 1)); out = set()
    states = (Fr(0),) if cap is None else (Fr(0), Fr(cap))
    for h in range(0, min(len(hinges), n - 2) + 1):
        nfix = n - 2 - h
        if nfix < 0:
            continue
        for H in itertools.combinations(hinges, h):
            for fixed in itertools.combinations(range(n), nfix):
                free = [t for t in range(n) if t not in fixed]
                rows = [[Fr(1)] * len(free), [Fr(t) for t in free]] + [[Fr(max(t - 2 * j, 0)) for t in free] for j in H]
                for vals in itertools.product(states, repeat=nfix):
                    fx = dict(zip(fixed, vals))
                    rhs = [k - sum(vals), k * (k - 1) - sum(t * v for t, v in fx.items())]
                    rhs += [hinge_rhs(k, j) - sum(max(t - 2 * j, 0) * v for t, v in fx.items()) for j in H]
                    sol = _solve_square(rows, [Fr(r) for r in rhs])
                    if sol is None:
                        break      # the matrix does not depend on vals: singular for all of them
                    x = [Fr(0)] * n
                    for t, v in fx.items():
                        x[t] = v
                    for t, v in zip(free, sol):
                        x[t] = v
                    if in_Q(x, k, cap):
                        out.add(tuple(x))
    return out


def _rows(k, cap):
    n = 2 * k - 1; ineq = []
    for t in range(n):
        r = [0] * (n + 1); r[1 + t] = 1; ineq.append(r)
        if cap is not None:
            r = [cap] + [0] * n; r[1 + t] = -1; ineq.append(r)
    for j in range(1, k - 1):
        ineq.append([hinge_rhs(k, j)] + [-max(t - 2 * j, 0) for t in range(n)])
    eq = [[-k] + [1] * n, [-k * (k - 1)] + list(range(n))]
    return ineq, eq


def vertices_cdd(k, cap=1):
    """Vertices of Q(k) with pycddlib (floating point), each recovered exactly from its active constraints and checked
    to lie in Q(k). Raises ImportError without pycddlib."""
    import cdd
    n = 2 * k - 1; ineq, eq = _rows(k, cap)
    M = cdd.matrix_from_array(ineq + eq, rep_type=cdd.RepType.INEQUALITY, lin_set=set(range(len(ineq), len(ineq) + 2)))
    out = set()
    for g in cdd.copy_generators(cdd.polyhedron_from_matrix(M)).array:
        vf = g[1:]
        act = [r for r in ineq if abs(r[0] + sum(a * b for a, b in zip(r[1:], vf))) < 1e-7]
        rows = [[Fr(c) for c in r[1:]] for r in eq + act]; rhs = [Fr(-r[0]) for r in eq + act]
        x = _least_rows_solve(rows, rhs, n)
        if x is None or not in_Q(x, k, cap):
            raise AssertionError(f"cdd vertex not recovered exactly at k = {k}")
        out.add(tuple(x))
    return out


def _least_rows_solve(rows, rhs, n):
    M = [list(r) + [h] for r, h in zip(rows, rhs)]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]; M[r] = [v / pv for v in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    if r < n or any(M[i][n] != 0 for i in range(r, len(M))):
        return None
    x = [Fr(0)] * n
    for i, c in enumerate(piv):
        x[c] = M[i][n]
    return x


# ---------------- the proof ----------------
def v_of(u):
    k = len(u)
    return [min(u[i] + u[t - i] for i in range(k) if 0 <= t - i < k) for t in range(2 * k - 1)]


def convex_minorant(u):
    """Contact points of the lower convex envelope of (i, u(i)) and phi(t/2) for t = 0..2k-2."""
    k = len(u); hull = []
    for i in range(k):
        while len(hull) >= 2:
            a, b = hull[-2], hull[-1]
            if (u[b] - u[a]) * (i - a) >= (u[i] - u[a]) * (b - a):
                hull.pop()
            else:
                break
        hull.append(i)
    phi2 = [None] * (2 * k - 1)
    for a, b in zip(hull, hull[1:]):
        for t in range(2 * a, 2 * b + 1):
            phi2[t] = u[a] + (u[b] - u[a]) * (Fr(t, 2) - a) / (b - a)
    if k == 1:
        phi2[0] = Fr(u[0])
    return hull, phi2


def minorant_lemma(u):
    """Lemma 3.3 for u, with its pointwise bounds; returns (lhs, rhs, v, phi(t/2), phi(i)) or raises."""
    k = len(u); v = v_of(u); hull, phi2 = convex_minorant(u)
    phi = [phi2[2 * i] for i in range(k)]
    if not (all(phi[i] <= u[i] for i in range(k)) and all(phi[i] == u[i] for i in hull)):
        raise AssertionError("not a minorant")
    convex = all(phi[i - 1] - 2 * phi[i] + phi[i + 1] >= 0 for i in range(1, k - 1))
    if not convex:
        raise AssertionError("minorant not convex")
    w = [u[i] - phi[i] for i in range(k)]
    for a, b in zip(hull, hull[1:]):
        for t in range(2 * a, a + b + 1):
            if not v[t] <= 2 * phi2[t] + w[t - a]:
                raise AssertionError("pointwise bound (left)")
        for t in range(a + b, 2 * b + 1):
            if not v[t] <= 2 * phi2[t] + w[t - b]:
                raise AssertionError("pointwise bound (right)")
    lhs = sum(max(v[t] - 2 * phi2[t], 0) for t in range(2 * k - 1)); rhs = 2 * sum(w)
    if not lhs <= rhs:
        raise AssertionError("minorant lemma")
    return lhs, rhs, v, phi2, phi


def chain_check(x, u):
    """The three-line chain of the proof, in count coordinates (multiply the displayed chain by k):
    sum v x <= sum (2 phi(t/2) + (v - 2 phi(t/2))^+) x <= sum 2 phi(t/2) x + sum (v - 2 phi(t/2))^+
            <= 2 sum phi(i) + sum (v - 2 phi(t/2))^+ <= 2 sum u(i)."""
    k = len(u); n = 2 * k - 1
    lhs, rhs, v, phi2, phi = minorant_lemma(u)
    s1 = sum(v[t] * x[t] for t in range(n))
    s2 = sum(2 * phi2[t] * x[t] + max(v[t] - 2 * phi2[t], 0) * x[t] for t in range(n))
    s3 = sum(2 * phi2[t] * x[t] for t in range(n)) + lhs
    s4 = 2 * sum(phi) + lhs
    s5 = 2 * sum(u)
    return s1 <= s2 <= s3 <= s4 <= s5


def random_u(k, rng):
    mode = rng.randrange(4)
    if mode == 0:
        return [Fr(rng.randint(-5, 5)) for _ in range(k)]
    if mode == 1:
        return [Fr(rng.randint(-60, 60), rng.randint(1, 9)) for _ in range(k)]
    if mode == 2:
        J = sorted(rng.sample(range(1, k), rng.randint(0, k - 1))) if k > 1 else []
        return [Fr(sum(1 for j in J if i >= j) * rng.choice([1, -1])) for i in range(k)]
    return [Fr(rng.randint(-3, 3) * (i - rng.randint(0, k)) ** 2 + rng.randint(-4, 4)) for i in range(k)]


def random_Q_point(k, rng):
    """A vertex of Q(k): exact simplex with slack variables and a random rational objective."""
    n = 2 * k - 1; H = list(range(1, k - 1)); nv = 2 * n + len(H)
    A = [[1] * n + [0] * (n + len(H)), list(range(n)) + [0] * (n + len(H))]; b = [k, k * (k - 1)]
    for t in range(n):
        r = [0] * nv; r[t] = 1; r[n + t] = 1; A.append(r); b.append(1)
    for a, j in enumerate(H):
        r = [max(t - 2 * j, 0) for t in range(n)] + [0] * n + [0] * len(H); r[2 * n + a] = 1; A.append(r); b.append(hinge_rhs(k, j))
    c = [Fr(rng.randint(-50, 50), rng.randint(1, 7)) for _ in range(n)] + [0] * (n + len(H))
    st, x = simplex(A, b, c)
    if st != "optimal":
        raise AssertionError("Q(k) is empty or unbounded")
    return x[:n]


# ---------------- grid duality ----------------
def duality_necessity(Z, u):
    """The necessity computation of Lemma 3.2 for a doubly stochastic Z with X-ray x:
    sum_t v(t) x(t) = sum_ij Z_ij v(i+j) <= sum_ij Z_ij (u(i) + u(j)) = 2 sum_i u(i)."""
    k = len(u); v = v_of(u)
    x = [sum(Z[i][t - i] for i in range(k) if 0 <= t - i < k) for t in range(2 * k - 1)]
    a = sum(v[t] * x[t] for t in range(2 * k - 1))
    b = sum(Z[i][j] * v[i + j] for i in range(k) for j in range(k))
    c = sum(Z[i][j] * (u[i] + u[j]) for i in range(k) for j in range(k))
    return a == b <= c == 2 * sum(u)


def separating_potential(x, k):
    """The sufficiency direction of Lemma 3.2, constructively: for x >= 0 outside hat-P(k), an exact u with
    sum_t v(t) x(t) > 2 sum_i u(i). Solves max sum_t c(t) x(t) - 2 sum_i u(i) subject to c(t) <= u(i) + u(j) for
    i + j = t and the normalization (objective) <= 1, by the exact simplex method (free variables split into two
    nonnegative parts). Returns u (then v >= c, so u separates) or None if the optimum is <= 0."""
    n = 2 * k - 1; nu = 2 * k; nc = 2 * n; npairs = k * k
    nvar = nu + nc + npairs + 1
    A = []; b = []
    for i in range(k):
        for j in range(k):
            r = [Fr(0)] * nvar
            r[nu + 2 * (i + j)] += 1; r[nu + 2 * (i + j) + 1] -= 1          # c(i+j)
            r[2 * i] -= 1; r[2 * i + 1] += 1; r[2 * j] -= 1; r[2 * j + 1] += 1   # - u(i) - u(j)
            r[nu + nc + i * k + j] = 1                                         # slack
            A.append(r); b.append(Fr(0))
    obj = [Fr(0)] * nvar
    for i in range(k):
        obj[2 * i] -= 2; obj[2 * i + 1] += 2
    for t in range(n):
        obj[nu + 2 * t] += x[t]; obj[nu + 2 * t + 1] -= x[t]
    r = list(obj); r[-1] = Fr(1); A.append(r); b.append(Fr(1))
    st, z = simplex(A, b, obj)
    if st != "optimal":
        return None
    val = sum(o * zz for o, zz in zip(obj, z))
    if val <= 0:
        return None
    u = [z[2 * i] - z[2 * i + 1] for i in range(k)]
    v = v_of(u)
    if not sum(v[t] * x[t] for t in range(n)) > 2 * sum(u):
        raise AssertionError("separating potential does not separate")
    return u
