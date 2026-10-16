"""Displacement windows of permutations (Section 6).

Lemma S3.1 (block structure). tau permutation of {0..m-1} with i - tau(i) >= -1 for all i  <=>  tau is a product of
consecutive cycles (s -> s+1 -> ... -> s+L-1 -> s). A block of length L has displacements (-1)^(L-1) and +(L-1).
So with window {-1, 0, 1, ..., w-2} the blocks have length 1..w-1, and the displacement histogram of tau is
x_{-1} = sum_L (L-1) n_L,  x_{L-1} = n_L (L >= 2),  x_0 = n_1,  with sum_L L n_L = m.

For w = 4 (window {-1,0,1,2}): Q(m) = conv{ (p + 2q, f, p, q) : f + 2p + 3q = m, f,p,q in Z>=0 }.
Candidate vertex rule (Lemma S3.2): the vertices of Q(m) are among
  O = (m,0,0), A = (m mod 2, floor(m/2), 0), B = (m mod 3, 0, floor(m/3)),
  E_p, E_q = the solutions of 2p + 3q = m with q minimal / maximal (if any),   written as (f,p,q).
"""
import itertools
import math
from fractions import Fraction as Fr


def _nint(x):
    """Nearest integer (halves upward)."""
    return math.floor(x + 0.5)


def lattice_points(m):
    return [(m - 2 * p - 3 * q, p, q) for q in range(m // 3 + 1) for p in range((m - 3 * q) // 2 + 1)]

def hull2d(points):
    """Exact convex hull (vertices only) of points given as (p, q) integer pairs; Andrew monotone chain."""
    pts = sorted(set(points))
    if len(pts) <= 2: return pts
    cross = lambda o, a, b: (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0: lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0: upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def vertices_exact(m):
    """Vertices of Q(m) as (f, p, q); the plane f + 2p + 3q = m is parametrised by (p, q)."""
    V = hull2d([(p, q) for (_, p, q) in lattice_points(m)])
    return sorted((m - 2 * p - 3 * q, p, q) for (p, q) in V)

def candidate_vertices(m):
    C = {(m, 0, 0), (m % 2, m // 2, 0), (m % 3, 0, m // 3)}
    sols = [(p, q) for q in range(m // 3 + 1) for p in [(m - 3 * q) // 2] if 2 * p + 3 * q == m]
    if sols:
        C.add((0, *min(sols, key=lambda s: s[1]))); C.add((0, *max(sols, key=lambda s: s[1])))
    return C

def to_histogram(fpq):
    f, p, q = fpq
    return (p + 2 * q, f, p, q)      # (x_{-1}, x_0, x_1, x_2)


def Q_law_vertices(m, mirror=False):
    """Vertices of Q(m)/m as laws on offsets (-1, 0, 1, 2) relative to the fiber mean (mirror: offsets (-2,-1,0,1))."""
    out = []
    for (f, p, q) in vertices_exact(m):
        x = [Fr(t, m) for t in to_histogram((f, p, q))]          # (x_-1, x_0, x_1, x_2)
        out.append(tuple(x[::-1]) if mirror else tuple(x))        # mirror: offsets (-2,-1,0,1) in increasing order
    return out

def fiber_on_window(m, mean_pos):
    """Law set of P_m restricted to a 4-window, where mean_pos = index of the fiber mean inside the window (0..3).
    Returns vertex list of laws on the 4 window points."""
    if mean_pos in (0, 3):
        e = [Fr(0)] * 4; e[mean_pos] = Fr(1); return [tuple(e)]
    if mean_pos == 1:                       # window offsets (-1, 0, 1, 2)
        return Q_law_vertices(m)
    return Q_law_vertices(m, mirror=True)   # mean_pos == 2: offsets (-2, -1, 0, 1)

def four_atom_prediction(b, a):
    """Four-atom window corollary: list of (weight, vertex-list) pairs whose Minkowski combination is the attainable set of
    laws on the window {a, a+b, a+2b, a+3b}; None if Theorem 2.2(i) says no lattice-supported sum exists."""
    L = 1 / b
    n = math.floor(L + 1e-12); theta = L - n
    if theta < 1e-12: theta = 0.0
    c = (a / b) % 1.0
    if abs(c - 1) < 1e-9: c = 0.0
    if theta > 0 and abs(c - theta) > 1e-9:
        return None
    A = _nint(a / b - c)                       # window positions A..A+3 in coordinates S/b - c
    if theta == 0:
        k = n
        parts = []
        if c > 1e-12:  parts.append((c, k, (k - 1) - A))          # fiber 1: P_k, mean k-1, no shift
        if c < 1 - 1e-12: parts.append((1 - c, k, k - A))         # fiber 2: P_k + 1, mean k
    else:
        parts = [(theta * (n + 1) / L, n + 1, n - A), ((1 - theta) * n / L, n, n - A)]
    out = []
    for wt, m, pos in parts:
        if not (0 <= pos <= 3): return None
        out.append((wt, fiber_on_window(m, pos)))
    return out

def support_of_prediction(pred, w):
    return sum(wt * max(sum(float(x) * wi for x, wi in zip(v, w)) for v in V) for wt, V in pred)


# ---------------- exact vertex counts of displacement-histogram polytopes ----------------
def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def hull_vertex_count_3d(points):
    """Number of vertices of the convex hull of integer points in Z^3 (exact integer arithmetic).
    A point is a vertex iff it is a vertex of the 2-D hull of some facet (or of the whole set when the affine
    dimension is below 3)."""
    P = sorted(set(tuple(int(x) for x in p) for p in points))
    n = len(P)
    if n <= 2:
        return n
    base = P[0]
    D = [_sub(p, base) for p in P[1:]]
    a = next((d for d in D if any(d)), None)
    b = next((d for d in D if any(_cross(a, d))), None)
    if b is None:                                   # collinear
        return 2
    nrm = _cross(a, b)
    if all(_dot(nrm, d) == 0 for d in D):           # coplanar
        drop = next(i for i in range(3) if nrm[i]); keep = [i for i in range(3) if i != drop]
        return len(hull2d([(p[keep[0]], p[keep[1]]) for p in P]))
    isv = [False] * n
    for i, j, l in itertools.combinations(range(n), 3):
        nv = _cross(_sub(P[j], P[i]), _sub(P[l], P[i]))
        if not any(nv):
            continue
        s = [_dot(_sub(p, P[i]), nv) for p in P]
        if all(x >= 0 for x in s) or all(x <= 0 for x in s):
            idx = [c for c in range(n) if s[c] == 0]
            drop = next(c for c in range(3) if nv[c]); keep = [c for c in range(3) if c != drop]
            hv = set(hull2d([(P[c][keep[0]], P[c][keep[1]]) for c in idx]))
            for c in idx:
                if (P[c][keep[0]], P[c][keep[1]]) in hv:
                    isv[c] = True
    return sum(isv)


def centered_window_histograms(m, width=2):
    """Displacement histograms (x_{-w}, ..., x_w) of all permutations tau of {0..m-1} with |i - tau(i)| <= w."""
    H = set()

    def rec(i, used, hist):
        if i == m:
            H.add(tuple(hist)); return
        for j in range(max(0, i - width), min(m, i + width + 1)):
            if not used >> j & 1:
                h = list(hist); h[i - j + width] += 1
                rec(i + 1, used | (1 << j), h)
    rec(0, 0, [0] * (2 * width + 1))
    return H


def five_point_window_vertex_count(m):
    """Vertices of the polytope of displacement histograms with displacements in {-2..2}. The histogram lies in the
    3-dimensional affine space sum x = m, sum d x_d = 0; coordinates (x_{-2}, x_{-1}, x_1) determine it."""
    pts = [(h[0], h[1], h[3]) for h in centered_window_histograms(m, 2)]
    return hull_vertex_count_3d(pts)
