"""Staircase inequalities for laws with density at most c (c = 2 in the main case): exact dual certificates.

All quantities are exact (fractions.Fraction). Notation as in Lemma 5.1 of the paper: for a staircase
u = sum_q 1{x >= eta_q} with distinct thresholds 0 < eta_1 < ... < eta_r < 1, N = u box u is the step function
N(s) = #{l : xi_l <= s} with the max-plus jumps xi_1 <= ... <= xi_2r (`jumps`), and a convex phi on [0, 2] certifies
the staircase inequality at cap c if

    J_c(phi) = (1/2) int_0^2 phi + c int_0^2 (N - phi)^+  <=  2 int_0^1 u = 2 sum_q (1 - eta_q).

`objective` evaluates J_c(phi) for a piecewise linear phi by exact integration, splitting [0, 2] at every jump of N
and every kink of phi and splitting each linear piece of N - phi at its zero.

Certificates (each phi is a nonnegative combination of hinges (s - t)^+, so convex):
  two_step(a, b)              two thresholds 0 < a < b < 1, cap 2, three cases by rho = (b - a)/(1 - b);
                              margin >= min(b - a, 1 - b)/21
  nonincreasing_gaps(g, c)    gaps g_1 >= ... >= g_r (the last gap is 1 - eta_r), any cap c; closed-form margin
  equally_spaced(r, d, c)     thresholds 1 - (r - j + 1) d, j = 1..r; closed-form margin
Sharpness: `equally_spaced_sharpness(r, c)` is a law with density c on 2r intervals, mass and mean one,
S <=cx 2U, violating the equally spaced staircase by (c - r - 1)/(2c(2r + 1)) (positive for c > r + 1).
`coupling_profiles` is an explicit piecewise affine coupling of two uniforms whose sum has density 2 on
[1/2, 3/4] u [5/4, 3/2]. `hinge_sum` is the sum of the single-step certificates of Proposition 5.2 for a list of
thresholds (it fails for the five thresholds 1/8, 1/4, 1/2, 5/8, 3/4). `verify_catalogue_record` checks a stored
certificate given by node values on the grid 1/(Km), with N derived independently from Minkowski sums of the
constant-level intervals of u.
"""
from fractions import Fraction as Fr


# ---------------- exact integration ----------------
def pos_integral(a, b, y0, y1):
    """int_a^b of (the linear function from y0 at a to y1 at b)^+."""
    if y0 >= 0 and y1 >= 0:
        return (b - a) * (y0 + y1) / 2
    if y0 <= 0 and y1 <= 0:
        return Fr(0)
    hi, lo = max(y0, y1), min(y0, y1)
    return (b - a) * hi * hi / (2 * (hi - lo))


def jumps(eta):
    """Max-plus jumps xi_1..xi_2r of N = u box u for distinct sorted thresholds eta (Lemma 4.1)."""
    r = len(eta); e = list(eta) + [Fr(1)]
    return [max(e[q - 1] + e[l - q] for q in range(1, r + 2) if 1 <= l + 1 - q <= r + 1) for l in range(1, 2 * r + 1)]


def hinge_phi(hinges):
    """phi(s) = sum slope * (s - t)^+ for hinges [(t, slope), ...]."""
    return lambda s: sum(sl * max(Fr(0), s - t) for t, sl in hinges)


def objective(xi, phi, kinks, c):
    """J_c(phi) on [0, 2] for N(s) = #{l : xi_l <= s} and a phi linear between consecutive points of kinks."""
    nodes = sorted({Fr(0), Fr(2), *(x for x in xi if 0 < x < 2), *(t for t in kinks if 0 < t < 2)})
    J = Fr(0)
    for a, b in zip(nodes, nodes[1:]):
        n = sum(1 for x in xi if x <= a)
        pa, pb = phi(a), phi(b)
        J += (b - a) * (pa + pb) / 4 + c * pos_integral(a, b, n - pa, n - pb)
    return J


def budget(eta):
    return 2 * sum(1 - e for e in eta)


# ---------------- two thresholds, cap 2 ----------------
def two_step(a, b):
    """Hinges, case label and predicted margin of the two-step certificate at cap 2."""
    a, b = Fr(a), Fr(b)
    assert 0 < a < b < 1
    d, e = b - a, 1 - b
    rho = d / e
    if rho <= Fr(3, 4):
        q = 2 / (rho + 1)
        alpha = Fr(8, 9) * (4 + q)
        hinges = [(2 - e * alpha / q, q / e)]
        case = "A"
        margin = e * (18 * rho * (rho + 3) - max(Fr(0), 13 * rho - 3) ** 2) / (81 * (rho + 1))
    elif rho <= 1:
        alpha = (64 + 16 * rho) / 17
        hinges = [(2 - e * alpha, 1 / e)]
        case = "B"
        margin = e * (-38 * rho ** 2 + 70 * rho - 30) / 17
    else:
        hinges = [(2 - 2 * e - Fr(8, 3) * d, 1 / d), (2 * b, 1 / e - 1 / d)]
        case = "C"
        margin = e / 9
    return hinges, case, margin


def two_step_check(a, b):
    """Exact J_2 of the two-step certificate; the actual margin is at least the predicted one (with equality when
    every hinge lies in [0, 2]) and at least min(b - a, 1 - b)/21."""
    a, b = Fr(a), Fr(b)
    hinges, case, predicted = two_step(a, b)
    xi = jumps([a, b])
    assert xi == [2 * a, a + b, max(1 + a, 2 * b), 1 + b]
    J = objective(xi, hinge_phi(hinges), [t for t, _ in hinges], Fr(2))
    actual = budget([a, b]) - J
    ok = (all(sl >= 0 for _, sl in hinges) and actual >= predicted > 0
          and (actual == predicted or any(t < 0 for t, _ in hinges))
          and actual >= min(b - a, 1 - b) / 21)
    return {"case": case, "hinges": hinges, "objective": J, "predicted_margin": predicted, "margin": actual, "ok": ok}


# ---------------- nonincreasing gaps, any cap ----------------
def nonincreasing_gaps(gaps, c):
    """Certificate for thresholds eta_0 = 1 - D, eta_j = eta_{j-1} + g_j (j < r), with gaps g_1 >= ... >= g_r and
    D = sum g < 1. Returns (margin = budget - J_c(phi), closed-form margin, t)."""
    gaps = [Fr(g) for g in gaps]; c = Fr(c)
    r = len(gaps); D = sum(gaps); d1 = gaps[0]
    assert 0 < D < 1 and all(gaps[j] >= gaps[j + 1] > 0 for j in range(r - 1))
    eta = [1 - D]
    for g in gaps[:-1]:
        eta.append(eta[-1] + g)
    xi = jumps(eta)
    balanced = []
    for j in range(r):
        balanced += [2 * eta[j], eta[j] + (eta[j + 1] if j + 1 < r else Fr(1))]
    assert xi == balanced, "jumps are not the balanced pair sums"
    alpha = 2 * D * (2 * c - 1) / (4 * c * D + d1)
    t = 2 * eta[0] - alpha * d1
    hinges = [(t, 1 / d1)] + [(2 * eta[j], 1 / gaps[j] - 1 / gaps[j - 1]) for j in range(1, r)]
    J = objective(xi, hinge_phi(hinges), [h for h, _ in hinges], c)
    closed = D * (D + (1 - c) * d1) / (4 * c * D + d1) + max(Fr(0), -t) ** 2 / (4 * d1)
    return budget(eta) - J, closed, t


# ---------------- equally spaced thresholds ----------------
def equally_spaced_eta(r, d):
    return [1 - (r - j + 1) * Fr(d) for j in range(1, r + 1)]


def equally_spaced(r, d, c):
    """(margin, closed-form margin, J, closed-form J) of the certificate phi = (s - t)^+/d at cap c."""
    d, c = Fr(d), Fr(c)
    eta = equally_spaced_eta(r, d); xi = jumps(eta)
    assert xi == [2 - 2 * r * d + (l - 1) * d for l in range(1, 2 * r + 1)]
    alpha = 2 * r * (2 * c - 1) / (4 * c * r + 1)
    t = 2 - (2 * r + alpha) * d
    J = objective(xi, hinge_phi([(t, 1 / d)]), [t], c)
    tm = max(Fr(0), -t)
    J_closed = c * r * (2 * r + 1) ** 2 * d / (4 * c * r + 1) - tm * tm / (4 * d)
    margin_closed = r * (r + 1 - c) * d / (4 * c * r + 1) + tm * tm / (4 * d)
    return budget(eta) - J, margin_closed, J, J_closed


def stop_loss(intervals, c, x):
    """E (S - x)^+ for density c on the given disjoint intervals."""
    out = Fr(0)
    for a, b in intervals:
        if x <= a:
            out += c * (b - a) * ((a + b) / 2 - x)
        elif x < b:
            out += c * (b - x) ** 2 / 2
    return out


def min_convex_order_gap(intervals, c):
    """min over [0, 2] of D(x) = (2 - x)^2/4 - E(S - x)^+ (piecewise quadratic: endpoints and stationary points)."""
    nodes = sorted({Fr(0), Fr(2), *(x for iv in intervals for x in iv)})
    best = None
    for a, b in zip(nodes, nodes[1:]):
        mid = (a + b) / 2
        dens = c if any(lo < mid < hi for lo, hi in intervals) else Fr(0)
        cdf_a = sum(c * max(Fr(0), min(a, hi) - lo) for lo, hi in intervals)
        deriv = a / 2 - cdf_a                   # D'(a+) = -(2 - a)/2 + P(S > a) = a/2 - P(S <= a)
        second = Fr(1, 2) - dens
        pts = [a, b]
        if second != 0:
            x0 = a - deriv / second
            if a < x0 < b:
                pts.append(x0)
        for x in pts:
            val = (2 - x) ** 2 / 4 - stop_loss(intervals, c, x)
            best = val if best is None or val < best else best
    return best


def equally_spaced_sharpness(r, c):
    """The law with density c on [x, x + delta] at the 2r jumps x of the equally spaced staircase with spacing
    d = (4cr + 1)/(2cr(2r + 1)), delta = 1/(2cr). Returns (eta, intervals, violation, min convex-order gap)."""
    c = Fr(c)
    d = (4 * c * r + 1) / (2 * c * r * (2 * r + 1)); delta = 1 / (2 * c * r)
    eta = equally_spaced_eta(r, d)
    xi = [2 - 2 * r * d + (l - 1) * d for l in range(1, 2 * r + 1)]
    assert xi == jumps(eta)
    intervals = [(x, x + delta) for x in xi]
    assert 0 < eta[0] and eta[-1] < 1
    assert intervals[0][0] > 0 and intervals[-1][1] < 2
    assert all(b < a2 for (_, b), (a2, _) in zip(intervals, intervals[1:]))
    assert c * sum(b - a for a, b in intervals) == 1
    assert c * sum((b * b - a * a) / 2 for a, b in intervals) == 1
    EN = sum(c * (b - a) * sum(1 for z in xi if z <= a) for a, b in intervals)
    return eta, intervals, EN - budget(eta), min_convex_order_gap(intervals, c)


def perturbed_margin(r, d, eps_list):
    """Margin at cap 2 of the equally spaced certificate phi (unchanged) for the thresholds eta_j + eps_j."""
    d = Fr(d); c = Fr(2)
    eta = equally_spaced_eta(r, d)
    alpha = 2 * r * (2 * c - 1) / (4 * c * r + 1)
    t = 2 - (2 * r + alpha) * d
    new = [e + x for e, x in zip(eta, eps_list)]
    assert all(0 < x < 1 for x in new) and all(x < y for x, y in zip(new, new[1:]))
    return budget(new) - objective(jumps(new), hinge_phi([(t, 1 / d)]), [t], c)


# ---------------- an explicit coupling ----------------
# Branches (x-interval, conditional probability, intercept, slope): given X = x in the interval, Y = intercept +
# slope * x with the stated probability.
COUPLING = [(Fr(0), Fr(1, 4), Fr(1), Fr(3, 4), Fr(-2)),
            (Fr(1, 4), Fr(3, 4), Fr(1, 2), Fr(3, 8), Fr(-1, 2)),
            (Fr(1, 4), Fr(3, 4), Fr(1, 2), Fr(9, 8), Fr(-1, 2)),
            (Fr(3, 4), Fr(1), Fr(1), Fr(9, 4), Fr(-2))]


def coupling_profiles(branches=COUPLING):
    """Piecewise constant densities of X, Y and S = X + Y, as lists of (lo, hi, density)."""
    def pushed(coord):
        out = []
        for lo, hi, p, b0, b1 in branches:
            if coord == "X":
                b0, b1 = Fr(0), Fr(1)
            elif coord == "S":
                b1 = b1 + 1
            a, b = sorted((b0 + b1 * lo, b0 + b1 * hi))
            out.append((a, b, p / abs(b1)))
        return out

    def profile(iv):
        pts = sorted({v for a, b, _ in iv for v in (a, b)})
        return [(a, b, sum(h for lo, hi, h in iv if lo <= (a + b) / 2 < hi)) for a, b in zip(pts, pts[1:])]
    return {k: profile(pushed(k)) for k in ("X", "Y", "S")}


def profile_moment(profile, power):
    return sum(h * (b ** (power + 1) - a ** (power + 1)) / (power + 1) for a, b, h in profile)


def profile_tail(profile, x):
    return sum(h * max(Fr(0), b - max(a, x)) for a, b, h in profile)


def profile_min_convex_order_gap(profile):
    """min over [0, 2] of (2 - x)^2/4 - E(S - x)^+ for a piecewise constant density on [0, 2]."""
    def sl(x):
        return sum(h * max(Fr(0), b - max(a, x)) * (b + max(a, x) - 2 * x) / 2 for a, b, h in profile if b > x)
    nodes = sorted({Fr(0), Fr(2), *(v for a, b, _ in profile for v in (a, b))})
    best = None
    for a, b in zip(nodes, nodes[1:]):
        dens = sum(h for lo, hi, h in profile if lo <= a and b <= hi)
        cdf_a = sum(h * max(Fr(0), min(a, hi) - lo) for lo, hi, h in profile)
        deriv = a / 2 - cdf_a; second = Fr(1, 2) - dens
        pts = [a, b]
        if second != 0:
            x0 = a - deriv / second
            if a < x0 < b:
                pts.append(x0)
        for x in pts:
            val = (2 - x) ** 2 / 4 - sl(x)
            best = val if best is None or val < best else best
    return best


# ---------------- sum of single-step certificates ----------------
def hinge_sum(eta, c=Fr(2)):
    """(J_c(phi), budget) for phi = sum_q (s - t_q)^+/(1 - eta_q), t_q = (8 eta_q - 2)/3 (Proposition 5.2)."""
    eta = [Fr(e) for e in eta]
    hinges = [((8 * e - 2) / 3, 1 / (1 - e)) for e in eta]
    return objective(jumps(eta), hinge_phi(hinges), [t for t, _ in hinges], c), budget(eta)


# ---------------- stored certificates on a grid ----------------
def min_convolution_cells(J, K, m):
    """N on the cells [i/(Km), (i+1)/(Km)) from the Minkowski sums of the constant-level intervals of u (independent
    of the max-plus formula)."""
    if not (K >= 2 and m >= 1 and list(J) == sorted(set(J)) and all(0 < j < K for j in J)):
        raise ValueError("invalid staircase thresholds")
    bounds = [0, *J, K]
    levels = [(bounds[p], bounds[p + 1], p) for p in range(len(bounds) - 1)]
    base = []
    for t in range(2 * K):
        s2 = 2 * t + 1
        vals = [p + q for a, b, p in levels for c, d, q in levels if 2 * (a + c) <= s2 < 2 * (b + d)]
        base.append(min(vals))
    return [v for v in base for _ in range(m)]


def verify_catalogue_record(rec, c=Fr(2)):
    """(certified, margin) for a record {J, K, m, phi: node values at i/(Km), margin}: phi convex and
    J_c(phi) <= budget, and the stored margin equals the recomputed one."""
    J, K, m = list(rec["J"]), int(rec["K"]), int(rec["m"])
    phi = [Fr(x) for x in rec["phi"]]
    if len(phi) != 2 * K * m + 1:
        raise ValueError("wrong number of nodes")
    sl = [b - a for a, b in zip(phi, phi[1:])]
    convex = all(a <= b for a, b in zip(sl, sl[1:]))
    N = min_convolution_cells(J, K, m)
    h = Fr(1, K * m)
    cost = sum((a + b) * h / 4 + c * pos_integral(Fr(0), h, n - a, n - b) for n, a, b in zip(N, phi, phi[1:]))
    margin = 2 * sum(Fr(K - j, K) for j in J) - cost
    if "margin" in rec and Fr(rec["margin"]) != margin:
        raise ValueError("stored margin differs from the recomputed margin")
    return convex and margin >= 0, margin
