#!/usr/bin/env python3
"""Corollary 3.2 (small atoms): the statement and every step of its proof, in exact rational arithmetic.

Count coordinates on the grid k: x(t) = k P(T = t), T = kS - 1. Q(k) is the set of x with 0 <= x(t) <= 1 on the
affine hull (sum x = k, sum t x(t) = k(k-1)) that satisfy the k-2 hinge inequalities I(h_j) of Proposition 3.3; by
Proposition 3.3(i) these say T <=cx 2U_k. The corollary says Q(k) = hat-P(k) intersected with the box.

  1. Statement, both directions (k <= 5 by default, k <= 7 with --full):
     necessity: every permutation X-ray satisfies T <=cx 2U_k at every integer hinge (so hat-P(k), their convex hull,
     does);
     sufficiency: every vertex of Q(k) lies in hat-P(k), certified by an exact doubly stochastic matrix with the
     prescribed anti-diagonal sums (found by a linear program and verified in rational arithmetic).
     Vertices: exhaustive exact enumeration of bases (complete by construction) for k <= --exact-kmax (default 5;
     6 with --full); beyond that, and as a cross-check where both run, pycddlib (floating point) with every vertex
     recovered exactly from its active constraints and checked to lie in Q(k). For k > --exact-kmax the vertex list
     relies on the completeness of the floating-point enumeration. Every vertex also satisfies T <=cx 2U_k at every
     integer hinge, odd ones included (Proposition 3.3(i) on Q(k)). Vertex counts 2, 7, 44, 258, 1493, 8300 (k = 2..7).
  2. Negative control: without x <= 1 the same checker finds vertices outside hat-P(k) (k = 3, 4, 5: 1 of 5, 4 of 13,
     15 of 33), e.g. (0, 3/2, 0, 3/2, 0) at k = 3; so the bound on the atoms is needed and the checker detects failures.
     Lemma 3.2, sufficiency: each vertex outside hat-P(k) gets an exact potential u with
     sum v(t) x(t) > 2 sum u(i) (exact LP), and the vertices inside get none.
  3. Lemma 3.3 for random potentials u (integer, rational, step-shaped, convex plus noise; k = 1..40) with the
     pointwise bounds on every piece of the convex minorant, and for every signed staircase potential with k <= 10
     (k <= 14 with --full; 32,766 potentials).
  4. The chain of the proof for k = 2..8 (k = 2..12 with --full): for each k, 12 random vertices of Q(k) (exact simplex
     with random objectives) and 6 random rational mixtures of them; for each point the four inequalities of the chain
     (Lemma 3.2 reduces membership to them) for 200 random u, an exact doubly stochastic matrix (the
     statement itself on the same instances), and the necessity computation of Lemma 3.2 with that matrix
     for 20 random u.
Default (about 20 s on one core): items 1 (k <= 5, exact enumeration and pycddlib if installed), 2, 3 (4,000 random u,
staircases k <= 10), 4 (k <= 8).
--full (about 3 minutes on one core): items 1 (k <= 7; exact enumeration k <= 6, pycddlib needed for k = 7), 3 (16,000
random u, staircases k <= 14), 4 (k <= 12).
"""
import argparse
import itertools
import random
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, report

from dsxray import small_atoms as SA


def have_cdd():
    try:
        import cdd  # noqa: F401
        return True
    except ImportError:
        return False


def statement(k, exact_kmax, use_cdd):
    nec = all(SA.satisfies_cx(SA.xray(s), k) for s in itertools.permutations(range(k)))
    V = None; src = []
    if k <= exact_kmax:
        V = SA.vertices_exact(k); src.append("exact")
    if use_cdd:
        C = SA.vertices_cdd(k); src.append("cdd")
        if V is not None and C != V:
            return {"necessity": nec, "vertices": len(V), "cdd_vertices": len(C), "agree": False, "outside": None}
        V = C if V is None else V
    if V is None:
        return None
    cx_all = all(SA.satisfies_cx(list(x), k) for x in V)
    outside = sum(1 for x in V if not SA.check_Z(SA.find_Z_fast(list(x), k) or [[None]], list(x), k))
    nonint = sum(1 for x in V if any(c.denominator != 1 for c in x))
    return {"necessity": nec, "vertices": len(V), "nonintegral": nonint, "outside": outside, "sources": src, "agree": True,
            "even_hinges_give_all_integer_hinges": cx_all}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/small_atoms")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--kmax", type=int, default=None, help="statement range (default 5; 7 with --full)")
    ap.add_argument("--exact-kmax", type=int, default=None, help="exact vertex enumeration up to this k (5; 6 with --full)")
    ap.add_argument("--chain-kmax", type=int, default=None, help="proof chain range (default 8; 12 with --full)")
    ap.add_argument("--seed", type=int, default=6)
    a = ap.parse_args()
    kmax = a.kmax or (7 if a.full else 5)
    exact_kmax = a.exact_kmax or (6 if a.full else 5)
    chain_kmax = a.chain_kmax or (12 if a.full else 8)
    n_random_u = 400 if a.full else 100
    stair_kmax = 14 if a.full else 10
    cdd_ok = have_cdd()
    res = {}; checks = {}
    expected = {2: 2, 3: 7, 4: 44, 5: 258, 6: 1493, 7: 8300}

    # 1. statement
    t = Timer(); rows = {}; ok = True; skipped = []
    for k in range(2, kmax + 1):
        r = statement(k, exact_kmax, cdd_ok)
        if r is None:
            skipped.append(k); print(f"[SKIP] k = {k}: needs pycddlib (or --exact-kmax {k})", flush=True); continue
        rows[k] = r
        ok &= r["necessity"] and r["agree"] and r.get("even_hinges_give_all_integer_hinges", False) and r["outside"] == 0 and (k not in expected or r["vertices"] == expected[k])
        print(f"  k = {k}: {r['vertices']} vertices of Q(k) ({r.get('nonintegral')} nonintegral, {'+'.join(r['sources'])}), "
              f"{r['outside']} outside hat-P(k)", flush=True)
    res["statement"] = {"rows": rows, "skipped": skipped, "pycddlib": cdd_ok}
    checks["statement"] = report(f"statement for k <= {max(rows)}: permutation X-rays satisfy T <=cx 2U_k; every vertex of "
                                 "Q(k) has an exact doubly stochastic matrix", ok and bool(rows), f"({t()} s)")

    # 2. negative control
    t = Timer(); neg = {}
    for k in (3, 4, 5):
        V = SA.vertices_exact(k, cap=None)
        bad = sorted(x for x in V if SA.find_Z(list(x), k) is None)
        neg[k] = {"vertices": len(V), "outside": len(bad), "example": [str(c) for c in bad[0]] if bad else None}
        for x in V:
            u = SA.separating_potential(list(x), k); Z = SA.find_Z(list(x), k)
            neg[k].setdefault("duality_consistent", True)
            neg[k]["duality_consistent"] &= (u is None) == (Z is not None)
    res["negative_control"] = neg
    ok = [(neg[k]["vertices"], neg[k]["outside"]) for k in (3, 4, 5)] == [(5, 1), (13, 4), (33, 15)]
    ok &= neg[3]["example"] == ["0", "3/2", "0", "3/2", "0"] and all(neg[k]["duality_consistent"] for k in neg)
    checks["negative_control"] = report("without x <= 1: 1 of 5, 4 of 13, 15 of 33 vertices outside hat-P(k) (k = 3, 4, 5); "
                                        "each has an exact separating potential, the others a doubly stochastic matrix",
                                        ok, f"({t()} s)")

    # 3. minorant lemma
    rng = random.Random(a.seed)
    t = Timer(); cnt = tight = 0; ok = True
    try:
        for k in range(1, 41):
            for _ in range(n_random_u):
                lhs, rhs, *_ = SA.minorant_lemma(SA.random_u(k, rng)); cnt += 1; tight += lhs == rhs
    except AssertionError as e:
        ok = False; print("   ", e)
    res["minorant_random"] = {"potentials": cnt, "equality": tight}
    checks["minorant_random"] = report(f"Lemma 3.3 and its pointwise bounds for {cnt} random u, k = 1..40",
                                       ok, f"equality in {tight} ({t()} s)")
    t = Timer(); tot = 0; ok = True
    try:
        for k in range(1, stair_kmax + 1):
            for mask in range(1 << (k - 1)):
                J = [j for j in range(1, k) if mask >> (j - 1) & 1]
                for s in (1, -1):
                    SA.minorant_lemma([Fr(s * sum(1 for j in J if i >= j)) for i in range(k)]); tot += 1
    except AssertionError as e:
        ok = False; print("   ", e)
    res["minorant_staircases"] = {"kmax": stair_kmax, "potentials": tot}
    checks["minorant_staircases"] = report(f"Lemma 3.3 for all {tot} signed staircases, k <= {stair_kmax}",
                                           ok and tot == 2 * (2 ** stair_kmax - 1), f"({t()} s)")

    # 4. the chain and the statement on random points of Q(k)
    t = Timer(); chain = {}; ok = True
    for k in range(2, chain_kmax + 1):
        verts = [SA.random_Q_point(k, rng) for _ in range(12)]
        pts = list(verts)
        for _ in range(6):
            lam = [Fr(rng.randint(1, 9)) for _ in verts]; s = sum(lam)
            pts.append([sum(l * vx[q] for l, vx in zip(lam, verts)) / s for q in range(2 * k - 1)])
        npairs = 0; maxden = 1
        for x in pts:
            ok &= SA.in_Q(x, k) and SA.satisfies_cx(x, k)
            for _ in range(200):
                ok &= SA.chain_check(x, SA.random_u(k, rng)); npairs += 1
            Z = SA.find_Z_fast(x, k)
            ok &= Z is not None and SA.check_Z(Z, x, k)
            ok &= all(SA.duality_necessity(Z, SA.random_u(k, rng)) for _ in range(20))
            maxden = max(maxden, max(c.denominator for c in x))
        chain[k] = {"points": len(pts), "pairs": npairs, "max_denominator": maxden}
        print(f"  k = {k}: {len(pts)} points of Q(k) (max denominator {maxden}), chain on {npairs} (x, u) pairs, "
              f"exact doubly stochastic matrix for every point", flush=True)
    res["chain"] = chain
    checks["chain"] = report(f"the chain of the proof and the statement on random points of Q(k), k = 2..{chain_kmax}",
                             ok, f"({t()} s)")
    res["checks"] = checks
    finish(res, out_dir(a.out), "small_atoms")


if __name__ == "__main__":
    main()
