#!/usr/bin/env python3
"""Remark S4.1: the atom bound 1/k of Corollary 3.2 cannot be raised to c/k with c = 101/100 (exact witnesses).

For each k in {9, 10, 11} a count vector x on {0, ..., 2k-2} (the law with masses x(t)/k at (t+1)/k) and a facet ray u
of hat-P(k) are given (WITNESSES below). Checked in exact rational arithmetic:
  - 0 <= x(t) <= 101/100 and x lies on the affine hull: sum x = k, sum t x(t) = k(k-1);
  - x satisfies all 2(k-2) hinge inequalities I(h_j), I(-h_j) of Proposition 3.3;
  - directly, without the hinges: T <=cx 2U_k (E(T-m)^+ <= E(2U_k-m)^+ for every integer m) and the folded condition
    U_k <=cx (law(T_flat) + law(kappa_k))/2 with kappa_k = (k-1) 1{T >= k}, T_flat = T - kappa_k;
  - u defines a facet of hat-P(k) (certified facet test) and x violates I(u) by 6/25, 9/25, 3/10 for k = 9, 10, 11
    (in count coordinates: sum_t v(t) x(t) - 2 sum_i u(i)).
Search (FLOAT; skipped with --no-search): the witnesses were found by a floating-point linear program over all facet
rays of hat-P(k) (data/facets), maximizing the violation subject to the bounds, the affine hull and the hinge
inequalities; every vertex found is recovered exactly from its active constraints and rechecked as above, and the
largest exactly verified violation equals the one listed. Default and --full are the same (about 10 s).
"""
import argparse
from fractions import Fraction as Fr

from common import Timer, finish, load_facet_rows, out_dir, report

from dsxray.facet_test import classify
from dsxray.tail_bounds import lattice_witness_check, lattice_witness_search, v_of

C = Fr(101, 100)
WITNESSES = {
    9: ([0, -3, -6, -1, -4, 1, -2, 3, 0],
        ["99/200", "101/100", "1/40", "47/50", "0", "101/100", "0", "101/100", "0", "101/100", "0", "101/100", "0",
         "101/100", "0", "101/100", "47/100"], Fr(6, 25)),
    10: ([0, -4, 1, -3, 2, -2, 3, -1, 4, 0],
         ["1", "0", "0", "101/100", "0", "101/100", "0", "101/100", "92/125", "101/100", "0", "101/100", "0",
          "101/100", "0", "101/100", "0", "101/100", "23/125"], Fr(9, 25)),
    11: ([0, -3, -6, -9, -12, -5, -8, -1, -4, 3, 0],
         ["0", "687/700", "101/100", "101/100", "0", "0", "101/100", "0", "0", "101/100", "0", "101/100", "313/700",
          "101/100", "0", "101/100", "337/700", "101/100", "0", "101/100", "0"], Fr(3, 10)),
}


def convex_order_direct(k, x):
    """T <=cx 2U_k for the law x(t)/k on {0..2k-2}: compare stop-loss at every integer m."""
    return all(sum(a * max(0, t - m) for t, a in enumerate(x)) / k <= Fr(sum(max(0, 2 * i - m) for i in range(k)), k)
               for m in range(-1, 2 * k + 1))


def folded_order_direct(k, x):
    """U_k <=cx (law(T_flat) + law(kappa_k)) / 2 (equal means; stop-loss at every integer m)."""
    M = {}
    for t, a in enumerate(x):
        p = a / k; kap = (k - 1) if t >= k else 0
        M[t - kap] = M.get(t - kap, 0) + p / 2; M[kap] = M.get(kap, 0) + p / 2
    mean_ok = sum(y * p for y, p in M.items()) == Fr(k - 1, 2)
    return mean_ok and all(Fr(sum(max(0, i - m) for i in range(k)), k) <= sum(p * max(0, y - m) for y, p in M.items())
                           for m in range(-1, k + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/lattice_tightness")
    ap.add_argument("--no-search", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}
    for k, (u, xs, target) in WITNESSES.items():
        t = Timer(); x = [Fr(s) for s in xs]
        hull, atoms, hinges, viol = lattice_witness_check(k, x, u, C)
        cx = convex_order_direct(k, x); fold = folded_order_direct(k, x); f = classify(u)["facet"]
        res[f"k{k}"] = {"affine_hull": hull, "atoms_at_most_101/100": atoms, "all_hinges": hinges,
                        "convex_order_direct": cx, "folded_order_direct": fold, "u_is_facet": f,
                        "violation": str(viol), "u": u, "x": xs}
        checks[f"witness_k{k}"] = report(f"k = {k}: atoms <= 101/100, hinges, convex and folded order, facet I(u) "
                                         f"violated by {target}", hull and atoms and hinges and cx and fold and f
                                         and viol == target and len(x) == 2 * k - 1, f"({t()} s)")
    if not a.no_search:
        found = {}
        for k in (9, 10, 11):
            t = Timer(); rays = [r["u"] for r in load_facet_rows(k) if r["facet"]]
            best = lattice_witness_search(k, C, rays)
            found[k] = None if best is None else {"violation": str(best[0]), "u": best[1], "x": [str(v) for v in best[2]]}
            checks[f"search_k{k}"] = report(f"FLOAT search, k = {k}: an exactly verified witness exists", best is not None
                                            and best[0] == WITNESSES[k][2], f"largest violation {best[0] if best else None} ({t()} s)")
        res["search"] = found
    res["checks"] = checks
    finish(res, out_dir(a.out), "lattice_tightness")


if __name__ == "__main__":
    main()
