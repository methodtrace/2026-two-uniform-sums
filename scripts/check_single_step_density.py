#!/usr/bin/env python3
"""Proposition 5.2 (single steps at density two) and the exact certificates of Lemma 5.1 after Problem 5.1.

  1. Proposition 5.2 for eta = j/240, j = 1..239, in rational arithmetic: with t = (8 eta - 2)/3 and
     phi(s) = (s - t)^+/(1 - eta), the quantity (1/2) int phi + 2 int (N - phi)^+ of Lemma 5.1 equals the closed form
     2(1 - eta) - (t^-)^2/(4(1 - eta)); and the intermediate facts of the proof: t <= 2 eta; phi = 2/3, 1, 5/3, 2 at
     2 eta, (5 eta + 1)/3, 1 + eta, (2 eta + 4)/3, in this order in [2 eta, 2]; the two triangles give
     2 int (N - phi)^+ = 2(1 - eta)/9; (1/2) int phi = 16(1 - eta)/9 for t >= 0 and (1 - t)/(1 - eta) for t < 0;
     t < 0 iff eta < 1/4.
  2. The paragraph after Problem 5.1: for every staircase with at least two distinct thresholds in {1/K..(K-1)/K}
     (120, 2,036 and 32,752 staircases for K = 8, 12, 16) a convex piecewise linear phi with rational slopes satisfies
     the hypothesis of Lemma 5.1 with c = 2. The candidate phi comes from a floating-point linear program (numpy,
     scipy) on the grid 1/(Km) (m = 4, 4, 3); it is rounded to a convex rational function and the hypothesis is
     verified in exact rational arithmetic. The thresholds xi of Lemma 4.1 used for N = u box u are compared with
     the independent implementation in tailsums.py on every staircase with K <= 8 and with direct minimization on
     a sample.
  3. FLOAT, numerical evidence only (the floating-point searches after Problem 5.1): coordinate ascent over staircases
     with r <= --float-rmax thresholds on the grid 1/--float-grid (repetitions and thresholds equal to 1 allowed)
     against the largest violation over densities <= 2 on cells of width 1/--float-grid with S <=cx 2U (linear
     programs, scipy HiGHS; convex order imposed at cell endpoints and midpoints). Expected: the largest violation
     is 0 (up to 1e-9), and every start that reaches it ends at a multiple of a single step (one distinct threshold
     below 1).
Default: K = 8; float search r <= 4 on the grid 1/20 (about 5 s). --full: K = 8, 12, 16 (K = 16: about 5
core-minutes, --workers splits the list); float search r <= 8 on the grid 1/40 (about 30 s).
"""
import argparse
import random
from fractions import Fraction as Fr

from common import Timer, finish, out_dir, pmap, report

from dsxray import density_certificates as DC
from dsxray import tailsums as TS

GRIDS = {8: (4, 120), 12: (4, 2036), 16: (3, 32752)}


def _cert(args):
    J, K, m = args
    ok, margin = DC.certify_staircase(J, K, m)
    return ok, margin


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/single_step_density")
    ap.add_argument("--grids", default="8")
    ap.add_argument("--float-rmax", type=int, default=4)
    ap.add_argument("--float-grid", type=int, default=20)
    ap.add_argument("--float-seeds", type=int, default=2)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    if a.full:
        a.grids = "8,12,16"; a.float_rmax = max(a.float_rmax, 8); a.float_grid = max(a.float_grid, 40)
        a.float_seeds = max(a.float_seeds, 3)
    out = out_dir(a.out); res = {}; checks = {}

    t = Timer(); etas = [Fr(j, 240) for j in range(1, 240)]
    closed = all(DC.single_step_J(e) == DC.single_step_closed_form(e) for e in etas)
    facts = [e for e in etas if not all(DC.single_step_facts(e).values())]
    eq = sum(1 for e in etas if DC.single_step_J(e) == 2 * (1 - e))
    res["prop_singlestep"] = {"etas": len(etas), "closed_form": closed, "facts_failures": [str(e) for e in facts],
                              "equality_cases": eq}
    checks["prop_singlestep"] = report("Proposition 5.2: J(phi) equals the closed form, and the facts of the proof, eta = j/240", closed and not facts,
                                       f"equality 2(1 - eta) for {eq} values (eta >= 1/4) ({t()} s)")

    t = Timer(); rng = random.Random(3); ok = True
    for J in DC.staircases(8, 1):
        ok &= [Fr(x, 8) for x in DC.xi_int(J, 8)] == TS.xi([Fr(j, 8) for j in J])
    for _ in range(40):
        K = rng.choice([8, 12, 16]); r = rng.randint(1, 5); J = tuple(sorted(rng.sample(range(1, K), r)))
        eta = [Fr(j, K) for j in J]; m = 2; N = DC.N_on_cells(J, K, m)
        for i in range(0, 2 * K * m, 3):
            ok &= TS.box_brute(eta, Fr(2 * i + 1, 2 * K * m)) == N[i]
    checks["thresholds"] = report("xi_l of Lemma 4.1 equal to tailsums.xi (K = 8) and N to direct minimization (sample)", ok, f"({t()} s)")

    for K in [int(x) for x in a.grids.split(",")]:
        t = Timer(); m, expected = GRIDS[K]
        Js = list(DC.staircases(K))
        outs = pmap(_cert, [(J, K, m) for J in Js], a.workers, chunksize=16)
        bad = [J for J, (good, _) in zip(Js, outs) if not good]
        worst = min(range(len(Js)), key=lambda i: outs[i][1])
        res[f"certificates_K{K}"] = {"staircases": len(Js), "grid_m": m, "not_certified": bad[:20],
                                    "min_margin": str(outs[worst][1]), "min_margin_at": Js[worst]}
        checks[f"certificates_K{K}"] = report(f"Lemma 5.1 certificates, c = 2, K = {K}: every staircase with >= 2 distinct thresholds",
                                              not bad and len(Js) == expected,
                                              f"{len(Js)} staircases, min exact margin {outs[worst][1]} at J = {Js[worst]} ({t()} s)")
    if a.float_rmax:
        import numpy as np
        t = Timer(); rng = np.random.default_rng(5); runs = {}; ok = True; best = -9.0
        for r in range(1, a.float_rmax + 1):
            o = DC.float_search(2.0, a.float_grid, a.float_grid, r, a.float_seeds, rng)
            runs[r] = [(float(f"{v:.12g}"), et) for v, et in o]
            for v, et in o:
                best = max(best, v)
                if v > -1e-9 and len({x for x in et if x < a.float_grid}) != 1:
                    ok = False
        ok &= best <= 1e-9
        res["float_search"] = {"grid": a.float_grid, "runs": runs, "max_violation": best}
        checks["float_search"] = report(f"FLOAT: staircases with <= {a.float_rmax} thresholds, grid 1/{a.float_grid}, c = 2: max violation {best:.2e}, attained only by single steps",
                                        ok, f"({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "single_step_density")


if __name__ == "__main__":
    main()
