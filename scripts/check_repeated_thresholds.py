#!/usr/bin/env python3
"""Example 4.2: on a fixed grid, repeated thresholds are needed (exact arithmetic throughout).

F_W has masses (1/12, 1/3, 1/6, 1/6, 1/4) on {1/3, 1/2, 1, 4/3, 5/3}; mass and mean are 1.
  (a) For the thresholds (1/6, 1/3, 2/3, 2/3): (xi_1..xi_8) = (1/3, 1/2, 5/6, 1, 4/3, 4/3, 5/3, 5/3) (closed formula
      of Lemma 4.1 and brute-force (u box u)), left side 53/12, right side 52/12. On the grid k = 6 this inequality is
      I(u) for u = (0, 1, 2, 2, 4, 4), divided by k; u is a ray of height two, its inequality is a facet of hat-P(6)
      (certified facet test), and by brute force over all 720 permutations it has 25 tight permutations whose X-rays
      span an affine space of dimension 8 = 2k - 4.
  (b) F_W satisfies I(u_J) for every step set J of {1, ..., K-1}, K in {6, 12, 18, 24}: all 2^(K-1) step sets
      (8,388,608 for K = 24), exact integer arithmetic. The C helper tools/c/grid_steps.c is used when a C compiler is
      available, and numpy (exact int64, chunked) for K <= --numpy-kmax (default 24); when both run they must agree.
  (c) J = {6, 12, 23, 24} on the grid K = 36 (thresholds (1/6, 1/3, 23/36, 2/3)) excludes F_W with violation 1/36,
      computed as a tail sum and as I(u_J)/K.
Default and --full are the same (about 10 seconds).
"""
import argparse
import os
import shutil
import subprocess
import json
import itertools
import tempfile
from fractions import Fraction as Fr

import numpy as np

from common import ROOT, Timer, finish, out_dir, report

from dsxray import tailsums as TS
from dsxray.facet_test import classify
from dsxray.padding import pattern

LAW = [(Fr(1, 3), Fr(1, 12)), (Fr(1, 2), Fr(1, 3)), (Fr(1), Fr(1, 6)), (Fr(4, 3), Fr(1, 6)), (Fr(5, 3), Fr(1, 4))]
ETA_A = [Fr(1, 6), Fr(1, 3), Fr(2, 3), Fr(2, 3)]
XI_A = [Fr(1, 3), Fr(1, 2), Fr(5, 6), Fr(1), Fr(4, 3), Fr(4, 3), Fr(5, 3), Fr(5, 3)]
U_A = [0, 1, 2, 2, 4, 4]
GRIDS = (6, 12, 18, 24)
J_C, K_C = (6, 12, 23, 24), 36


def grid_counts(K):
    """Integer form of F_W on the grid 1/K: (t, w) with mass w/12 at (t+1)/K."""
    return [(int(K * a) - 1, int(12 * m)) for a, m in LAW], 12


def numpy_all_step_sets(K):
    tw, D = grid_counts(K); n = 1 << (K - 1); chunk = 1 << 16; viol = 0; best = None
    for start in range(0, n, chunk):
        J = np.arange(start, min(n, start + chunk), dtype=np.int64)
        bits = ((J[:, None] >> np.arange(K - 1)[None, :]) & 1).astype(np.int64)
        U = np.concatenate([np.zeros((len(J), 1), dtype=np.int64), np.cumsum(bits, axis=1)], axis=1)
        lhs = np.zeros(len(J), dtype=np.int64)
        for t, w in tw:
            lo, hi = max(0, t - K + 1), min(K - 1, t)
            mn = U[:, lo] + U[:, t - lo]
            for i in range(lo + 1, hi + 1):
                mn = np.minimum(mn, U[:, i] + U[:, t - i])
            lhs += w * mn
        diff = lhs * K - 2 * D * U.sum(axis=1)
        viol += int((diff > 0).sum()); m = int(diff.max()); best = m if best is None else max(best, m)
    return {"K": K, "step_sets": n, "violated": viol, "max_lhs_minus_rhs": best}


def c_all_step_sets(K, exe):
    tw, D = grid_counts(K)
    args = [exe, str(K), str(D)] + [str(x) for p in tw for x in p]
    return json.loads(subprocess.run(args, check=True, capture_output=True, text=True).stdout)


def build_c(workdir):
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if not cc:
        return None
    exe = os.path.join(workdir, "grid_steps")
    try:   # a compiler stub that cannot compile (macOS without the command-line tools) counts as no compiler
        subprocess.run([cc, "-O2", "-o", exe, os.path.join(ROOT, "tools", "c", "grid_steps.c")], check=True,
                       capture_output=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    return exe if os.path.exists(exe) else None


def tight_permutation_xrays(u):
    k = len(u); v = TS.grid_v(u); out = []
    for s in itertools.permutations(range(k)):
        if all(u[i] + u[s[i]] == v[i + s[i]] for i in range(k)):
            x = [0] * (2 * k - 1)
            for i in range(k):
                x[i + s[i]] += 1
            out.append(x)
    return out


def rank(rows):
    M = [[Fr(x) for x in r] for r in rows]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/repeated_thresholds")
    ap.add_argument("--numpy-kmax", type=int, default=24)
    ap.add_argument("--no-c", action="store_true")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    res = {}; checks = {}
    mass = sum(m for _, m in LAW); mean = sum(x * m for x, m in LAW)
    checks["law"] = report("F_W has mass 1 and mean 1", mass == 1 and mean == 1)

    t = Timer()
    X = TS.xi(ETA_A)
    lhs = sum(m * sum(1 for x in X if x <= s) for s, m in LAW); rhs = 2 * sum(1 - e for e in ETA_A)
    brute = TS.violation_brute(ETA_A, LAW)
    k = 6; x = [Fr(0)] * (2 * k - 1)
    for s, m in LAW:
        x[int(k * s) - 1] += k * m
    v = TS.grid_v(U_A); iu = (sum(p * q for p, q in zip(v, x)) - 2 * sum(U_A)) / k
    c = classify(U_A); xr = tight_permutation_xrays(U_A)
    dim = rank([[p - q for p, q in zip(r, xr[0])] for r in xr[1:]])
    res["part_a"] = {"xi": [str(z) for z in X], "lhs": str(lhs), "rhs": str(rhs), "I_u_over_k": str(iu),
                     "height": max(pattern(U_A)), "ray": c["ray"], "facet": c["facet"], "certificate": c["cert"],
                     "tight_permutations": len(xr), "affine_dimension": dim}
    ok = X == XI_A and lhs == Fr(53, 12) and rhs == Fr(52, 12) and brute == lhs - rhs == iu
    ok &= max(pattern(U_A)) == 2 and c["ray"] and c["facet"] and len(xr) == 25 and dim == 8 == 2 * k - 4
    checks["part_a"] = report("(a) xi, 53/12 > 52/12, I(0,1,2,2,4,4) on k = 6: ray of height 2, facet, 25 tight "
                              "permutations spanning dimension 8", ok, f"({t()} s)")

    t = Timer(); ok = True; rows = {}
    with tempfile.TemporaryDirectory() as tmp:
        exe = None if a.no_c else build_c(tmp)
        for K in GRIDS:
            r_np = numpy_all_step_sets(K) if (K <= a.numpy_kmax or exe is None) else None
            r_c = c_all_step_sets(K, exe) if exe else None
            if r_np and r_c:
                ok &= all(r_np[f] == r_c[f] for f in ("step_sets", "violated", "max_lhs_minus_rhs"))
            r = r_c or r_np
            rows[K] = {"step_sets": r["step_sets"], "violated": r["violated"], "max_lhs_minus_rhs": r["max_lhs_minus_rhs"],
                       "engines": [e for e, rr in (("numpy", r_np), ("C", r_c)) if rr]}
            ok &= r["step_sets"] == 2 ** (K - 1) and r["violated"] == 0
    res["part_b"] = rows
    checks["part_b"] = report("(b) no step set J on the grids K = 6, 12, 18, 24 excludes F_W (8,388,608 at K = 24)",
                              ok and rows[24]["step_sets"] == 8388608, f"({t()} s)")

    eta_c = [Fr(j, K_C) for j in J_C]
    viol = TS.violation(eta_c, LAW); viol_b = TS.violation_brute(eta_c, LAW)
    uJ = [sum(i >= j for j in J_C) for i in range(K_C)]; vJ = TS.grid_v(uJ)
    xc = [Fr(0)] * (2 * K_C - 1)
    for s, m in LAW:
        xc[int(K_C * s) - 1] += K_C * m
    iuJ = (sum(p * q for p, q in zip(vJ, xc)) - 2 * sum(uJ)) / K_C
    res["part_c"] = {"thresholds": [str(e) for e in eta_c], "violation": str(viol), "I_uJ_over_K": str(iuJ)}
    checks["part_c"] = report("(c) J = {6, 12, 23, 24} on K = 36 excludes F_W by 1/36",
                              viol == viol_b == iuJ == Fr(1, 36) and eta_c == [Fr(1, 6), Fr(1, 3), Fr(23, 36), Fr(2, 3)])
    res["checks"] = checks
    finish(res, out_dir(a.out), "repeated_thresholds")


if __name__ == "__main__":
    main()
