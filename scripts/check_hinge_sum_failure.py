#!/usr/bin/env python3
"""The sum of the single-step certificates does not certify multistep staircases in general.

Claim: for the staircase with thresholds 1/8, 1/4, 1/2, 5/8, 3/4, the function phi = sum_q (s - t_q)^+/(1 - eta_q),
t_q = (8 eta_q - 2)/3 (the sum of the certificates of Proposition 5.2 for its single steps) has
(1/2) int phi + 2 int (N - phi)^+ = 262105/39312, which exceeds the budget 2 sum (1 - eta_q) = 11/2 by 45889/39312.
This is a failure of that construction of a certificate only; the staircase inequality itself is certified for this
staircase by the catalogue on the grid 1/8 (check_density_catalogues.py).

Checks, exact rational arithmetic: the jumps of N (1/4, 3/8, 5/8, 3/4, 1, 9/8, 5/4, 3/2, 13/8, 7/4), the value and
the excess; that the excess is the same when N is recomputed by direct minimization (tailsums.box_brute) on every
interval between breakpoints; and that the deposited certificate for this staircase in the grid-1/8 catalogue has a
nonnegative margin. Seconds.
"""
import argparse
import os
from fractions import Fraction as Fr

from common import DATA, Timer, finish, out_dir, read_jsonl, report

from dsxray import density_two as D2
from dsxray import tailsums as TS

ETA = [Fr(1, 8), Fr(1, 4), Fr(1, 2), Fr(5, 8), Fr(3, 4)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="runs/hinge_sum_failure")
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    out = out_dir(a.out); res = {}; checks = {}
    t = Timer()
    xi = D2.jumps(ETA)
    J, bud = D2.hinge_sum(ETA)
    ok = xi == [Fr(x, 8) for x in (2, 3, 5, 6, 8, 9, 10, 12, 13, 14)] and J == Fr(262105, 39312) and bud == Fr(11, 2)
    ok &= J - bud == Fr(45889, 39312)
    # the same integral with N from direct minimization
    hinges = [((8 * e - 2) / 3, 1 / (1 - e)) for e in ETA]
    phi = D2.hinge_phi(hinges)
    nodes = sorted({Fr(0), Fr(2), *xi, *(h for h, _ in hinges if 0 < h < 2)})
    J2 = Fr(0)
    for lo, hi in zip(nodes, nodes[1:]):
        n = TS.box_brute(ETA, (lo + hi) / 2)
        J2 += (hi - lo) * (phi(lo) + phi(hi)) / 4 + 2 * D2.pos_integral(lo, hi, n - phi(lo), n - phi(hi))
    ok &= J2 == J
    res["hinge_sum"] = {"thresholds": [str(e) for e in ETA], "jumps": [str(x) for x in xi], "value": str(J),
                        "budget": str(bud), "excess": str(J - bud)}
    checks["hinge_sum_excess"] = report("sum of single-step certificates exceeds the budget by 45889/39312", ok, f"({t()} s)")
    t = Timer(); found = None
    for rec in read_jsonl(os.path.join(DATA, "density", "density_two_certificates_K08.jsonl.gz")):
        if rec["J"] == [1, 2, 4, 5, 6]:
            found = D2.verify_catalogue_record(rec)
    res["catalogue_margin"] = str(found[1]) if found else None
    checks["staircase_certified"] = report("the staircase itself is certified by the grid-1/8 catalogue", bool(found and found[0]),
                                           f"margin {found[1] if found else None} ({t()} s)")
    finish({"args": vars(a), "results": res, "checks": checks}, out, "hinge_sum_failure")


if __name__ == "__main__":
    main()
