"""Density-two certificates: two-step, nonincreasing gaps, equally spaced, coupling, hinge sum, catalogue verifier."""
import os
from fractions import Fraction as Fr

import pytest

from dsxray import density_certificates as DC
from dsxray import density_two as D2
from dsxray import tailsums as TS

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_jumps_match_tailsums_and_brute_force():
    for eta in ([Fr(1, 5)], [Fr(1, 7), Fr(3, 7)], [Fr(1, 9), Fr(4, 9), Fr(5, 9), Fr(8, 9)]):
        xi = D2.jumps(eta)
        assert xi == TS.xi(eta)
        for j in range(1, 80, 3):
            s = Fr(j, 40)
            assert TS.box_brute(eta, s) == sum(1 for x in xi if x <= s)


def test_two_step_cases_and_margin():
    seen = set()
    for i in range(1, 24):
        for j in range(i + 1, 24):
            r = D2.two_step_check(Fr(i, 24), Fr(j, 24))
            assert r["ok"]
            seen.add(r["case"])
    assert seen == {"A", "B", "C"}
    assert D2.two_step(Fr(1, 4), Fr(1, 2) + Fr(1, 8))[2] > 0


def test_nonincreasing_gaps_closed_form_and_single_step():
    m, closed, _ = D2.nonincreasing_gaps([Fr(1, 5), Fr(1, 5), Fr(1, 10)], 2)
    assert m == closed > 0
    eta = Fr(3, 10)
    m, _, t = D2.nonincreasing_gaps([1 - eta], 2)
    assert t == (8 * eta - 2) / 3 and m == 2 * (1 - eta) - DC.single_step_J(eta)


def test_equally_spaced_certificate_and_sharpness():
    for r in (1, 2, 5):
        m, mc, J, Jc = D2.equally_spaced(r, Fr(1, 2 * r), r + 1)
        assert m == mc >= 0 and J == Jc
        eta, iv, viol, gap = D2.equally_spaced_sharpness(r, r + 2)
        assert viol == Fr(1, 2 * (r + 2) * (2 * r + 1)) and gap >= 0


def test_coupling():
    P = D2.coupling_profiles()
    assert all(h == 1 for _, _, h in P["X"] + P["Y"])
    assert D2.profile_moment(P["S"], 1) == 1 and D2.profile_min_convex_order_gap(P["S"]) == 0


def test_hinge_sum_excess():
    J, bud = D2.hinge_sum([Fr(1, 8), Fr(1, 4), Fr(1, 2), Fr(5, 8), Fr(3, 4)])
    assert J - bud == Fr(45889, 39312)


def test_catalogue_verifier_accepts_and_rejects():
    phi = [max(Fr(0), 2 * (Fr(i, 6) - Fr(2, 3))) for i in range(13)]
    ok, margin = D2.verify_catalogue_record({"J": [1], "K": 2, "m": 3, "phi": [str(x) for x in phi]})
    assert ok and margin == 0
    ok, margin = D2.verify_catalogue_record({"J": [1], "K": 2, "m": 1, "phi": ["0"] * 5})
    assert not ok and margin == -2
    with pytest.raises(ValueError):
        D2.verify_catalogue_record({"J": [1], "K": 2, "m": 1, "phi": ["0"] * 5, "margin": "0"})
    with pytest.raises(ValueError):
        D2.verify_catalogue_record({"J": [1, 1], "K": 2, "m": 1, "phi": ["0"] * 5})


def test_deposited_catalogue_record():
    import gzip
    import json
    with gzip.open(os.path.join(ROOT, "data", "density", "density_two_certificates_K08.jsonl.gz"), "rt") as f:
        rec = json.loads(f.readline())
    ok, margin = D2.verify_catalogue_record(rec)
    assert ok and margin > 0
    assert D2.min_convolution_cells(rec["J"], 8, 4) == DC.N_on_cells(tuple(rec["J"]), 8, 4)
    bad = dict(rec, phi=list(rec["phi"]))
    bad["phi"][5] = str(Fr(bad["phi"][5]) + 1)          # breaks convexity
    assert not D2.verify_catalogue_record({k: v for k, v in bad.items() if k != "margin"})[0]
