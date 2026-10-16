"""Proposition 5.2 and the Lemma 5.1 certificates (exact verification of LP candidates)."""
from fractions import Fraction as Fr

from dsxray import density_certificates as DC
from dsxray import tailsums as TS


def test_single_step_closed_form_and_facts():
    for j in (1, 7, 30, 59, 60, 61, 120, 200, 239):
        e = Fr(j, 240)
        assert DC.single_step_J(e) == DC.single_step_closed_form(e)
        assert all(DC.single_step_facts(e).values())


def test_xi_matches_tailsums():
    for J in DC.staircases(7, 1):
        assert [Fr(x, 7) for x in DC.xi_int(J, 7)] == TS.xi([Fr(j, 7) for j in J])


def test_certificates_K6():
    for J in DC.staircases(6):
        ok, margin = DC.certify_staircase(J, 6, 4)
        assert ok and margin >= 0


def test_certificate_value_rejects_bad_phi():
    J, K, m = (2, 4), 6, 2
    phi = [Fr(0)] * (2 * K * m + 1)          # phi = 0: J(phi) = 2 int N, too large
    tot, bud = DC.certificate_value(phi, J, K, m)
    assert tot > bud
