# -*- coding: utf-8 -*-
"""
Regression tests for the FullThermoPkg compartment-potential unit bug fix.

Background:
    build_gibbs_constraint assembles a Gibbs constraint whose membrane term is:
        charge * Faraday(kC/mol) * compartment_potential(mV) / 1000

    Before the fix, the /1000 divisor was missing.  The shipped default of
    c0 = -160 mV was therefore treated as -160 V, inflating the membrane
    contribution by ~1000x and driving TMFA infeasible.

These tests verify:
    1. The constant embedded in the potc constraint reflects mV-scaled arithmetic
       (charge * 96.485 kC/mol * 0.16 V = charge * 15.44 kJ/mol for c0 = -160 mV).
    2. The default compartment_potential dict still uses -160 (mV) for c0.
    3. The warns-only guard fires for > 300 mV and > 1000 mV inputs.
"""

import logging
import math
import pytest

import cobra

from modelseedpy.fbapkg.fullthermopkg import FullThermoPkg


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_minimal_model():
    """Return a tiny two-metabolite / one-reaction COBRA model.

    Metabolites carry ModelSEED-compatible IDs and store deltaG in notes so
    build_gibbs_constraint can build potc constraints without hitting the DB.

    The charged metabolite (cpd_charged_c0, charge=-2) lives in compartment
    c0, which receives the shipped default potential of -160 mV.  The neutral
    metabolite (cpd_neutral_c0, charge=0) acts as a substrate so the model
    has at least one reaction and a non-trivial LP.
    """
    model = cobra.Model("test_minimal_tmfa")

    # Charged metabolite in c0 — the membrane term should be non-zero.
    met_charged = cobra.Metabolite(
        "cpd99001_c0",
        name="charged_cpd",
        compartment="c0",
        charge=-2,
        formula="C",
    )
    met_charged.notes = {"deltag": -200.0}  # kcal/mol (pre-conversion), arbitrary

    # Neutral metabolite in c0 — membrane term should be zero regardless of fix.
    met_neutral = cobra.Metabolite(
        "cpd99002_c0",
        name="neutral_cpd",
        compartment="c0",
        charge=0,
        formula="C",
    )
    met_neutral.notes = {"deltag": -100.0}

    # Neutral metabolite in e0 — zero-reference compartment.
    met_ext = cobra.Metabolite(
        "cpd99003_e0",
        name="ext_cpd",
        compartment="e0",
        charge=0,
        formula="C",
    )
    met_ext.notes = {"deltag": -50.0}

    # Single conversion reaction: neutral_e -> charged_c + neutral_c
    rxn = cobra.Reaction("rxn_test_c0")
    rxn.name = "test reaction"
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    rxn.add_metabolites({
        met_ext: -1,
        met_charged: 1,
        met_neutral: 1,
    })

    # Exchange reaction so the LP is not degenerate
    ex_rxn = cobra.Reaction("EX_cpd99003_e0")
    ex_rxn.lower_bound = -1000
    ex_rxn.upper_bound = 1000
    ex_rxn.add_metabolites({met_ext: -1})

    # Simple demand on the charged metabolite as objective proxy
    obj_rxn = cobra.Reaction("DM_cpd99001_c0")
    obj_rxn.lower_bound = 0
    obj_rxn.upper_bound = 1000
    obj_rxn.add_metabolites({met_charged: -1})

    model.add_reactions([rxn, ex_rxn, obj_rxn])
    model.objective = "DM_cpd99001_c0"

    return model


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestDefaultCompartmentPotential:
    """Verify that the shipped default keeps -160 mV for c0."""

    def test_default_c0_is_negative_160_mv(self):
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        pkg.build_package({})
        cp = pkg.parameters["compartment_potential"]
        assert "c0" in cp, "c0 must be present in compartment_potential defaults"
        assert cp["c0"] == -160, (
            f"Default c0 compartment_potential must be -160 mV, got {cp['c0']}"
        )

    def test_default_e0_is_zero(self):
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        pkg.build_package({})
        cp = pkg.parameters["compartment_potential"]
        assert cp.get("e0", None) == 0, "e0 (zero reference) must be 0 mV"


class TestMembraneTermScaling:
    """The potc constraint constant must reflect mV-scaled membrane arithmetic."""

    Faraday = 96.485  # kC/mol
    COMP_POT_MV = -160  # mV default for c0

    def _get_potc_constant(self, model, met_id):
        """Extract lb (== ub == constant) from the potc constraint for met_id."""
        pkg = FullThermoPkg(model)
        pkg.build_package({})
        potc = pkg.constraints["potc"]
        if met_id not in potc:
            return None
        con = potc[met_id]
        assert con.lb == con.ub, "potc constraint must be an equality (lb == ub)"
        return con.lb

    def test_charged_met_constant_is_mv_scaled(self):
        """charge * Faraday * (mV/1000) must match the embedded constant."""
        model = _make_minimal_model()
        # cpd99001_c0 has charge=-2, deltag=-200 kcal/mol, compartment c0 (-160 mV)
        constant = self._get_potc_constant(model, "cpd99001_c0")
        assert constant is not None, "potc constraint for cpd99001_c0 not built"

        charge = -2
        deltag_kcal = -200.0
        deltag_kj = deltag_kcal * 4.184  # -836.8 kJ/mol
        membrane_kj = charge * self.Faraday * (self.COMP_POT_MV / 1000)
        # = -2 * 96.485 * (-0.16) = +30.875 kJ/mol
        expected = deltag_kj + membrane_kj

        assert math.isclose(constant, expected, rel_tol=1e-6), (
            f"potc constant for charged c0 metabolite: expected {expected:.4f} kJ/mol "
            f"(mV-scaled), got {constant:.4f} kJ/mol. "
            f"A ~1000x discrepancy indicates the /1000 divisor is missing."
        )

    def test_membrane_contribution_magnitude_is_not_1000x_inflated(self):
        """The membrane term magnitude must be on the order of kJ/mol, not MJ/mol."""
        model = _make_minimal_model()
        constant = self._get_potc_constant(model, "cpd99001_c0")
        assert constant is not None

        charge = -2
        expected_membrane_kj = charge * self.Faraday * (self.COMP_POT_MV / 1000)
        # +30.875 kJ/mol

        # Pre-fix (missing /1000) would give: -2 * 96.485 * (-160) = +30875.2 kJ/mol
        # We assert the constant is NOT inflated by ~1000x.
        inflated_membrane = charge * self.Faraday * self.COMP_POT_MV
        deltag_kj = -200.0 * 4.184
        inflated_constant = deltag_kj + inflated_membrane

        assert not math.isclose(constant, inflated_constant, rel_tol=0.01), (
            "potc constant matches the pre-fix inflated value — /1000 divisor is missing"
        )

    def test_neutral_met_has_zero_membrane_term(self):
        """A charge=0 metabolite must have no membrane contribution."""
        model = _make_minimal_model()
        constant = self._get_potc_constant(model, "cpd99002_c0")
        assert constant is not None, "potc constraint for cpd99002_c0 not built"

        deltag_kj = -100.0 * 4.184
        # membrane term = 0 * anything = 0
        expected = deltag_kj
        assert math.isclose(constant, expected, rel_tol=1e-6), (
            f"Neutral metabolite potc constant: expected {expected:.4f}, got {constant:.4f}"
        )

    def test_zero_reference_compartment_has_zero_membrane_term(self):
        """Metabolites in e0 (potential=0) must have no membrane contribution."""
        model = _make_minimal_model()
        # cpd99003_e0 is in e0 which has compartment_potential=0
        constant = self._get_potc_constant(model, "cpd99003_e0")
        assert constant is not None, "potc constraint for cpd99003_e0 not built"

        deltag_kj = -50.0 * 4.184
        expected = deltag_kj + 0  # e0 potential = 0
        assert math.isclose(constant, expected, rel_tol=1e-6), (
            f"e0 metabolite potc constant: expected {expected:.4f}, got {constant:.4f}"
        )


class TestWarnsOnlyGuard:
    """The guard must warn (not raise) for large compartment_potential values."""

    def test_warns_above_300_mv(self, caplog):
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        with caplog.at_level(logging.WARNING):
            pkg.build_package({"compartment_potential": {"c0": 400, "e0": 0}})
        # At least one warning about c0 being unusually large
        msgs = [r.message for r in caplog.records if "c0" in r.message and "400" in r.message]
        assert msgs, (
            "Expected a WARNING for compartment_potential=400 mV (> 300 mV) but none was emitted"
        )

    def test_warns_above_1000_mv(self, caplog):
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        with caplog.at_level(logging.WARNING):
            pkg.build_package({"compartment_potential": {"c0": 1500, "e0": 0}})
        msgs = [r.message for r in caplog.records if "c0" in r.message and "1500" in r.message]
        assert msgs, (
            "Expected a WARNING for compartment_potential=1500 mV (> 1000 mV) but none was emitted"
        )

    def test_no_warning_for_default_mv(self, caplog):
        """The default -160 mV must not trigger any warning."""
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        with caplog.at_level(logging.WARNING):
            pkg.build_package({})  # uses defaults: c0=-160, e0=0
        suspect = [
            r for r in caplog.records
            if r.levelno >= logging.WARNING and "compartment_potential" in r.message
        ]
        assert not suspect, (
            f"Unexpected warning(s) for default compartment_potential: {[r.message for r in suspect]}"
        )

    def test_does_not_raise_for_large_value(self):
        """Even egregiously large values must not raise an exception."""
        model = _make_minimal_model()
        pkg = FullThermoPkg(model)
        # Should not raise — only warn
        pkg.build_package({"compartment_potential": {"c0": 5000, "e0": 0}})
