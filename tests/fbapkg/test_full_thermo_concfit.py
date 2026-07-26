# -*- coding: utf-8 -*-
"""
Tests for the weighted L1/L2 TMFA concentration-fit objective builders added
to FullThermoPkg (agent-io/prds/tmfa-weighted-concfit/fullprompt.md, Modules
1-3):

    - build_concfit_objective(target_concentrations, weights=None)   [L1, MILP]
    - build_concfit_objective_l2(target_concentrations, weights=None) [L2, MIQP]
    - build_package(..., concfit_objective_type="l1"|"l2", concfit_weights=...)

Follows the fbapkg test pattern used elsewhere in this directory (build a
minimal hand-rolled cobra model, run build_package, assert on the resulting
objective/solution) rather than exercising internals.

Some tests require a QP-capable optlang backend (gurobi or cplex) and are
skipped if neither is importable/usable in this environment.
"""

import math

import cobra
import pytest

from modelseedpy.fbapkg.fullthermopkg import FullThermoPkg


# ---------------------------------------------------------------------------
# QP-backend availability probe (skip L2-positive tests if unavailable)
# ---------------------------------------------------------------------------

def _gurobi_available():
    try:
        model = cobra.Model("probe")
        model.solver = "gurobi"
        return True
    except Exception:
        return False


GUROBI_AVAILABLE = _gurobi_available()
requires_gurobi = pytest.mark.skipif(
    not GUROBI_AVAILABLE, reason="gurobi optlang backend not available in this environment"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_two_metabolite_model():
    """A minimal two-metabolite, zero-reaction cobra model.

    No deltaG/notes are attached — build_gibbs_constraint gracefully no-ops
    for metabolites lacking a deltaG value (see fullthermopkg.py), so this is
    sufficient to exercise the logconc/pconcfit/nconcfit machinery that the
    concfit objective builders operate on, independent of the rest of the
    TMFA formulation.
    """
    model = cobra.Model("test_concfit")
    met_a = cobra.Metabolite("metA_c0", name="metA", compartment="c0", charge=0, formula="C")
    met_b = cobra.Metabolite("metB_c0", name="metB", compartment="c0", charge=0, formula="C")
    model.add_metabolites([met_a, met_b])
    return model


def _build_pkg(solver="glpk"):
    model = _make_two_metabolite_model()
    model.solver = solver
    pkg = FullThermoPkg(model)
    pkg.build_package(
        {
            # Wide concentration bounds (M) so the hand-checked ln(M) targets
            # used below are always within the logconc variable's bounds.
            "custom_concentrations": {
                "metA_c0": [1e-10, 10],
                "metB_c0": [1e-10, 10],
            },
        }
    )
    return pkg


def _add_coupling_constraint(pkg, total):
    """Manually tie logconc_A + logconc_B = total.

    This hand-added constraint (not part of TMFA itself) is the couplings
    device the tests use to force a genuine trade-off between the two
    fitted metabolites' targets, so weighted-fit behavior is hand-checkable.
    """
    log_a = pkg.variables["logconc"]["metA_c0"]
    log_b = pkg.variables["logconc"]["metB_c0"]
    coupling = pkg.model.problem.Constraint(
        log_a + log_b, lb=total, ub=total, name="test_coupling"
    )
    pkg.model.add_cons_vars(coupling)


# ---------------------------------------------------------------------------
# Module 1 — weighted L1
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    """AC1: weights=None must reproduce the pre-change unweighted objective —
    same terms, unit (1.0) coefficients."""

    def test_weights_none_gives_unit_coefficients(self):
        pkg = _build_pkg("glpk")
        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        pkg.build_concfit_objective(targets, weights=None)

        coefs = pkg.model.objective.get_linear_coefficients(
            [
                pkg.variables["pconcfit"]["metA_c0"],
                pkg.variables["nconcfit"]["metA_c0"],
                pkg.variables["pconcfit"]["metB_c0"],
                pkg.variables["nconcfit"]["metB_c0"],
            ]
        )
        for var, coef in coefs.items():
            assert coef == 1.0, f"expected unit coefficient, got {coef} for {var}"

    def test_weights_none_and_omitted_are_equivalent(self):
        """Calling with weights=None and calling with no weights arg at all
        must produce the same objective coefficients."""
        pkg1 = _build_pkg("glpk")
        pkg2 = _build_pkg("glpk")
        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        pkg1.build_concfit_objective(targets, weights=None)
        pkg2.build_concfit_objective(targets)

        for pkg in (pkg1, pkg2):
            coefs = pkg.model.objective.get_linear_coefficients(
                [pkg.variables["pconcfit"]["metA_c0"], pkg.variables["nconcfit"]["metA_c0"]]
            )
            assert all(c == 1.0 for c in coefs.values())

    def test_existing_build_package_callers_unchanged(self):
        """build_package callers that never pass concfit_objective_type /
        concfit_weights still get l1 + weights=None behavior (AC9)."""
        pkg = _build_pkg("glpk")
        assert pkg.parameters["concfit_objective_type"] == "l1"
        assert pkg.parameters["concfit_weights"] is None


class TestWeightedL1:
    """AC3, AC10: weight validation, weight-0 drop, and objective replacement."""

    def test_nonuniform_weights_reflected_in_coefficients(self):
        pkg = _build_pkg("glpk")
        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        pkg.build_concfit_objective(targets, weights={"metA_c0": 5.0, "metB_c0": 2.0})
        coefs = pkg.model.objective.get_linear_coefficients(
            [
                pkg.variables["pconcfit"]["metA_c0"],
                pkg.variables["nconcfit"]["metA_c0"],
                pkg.variables["pconcfit"]["metB_c0"],
                pkg.variables["nconcfit"]["metB_c0"],
            ]
        )
        assert coefs[pkg.variables["pconcfit"]["metA_c0"]] == 5.0
        assert coefs[pkg.variables["nconcfit"]["metA_c0"]] == 5.0
        assert coefs[pkg.variables["pconcfit"]["metB_c0"]] == 2.0
        assert coefs[pkg.variables["nconcfit"]["metB_c0"]] == 2.0

    def test_negative_weight_raises_value_error(self):
        pkg = _build_pkg("glpk")
        targets = {"metA_c0": -2.0}
        with pytest.raises(ValueError):
            pkg.build_concfit_objective(targets, weights={"metA_c0": -1.0})

    def test_nan_or_inf_weight_raises_value_error(self):
        pkg = _build_pkg("glpk")
        targets = {"metA_c0": -2.0}
        with pytest.raises(ValueError):
            pkg.build_concfit_objective(targets, weights={"metA_c0": float("nan")})
        with pytest.raises(ValueError):
            pkg.build_concfit_objective(targets, weights={"metA_c0": float("inf")})

    def test_zero_weight_drops_term_from_objective(self):
        pkg = _build_pkg("glpk")
        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        pkg.build_concfit_objective(targets, weights={"metA_c0": 0.0, "metB_c0": 1.0})
        coefs = pkg.model.objective.get_linear_coefficients(
            [
                pkg.variables["pconcfit"]["metA_c0"],
                pkg.variables["nconcfit"]["metA_c0"],
                pkg.variables["pconcfit"]["metB_c0"],
                pkg.variables["nconcfit"]["metB_c0"],
            ]
        )
        assert coefs[pkg.variables["pconcfit"]["metA_c0"]] == 0.0
        assert coefs[pkg.variables["nconcfit"]["metA_c0"]] == 0.0
        assert coefs[pkg.variables["pconcfit"]["metB_c0"]] == 1.0
        assert coefs[pkg.variables["nconcfit"]["metB_c0"]] == 1.0

    def test_msid_only_key_contributes_no_term(self):
        """AC2: a bare-MSID key (no compartment suffix) does not match any
        metabolite.id in this model (metA_c0/metB_c0), so it is silently
        dropped rather than fitted."""
        pkg = _build_pkg("glpk")
        pkg.build_concfit_objective({"metA": -2.0}, weights=None)
        assert pkg.variables["pconcfit"] == {}
        assert pkg.variables["nconcfit"] == {}
        assert pkg.model.objective.expression == 0

    def test_weighted_fit_prefers_high_weight_target_when_conflicting(self):
        """Hand-checkable conflict case (PRD Testing Decisions / AC3):
        couple logconc_A + logconc_B = K via a manual constraint so the two
        targets cannot both be hit exactly, then verify the higher-weight
        metabolite's target wins the trade-off.

        Analytically: minimize w_A|x-a| + w_B|y-b| s.t. x+y=K.
        Substituting y=K-x gives a piecewise-linear convex 1-D problem in x
        whose minimizer is x=a whenever w_A > w_B (the higher-weight point
        dominates a two-point weighted-median problem).
        """
        pkg = _build_pkg("glpk")
        _add_coupling_constraint(pkg, total=-6.5)

        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        weights = {"metA_c0": 10.0, "metB_c0": 1.0}
        pkg.build_concfit_objective(targets, weights)

        status = pkg.model.solver.optimize()
        assert status == "optimal"

        log_a = pkg.variables["logconc"]["metA_c0"].primal
        log_b = pkg.variables["logconc"]["metB_c0"].primal

        assert math.isclose(log_a, targets["metA_c0"], abs_tol=1e-6), (
            f"expected high-weight metabolite A to hit its target exactly, got logconc_A={log_a}"
        )
        assert math.isclose(log_b, -6.5 - targets["metA_c0"], abs_tol=1e-6)
        # A's residual is ~0 while B absorbs the whole conflict.
        assert abs(log_a - targets["metA_c0"]) < abs(log_b - targets["metB_c0"])


class TestObjectiveReplacement:
    """AC10: build_concfit_objective* replace the model objective entirely —
    they do not combine with a prior min-error (or other) objective."""

    def test_concfit_replaces_min_error_objective(self):
        pkg = _build_pkg("glpk")
        pkg.build_min_error_objective()
        min_error_expr_str = str(pkg.model.objective.expression)
        assert "pdgerr" in min_error_expr_str or "ndgerr" in min_error_expr_str or min_error_expr_str == "0"

        pkg.build_concfit_objective({"metA_c0": -2.0}, weights=None)
        new_expr_str = str(pkg.model.objective.expression)
        assert "pdgerr" not in new_expr_str
        assert "ndgerr" not in new_expr_str
        assert "pconcfit" in new_expr_str or "nconcfit" in new_expr_str


# ---------------------------------------------------------------------------
# Module 2 — weighted L2
# ---------------------------------------------------------------------------

class TestL2Guard:
    """AC5: the QP-capability guard must raise NotImplementedError (not a
    deep-stack solver ValueError) when the backend is not gurobi/cplex."""

    def test_l2_raises_notimplementederror_on_glpk(self):
        pkg = _build_pkg("glpk")
        with pytest.raises(NotImplementedError, match="QP-capable"):
            pkg.build_concfit_objective_l2({"metA_c0": -2.0}, weights=None)

    def test_build_package_l2_raises_notimplementederror_on_glpk(self):
        model = _make_two_metabolite_model()
        model.solver = "glpk"
        pkg = FullThermoPkg(model)
        with pytest.raises(NotImplementedError, match="QP-capable"):
            pkg.build_package(
                {
                    "custom_concentrations": {"metA_c0": [1e-10, 10], "metB_c0": [1e-10, 10]},
                    "set_concfit_objective": True,
                    "target_concentrations": {"metA_c0": -2.0},
                    "concfit_objective_type": "l2",
                }
            )

    def test_guard_checked_before_weight_validation(self):
        """A negative weight combined with a GLPK backend should still
        surface the (more actionable) NotImplementedError, per the 'early'
        guard specified in the PRD."""
        pkg = _build_pkg("glpk")
        with pytest.raises(NotImplementedError):
            pkg.build_concfit_objective_l2({"metA_c0": -2.0}, weights={"metA_c0": -1.0})


@requires_gurobi
class TestL2Quadratic:
    """AC4: on a QP-capable backend, build_concfit_objective_l2 must produce
    a genuinely quadratic objective (squared logconc terms), and its optimum
    on a hand-solvable coupled instance must equal the analytic weighted
    least-squares point."""

    def test_l2_objective_is_quadratic(self):
        pkg = _build_pkg("gurobi")
        pkg.build_concfit_objective_l2(
            {"metA_c0": -2.0, "metB_c0": -3.0}, weights={"metA_c0": 2.0}
        )
        expr_str = str(pkg.model.objective.expression)
        assert "**2" in expr_str, f"expected a squared term in the objective, got: {expr_str}"

    def test_l2_no_slack_variables_created(self):
        pkg = _build_pkg("gurobi")
        pkg.build_concfit_objective_l2({"metA_c0": -2.0, "metB_c0": -3.0}, weights=None)
        assert pkg.variables["pconcfit"] == {}
        assert pkg.variables["nconcfit"] == {}
        assert pkg.constraints["concfit"] == {}

    def test_l2_analytic_weighted_least_squares_optimum(self):
        pkg = _build_pkg("gurobi")
        total = -6.5
        _add_coupling_constraint(pkg, total)

        w_a, w_b = 10.0, 1.0
        target_a, target_b = -2.0, -3.0
        pkg.build_concfit_objective_l2(
            {"metA_c0": target_a, "metB_c0": target_b},
            weights={"metA_c0": w_a, "metB_c0": w_b},
        )
        status = pkg.model.solver.optimize()
        assert status == "optimal"

        # Analytic weighted-least-squares solution for
        # minimize w_a*(x-target_a)^2 + w_b*(y-target_b)^2 s.t. x + y = total
        # (substitute y = total - x and set the derivative to zero).
        expected_x = (w_a * target_a + w_b * (total - target_b)) / (w_a + w_b)
        expected_y = total - expected_x

        log_a = pkg.variables["logconc"]["metA_c0"].primal
        log_b = pkg.variables["logconc"]["metB_c0"].primal
        assert math.isclose(log_a, expected_x, rel_tol=1e-4)
        assert math.isclose(log_b, expected_y, rel_tol=1e-4)

    def test_l2_after_l1_removes_stale_slacks_and_constraints(self):
        """AC7: building L2 after L1 for the same fitted set must remove
        the L1 pconcfit/nconcfit variables and concfit constraints — no
        stale linear-fit machinery left behind."""
        pkg = _build_pkg("gurobi")
        targets = {"metA_c0": -2.0, "metB_c0": -3.0}
        pkg.build_concfit_objective(targets, weights=None)
        assert set(pkg.variables["pconcfit"].keys()) == {"metA_c0", "metB_c0"}
        assert set(pkg.constraints["concfit"].keys()) == {"metA_c0", "metB_c0"}

        pkg.build_concfit_objective_l2(targets, weights={"metA_c0": 2.0})

        assert pkg.variables["pconcfit"] == {}
        assert pkg.variables["nconcfit"] == {}
        assert pkg.constraints["concfit"] == {}
        # And the underlying optlang model no longer carries those variables/constraints.
        remaining_names = [v.name for v in pkg.model.solver.variables]
        assert "metA_c0_pconcfit" not in remaining_names
        assert "metA_c0_nconcfit" not in remaining_names
        assert "metB_c0_pconcfit" not in remaining_names
        assert "metB_c0_nconcfit" not in remaining_names

    def test_l2_replaces_objective_entirely(self):
        # Note: build_min_error_objective's literal `1 * var` terms hit a
        # pre-existing, unrelated optlang/gurobi quirk (bare-Symbol addends
        # collapse and fail Gurobi's quadratic-expression parser) -- out of
        # scope for this PRD. Use the (already-fixed-for-gurobi) L1 concfit
        # objective as the "prior objective" instead, which is enough to
        # verify L2 fully replaces rather than combines with it.
        pkg = _build_pkg("gurobi")
        pkg.build_concfit_objective({"metA_c0": -2.0, "metB_c0": -3.0}, weights=None)
        pkg.build_concfit_objective_l2({"metA_c0": -2.0}, weights=None)
        expr_str = str(pkg.model.objective.expression)
        assert "pconcfit" not in expr_str
        assert "nconcfit" not in expr_str
        assert "**2" in expr_str


# ---------------------------------------------------------------------------
# Module 3 — build_package dispatch
# ---------------------------------------------------------------------------

class TestBuildPackageDispatch:
    def test_defaults_are_l1_and_none(self):
        pkg = _build_pkg("glpk")
        assert pkg.parameters["concfit_objective_type"] == "l1"
        assert pkg.parameters["concfit_weights"] is None

    def test_missing_target_concentrations_raises_value_error(self):
        """AC8: ValueError, not KeyError."""
        model = _make_two_metabolite_model()
        model.solver = "glpk"
        pkg = FullThermoPkg(model)
        with pytest.raises(ValueError):
            pkg.build_package(
                {
                    "custom_concentrations": {"metA_c0": [1e-10, 10], "metB_c0": [1e-10, 10]},
                    "set_concfit_objective": True,
                }
            )

    def test_dispatch_l1_by_default(self):
        model = _make_two_metabolite_model()
        model.solver = "glpk"
        pkg = FullThermoPkg(model)
        pkg.build_package(
            {
                "custom_concentrations": {"metA_c0": [1e-10, 10], "metB_c0": [1e-10, 10]},
                "set_concfit_objective": True,
                "target_concentrations": {"metA_c0": -2.0, "metB_c0": -3.0},
            }
        )
        assert set(pkg.variables["pconcfit"].keys()) == {"metA_c0", "metB_c0"}

    @requires_gurobi
    def test_dispatch_l2_with_weights(self):
        model = _make_two_metabolite_model()
        model.solver = "gurobi"
        pkg = FullThermoPkg(model)
        pkg.build_package(
            {
                "custom_concentrations": {"metA_c0": [1e-10, 10], "metB_c0": [1e-10, 10]},
                "set_concfit_objective": True,
                "target_concentrations": {"metA_c0": -2.0, "metB_c0": -3.0},
                "concfit_objective_type": "l2",
                "concfit_weights": {"metA_c0": 5.0},
            }
        )
        assert pkg.variables["pconcfit"] == {}
        expr_str = str(pkg.model.objective.expression)
        assert "**2" in expr_str

    def test_invalid_concfit_objective_type_raises_value_error(self):
        model = _make_two_metabolite_model()
        model.solver = "glpk"
        pkg = FullThermoPkg(model)
        with pytest.raises(ValueError):
            pkg.build_package(
                {
                    "custom_concentrations": {"metA_c0": [1e-10, 10], "metB_c0": [1e-10, 10]},
                    "set_concfit_objective": True,
                    "target_concentrations": {"metA_c0": -2.0},
                    "concfit_objective_type": "bogus",
                }
            )
