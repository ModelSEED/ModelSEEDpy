# -*- coding: utf-8 -*-
import pytest
from cobra.core import Metabolite, Model, Reaction

from modelseedpy.fbapkg.objectivepkg import (
    ObjectiveData,
    ObjectivePkg,
    ObjectiveTerm,
)


@pytest.fixture
def model():
    """Two-reaction model: enough to exercise forward, reverse and unsigned terms."""
    m = Model("objtest")
    a = Metabolite("cpd00001_c0", compartment="c0")
    rxn1 = Reaction("rxn00001_c0", lower_bound=-1000, upper_bound=1000)
    rxn1.add_metabolites({a: -1})
    rxn2 = Reaction("rxn00002_c0", lower_bound=-1000, upper_bound=1000)
    rxn2.add_metabolites({a: 1})
    m.add_reactions([rxn1, rxn2])
    return m


@pytest.mark.parametrize(
    "term_string,coefficient,direction,variable",
    [
        ("rxn00001_c0", 1, "", "rxn00001_c0"),
        ("+rxn00001_c0", 1, "+", "rxn00001_c0"),
        ("-rxn00001_c0", 1, "-", "rxn00001_c0"),
        ("(1)+rxn00001_c0", 1, "+", "rxn00001_c0"),
        # Multi-digit and decimal coefficients: these were truncated to a single
        # character, so (12) silently parsed as 1 and (2.5) as 2.
        ("(12)+rxn00001_c0", 12, "+", "rxn00001_c0"),
        ("(2.5)-rxn00002_c0", 2.5, "-", "rxn00002_c0"),
        ("(100)rxn00002_c0", 100, "", "rxn00002_c0"),
    ],
)
def test_objective_term_from_string(term_string, coefficient, direction, variable):
    """ObjectiveTerm.from_string must preserve the whole coefficient, not its first digit."""
    term = ObjectiveTerm.from_string(term_string)
    assert term.coefficient == coefficient
    assert term.direction == direction
    assert term.variable == variable


@pytest.mark.parametrize(
    "objective_string",
    [
        "MAX{(1)+rxn00001_c0}",
        "MIN{(1)-rxn00002_c0}",
        "MAX{(12)+rxn00001_c0|(2.5)-rxn00002_c0}",
    ],
)
def test_objective_data_string_round_trip(objective_string):
    """from_string -> to_string must be lossless, including the sense and every coefficient."""
    assert ObjectiveData.from_string(objective_string).to_string() == objective_string


def test_objective_data_sense():
    assert ObjectiveData.from_string("MAX{(1)+rxn00001_c0}").sense == "max"
    assert ObjectiveData.from_string("MIN{(1)+rxn00001_c0}").sense == "min"


def test_to_cobrapy_objective_directions(model):
    """A signed term weights one variable; an unsigned term weights both, reverse negated."""
    data = ObjectiveData.from_string("MAX{(2)+rxn00001_c0|(3)rxn00002_c0}")
    objective = data.to_cobrapy_objective(model)

    assert objective.direction == "max"
    coefficients = objective.get_linear_coefficients(
        [
            model.reactions.rxn00001_c0.forward_variable,
            model.reactions.rxn00001_c0.reverse_variable,
            model.reactions.rxn00002_c0.forward_variable,
            model.reactions.rxn00002_c0.reverse_variable,
        ]
    )
    assert coefficients[model.reactions.rxn00001_c0.forward_variable] == 2
    assert coefficients[model.reactions.rxn00001_c0.reverse_variable] == 0
    assert coefficients[model.reactions.rxn00002_c0.forward_variable] == 3
    assert coefficients[model.reactions.rxn00002_c0.reverse_variable] == -3


def test_to_cobrapy_objective_warns_on_missing_reaction(model, caplog):
    """A term naming a reaction outside the model is reported, not silently dropped."""
    ObjectiveData.from_string("MAX{(1)+rxn99999_c0}").to_cobrapy_objective(model)
    assert "rxn99999_c0" in caplog.text


def test_build_package_returns_derived_name(model):
    """With no explicit name the derived name is returned, so the caller can restore by it."""
    pkg = ObjectivePkg(model)
    name = pkg.build_package("MAX{(1)+rxn00001_c0}")

    assert name == "MAX{(1)+rxn00001_c0}"
    assert name in pkg.objective_data_cache
    assert pkg.model.objective.direction == "max"


def test_build_package_explicit_name_and_restore(model):
    """restore_objective must rebuild a cached objective and keep the cached original intact."""
    pkg = ObjectivePkg(model)
    original = model.objective
    pkg.build_package("MAX{(1)+rxn00001_c0}", objective_name="growth")
    pkg.build_package("MIN{(1)+rxn00002_c0}", objective_name="maintenance")

    assert pkg.original_model_objective is not None
    pkg.restore_objective("growth")

    assert pkg.model.objective.direction == "max"
    assert set(pkg.objective_data_cache) == {"growth", "maintenance"}
    assert original is not None


def test_build_package_without_setting_objective(model):
    """set_objective=False caches the objective without applying it to the model."""
    pkg = ObjectivePkg(model)
    before = model.objective.direction
    pkg.build_package("MIN{(1)+rxn00001_c0}", objective_name="dry_run", set_objective=False)

    assert "dry_run" in pkg.objective_data_cache
    assert model.objective.direction == before


def test_build_package_rejects_bad_input(model):
    with pytest.raises(TypeError):
        ObjectivePkg(model).build_package(["not", "an", "objective"])
