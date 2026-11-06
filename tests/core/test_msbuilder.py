# -*- coding: utf-8 -*-
from modelseedpy.core.msmodel import get_direction_from_constraints
from modelseedpy.core.msbuilder import MSBuilder
from tests.test_data.mock_data import mock_template, mock_genome_rast, mock_model


def test_get_direction_from_constraints1():
    assert get_direction_from_constraints(-10, 10) == "="


def test_get_direction_from_constraints2():
    assert get_direction_from_constraints(0, 10) == ">"


def test_get_direction_from_constraints3():
    assert get_direction_from_constraints(5, 10) == ">"


def test_get_direction_from_constraints4():
    assert get_direction_from_constraints(-10, 0) == "<"


def test_get_direction_from_constraints5():
    assert get_direction_from_constraints(-10, -5) == "<"


def test_get_direction_from_constraints6():
    assert get_direction_from_constraints(0, 0) == "?"


def test_apply_gapfill_reversals():
    from cobra.core import Model, Metabolite, Reaction

    model = Model("test")
    a = Metabolite("a_c0", compartment="c0")
    b = Metabolite("b_c0", compartment="c0")
    for rxn_id, lb, ub in [("r1_c0", 0, 1000), ("r2_c0", -5, 0), ("r3_c0", 0, 2000)]:
        rxn = Reaction(rxn_id, lower_bound=lb, upper_bound=ub)
        rxn.add_metabolites({a: -1, b: 1})
        model.add_reactions([rxn])

    applied = MSBuilder.apply_gapfill_reversals(
        model, {"r1_c0": "<", "r2_c0": ">", "r3_c0": ">", "r4_c0": "<"}
    )

    assert applied == ["r1_c0", "r2_c0", "r3_c0"]
    assert model.reactions.r1_c0.bounds == (-1000, 1000)
    assert model.reactions.r2_c0.bounds == (-5, 1000)
    assert model.reactions.r3_c0.bounds == (0, 2000)


def test_some_gpr():
    # using this for later
    reaction_complex_role_gene_mapping = {
        "cpx1": {
            "role1": ["role1function", True, False, {"g1", "g2"}],
            "role2": ["role2function", True, False, {"g3"}],
        },
        "cpx2": {
            "role1": ["role1function", True, False, {"g1", "g2"}],
            "role3": ["role2function", True, False, {"g4", "g5"}],
        },
    }

    # reaction_complex_role_gene_mapping -> GPR_func -> expected

    expected = {
        frozenset({"g1", "g3"}),
        frozenset({"g2", "g3"}),
        frozenset({"g1", "g4"}),
        frozenset({"g1", "g5"}),
        frozenset({"g2", "g4"}),
        frozenset({"g2", "g5"}),
    }


def test_build():
    template = mock_template()
    genome = mock_genome_rast()
    # builder = MSBuilder(genome, template)
    # model = builder.build('test_model', '0', True, False)

    expect = mock_model()

    pass
