# -*- coding: utf-8 -*-
import pytest

from modelseedpy.core.mstemplate import (
    MSTemplate,
    MSTemplateMetabolite,
    MSTemplateReaction,
    MSTemplateSpecies,
)
from modelseedpy.core.mstemplate import (
    NewModelTemplateRole,
    NewModelTemplateComplex,
    MSTemplateCompartment,
)


@pytest.fixture
def empty_template():
    return MSTemplate("test", "test name", "test")


def test_empty_template():
    template = MSTemplate("test", "test name", "test")
    assert template.id == "test"
    assert template.name == "test name"
    assert len(template.roles) == 0
    assert len(template.complexes) == 0
    assert len(template.compounds) == 0
    assert len(template.compcompounds) == 0
    assert len(template.reactions) == 0


def test_template_add_role(empty_template):
    role = NewModelTemplateRole("role1", "metabolic function")
    empty_template.add_roles([role])
    assert len(empty_template.roles) == 1


def test_template_add_role_mult(empty_template):
    role_a = NewModelTemplateRole("roleA", "metabolic function A")
    role_b = NewModelTemplateRole("roleB", "metabolic function B")
    role_c = NewModelTemplateRole("roleC", "metabolic function C")
    empty_template.add_roles([role_a, role_b, role_c])
    assert len(empty_template.roles) == 3


def test_template_add_simple_complex(empty_template):
    role = NewModelTemplateRole("role1", "metabolic function")
    empty_template.add_roles([role])

    seed_complex = NewModelTemplateComplex("complex1", "example complex")

    seed_complex.add_role(empty_template.roles.role1)

    empty_template.add_complexes([seed_complex])

    assert len(empty_template.complexes) == 1


def test_template_add_simple_metabolite(empty_template):
    cpd_apple = MSTemplateMetabolite("apple", "C100", "just a apple")
    empty_template.add_compounds([cpd_apple])

    assert len(empty_template.compounds) == 1


def test_template_add_simple_metabolite_species(empty_template):
    cpd_apple = MSTemplateMetabolite("apple", "C100", "just a apple")
    empty_template.add_compounds([cpd_apple])

    comp_cpd_apple = MSTemplateSpecies("apple_k", 0, "k", "apple")
    empty_template.add_comp_compounds([comp_cpd_apple])

    assert len(empty_template.compounds) == 1
    assert len(empty_template.compcompounds) == 1
    assert empty_template.compcompounds.apple_k.compound
    assert empty_template.compcompounds.apple_k.compound.name == "just a apple"
    assert len(empty_template.compounds.apple.species) == 1


def test_template_add_compartment(empty_template):
    empty_template.compartments += [MSTemplateCompartment("w", "world", 4)]

    assert len(empty_template.compartments) == 1


def test_template_add_reaction(empty_template):
    cpd_apple = MSTemplateMetabolite("apple", "C100", "just a apple")
    cpd_apple_pie = MSTemplateMetabolite("appie", "C1000", "apple pie (10 apples)")
    empty_template.add_compounds([cpd_apple, cpd_apple_pie])

    comp_cpd_apple = MSTemplateSpecies("apple_k", 0, "k", "apple")
    comp_cpd_apple_pie = MSTemplateSpecies("appie_k", 0, "k", "appie")
    empty_template.add_comp_compounds([comp_cpd_apple, comp_cpd_apple_pie])

    rxn_make_pie = MSTemplateReaction(
        "rxn_pie_k", "rxn00000", "make pie", "pie", 0, 1000
    )
    rxn_make_pie.add_metabolites(
        {
            empty_template.compcompounds.apple_k: -10,
            empty_template.compcompounds.appie_k: 1,
        }
    )

    empty_template.add_reactions([rxn_make_pie])

    assert len(empty_template.reactions) == 1
    assert empty_template.reactions.rxn_pie_k.check_mass_balance() == {}
