# -*- coding: utf-8 -*-

# cobra_config = cobra.Configuration()

# from cobra import Model, Reaction, Metabolite
import cobra.test
import os
from os.path import join

data_dir = cobra.test.data_dir


def test_edit_reaction2():
    # change a forward reaction to a reversible reaction
    rxn = "HXPRT"
    editor.edit_reaction(model, rxn, "<=>")
    assert model.reactions.get_by_id(rxn).reversibility == True
    assert model.reactions.get_by_id(rxn).lower_bound == -1000
    assert model.reactions.get_by_id(rxn).upper_bound == 1000


def test_edit_reaction3():
    # change a forward reaction to a reverse reaction
    rxn = "CYTDK2"
    editor.edit_reaction(model, rxn, "<=")
    assert model.reactions.get_by_id(rxn).reversibility == False
    assert model.reactions.get_by_id(rxn).lower_bound == -1000
    assert model.reactions.get_by_id(rxn).upper_bound == 0


# may throw exception because error handling was difficult, uncomment at your own risk
# def test_edit_reaction4():
# invalid gpr
# with pytest.raises(SyntaxError):
# editor.edit_reaction(model,rxn1,"=>","(b0001 and b0002 or b1010")

# testing for edit_biomass_compound
def test_edit_biomass_compound1():
    # properly change coefficient
    model2 = model
    editor.edit_biomass_compound(model2, "CYTDK2", "cmp_c", 2)
    assert model2.reactions.get_by_id("CYTDK2").get_coefficient("cmp_c") == 2


def test_edit_biomass_compound2():
    # nonexistent reaction
    model3 = model
    with pytest.raises(Exception):
        editor.edit_biomass_compound(model3, "CYTDK", "cmp_c", 2)


def test_edit_biomass_compound3():
    # nonexistent metabolite
    model4 = model
    with pytest.raises(Exception):
        editor.edit_biomass_compound(model4, "CYTDK2", "cmp_", 2)


# test molecular weight test
def test_molecular_weight_1():
    model = cobra.io.load_json_model("iML1515.json")
    assert editor.compute_molecular_weight(model, "octapb_c") == 127.2041


# tests for adding a custom reaction
def test_add_custom_reaction_1():
    test = editor.build_from_palsson_string("octapb + cysi__L[e] => (2)dhap + prbatp")
    model = cobra.io.load_json_model("iML1515.json")
    editor.add_custom_reaction(model, "test_id", test)
    assert model.reactions.get_by_id("test_id").lower_bound == 0
    assert model.reactions.get_by_id("test_id").upper_bound == 1000


def test_add_custom_reaction_2():
    test = editor.build_from_palsson_string("octapb + cysi__L[e] <=> (2)dhap + prbatp")
    model = cobra.io.load_json_model("iML1515.json")
    editor.add_custom_reaction(model, "test_id", test)
    assert model.reactions.get_by_id("test_id").lower_bound == -1000
    assert model.reactions.get_by_id("test_id").upper_bound == 1000


def test_add_custom_reaction_3():
    test = editor.build_from_palsson_string("octapb + cysi__L[e] <= (2)dhap + prbatp")
    model = cobra.io.load_json_model("iML1515.json")
    editor.add_custom_reaction(model, "test_id", test)
    assert model.reactions.get_by_id("test_id").lower_bound == -1000
    assert model.reactions.get_by_id("test_id").upper_bound == 0


# testing for copy_model_reactions
def test_copy_model_reactions1():
    model2 = model
    lst = ["CYTDK2"]
    editor.remove_reactions(model2, lst)
    assert len(model2.reactions) == 2712 - len(lst)
    editor.copy_model_reactions(model2, model, lst)
    assert len(model.reactions) == 2712
    assert len(model2.reactions) == 2712


def test_copy_model_reactions2():
    # copying an id that already exists in the duplicate
    model3 = model
    lst = ["CYTDK2"]
    editor.copy_model_reactions(model3, model, lst)
    assert len(model3.reactions) == 2712


def test_copy_model_reactions3(self):
    # testing on copying a id that does not exist in the origional
    model4 = model
    lst = ["C"]
    with pytest.raises(Exception):
        editor.copy_model_reactions(model4, model, lst)
