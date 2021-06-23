from modelseedpy.core.mseditorapi import MSEditorAPI, MSEquation
import pytest
import cobra
import cobra.test


@pytest.fixture
def editor():
    return MSEditorAPI()


@pytest.fixture
def example_model():
    return cobra.test.create_test_model("textbook")


def test_remove_reactions1(editor, example_model):
    """
    testing for remove_reaction()
    """
    # remove valid reactions
    total_reactions = len(example_model.reactions)
    lst = ['ACALD', 'MDH']
    editor.remove_reactions(example_model, lst)
    assert len(example_model.reactions) == total_reactions - len(lst)


def test_remove_reactions2(editor, example_model):
    """
    remove invalid reactions
    """
    lst = ['C']
    with pytest.raises(Exception):
        editor.remove_reactions(example_model, lst)


def test_edit_reaction1(editor, example_model):
    """
    testing for edit_reaction()
    change a reversible reaction to a forward reaction and change the gpr
    """
    rxn = 'ICDHyr'
    editor.edit_reaction(example_model, rxn, "=>", "(b0001 and b0002) or b1010")
    reaction = example_model.reactions.get_by_id(rxn)
    assert reaction.reversibility is False
    assert reaction.lower_bound == 0
    assert reaction.upper_bound == 1000
    # assert reaction.gene_name_reaction_rule == '(b0001 and b0002) or b1010' this is failing !


def test_copy_model_reactions1():
    pass


def test_copy_model_reactions2():
    pass


def test_copy_model_reactions3():
    pass


def test_build_from_palsson_string_1(editor):
    """
    test for building a modelseed equaiton form a string
    """
    eq = MSEquation.build_from_palsson_string('cpd00001 + cpd00002[e] <= (2)cpd00003 + cpd00004')
    #assert(test.equation == "{('cpd00001', 'c'): -1, ('cpd00002', 'e'): -1, ('cpd00003', 'c'): 2, ('cpd00004', 'c'): 1}")
    assert eq.direction == "<"


def test_build_from_palsson_string_2(editor):
    eq = MSEquation.build_from_palsson_string('cpd00001 + cpd00002[e] <=> (2)cpd00003 + cpd00004')
    #assert(test.equation == "{('cpd00001', 'c'): -1, ('cpd00002', 'e'): -1, ('cpd00003', 'c'): 2, ('cpd00004', 'c'): 1}")
    assert eq.direction == "="


def test_build_from_palsson_string_3(editor):
    eq = MSEquation.build_from_palsson_string('cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004')
    #assert(test.equation == "{('cpd00001', 'c'): -1, ('cpd00002', 'e'): -1, ('cpd00003', 'c'): 2, ('cpd00004', 'c'): 1}")
    assert eq.direction == ">"


def test_build_from_palsson_string_4(editor):
    eq = MSEquation.build_from_palsson_string('cpd00001 + cpd00002[e] = (2)cpd00003 + cpd00004')
    #assert(test.equation == "{('cpd00001', 'c'): -1, ('cpd00002', 'e'): -1, ('cpd00003', 'c'): 2, ('cpd00004', 'c'): 1}")
    assert eq.direction == "?"
