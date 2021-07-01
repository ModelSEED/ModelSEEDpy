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
    
    
@pytest.fixture
def example_model2():
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
    assert len(reaction.genes) == 3

def test_edit_reaction2(editor, example_model):
    """
    testing for edit_reaction()
    change a forward reaction to a reversible reaction
    """
    rxn = 'CS'
    editor.edit_reaction(example_model, rxn, "<=>")
    reaction = example_model.reactions.get_by_id(rxn)
    assert reaction.reversibility is True
    assert reaction.lower_bound == -1000
    assert reaction.upper_bound == 1000

def test_edit_reaction3(editor, example_model):
    """
    testing for edit_reaction()
    change a forward reaction to a reverse reaction
    """
    rxn = 'CYTBD'
    editor.edit_reaction(example_model, rxn, "<=")
    reaction = example_model.reactions.get_by_id(rxn)
    assert reaction.reversibility is False
    assert reaction.lower_bound == -1000
    assert reaction.upper_bound == 0



def test_copy_model_reactions1(editor,example_model,example_model2):
    """
    testing for copy_model_reaction()
    copying reactions from one list to another
    """
    lst = ['ICDHyr','CYTBD']
    editor.remove_reactions(example_model2,lst)
    assert len(example_model2.reactions) == 95 - len(lst)
    editor.copy_model_reactions(example_model2, example_model, lst)
    assert len(example_model2.reactions) == 95

def test_copy_model_reactions2(editor,example_model,example_model2):
    """
    testing for copy_model_reaction()
    copying a reaction that already in the original model
    """
    lst = ['ICDHyr','CYTBD']
    editor.copy_model_reactions(example_model2, example_model, lst)
    assert len(example_model2.reactions) == 95

def test_copy_model_reactions3(editor,example_model,example_model2):
    """
    testing for copy_model_reaction()
    copying a reaction that doesn't exist in the original model
    """
    lst = ['C','ICDHyr']
    with pytest.raises(Exception):
        editor.copy_model_reactions(example_model2, example_model, lst)
    assert len(example_model2.reactions) == 95
    
    
    
    
def test_copy_all_model_reactions1(editor,example_model,example_model2):
    pass
    """
    testing for copy_all_model_reactions()
    copying all reactions from a source model that don't already exist in the receiving model
    """
    lst = ['ICDHyr']
    editor.remove_reactions(example_model2,lst)
    editor.copy_all_model_reactions(example_model2,example_model)
    assert len(example_model2.reactions) == 95



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
