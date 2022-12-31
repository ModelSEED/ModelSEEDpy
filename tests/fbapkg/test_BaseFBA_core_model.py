# -*- coding: utf-8 -*-
# import cobrakbase
# FIXME: COMMENTING ALL OF THIS will provide a model later
"""

import optlang
import os
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'KCYWCM5Q3ZFAGQKDG3ESMJXGHQPK7UMN'
kbase = cobrakbase.KBaseAPI(token)

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
base = BaseFBAPkg(model = model, name = 'test_model', variable_types = {'concentration': 'reaction'}, constraint_types = {'concentration': 'reaction'})

# ------------------------ test the BaseFBA Package ---------------------------------------

def test_init():

    # assert results of the model
    assert type(base.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(base.name) is str
    assert type(base.variable_types) is dict
    assert type(base.constraint_types) is dict


def test_validate_parameters():
    # define arbitrary argument content for the function test
    params = {'a':2, 'b':4, 'c':3}
    required = ['a', 'b', 'c']
    defaults = {'a':1, 'b':1, 'c':1, 'd':1}

    # execute the function and assert results of the function
    base.validate_parameters(params, required, defaults)
    assert base.parameters == {'a': 2, 'b': 4, 'c': 3, 'd': 1}


def test_build_reaction_variable():
    # define the model instance of the test reaction
    test_reaction = 'rxn05297_c0'
    ub = 133
    lb = 10
    var_type = 'continuous'
    variable_type = 'concentration'
    model_reaction = model.reactions.get_by_id(test_reaction)

    # execute the function
    built_variable = base.build_variable(type = variable_type, lower_bound = lb, upper_bound = ub, vartype = var_type, object = model_reaction)

    # assert results of the function
    assert built_variable.ub == ub
    assert built_variable.lb == lb
    assert built_variable.type == var_type
    assert type(built_variable) is optlang.cplex_interface.Variable
    assert len(base.variables['concentration']) == 1
    assert built_variable.name == '{}_{}'.format(test_reaction, variable_type)


def test_build_reaction_constraint():
    # define arbitrary argument content for the function test
    test_reaction = 'rxn05297_c0'
    ub = 133
    lb = 10
    constraint_type = 'concentration'
    model_reaction = model.reactions.get_by_id(test_reaction)

    # execute the function
    built_constraint = base.build_constraint(type = constraint_type, lower_bound = lb, upper_bound = ub, object = model_reaction)

    # assert results of the function
    assert built_constraint.ub == ub
    assert built_constraint.lb == lb
    assert len(base.constraints['concentration']) == 1
    assert type(built_constraint) is optlang.cplex_interface.Constraint
    assert built_constraint.name == '{}_{}'.format(test_reaction, constraint_type)


def test_all_variables():
    # define the initial conditions of the model
    variables_quantity = 0
    for var in base.variables:
        variables_quantity += 1

    # execute the function
    instance_variables = len(base.all_variables())

    # assert results of the function
    assert variables_quantity == instance_variables


def test_all_constraints():
    # define the initial conditions of the model
    constraints_quantity = 0
    for var in base.constraints:
        constraints_quantity += 1

    # execute the function
    instance_constraints = len(base.all_constraints())

    # assert results of the function
    assert constraints_quantity == instance_constraints


def clear():
    # define the initial conditions of the model
    variables_quantity = 0
    constraints_quantity = 0
    for type in base.variables:
        for object in base.variables[type]:
            variables_quantity += 1
    for type in base.constraints:
        for object in base.constraints[type]:
            constraints_quantity += 1

    # execute the function
    instance_constraints = base.clear()

    # evaluate the reaction names between ModelSEED and COBRA
    assert variables_quantity > len(base.variables)
    assert len(base.variables) == 0
    assert constraints_quantity > len(base.constraints)
    assert len(base.constraints) == 0

"""
