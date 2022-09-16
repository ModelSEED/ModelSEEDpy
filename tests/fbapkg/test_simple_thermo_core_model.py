# -*- coding: utf-8 -*-
# import cobrakbase
# FIXME: COMMENTING ALL OF THIS will provide a model later

"""
import os
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'KCYWCM5Q3ZFAGQKDG3ESMJXGHQPK7UMN'
kbase = cobrakbase.KBaseAPI(token)

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg
simple_thermo = SimpleThermoPkg(model = model)

import re

# ------------------------ test the Simple Thermo Package ---------------------------------------

def test_init():

    # assert results of the model
    assert type(simple_thermo.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(simple_thermo.name) is str
    assert type(simple_thermo.variable_types) is dict
    assert type(simple_thermo.constraint_types) is dict
    assert 'potential' in list(simple_thermo.variables.keys())
    assert 'thermo' in list(simple_thermo.constraints.keys())


def test_build_package():

    # execute the function
    simple_thermo.build_package(simple_thermo.parameters)

    # assert results from creating the potential variables and constraints
    added_parameters = ['filter', 'min_potential', 'max_potential']
    for param in added_parameters:
        assert param in simple_thermo.parameters
    assert 'potential' in list(simple_thermo.variables.keys())
    assert 'dgbinF' in list(simple_thermo.variables.keys())
    assert 'dgbinR' in list(simple_thermo.variables.keys())

    for metabolite in simple_thermo.model.metabolites:
        simple_thermo_var = simple_thermo.variables["potential"][metabolite.id]
        assert simple_thermo_var.type == 'continuous'
        assert simple_thermo_var.lb == 0
        assert simple_thermo_var.ub == 1000


    # assert results from creating the revbin variables and constraints
    constraint_types = {"F": [None, 0],
                        'R': [None, 1000]}

    for reaction in simple_thermo.model.reactions:
        revbin_var = simple_thermo.pkgmgr.getpkg("RevBinPkg").variables["revbin"][reaction.id]
        assert revbin_var
        assert revbin_var.lb == 0
        assert revbin_var.ub == 1
        assert revbin_var.type == 'binary'

        for type in constraint_types:
            revbin_cons = simple_thermo.pkgmgr.getpkg("RevBinPkg").constraints['revbin{}'.format(type)][reaction.id]
            assert revbin_cons
            assert revbin_cons.lb == constraint_types[type][0]
            assert revbin_cons.ub == constraint_types[type][1]
            assert revbin_cons.name == '{}_revbin{}'.format(reaction.id, type)


def test_dgbin():

    # execute the funtion
    simple_thermo.build_package({'dgbin': True})

    # assert results from the function
    for reaction in simple_thermo.model.reactions:
        if re.search('^EX_', reaction.id) == None and re.search('^SK', reaction.id) == None and re.search('^DM_', reaction.id) == None:

            dgbinF_var = simple_thermo.variables['dgbinF'][reaction.id]
            assert dgbinF_var
            assert dgbinF_var.lb == 0
            assert dgbinF_var.ub == 1
            assert dgbinF_var.type == "binary"

            dgbinR_var = simple_thermo.variables['dgbinR'][reaction.id]
            assert dgbinR_var
            assert dgbinR_var.lb == 0
            assert dgbinR_var.ub == 1
            assert dgbinR_var.type == "binary"

"""
