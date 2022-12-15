# -*- coding: utf-8 -*-
# import the COBRA model
# FIXME: COMMENTING ALL OF THIS will provide a model later

"""
import cobra
import optlang

# import cobrakbase
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
from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg
from modelseedpy.fbapkg.fullthermopkg import FullThermoPkg

# import the modelseed database and load the model
from modelseedpy.core.fbahelper import FBAHelper
modelseed_db_path = '..\..\..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
full_thermo = FullThermoPkg(model = model)
ms_api = FBAHelper.get_modelseed_db_api(modelseed_db_path)

# ------------------------ test the Full Thermo Package ---------------------------------------

def test_init():

    # assert results of the model
    assert type(full_thermo.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(full_thermo.name) is str
    assert type(full_thermo.variable_types) is dict
    assert type(full_thermo.constraint_types) is dict
    assert 'logconc' in list(full_thermo.variables.keys())
    assert 'dgerr' in list(full_thermo.variables.keys())
    assert 'potc' in list(full_thermo.constraints.keys())


def test_build_package():

    # execute the function
    full_thermo.build_package(parameters = {'modelseed_path':modelseed_db_path})

    # assert the addition of added parameters
    added_parameters = ["default_max_conc", "default_min_conc", "default_max_error", "custom_concentrations",'modelseed_api',
                        "custom_deltaG_error", "compartment_potential", "temperature", "filter",'modelseed_path']
    for param in added_parameters:
        assert param in full_thermo.parameters

    # assert results of the function
    added_parameters = ['combined_custom_concentrations', 'combined_custom_deltaG_error', 'combined_custom_comp_pot']

    for param in added_parameters:
        assert param in full_thermo.parameters

    for metabolite in full_thermo.model.metabolites:
        msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
        mscpd = full_thermo.parameters["modelseed_api"].get_seed_compound(msid)
        if (msid and mscpd) != None and mscpd.deltag != 10000000 and msid != "cpd00001":
            assert full_thermo.variables['dgerr'][metabolite.id]
            assert full_thermo.variables['logconc'][metabolite.id]
            assert full_thermo.constraints['potc'][metabolite.id]


    from scipy.constants import physical_constants, calorie, kilo, R
    Faraday = physical_constants['Faraday constant'][0]#C/mol
    from numpy import float64
    import re

    # define arbitrary initial conditions to test
    temperature = 298

    # assert results of the function
    for metabolite in full_thermo.model.metabolites:

        # evaluate the thermodynamic calculations
        compartment_potential = 0
        if metabolite.compartment in full_thermo.parameters["combined_custom_comp_pot"]:
            compartment_potential = full_thermo.parameters["combined_custom_comp_pot"][metabolite.compartment]

        msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
        mscpd = ms_api.get_seed_compound(msid)
        if mscpd != None:
            constant = mscpd.deltag/calorie + metabolite.charge*Faraday*compartment_potential/kilo

            assert type(constant) is float

        # execute and assert the results of the function
        if msid != "cpd00001":
            kind = 'logconc'
            full_thermo_conc_var = full_thermo.variables[kind][metabolite.id]
            assert full_thermo_conc_var.type == 'continuous'
            assert type(full_thermo_conc_var.lb) is float64
            assert type(full_thermo_conc_var.ub) is float64

        kind = 'dgerr'
        full_thermo_dgerr_var = full_thermo.variables[kind][metabolite.id]
        assert full_thermo_dgerr_var.type == 'continuous'
        assert full_thermo_dgerr_var.lb == -full_thermo_dgerr_var.ub

"""
