# import cobrakbase
import os
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = '5QQKGJK7BYX7HF7M2TFI3EVJXC7NE67T'
kbase = cobrakbase.KBaseAPI(token)

# define the example individual model and associated API media package
model = kbase.get_from_ws('E_iAH991V2', 93832)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg
simple_thermo = SimpleThermoPkg(model = model)  

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
    
    for metabolite in simple_thermo.model.metabolites:
        simple_thermo_var = simple_thermo.variables["potential"][metabolite.id]
        assert simple_thermo_var.type == 'continuous'
        assert simple_thermo_var.lb == 0
        assert simple_thermo_var.ub == 1000
        
        
    # assert results from creating the revbin variables and constraints
    constraint_types = {"F": [None, 0],
                        'R': [None, 1000]}
    
    for reaction in simple_thermo.model.reactions:
        revbin_var = simple_thermo.variables["revbin"][reaction.id]
        assert revbin_var
        assert revbin_var.lb == 0
        assert revbin_var.ub == 1
        assert revbin_var.type == 'binary'
        
        for type in constraint_types:
            revbin_cons = simple_thermo.constraints['revbin{}'.format(type)][reaction.id]
            assert revbin_cons
            assert revbin_cons.lb == constraint_types[type][0]
            assert revbin_cons.ub == constraint_types[type][1]    
            assert revbin_cons.name == '{}_revbin{}'.format(reaction.id, type)

            
def test_build_constraint():
    
    # define arbitrary initial conditions to test
    test_reaction = '3OAR140_c0'
    model_reaction = simple_thermo.model.reactions.get_by_id(test_reaction)
    
    # execute the function    
    built_constraint = simple_thermo.build_constraint(object = model_reaction)
      
    # assert results of the function
    import re
    stoichiometric_list = []
    for var in built_constraint.variables:
        if not re.search('_revbin', str(var)):
            stoich = float(built_constraint.expression.coeff(var))
            stoichiometric_list.append(stoich)
        
    stoichiometric_list2 = []
    for metabolite in model_reaction.metabolites:
        stoichiometry = float(model_reaction.metabolites[metabolite])
        stoichiometric_list2.append(stoichiometry)

    assert set(stoichiometric_list) == set(stoichiometric_list2)