# import the packages
import os
import cobra
import numpy
# FIXME: COMMENTING ALL OF THIS will provide a model later
"""
from glob import glob
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'xx'
kbase = cobrakbase.KBaseAPI(token)
import re

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
import modelseedpy
from modelseedpy.core.msgapfill import MSGapfill
gapfill = MSGapfill(model)

def test_init():
    assert type(gapfill.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    assert type(gapfill.blacklist) is list
    assert type(gapfill.solutions) is dict
    
def test_run_gapfilling_and_integrate_gapfill_solution():
    solutions = gapfill.run_gapfilling()
    
    # test that the objective expression is correctly set
    if solutions is not None:
        assert type(solutions) is dict
        
        # verify the integrate_gapfill_solution function
        model_2 = gapfill.integrate_gapfill_solution(solutions)
        assert type(model_2) is cobrakbase.core.kbasefba.fbamodel.FBAModel
        
        for reaction in solutions['reversed']:
            if solution["reversed"][reaction] == ">":
                assert reaction.upper_bound == 100
            else:
                assert reaction.lower_bound == -100
                
        for reaction in solutions['new']:
            if solution["new"][reaction] == ">":
                assert reaction.upper_bound == 100
                assert reaction.lower_bound == 0 
            else:
                assert reaction.upper_bound == 0
                assert reaction.lower_bound == -100

def test_gapfill():
    pass   
"""