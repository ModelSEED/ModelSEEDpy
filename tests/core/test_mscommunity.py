# import the packages
import os
os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
import cobrakbase
token = 'TKWQA762H2SMAYRES3BFUP37LKGNGQSM'
kbase = cobrakbase.KBaseAPI(token)
import re

# define the example individual model and associated API media package
model = kbase.get_from_ws('e_coli_core.kb', 95098)
model.solver = 'optlang-cplex'

# import the modelseedpy packages
from modelseedpy.core.mscommunity import MSCommunity
ftp_path = '../../../ModelSEEDDatabase'
cfba = MSCommunity(model)

def test_init():
    assert type(cfba.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    print(cfba.pkgmgr)
    assert cfba.pkgmgr is dict
    assert cfba.metabolite_uptake is dict
    
def test_set_objective():
    cfba.set_objective()
    assert cfba.model.objective.direction == 'max'
    assert cfba.model.objective.expression == cfba.model.reactions.get_by_id('bio1').flux_expression
    
def test_gapfill():
    solution = cfba.gapfill()
    # test that the objective expression is correctly set
    assert cfba.model.objective.direction == 'max'
    assert cfba.model.objective.expression == cfba.model.reactions.get_by_id('bio1').flux_expression
    
    # test that the gapfilling yields a solution
    assert solution is not None
    
def test_run():
    
    
    
    
def test_compute_interactions():
    
    
    
def test_constrain():
    
    
    
    
def test_visualize():
    
    
    
def test_community_fba():
    
    
    