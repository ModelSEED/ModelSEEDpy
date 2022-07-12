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
from modelseedpy.core.mscommunity import MSCommunity
ftp_path = '../../../ModelSEEDDatabase'
cfba = MSCommunity(model)

def test_init():
    assert type(cfba.model) is cobrakbase.core.kbasefba.fbamodel.FBAModel
    print(cfba.pkgmgr)
    assert type(cfba.pkgmgr) is modelseedpy.fbapkg.mspackagemanager.MSPackageManager
    
def test_gapfill():
    solution = cfba.gapfill()
    
    # test that the objective expression is correctly set
    assert cfba.model.objective.direction == 'max'
    assert cfba.model.objective.expression == cfba.model.reactions.get_by_id('bio1').flux_expression
    
    # test that the gapfilling yields a solution
    if solution is not None:
        assert type(solution) is cobrakbase.core.kbasefba.fbamodel.FBAModel   
    
def test_drain_fluxes():
    drains = cfba.drain_fluxes()
    
    # verify that drain reactions were added
    reactions = []
    for reaction in cfba.model.reactions:
        reactions.append(reaction.id)
    reactions_string = ' '.join(reactions)
    match = re.search('DM_', reactions_string)
    if match is not None:
        assert type(match.group()) is str 
    
    # verify that the returned dictionary contains appropriate content
    assert type(drains) is dict
    for flux in drains:
        assert flux.isnumeric()
        assert type(drains[flux]) is float
    
def test_run_and_compute_interactions():
    cfba.constrain()
    print('suffix:', cfba.suffix)
    solution = cfba.run(print_lp = True)
    
    # verify that the LP file was created
    assert len(glob('*.lp')) >= 1
    # verify that the object is properly set
    assert cfba.model.objective.direction == 'max'
    assert cfba.model.objective.expression == cfba.model.reactions.get_by_id('bio1').flux_expression
    # verify that the optimize function executed
    if solution is not None:
        assert type(solution) is cobra.core.solution.Solution
    
    cfba.compute_interactions(solution)
    
    # verify the population of the metabolite_uptake dictionary
    assert type(cfba.metabolite_uptake) is dict
    for (id, index), rate in cfba.metabolite_uptake:
        assert type(id) is string
        assert index.isnumeric()
        assert type(rate) is float
    
    # verify the creation of the consumption and production matrices 
    assert type(cfba.consumption) is numpy.ndarray
    assert type(cfba.production) is numpy.ndarray
    
def test_constrain():
    # verify that the appropriate constraints could be implemented    
    pkg = cfba.pkgmgr.getpkg('KBaseMediaPkg')
    assert type(pkg) is modelseedpy.fbapkg.kbasemediapkg.KBaseMediaPkg
    
    pkg = cfba.pkgmgr.getpkg('ElementUptakePkg')
    assert type(pkg) is modelseedpy.fbapkg.elementuptakepkg.ElementUptakePkg
    
    pkg = cfba.pkgmgr.getpkg('CommKineticPkg')
    assert type(pkg) is modelseedpy.fbapkg.commkineticpkg.CommKineticPkg
    
    pkg = cfba.pkgmgr.getpkg('FullThermoPkg')
    assert type(pkg) is modelseedpy.fbapkg.fullthermopkg.FullThermoPkg

def test_community_fba():
    pass   
"""