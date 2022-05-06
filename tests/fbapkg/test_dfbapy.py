# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 12:03:59 2022

@author: Andrew Freiburger
"""
import matplotlib, pandas
import shutil
import dfbapy
import os

def isnumber(string):
    try:
        float(string)
        return True
    except:
        try:
            int(string)
            return True
        except:
            return False

model_path = os.path.join(os.path.dirname(__file__),'iSB619.xml')
kinetics_path = os.path.join(os.path.dirname(__file__),'model_kinetics.json')

def test_init():
    dfba = dfbapy.dFBA(bigg_model_path, 'glpk')
    
    # assert qualities of the content
    for dic in [
            dfba.bigg_metabolites_ids,
            dfba.bigg_metabolites_names,
            dfba.parameters,
            dfba.variables,
            dfba.variables['concentrations'],
            dfba.variables['time_series'],
            ]:
        assert type(dic) is dict
    for lis in [dfba.model_ids,dfba.model_names]:
        assert type(lis) is list
    for string in [dfba.parameters['bigg_model_name']]:
        assert type(string) is str
    for bol in [dfba.verbose,dfba.printing]:
        assert type(bol) is bool
        
def test_simulate():
    # load the model
    dfba = dfbapy.dFBA(bigg_model_path, 'glpk')
    dfba.simulate(kinetics_path, export_name = 'test_dfba')
    
    # assert qualities of the simulation
    for dic in [dfba.kinetics_data,dfba.defined_reactions]:
        assert type(dic) is dict
    for lis in [dfba.constrained,dfba.solutions]:
        assert type(lis) is list
    for string in [dfba.col]:
        assert type(string) is str
    for df in [dfba.concentrations,dfba.fluxes]:
        assert type(df) is pandas.DataFrame
    for quant in [
            dfba.minimum,
            dfba.parameters['timesteps'],
            dfba.timestep_value,
            dfba.total_time,
            dfba.variables['elapsed_time'],
            dfba.parameters['pH'],
            dfba.parameters['temperature']
            ]:
        assert isnumber(quant)
    for fig in [dfba.figure]:
        assert type(dfba.figure) is matplotlib.figure.Figure
    for st in [dfba.changed,dfba.unchanged]:
        assert type(st) is set
        
    assert os.path.exists(os.path.join(os.getcwd(), 'test_dfba'))
    shutil.rmtree(dfba.simulation_path)