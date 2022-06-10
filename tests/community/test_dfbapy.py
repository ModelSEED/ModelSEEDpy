# # -*- coding: utf-8 -*-

# import matplotlib, pandas
# import shutil
# from modelseedpy.fbapkg import dfbapkg

# # define the environment path 
# import os
# local_cobrakbase_path = os.path.join('C:', 'Users', 'Andrew Freiburger','Documents','Argonne','cobrakbase')
# os.environ["HOME"] = local_cobrakbase_path

# # import the models
# import cobrakbase
# token = 'WE6CHYRDTJSGOHFIDGPE7WYFT6PRPXJL'
# kbase_api = cobrakbase.KBaseAPI(token)
# model1 = kbase_api.get_from_ws("iML1515",76994)
# # model2 = kbase_api.get_from_ws("iSB1139.kb.gf",30650)

# # define in the initial conditions
# modelseed_db_path = os.path.join('..', '..', '..', 'ModelSEEDDatabase')
# initial_concentrations = {"ATP_c0":200, "ADP_c0":100}
# total_time = 5
# timestep = 1

# # execute the package
# dfba = dfbapkg.dFBA(model1, modelseed_db_path, warnings = False, verbose = False)
# conc, fluxes = dfba.simulate('model_kinetics.json', initial_concentrations, total_time, timestep, labeled_plots = False, export_name = 'test_dfba')

# def isnumber(string):
#     try:
#         float(string)
#         return True
#     except:
#         try:
#             int(string)
#             return True
#         except:
#             return False

# def test_init():
#     # assert qualities of the content
#     for dic in [
#             dfba.compound_names,
#             dfba.parameters,
#             dfba.variables,
#             dfba.variables['concentrations'],
#             dfba.variables['time_series'],
#             ]:
#         assert type(dic) is dict
#     assert type(dfba.met_names) is list
#     for bol in [dfba.verbose,dfba.printing, dfba.warnings, dfba.jupyter]:
#         assert type(bol) is bool
        
# def test_lp():
#     dfba.lp('test')
#     assert os.path.exists('test.lp')
#     shutil.rmtree('test.lp')
        
# def test_simulate():
#     for dic in [dfba.kinetics_data,dfba.defined_reactions]:
#         assert type(dic) is dict
#     for lis in [dfba.constrained,dfba.solutions]:
#         assert type(lis) is list
#     for string in [dfba.col]:
#         assert type(string) is str
#     for df in [dfba.concentrations,dfba.fluxes]:
#         assert type(df) is pandas.DataFrame
#     for quant in [
#             dfba.minimum,
#             dfba.parameters['timesteps'],
#             dfba.timestep_value,
#             dfba.total_time,
#             dfba.variables['elapsed_time'],
#             dfba.parameters['pH'],
#             dfba.parameters['temperature'],
#             dfba.cellular_dry_mass_fg,
#             dfba.cellular_fL
#             ]:
#         assert isnumber(quant)
#     for fig in [dfba.figure]:
#         assert type(dfba.figure) is matplotlib.figure.Figure
#     for st in [dfba.changed,dfba.unchanged]:
#         assert type(st) is set
        
#     assert os.path.exists(os.path.join(os.getcwd(), 'test_dfba'))
#     shutil.rmtree(dfba.simulation_path)