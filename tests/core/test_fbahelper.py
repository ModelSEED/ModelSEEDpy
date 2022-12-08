# -*- coding: utf-8 -*-
# # import the packages
# import os
# import cobra
# import numpy
# from glob import glob
# os.environ["HOME"] = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\cobrakbase'
# import cobrakbase
# token = 'xx'
# kbase = cobrakbase.KBaseAPI(token)
# import re

# # define the example individual model
# model = kbase.get_from_ws('e_coli_core.kb', 95098)
# media = kbase.get_from_ws('e_coli_core.media', 95098)
# model.solver = 'optlang-cplex'

# # import the modelseedpy packages
# import modelseedpy
# from modelseedpy import biochem
# from modelseedpy.core.fbahelper import FBAHelper
# from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
# pkgmgr = MSPackageManager.get_pkg_mgr(model)


# def test_drain_functions():
#     # define an arbitrary auto_sink
#     auto_sink_compound = 'cpd00169'
#     met_id = auto_sink_compound + '_c0'
#     cobra_id = model.metabolites.get_by_id(met_id)

#     drain_reactions = FBAHelper.add_autodrain_reactions_to_community_model(model, [auto_sink_compound])
#     model.add_reactions(drain_reactions)
#     for reaction in model.reactions:
#         if re.search('DM_', reaction.id):
#             assert reaction.id == f'DM_{met_id}'
#             assert reaction.metabolites[cobra_id] == -1
#             assert reaction.lower_bound == 0
#             assert reaction.upper_bound == 100
#             assert reaction.annotation['sbo'] == 'SBO:0000627'

# def test_test_condition_list():
#     pkgmgr = MSPackageManager.get_pkg_mgr(model)

#     condition_list = [{'objective':'bio1', 'media':media, 'threshold':0.9,'is_max_threshold':False},{'objective':'bio1','media':media,'threshold':1,'is_max_threshold':True}]

#     passage = FBAHelper.test_condition_list(model, condition_list, pkgmgr)

#     assert passage is True

# def test_reaction_expansion_test():
#     reaction_list = []
#     count = 0
#     for reaction in model.reactions:
#         if count < len(model.reactions)/2:
#             reaction_list.append({0:reaction, 1:'>'})
#         else:
#             reaction_list.append({0:reaction, 1:'<'})

#     condition_list = [{'objective':'bio1', 'media':media, 'threshold':0.9,'is_max_threshold':False},{'objective':'bio1','media':media,'threshold':1,'is_max_threshold':True}]

#     filtered_list = FBAHelper.reaction_expansion_test(model, reaction_list, condition_list, pkgmgr)

#     assert type(filtered_list) is list

# def test_set_reaction_bounds_from_direction():
#     count = 0
#     for reaction in model.reactions:
#         direction = '<'
#         if count > len(model.reactions)/2:
#             direction = '>'
#         FBAHelper.set_reaction_bounds_from_direction(reaction, direction)

#         if direction == '<':
#             assert reaction.lower_bound == -100
#             assert reaction.upper_bound == 0
#         elif direction == '>':
#             assert reaction.lower_bound == 0
#             assert reaction.upper_bound == 100

# def test_set_objective_from_target_reaction():
#     target_reaction = 'bio1'
#     reaction = FBAHelper.set_objective_from_target_reaction(model, target_reaction)

#     assert type(reaction) is cobrakbase.core.kbasefba.fbamodel_biomass.Biomass
#     assert str(model.objective.expression) == '1.0*bio1 - 1.0*bio1_reverse_b18f7'
#     assert model.objective.direction == 'max'

#     reaction = FBAHelper.set_objective_from_target_reaction(model, target_reaction, minimize = True)
#     assert model.objective.direction == 'min'

# def test_compute_flux_values_from_variables():
#     flux_values = FBAHelper.compute_flux_values_from_variables(model)

#     assert type(flux_values) is dict

# def test_modelseed_id_from_cobra_metabolite():
#     for met in model.metabolites:
#         returned = FBAHelper.modelseed_id_from_cobra_metabolite(met)
#         match = re.search('^(cpd\d+)', str(returned))

#         assert match is not None
#         assert type(returned) is str

# def test_modelseed_id_from_cobra_reaction():
#     for rxn in model.reactions:
#         returned = FBAHelper.modelseed_id_from_cobra_reaction(rxn)
#         match = re.search('^(rxn\d+)', str(returned))

#         if match:
#             assert match is not None
#             assert type(returned) is str
#         else:
#             assert match is None

# def test_metabolite_mw():
#     pyruvate = 'cpd00020_c0'
#     pyruvate_mass = 87.05  # https://pubchem.ncbi.nlm.nih.gov/compound/pyruvate

#     for met in model.metabolites:
#         if met.id == pyruvate:
#             pyruvate = met

#     calculated_mass = FBAHelper.metabolite_mw(pyruvate)

#     assert pyruvate_mass == round(calculated_mass, 2)

# def test_elemental_mass():
#     elementmass = FBAHelper.elemental_mass()

#     assert type(elementmass) is dict

# def test_get_modelseed_db_api():
#     msdb_path = 'C:\\Users\\Andrew Freiburger\\Dropbox\\My PC (DESKTOP-M302P50)\\Documents\\UVic Civil Engineering\\Internships\\Agronne\\TFA\\ModelSEEDDatabase'
#     db_api = FBAHelper.get_modelseed_db_api(msdb_path)

#     assert type(db_api) is modelseedpy.biochem.modelseed_biochem.ModelSEEDBiochem

# def test_is_ex():
#     for reaction in model.reactions:
#         result = FBAHelper.is_ex(reaction)

#         if len(reaction.id) > 3 and reaction.id[0:3] in ["EX_", "DM_", "SK_"]:
#             assert result is True
#         else:
#             assert result is False

# def test_is_biomass():
#     for reaction in model.reactions:
#         result = FBAHelper.is_biomass(reaction)

#         if reaction.id[0:3] == 'bio':
#             assert result is True
#         else:
#             assert result is False
