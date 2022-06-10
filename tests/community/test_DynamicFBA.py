# # ------------------------ Define an arbitrary model --------------------------------------

# import optlang
# from cobra import Metabolite, Reaction, Model

# # define metabolites
# ACP_c = Metabolite(
#     'ACP_c',
#     formula='C11H21N2O7PRS',
#     name='acyl-carrier-protein',
#     compartment='c')
# omrsACP_c = Metabolite(
#     'M3omrsACP_c',
#     formula='C25H45N2O9PRS',
#     name='3-Oxotetradecanoyl-acyl-carrier-protein',
#     compartment='c')
# co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')
# malACP_c = Metabolite(
#     'malACP_c',
#     formula='C14H22N2O10PRS',
#     name='Malonyl-acyl-carrier-protein',
#     compartment='c')
# h_c = Metabolite('h_c', formula='H', name='H', compartment='c')
# ddcaACP_c = Metabolite(
#     'ddcaACP_c',
#     formula='C23H43N2O8PRS',
#     name='Dodecanoyl-ACP-n-C120ACP',
#     compartment='c')

# # define reactions, with the metabolites
# reaction = Reaction('R_3OAS140')
# reaction.name = 'R_3OAS140'
# reaction.add_metabolites({
#     malACP_c: -1.0,
#     h_c: -1.0,
#     ddcaACP_c: -1.0,
#     co2_c: 1.0,
#     ACP_c: 1.0,
#     omrsACP_c: 1.0
# })

# reaction2 = Reaction('rxn2')
# reaction2.add_metabolites({
#     malACP_c: -1.0,
#     h_c: 1.0,
#     ddcaACP_c: -1.0,
#     co2_c: -1.0,
#     ACP_c: 1.0
# })
# reaction2.name = 'rxn2'

# # define the model, with the reactions
# model = Model('test')
# model.add_reactions([reaction, reaction2])
# model.objective = reaction.flux_expression

# # ------------------ define arbitrary initial conditions and kinetics data -------------------------------

# source = 'custom'
# initial_concentrations = {"Dodecanoyl-ACP-n-C120ACP":200, 'Malonyl-acyl-carrier-protein':22, 'CO2': 12, 'H': 14}
# reaction_kinetics = {
#     'R_3OAS140': {
#         'SubstitutedRateLaw': '(68.0*A*B)/(50.0*0.34*C+360.0*B+0.34*A+A*B*C)',
#         'parameters': {
#             'A': 'Dodecanoyl-ACP-n-C120ACP',
#             'B': 'H',
#             'C': 'malACP_c'
#         }
#     },
#     'rxn2': {
#         'SubstitutedRateLaw': '(A*B)/(50.0*0.34*C+3*B+0.34*A+C)',
#         'parameters': {
#             'A': 'Dodecanoyl-ACP-n-C120ACP',
#             'B': 'CO2',
#             'C': 'Malonyl-acyl-carrier-protein'
#         }
#     }
# }

# # ------------------------ test the dynamicFBA Package ---------------------------------------

# # import statements
# from numpy import float64 
# from math import inf
# from numpy import nan   
# import pandas
# import cobra
# import json
# import re

# empty = ['', '-', '', None, nan, 'NaN']
   
# # add the units of logarithm to the Magnesium concentration
# def isnumber(string):
#     string = str(string)
#     try:
#         string = string.strip()
#         if re.sub('([0-9\.\-])', '',string) == '':
#             return True
#     except:
#         print(string)
#         return False
    
# def average(num_1, num_2 = None):
#     if type(num_1) is list:
#         average = sum(num_1) / len(num_1)
#         return average
#     elif isnumber(num_1): 
#         if num_2 is not None:
#             numbers = [num_1, num_2]
#             average = sum(numbers) / len(numbers)
#             return average
#         else:
#             return num_1
#     else:
#         return None

# # import the modelseedpy packages
# from modelseedpy.fbapkg.dynamicfbapkg import dynamicFBAPkg
# dfba = dynamicFBAPkg(model, total_time = 5, timestep = 1, initial_concentrations = initial_concentrations, reaction_kinetics = reaction_kinetics, logging = True)

# def test_init():
#     # assert results of the model 
#     assert type(dfba.model) is cobra.core.model.Model
#     assert dfba.verbose is False
#     assert type(dfba.parameters) is dict
#     assert type(dfba.variables) is dict
#     assert ('timesteps' and 'reaction_kinetics' and 'initial_concentrations') in list(dfba.parameters.keys())
#     assert ('time_series' and 'concentrations') in list(dfba.variables.keys())
    
#     # ensure that the dataframes are generated
#     assert type(dfba.fluxes_df) is pandas.core.frame.DataFrame
#     assert type(dfba.concentrations_df) is pandas.core.frame.DataFrame
    
# def test_simulate():
#     # execute the function
# #     dfba.simulate(source = 'custom')
#     temperature = 25
#     p_h = 7
    
#     # ------------------------ execute the simulate mechanics -----------------------------------
#     solutions = []
#     for dfba.timestep in range(1,dfba.parameters['timesteps']+1):
#         if dfba.timestep == 1:
#             dfba.calculate_kinetics(dfba.parameters['initial_concentrations'], source, temperature, p_h)
#         else:
#             dfba.calculate_kinetics(dfba.variables['concentrations'], source, temperature, p_h)

#         #Calcuate and parameterize fluxes from michaelis-menten kinetics
#         for reaction in dfba.model.reactions:
#             if any(rxn.lower() == reaction.name.lower() for rxn in dfba.parameters['calculated_rate_laws']):
#                 kinetic_flux = dfba.parameters['calculated_rate_laws'][reaction.name.lower()]
#             elif any(rxn.lower() == reaction.id.lower() for rxn in dfba.parameters['calculated_rate_laws']):
#                 kinetic_flux = dfba.parameters['calculated_rate_laws'][reaction.id.lower()]
#             else:
#                 if dfba.verbose:
#                     print(f'--> ERROR: The {reaction.name} reaction is not defined in the kinetics data.')    
#                 continue

#             if not isnumber(kinetic_flux):
#                 print(f'--> ERROR: The constant for {reaction.name} is erronenous.')
#                 continue
#             if kinetic_flux > reaction.lower_bound:
#                 reaction.upper_bound = kinetic_flux
#                 reaction.lower_bound = kinetic_flux
#             else:
#                 reaction.lower_bound = kinetic_flux
#                 reaction.upper_bound = kinetic_flux
                
#             assert isnumber(reaction.lower_bound)
#             assert isnumber(reaction.upper_bound)
#             assert reaction.lower_bound == reaction.upper_bound == kinetic_flux

#             print(f'enzyme {reaction.name}', 'kinetic_flux', kinetic_flux)

#         # execute the model and update concentrations
#         solution = dfba.model.optimize()
#         solutions.append(solution)
#         print(f'\nobjective value for timestep {dfba.timestep}: ', solution.objective_value, '\n\n')
#         dfba.update_concentrations(solution)
    
#     for solution in solutions:
#         assert type(solution) is cobra.core.solution.Solution
    
#     # evaluate the dataframes
#     assert type(dfba.fluxes_df) is pandas.core.frame.DataFrame    
#     for index, row in dfba.fluxes_df.iterrows():
#         for entry in row:
#             if type(entry) is float:
#                 assert type(entry) == float
#             elif type(entry) is float64:
#                 assert type(entry) == float64
#             else:
#                 assert type(entry) is None
       
#     assert type(dfba.concentrations_df) is pandas.core.frame.DataFrame
#     for index, row in dfba.concentrations_df.iterrows():
#         for entry in row:
#             assert type(entry) == (float64 or float)