# -*- coding: utf-8 -*-

# import statements
from modelseedpy import BaseFBAPkg
from math import inf
from numpy import nan   
import pandas
import cobra
import json
import re

empty = ['', '-', '', None, nan, 'NaN']
   
# add the units of logarithm to the Magnesium concentration
def isnumber(string):
    string = str(string)
    try:
        string = string.strip()
        if re.sub('([0-9]|\.)', '',string) == '':
            return True
    except:
        return False
    
def average(num_1, num_2 = None):
    if type(num_1) is list:
        average = sum(num_1) / len(num_1)
        return average
    elif isnumber(num_1): 
        if num_2 is not None:
            numbers = [num_1, num_2]
            average = sum(numbers) / len(numbers)
            return average
        else:
            return num_1
    else:
        return None
            
# define chemical concentrations
class dynamicFBAPkg(BaseFBAPkg):
    def __init__(self, model, total_time, timestep, initial_concentrations, reaction_kinetics, logging = True, verbose = True):
        self.model = model
        self.verbose = verbose
        self.parameters = {}
        self.variables = {}
#         self.constraints = {}
#         self.constraint_types = {}
        
#         self.constraints['kinetic'] = {}
#         self.constraint_types['kinetic'] = {}
        self.parameters['timesteps'] = round(total_time/timestep)
        self.parameters['reaction_kinetics'] = reaction_kinetics
        self.parameters['initial_concentrations'] = initial_concentrations
        self.variables['time_series'] = {}
        self.variables['concentrations'] = {}
        for metabolite in self.model.metabolites:
            self.variables['time_series'][metabolite.name] = []
            
        # define the dataframe for the time series content
        if logging:
            # define the column for each timestep
            columns = [f'{(time+1)*timestep} (min)' for time in range(self.parameters['timesteps'])]
            
            # create the fluxes dataframe
            indices = [enzyme for enzyme in reaction_kinetics]
            self.fluxes_df = pandas.DataFrame(index = indices, columns = columns)
            self.fluxes_df.index.name = 'Enzyme names'
            
            # create the concentrations dataframe
            indices = list(set(f'{met.name} ({met.compartment})' for met in self.model.metabolites))
            self.concentrations_df = pandas.DataFrame(index = indices, columns = columns)
            self.concentrations_df.index.name = 'Metabolite names (compartment)'
            

    def simulate(self, source = 'sabio', temperature = 25, p_h = 7, visualize = True):
        #Simulate for each timestep within total time frame
        solutions = []
        for self.timestep in range(1,self.parameters['timesteps']+1):
            if self.verbose:
                print('timestep', self.timestep)
            if self.timestep == 1:
                self.calculate_kinetics(self.parameters['initial_concentrations'], source, temperature, p_h)
            else:
                self.calculate_kinetics(self.variables['concentrations'], source, temperature, p_h)

            #Calcuate and parameterize fluxes from michaelis-menten kinetics
            for reaction in self.model.reactions:
                if any(rxn.lower() == reaction.name.lower() for rxn in self.parameters['calculated_rate_laws']):
                    kinetic_flux = self.parameters['calculated_rate_laws'][reaction.name.lower()]
                elif any(rxn.lower() == reaction.id.lower() for rxn in self.parameters['calculated_rate_laws']):
                    kinetic_flux = self.parameters['calculated_rate_laws'][reaction.id.lower()]
                else:
                    if self.verbose:
                        print(f'--> ERROR: The {reaction.name} reaction is not defined in the kinetics data.')    
                    continue
                        
                if not isnumber(kinetic_flux):
                    print(f'--> ERROR: The constant for {reaction.name} is erronenous.')
                    continue
#                 else:
#                     coef = {reaction.forward_variable:1, reaction.reverse_variable:1}
#                     BaseFBAPkg.build_constraint(self, "kinetic",kinetic_flux,kinetic_flux,coef,reaction)
                if kinetic_flux > reaction.lower_bound:
                    reaction.upper_bound = kinetic_flux
                    reaction.lower_bound = kinetic_flux
                else:
                    reaction.lower_bound = kinetic_flux
                    reaction.upper_bound = kinetic_flux
                    
                print(f'enzyme {reaction.name}', 'kinetic_flux', kinetic_flux)

            # execute the model and update concentrations
            solution = self.model.optimize()
            solutions.append(solution)
            print(f'\nobjective value for timestep {self.timestep}: ', solution.objective_value, '\n\n')
            self.update_concentrations(solution)
        
        if visualize:
            self.visualize()
        
        return self.fluxes_df, self.concentrations_df, solutions
    
    def calculate_kinetics(self,concentrations, source = 'sabio', temperature = 25, p_h = 7):
        self.parameters['calculated_rate_laws'] = {}

        for enzyme in self.parameters['reaction_kinetics']:
            entry = False
            A = 0
            B = 0
            C = 0
            S = 0
            
            if source == 'custom':
                remainder = re.sub('([0-9ABCS\-\/\(\)\+\.\*])', '', self.parameters['reaction_kinetics'][enzyme]["SubstitutedRateLaw"])
                if remainder == '':
                    if "A" in self.parameters['reaction_kinetics'][enzyme]['parameters']:   
                        if self.parameters['reaction_kinetics'][enzyme]['parameters']["A"] in concentrations: 
                            A = concentrations[self.parameters['reaction_kinetics'][enzyme]['parameters']["A"]]

                    if "B" in self.parameters['reaction_kinetics'][enzyme]['parameters']:   
                        if self.parameters['reaction_kinetics'][enzyme]['parameters']["B"] in concentrations:    
                            B = concentrations[self.parameters['reaction_kinetics'][enzyme]['parameters']["B"]]

                    if "C" in self.parameters['reaction_kinetics'][enzyme]['parameters']:   
                        if self.parameters['reaction_kinetics'][enzyme]['parameters']["C"] in concentrations:    
                            C = concentrations[self.parameters['reaction_kinetics'][enzyme]['parameters']["C"]]
                            
                    if "S" in self.parameters['reaction_kinetics'][enzyme]['parameters']:   
                        if self.parameters['reaction_kinetics'][enzyme]['parameters']["S"] in concentrations:    
                            S = concentrations[self.parameters['reaction_kinetics'][enzyme]['parameters']["S"]]

                    try:
                        fluxes = [eval(self.parameters['reaction_kinetics'][enzyme]["SubstitutedRateLaw"])]
                        if self.verbose:
                            print('enzyme', enzyme, 'rate', self.parameters['reaction_kinetics'][enzyme]["SubstitutedRateLaw"], 'fluxes', fluxes)
                        self.fluxes_df.at[enzyme, f'{self.timestep} (min)'] = average(fluxes)
                        self.parameters['calculated_rate_laws'][enzyme.lower()] = average(fluxes)
                    except:
                        pass
                else:
                    if self.verbose:
                        print('--> ERROR: The {} rate law has an unaccepted form.'.format(self.parameters['reaction_kinetics'][enzyme]["SubstitutedRateLaw"]))
                
            elif source == 'sabio':
                minimum = inf
                for reaction in self.parameters['reaction_kinetics'][enzyme]:           
                    for entry_id in self.parameters['reaction_kinetics'][enzyme][reaction]:
                        if "SubstitutedRateLaw" in self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]:
                            if self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]['SubstitutedRateLaw'] not in empty:
                                remainder = re.sub('([0-9ABCS\-\/\(\)\+\.\*])', '', self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["SubstitutedRateLaw"])
                                if remainder == '':
                                    A = 0
                                    B = 0
                                    C = 0


                                    if "A" in self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]:   
                                        if self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["A"]["species"] in concentrations: 
                                            A = concentrations[self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["A"]["species"]]

                                    if "B" in self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]:   
                                        if self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["B"]["species"] in concentrations:    
                                            B = concentrations[self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["B"]["species"]]

                                    if "C" in self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]:   
                                        if self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["C"]["species"] in concentrations:    
                                            C = concentrations[self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["C"]["species"]]

                                    if "S" in self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]:   
                                        if self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["S"]["species"] in concentrations:    
                                            S = concentrations[self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Parameters"]["S"]["species"]]

                                    try:
                                        flux = eval(self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["SubstitutedRateLaw"])
                                    except:
                                        flux = None
                                        pass

                                    if flux:

                                        # define the closest match of the data to the parameterized conditions
                                        if (temperature and self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Temperature"]) not in empty:
                                            temperature_deviation = abs(temperature - float(self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["Temperature"]))/temperature
                                        else:
                                            temperature_deviation = 0
                                        if (p_h and self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["pH"]) not in empty:
                                            ph_deviation = abs(p_h - float(self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["pH"]))/p_h
                                        else:
                                            ph_deviation = 0

                                        old_minimum = minimum
                                        deviation = average(temperature_deviation, ph_deviation)
                                        minimum = min(deviation, minimum)
            #                             print('minimum', minimum)
            #                             print('deviation', deviation)

                                        if old_minimum == minimum:
                                            fluxes.append(flux)
                                        elif deviation == minimum:
                                            fluxes = [flux]

                                        if self.verbose:
                                            print('enzyme:', enzyme, '\trate:', self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["SubstitutedRateLaw"], '\tfluxes:', fluxes)
                                        self.fluxes_df.at[enzyme, f'{self.timestep} (min)'] = average(fluxes)
                                else:
                                    if self.verbose:
                                        print('--> ERROR: The {} rate law has an unaccepted form.'.format(self.parameters['reaction_kinetics'][enzyme][reaction][entry_id]["SubstitutedRateLaw"]))
                
    def update_concentrations(self, solution):
        for reaction in self.model.reactions:
            print('upper_bound', reaction.upper_bound)
            print('lower_bound', reaction.lower_bound)
            print('flux', solution.fluxes[reaction.id])
            #print(f'{reaction.id} \tflux: {solution.fluxes[reaction.id]}') # \tExpression: {reaction.reaction}
            for metabolite in reaction.metabolites:
                delta_conc = reaction.metabolites[metabolite] * solution.fluxes[reaction.id]
#                 print('stoich', reaction.metabolites[metabolite])
                print('delta_conc', delta_conc)
                if self.timestep == 1:
                    initial_conc = 0
                    if metabolite.name in self.parameters['initial_concentrations']:
                        initial_conc = self.parameters['initial_concentrations'][metabolite.name]
                    self.variables['concentrations'][metabolite.name] = self.concentrations_df.at[f'{metabolite.name} ({metabolite.compartment})', f'{self.timestep} (min)'] = initial_conc + delta_conc
                else:
                    before = self.variables['concentrations'][metabolite.name]
                    self.variables['concentrations'][metabolite.name] += delta_conc
                    self.concentrations_df.at[f'{metabolite.name} ({metabolite.compartment})', f'{self.timestep} (min)'] = self.variables['concentrations'][metabolite.name]
                    after = self.variables['concentrations'][metabolite.name]
                    print('before', before)
                    print('after', after)
                    if self.verbose:
                        if before != after:
                            print(f'{metabolite.name} changed concentration, timestep {self.timestep}: {before}\t{after}')
                        
    def visualize(self):
        # print the reaction fluxes dataframe
        print('\n\nKinetically-constrained reaction fluxes\n', '='*50)
        print(self.fluxes_df.to_string())
        
        # print the metabolite concentrations dataframe
        print('\n\nMetabolite concentrations\n', '='*50)
        print(self.concentrations_df.to_string())