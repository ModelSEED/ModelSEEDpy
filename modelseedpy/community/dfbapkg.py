# -*- coding: utf-8 -*-

from scipy.constants import milli, hour, minute, day, femto
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
# from modelseedpy.core.fbahelper import FBAHelper
from collections import OrderedDict
from optlang.symbolics import Zero
from numpy import log10, nan, mean
from warnings import warn
from matplotlib import pyplot
from pprint import pprint
from datetime import date
from math import inf 
import pandas
# import cython
import json, re, os

   
def isnumber(string):
    try:
        float(string)
        remainder = re.sub('([0-9.\-eE])', '', str(string))
    except:
        try:
            int(string)
            remainder = re.sub('[0-9.-eE])', '', str(string))
        except:
            return False
    if remainder == '':
        return True

class dFBAPkg(BaseFBAPkg):
    def __init__(self, 
                 model,                   # Cobrakbase model
                 modelseed_db_path: str,  # local path to the ModelSEED Database
                 solver: str = 'glpk',    # specifies the LP solver
                 warnings: bool = True, verbose: bool = False, printing: bool = False, jupyter: bool = False
                 ):
        # define the parameter and variable dictionaries
        BaseFBAPkg.__init__(self, model, "BasedFBA", {"met": "metabolite"}, {"conc": 'metabolite'})
        self.pkgmgr.addpkgs(["FullThermoPkg"])
        # self.parameters["modelseed_api"] = FBAHelper.get_modelseed_db_api(self.parameters["modelseed_db_path"])
            
        # define simulation conditions
        self.warnings: bool = warnings; self.verbose: bool = verbose; self.printing: bool = printing; self.jupyter: bool = jupyter 
        self.model = model
        self.model.solver = solver

        # define a list of metabolite ids
        self.met_ids, self.met_names = OrderedDict(), []
        for met in self.model.metabolites:
            self.met_ids[met.id] = met.name
            self.met_names.append(met.name)
            
    # user functions
    def print_lp(self, filename=None):
        BaseFBAPkg.print_lp(self, filename.removesuffix('.lp'))
                            
    def simulate(self, 
                 kinetics_path: str = None,                             # the path of the kinetics data JSON file
                 initial_concentrations_M: dict = {},                   # a dictionary of the initial metabolic concentrations , which supplants concentrations from the defined kinetics data
                 total_time: float = 200, timestep: float = 20,  # total simulation time and the simulation timestep in minutes
                 export_name: str = None, export_directory: str = None, # the location to which simulation content will be exported
                 chemostat_L: float = None, feed_profile: dict = {},    # the volume (l) and feed profile for a chemostat simulation, where None ignores a chemostat
                 exchange_rate: float = None,                    # the flow rate (Molar/Liter) of addition to and removal from the chemostat system 
                 thermo_constraints: bool = False,                      # specifies whether thermodynamic constraints will be layered with the kinetic constraints
                 kinetics_data: dict = {},                              # A dictionary of custom kinetics data
                 temperature: float = 25, p_h: float = 7, # simulation conditions 
                 cellular_dry_mass_fg: float = 222,              # mass of the simulated cell in femtograms
                 cellular_fL: float = 1,                         # volume of the simulated cell in femtoliters
                 figure_title: str = 'Metabolic perturbation',          # title of the concentrations figure
                 included_metabolites: list = [],                       # A list of the metabolites that will be graphically displayed
                 labeled_plots: bool = True,                            # specifies whether plots will be individually labeled 
                 visualize: bool = True, export: bool = True            # specifies whether simulation content will be visualized or exported, respectively
                 ):
        # define the dataframe for the time series content
        self.cellular_dry_mass_fg: float = cellular_dry_mass_fg*femto; self.cellular_fL: float = cellular_fL*femto
        self.parameters['timesteps'] = int(total_time/timestep)    
        self.timestep_value: float = timestep; self.total_time: float = total_time
        self.constrained = OrderedDict()
        self.solutions = []
        self.minimum = inf
        
        # define experimental conditions
        self.parameters['pH'], self.parameters['temperature'] = p_h, temperature
        self.variables['elapsed_time'] = 0
        
        # define initial concentrations
        self.initial_concentrations_M = initial_concentrations_M
        self._initial_concentrations(kinetics_path, kinetics_data)
                
        # apply constraints and system specifications
        chemostat_requirements = [isnumber(chemostat_L), feed_profile != {}, isnumber(exchange_rate)]
        if any(chemostat_requirements) and not all(chemostat_requirements):
            warn(f'The chemostat_L ({chemostat_L}), feed_profile ({feed_profile}), and exchange_rate ({exchange_rate}) parameters must all be defined to simulate a chemostat.')
        if thermo_constraints:
            self.pkgmgr.addpkgs(["FullThermoPkg"])
            self.pkgmgr.getpkg("FullThermoPkg").build_package()
            
        # determine the reactions for which kinetics are predefined
        self.defined_reactions = {}
        for rxn in self.model.reactions:
            if rxn.name in kinetics_data:
                self.defined_reactions[rxn.name] = rxn

        # execute FBA for each timestep
        for self.timestep in range(1,self.parameters['timesteps']+1):
            # calculate custom fluxes, constrain the model, and update concentrations
            self._define_timestep()
            self._build_constraints()
            self._calculate_kinetics()                    
            self._execute_cobra()
            self._update_concentrations()
            if all(chemostat_requirements):
                self.chemical_moles[self.col] = self.concentrations[self.col]*milli * chemostat_L
                self._chemostat(feed_profile, exchange_rate, chemostat_L)
            
            self.variables['elapsed_time'] += self.timestep
            if self.printing:
                print(f'\nobjective value for timestep {self.timestep}: ', self.solutions[-1].objective_value)                
        
        # identify the chemicals that dynamically changed in concentrations
        self.changed, self.unchanged = set(), set()
        for met_name in self.met_names:
            first = self.concentrations.at[met_name, '0 min']
            final = self.concentrations.at[met_name, self.col]
            if first != final:
                self.changed.add(met_name)
            if first == final:
                self.unchanged.add(met_name)
                
        # visualize concentration changes over time
        if visualize:
            self._visualize(figure_title,included_metabolites,labeled_plots)
        if export:
            self._export(export_name, export_directory)
        
        # view calculations and results
        if self.verbose:
            print('\n\n', 'Changed metabolite  concentrations\n', '='*2*len('changed metabolites'), '\n', self.changed)
            print('\nConstrained reactions:', self.constrained.keys())     
        elif self.printing:
            if self.jupyter:
                pandas.set_option('max_rows', None)
                display(self.concentrations)
                display(self.fluxes)
            if self.unchanged == set():
                print('\nAll of the metabolites changed concentration over the simulation')
            else:
                print('\n\nUnchanged metabolite concentrations', '\n', '='*2*len('unchanged metabolites'), '\n', self.unchanged)       
                
        return self.concentrations, self.fluxes
            
    #utility functions
    def _initial_concentrations(self,
                                kinetics_path: str = None,  # the absolute path to a JSON file of kinetics data 
                                kinetics_data: dict = {},   # a dictionary of kinetics data, which supplants imported data from the kinetics_path
                                ):
        # define kinetics of the system
        self.kinetics_data = {}
        if kinetics_path:
            if not os.path.exists(kinetics_path):
                raise ValueError('The path {kinetics_data} is not a valid path')
            with open(kinetics_path) as data:
                self.kinetics_data = json.load(data)
        if kinetics_data != {}:
            for reaction in kinetics_data:
                self.kinetics_data[reaction] = kinetics_data[reaction]
        if self.kinetics_data == {}:
            raise ValueError('Kinetics data must be defined.')
                
        # define the DataFrames
        self.col = '0 min'     
        self.concentrations = pandas.DataFrame(index=set(self.met_names), columns=[self.col])
        self.chemical_moles = pandas.DataFrame(index=set(self.met_names), columns=[self.col])
        self.concentrations.index.name = self.chemical_moles.index.name = 'metabolite (\u0394mM)'
        
        self.fluxes = pandas.DataFrame(index = set(rxn.name for rxn in self.model.reactions), columns = [self.col])
        self.fluxes.index.name = 'reactions (mmol/g_(dw)/hr)'
        
        # parse the kinetics data 
        initial_concentrations = {}
        for met in self.met_names:
            self.concentrations.at[str(met), self.col] = float(0)
        for reaction_name in self.kinetics_data:
            for condition, datum in self.kinetics_data[reaction_name].items():
                for var in datum['initial_concentrations_M']:
                    met_id = datum['met_id'][var]
                    if met_id in self.met_ids:
                        name = self.met_ids[met_id]
                        if name in self.met_names:                                    
                            self.concentrations.at[name, self.col] += datum['initial_concentrations_M'][var]/milli
                            initial_concentrations[met_id] = self.concentrations.at[name, self.col]
                        else:
                            if self.warnings:
                                warn(f"KineticsError: The {name} reagent ({var}) in the {datum['substituted_rate_law']} rate law is not defined by the model.")
                    else:
                        if self.warnings:
                            warn(f"KineticsError: The {name} reagent ({var}) in the {datum['substituted_rate_law']} rate law is not recognized by the ModelSEED Database.")
                        
        # incorporate custom initial concentrations
        if isinstance(self.initial_concentrations_M, dict) and self.initial_concentrations_M != {}:
            for met_id in self.initial_concentrations_M:
                met_name = self.met_ids[met_id]
                if met_name not in self.concentrations.index:
                    if self.warnings:
                        warn(f'InitialConcError: The {met_id} ({met_name}) metabolite is not defined by the model.')
                else:
                    self.concentrations.at[met_name, self.col] = self.initial_concentrations_M[met_id]*milli
                    initial_concentrations[met_id] = self.concentrations.at[name, self.col]
        self.initial_concentrations_M = initial_concentrations
                    

    def _define_timestep(self,):
        self.col = f'{self.timestep*self.timestep_value} min'
        self.previous_col = f'{(self.timestep-1)*self.timestep_value} min'
        self.concentrations[self.col] = [float(0) for ind in self.concentrations.index]  #!!! 
        self.fluxes[self.col] = [nan for ind in self.fluxes.index]
        
    
    def _calculate_kinetics(self):        
        for reaction_name in self.kinetics_data:
            fluxes = []
            for source in self.kinetics_data[reaction_name]: 
                incalculable = False
                datum = self.kinetics_data[reaction_name][source]
                if "substituted_rate_law" in datum:     #!!! Statistics of aggregating each condition should be provided for provenance.
                    remainder = re.sub('([0-9A-Za-z/()e\-\+\.\*])', '', datum["substituted_rate_law"])
                    if remainder == '':
                        # define each variable concentration
                        conc_dict = {}
                        for var in datum['met_id']:
                            met_id = datum['met_id'][var]
                            if len(var) == 1:
                                conc_dict[var] = self.concentrations.at[self.met_ids[met_id], self.previous_col]*milli    # concentrations are mM
                                # warn(f'MetaboliteError: The {self.met_ids[met_id]} chemical is not recognized by the ModelSEED Database.')

                        if conc_dict != {}:
                            locals().update(conc_dict)
                            flux = eval(datum["substituted_rate_law"])
                            
                            # average or overwrite flux calculations based upon the alignment of the data conditions with the simulation conditions
                            add_or_write = 'a'
                            if 'metadata' in self.kinetics_data[reaction_name][source]:
                                add_or_write = self.__find_data_match(reaction_name, source)
                            if add_or_write == 'a':
                                fluxes.append(flux) 
                            elif add_or_write == 'w':
                                fluxes = [flux]
                        else:
                            if self.warnings:
                                warn(f'MetaboliteError: The {reaction_name} reaction possesses unpredictable chemicals.')
                    else:
                         if self.warnings:
                             warn('RateLawError: The {datum["substituted_rate_law"]} rate law contains unknown characters: {remainder}')
                else:    
                    if self.warnings:
                        warn(f'RateLawError: The {datum} datum lacks a rate law.')
                        
            flux = mean(fluxes)
            if isnumber(flux):
                if reaction_name in self.defined_reactions:
                    self.__set_constraints(reaction_name, flux)
                    self.fluxes.at[reaction_name, self.col] = flux 
                    if self.printing:
                        print('\n')
                else:
                    if self.warnings:
                        warn(f'ReactionError: The {reaction_name} reaction, with a flux of {flux}, is not described by the model.')
            else:
                 if self.warnings:
                     warn(f'FluxError: The {reaction_name} reaction flux {datum["substituted_rate_law"]} value {flux} is not numberic.')
                
    def _execute_cobra(self):
        # execute the COBRA model 
        solution = self.model.optimize()
        self.solutions.append(solution)
        for rxn in self.model.reactions:
            if not isnumber(self.fluxes.at[rxn.name, self.col]):
                self.fluxes.at[rxn.name, self.col] = solution.fluxes[rxn.id]
            
    def _update_concentrations(self,):
        for met in self.model.metabolites:      
            self.concentrations.at[self.met_ids[met.id], self.col] = 0
            for rxn in met.reactions:  # flux units: mmol/(g_(dry weight)*hour)
                stoich = rxn.metabolites[met]
                flux = self.fluxes.at[rxn.name, self.col] 
                delta_conc = stoich * (flux * self.timestep_value*(minute/hour) * self.cellular_dry_mass_fg/self.cellular_fL) 
                self.concentrations.at[self.met_ids[met.id], self.col] += delta_conc
                    
    def _visualize(self,
                   figure_title,          # defines the title of the concentrations figure
                   included_metabolites,  # specifies which metabolites will be included in the figure
                   labeled_plots          # specifies which plots will be labeled in the figure
                   ):
        # define the figure
        pyplot.rcParams['figure.figsize'] = (11, 7)
        pyplot.rcParams['figure.dpi'] = 150
        
        self.figure, ax = pyplot.subplots()
        ax.set_title(figure_title)
        ax.set_ylabel('Concentrations (mM)') 
        
        x_axis_scalar, unit = self.__x_axis_determination()
        ax.set_xlabel('Time '+unit)
        legend_list, times = [], [t*self.timestep_value*x_axis_scalar for t in range(self.parameters['timesteps']+1)]
        
        # determine the plotted metabolites and the scale of the figure axis
        bbox = (1,1)
        if included_metabolites == []:
            bbox = (1.7,1)
            for chem in self.changed:
                if max(self.concentrations.loc[[chem]].values[0].tolist()) > 1e-2:   # an arbitrary concentration threshold for plotting on the figure
                    included_metabolites.append(chem)
        
        log_axis = False
        minimum, maximum = inf, -inf
        printed_concentrations = {}
        for chem in self.changed:
            if chem in included_metabolites:
                concentrations = self.concentrations.loc[[chem]].values[0].tolist()  # molar
                
                # determine the concentration range
                max_conc = max([x if x > 1e-9 else 0 for x in concentrations])
                maximum = max(maximum, max_conc)
                min_conc = min([x if x > 1e-9 else 0 for x in concentrations])
                minimum = min(minimum, min_conc)
                
                # plot chemicals with perturbed concentrations
                ax.plot(times, concentrations)
                if len(chem) > 25:
                    chem = list(self.met_ids.keys())[self.met_names.index(chem)]
                if not concentrations[0] < 1e-9:
                    legend_list.append(chem)
                else:
                    legend_list.append(f'(rel) {chem}')
                    
                # design the proper location of the overlaid labels in the figure
                if labeled_plots:
                    for i, conc in enumerate(concentrations):
                        if conc > 1e-9:
                            x_value = i*self.timestep_value
                            vertical_adjustment = 0
                            if x_value in printed_concentrations:
                                vertical_adjustment = (maximum - minimum)*.05
                                if log_axis:
                                    vertical_adjustment = log10(maximum - minimum)/3
                            ax.text(x_value, conc+vertical_adjustment, f"{chem} - {round(conc, 4)}", ha="left")
                            printed_concentrations[x_value] = conc
                            break

        # finalize figure details
        if maximum > 10*minimum:
            log_axis = True
            ax.set_yscale('log')
        ax.set_xticks(times)
        ax.grid(True)
        ax.legend(legend_list, title = 'Changed chemicals', loc='upper right', bbox_to_anchor = bbox, title_fontsize = 'x-large', fontsize = 'large') 
        
        
    def _export(self, 
                export_name: str,      # the folder name to which the simulation content will be exported
                export_directory: str  # the directory within which the simulation folder will be created
                ):
        # define a unique simulation name 
        directory = os.getcwd()
        if export_directory is not None:
            directory = os.path.dirname(export_directory)
        if export_name is None:
            export_name = '-'.join([re.sub(' ', '_', str(x)) for x in [date.today(), 'dFBA', self.model.name, f'{self.total_time} min']])
            
        simulation_number = -1
        while os.path.exists(os.path.join(directory, export_name)):
            simulation_number += 1
            export_name = re.sub('(\-\d+$)', '', export_name)
            export_name = '-'.join([export_name, str(simulation_number)])
            
        self.parameters['simulation_path'] = self.simulation_path = os.path.join(directory, export_name)
        os.mkdir(self.simulation_path)
        
        # export simulation content
        self.fluxes.to_csv(os.path.join(self.simulation_path, 'fluxes.csv'))
        self.concentrations.to_csv(os.path.join(self.simulation_path, 'concentrations.csv'))
        
        times = self.fluxes.columns
        with open(os.path.join(self.simulation_path, 'objective_values.csv'), 'w') as obj_val:   
            obj_val.write('min,objective_value') 
            for sol in self.solutions:
                index = self.solutions.index(sol)
                time = re.sub('(\smin)', '', times[index])
                obj_val.write(f'\n{time},{sol.objective_value}')              
        
        # export the parameters
        parameters = {'parameter':[], 'value':[]}
        for parameter in self.parameters:
            parameters['parameter'].append(parameter)
            parameters['value'].append(self.parameters[parameter])
            
        parameters_table = pandas.DataFrame(parameters)
        parameters_table.to_csv(os.path.join(self.simulation_path, 'parameters.csv'))
        
        # export the figure
        self.figure.savefig(os.path.join(self.simulation_path, 'changed_concentrations.svg'))
        if self.verbose:
            if not self.jupyter:
                self.figure.show()    
                                
    def _build_constraints(self):
        # create a metabolite variable that prevents negative concentrations
        timestep_hr = self.timestep_value * (minute/hour)
        for met in self.model.metabolites: 
            if met.id in self.initial_concentrations_M:
                coef = {}
                for rxn in met.reactions:
                    stoich = timestep_hr*rxn.metabolites[met]  # The product of the reaction stoichiometry and the timestep, for the integration of the steady-state
                    coef[rxn.forward_variable] = stoich
                    coef[rxn.reverse_variable] = -stoich
                    if stoich < 0:
                        coef[rxn.forward_variable] = -stoich
                        coef[rxn.reverse_variable] = stoich
                    
                # build the metabolite constraint
                if met.id in self.constraints["conc"]:
                    self.model.remove_cons_vars(self.constraints["conc"][met.id])
                self.constraints["conc"][met.id] = self.model.problem.Constraint(
                    Zero, lb=-self.concentrations.at[met.name, self.previous_col], ub=None, name=f'{met.id}_conc'
                    )
                self.model.add_cons_vars(self.constraints["conc"][met.id])
                self.model.solver.update()
                
                self.constraints["conc"][met.id].set_linear_coefficients(coef)
                self.model.solver.update()
                
                # var = BaseFBAPkg.build_variable(self,"met",0,None,"continuous",met)  
                # BaseFBAPkg.build_constraint(self,"conc",0,None,{met:1},met)
                
    def _chemostat(self, feed_profile, # a dictionary of the chemicals and their concentrations in the influent feed
                  exchange_rate,       # the L/hr flow rate of the feed and extraction of the chemostat
                  chemostat_L          # the volume (l) of the chemostat 
                  ):
        L_changed = exchange_rate*self.timestep_value
        # chemostat addition
        for chem_id, conc in feed_profile.items():
            chem_name = self.met_ids[chem_id]
            self.chemical_moles.at[chem_name, self.col] += conc*L_changed
            self.concentrations.at[chem_name, self.col] = (self.chemical_moles.at[chem_name, self.col]/milli/chemostat_L)    # normalize to the chemostat volume
            
        # chemostat subtraction
        for met in self.model.metabolites:
            if met.compartment[0] == 'e':
                self.chemical_moles.at[met.name, self.col] -= (self.concentrations.at[met.name, self.col]*L_changed)  
                self.concentrations.at[met.name, self.col] = (self.chemical_moles.at[met.name, self.col]/milli/chemostat_L)  # normalize to the chemostat volume
        
        
    # nested functions
    def __find_data_match(self,
                          reaction_name: str, # specifies the name of the given reaction
                          source: str         # specifies which datum of the enzymatic data will be used, where multiple data entries are present
                          ):
        # identifies the datum whose experimental conditions most closely matches the simulation conditions
        temperature_deviation = ph_deviation = 0
        if isnumber(self.kinetics_data[reaction_name][source]['metadata']["Temperature"]):
            temperature_deviation = abs(self.parameters['temperature'] - float(self.kinetics_data[reaction_name][source]['metadata']["Temperature"]))/self.parameters['temperature']
        if isnumber(self.kinetics_data[reaction_name][source]['metadata']["pH"]):
            ph_deviation = abs(self.parameters['pH'] - float(self.kinetics_data[reaction_name][source]['metadata']["pH"]))/self.parameters['pH']

        # equally weight between temperature and pH deviation from the simulation conditions
        old_minimum = self.minimum
        deviation = mean(temperature_deviation, ph_deviation)
        self.minimum = min(deviation, self.minimum)

        if old_minimum == self.minimum:
            return 'a'                  # append to an existing list of data
        elif deviation == self.minimum:
            return 'w'                  # construct a new list of data
        
    def __set_constraints(self, 
                          reaction_name: str, flux: float   # specify the name and flux of the given reaction, respectively
                          ):           
        rxn = self.defined_reactions[reaction_name]
        rxn_name = re.sub(' ', '_', rxn.name) 
        if rxn_name in self.constrained:
            self.model.remove_cons_vars(self.constrained[rxn_name])
        self.model.solver.update()
        self.constrained[rxn_name] = self.model.problem.Constraint(rxn.flux_expression, lb=flux, ub=flux, name=f'{rxn_name}_kinetics')
        self.model.add_cons_vars(self.constrained[rxn_name])
        self.model.solver.update()
        if self.verbose:
            print(self.model.constraints[f'{rxn_name}_kinetics'])
            
        
    def __x_axis_determination(self,):
        scalar = minute
        time = self.total_time*scalar
        unit = 's'
        if time > 600:
            unit = 'min'
            scalar = 1
            if time > 7200:
                unit = 'hr'
                scalar = 1/hour
                if time > 2e5:
                    scalar = 1/day
                    unit = 'days'
        return scalar, unit