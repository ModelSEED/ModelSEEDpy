import logging
import itertools
import cobra
import re
import os
from numpy import zeros
from itertools import combinations
from optlang.symbolics import Zero
import networkx
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from to_precision import auto_notation
from matplotlib import pyplot

logger = logging.getLogger(__name__)

#Adding a few exception classes to handle different types of errors
class CommunityError(Exception):
    """Error in community FBA"""
    pass

class NoSolutionExists(Exception):
    """Error execution in gapfilling."""   
    def __str__(self):
        return 'The simulation lacks a solution, and cannot be gapfilled.'



class MSCommunity:

    def __init__(self,model):
        self.model = model
        self.gapfilling  = None
        self.constrained = False
        self.compartments = None
        self.production = None
        self.consumption = None
        self.solution = None
        self.suffix = ""
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)     
    
    def set_objective(self,target = "bio1",minimize = False):
        #Setting objective function
        sense = "max"
        if minimize:
            sense = "min"
        self.model.objective = self.model.problem.Objective(
            self.model.reactions.get_by_id(target).flux_expression, 
            direction=sense
        )
    
    def gapfill(self, media = None, target = "bio1", minimize = False, default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = []):
        self.set_objective(target,minimize)
        self.gapfilling = MSGapfill(self.model, default_gapfill_templates, default_gapfill_models, test_conditions, reaction_scores, blacklist)
        gfresults = self.gapfilling.run_gapfilling(media,target)
        if gfresults is None:
            print('\n--> ERROR: The simulation lacks a solution, and cannot be gapfilled.\n')
            #raise NoSolutionExists()
            return None
        else:
            return self.gapfilling.integrate_gapfill_solution(gfresults)
    
    def constrain(self,media = None,element_uptake_limit = None, kinetic_coeff = None, abundances = None, msdb_path_for_fullthermo = None, verbose = True):
        # applying media constraints
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        # applying uptake constraints
        element_contraint_name = ''
        if element_uptake_limit is not None:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
            element_contraint_name = 'eup'
        # applying kinetic constraints
        kinetic_contraint_name = ''
        if kinetic_coeff is not None:
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,abundances)
            kinetic_contraint_name = 'ckp'
        # applying FullThermo constraints
        thermo_contraint_name = ''
        if msdb_path_for_fullthermo is not None:
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_path':msdb_path_for_fullthermo}, verbose)
            thermo_contraint_name = 'ftp'
        self.suffix = '_'.join([self.model.id, thermo_contraint_name, kinetic_contraint_name, element_contraint_name])+".lp"
    
    def drain_fluxes(self, media = None, predict_abundances = True):
        biomass_drains = {}
        if predict_abundances:
            # parse the metabolites in the biomass reactions
            for reaction in self.model.reactions:
                if re.search('^bio\d+$', reaction.id):
                    for metabolite in reaction.metabolites:
                        # identify the biomass metabolites
                        msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
                        m = re.search('[a-z](\d+)', metabolite.compartment).group()
                        if msid == "cpd11416" and m != None:
                            # evaluate only cross-feeding
                            index = m[1]
                            if index != "0" and index not in biomass_drains:
                                print(f"Making biomass drain: {metabolite.id}")
                                drain_reaction = FBAHelper.add_drain_from_metabolite_id(self.model, metabolite.id,0,100,"DM_")
                                self.model.add_reactions([drain_reaction])
                                biomass_drains[index] = drain_reaction

        # print each optimal drain flux in the model
        for i in range(1,len(biomass_drains)+1):
            if f'DM_cpd11416_c{i}' not in self.model.reactions:
                print(f'--> ERROR: the reaction < DM_cpd11416_c{i} > is malformed.')
            FBAHelper.set_objective_from_target_reaction(self.model,f"DM_cpd11416_c{i}")
            self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            sol=self.model.optimize()
            print(f"species {i} drain-flux objective value: {sol.objective_value}")
            
            
        return biomass_drains
    
    def run(self,target = "bio1",minimize = False,pfba = True,print_lp = False,summary = False,objective_value=True):
        # conditionally print the LP file of the model
        if print_lp:
            count_iteration = 0
            file_name = re.sub('.lp', f'_{count_iteration}.lp', self.suffix)
            while os.path.exists(file_name):
                count_iteration += 1
                file_name = re.sub('(_\d).lp', f'_{count_iteration}.lp', file_name)
                print(count_iteration)
            with open(file_name, 'w') as out:
                out.write(str(self.model.solver))
                out.close()
                
        #Solving model
        self.set_objective(target,minimize)
        solution = self.model.optimize()
        if objective_value:
            print('\nModel objective value:', solution.objective_value)
        print('\n')
        if summary:
            print(self.model.summary())
        if pfba:
            solution = cobra.flux_analysis.pfba(self.model)
        #Checking that an optimal solution exists
        if solution.status != 'optimal':
            logger.warning("No solution found for the simulation.")
            solution = None
        return solution
        
    def compute_interactions(self,solution):
        # calculate the metabolic exchanges
        self.compartments = set()
        self.metabolite_uptake = {}
        for rxn in self.model.reactions:
            cmp = rxn.id.split("_").pop()
            comp_index = cmp[1:]
            if comp_index.isnumeric() and comp_index != '0':
                comp_index = int(comp_index)
                self.compartments.add(comp_index)
                for metabolite in rxn.metabolites:
                    if metabolite.compartment == "e0":
                        flux = solution.fluxes[rxn.id]
                        if flux != 0:
                            rate_law = rxn.metabolites[metabolite]*flux
                            self.metabolite_uptake[(metabolite.id, comp_index)] = rate_law
        # calculate cross feeding of extracellular metabolites 
        self.number_of_compartments = len(self.compartments)
        self.production = zeros((self.number_of_compartments, self.number_of_compartments)) 
        self.consumption = zeros((self.number_of_compartments, self.number_of_compartments))
        
        cross_all = []
        for rxn in self.model.reactions:
            for metabolite in rxn.metabolites:
                if metabolite.compartment == "e0":
                    # determine each directional flux rate 
                    rate_out = {compartment_number: rate for (metabolite_id, compartment_number), rate in self.metabolite_uptake.items() if metabolite_id == metabolite.id and rate > 0}
                    rate_in = {compartment_number: abs(rate) for (metabolite_id, compartment_number), rate in self.metabolite_uptake.items() if metabolite_id == metabolite.id and rate < 0}
                    # determine total directional flux rate 
                    total_in = sum(rate_in.values())
                    total_out = sum(rate_out.values())
                    max_total_rate = max(total_in, total_out)
                    # determine net flux 
                    net_flux = total_in - total_out
                    if net_flux > 0:
                        rate_out[None] = net_flux
                    if net_flux < 0:
                        rate_in[None] = abs(net_flux)
                    # establish the metabolites that partake in cross feeding 
                    for donor, rate_1 in rate_out.items():
                        if donor is not None:
                            donor_index = int(donor) - 1           
                            for receiver, rate_2 in rate_in.items():
                                if receiver is not None:
                                    receiver_index = int(receiver) - 1

                                    # assign calculated feeding rates to the production and consumption matrices
                                    rate = rate_1 * rate_2 / max_total_rate
                                    self.production[donor_index][receiver_index] += rate
                                    self.consumption[receiver_index][donor_index] += rate
        
    def visualize(self, graph = True, table = True): 
        # view the cross feeding matrices
        if table:
            from pandas import DataFrame as df

            species = [num+1 for num in range(len(self.production))]

            print('\nProduction matrix:')
            prod_df = df(self.production)
            prod_df.index = prod_df.columns = species
            prod_df.index.name = 'Donor'
            print(prod_df)

            print('\n\nConsumption matrix:')
            cons_df = df(self.consumption)
            cons_df.index = cons_df.columns = species
            cons_df.index.name = 'Receiver'
            print(cons_df)
            print('\n')
            
        # graph the community network
        if graph:
            graph = networkx.Graph()
            for num in self.compartments:
                graph.add_node(num)
            for com in combinations(self.compartments, 2):
                species_1 = int(com[0])-1
                species_2 = int(com[1])-1

                interaction_net_flux = auto_notation((self.production[species_1][species_2] - self.consumption[species_1][species_2]), 3)
                if species_1 < species_2:
                    graph.add_edge(com[0],com[1],flux = interaction_net_flux)
                elif species_1 > species_2:
                    graph.add_edge(com[0],com[1],flux = -interaction_net_flux)

            pos = networkx.circular_layout(graph)
            networkx.draw_networkx(graph,pos)
            labels = networkx.get_edge_attributes(graph,'flux')
            networkx.draw_networkx_edge_labels(graph,pos,edge_labels=labels)
            pyplot.show()
            return graph, pos, labels

    @staticmethod
    def community_fba(model,media = None,target = "bio1",minimize = False,pfba = True,print_lp = False,summary = False,element_uptake_limit = None, kinetic_coeff = None, abundances = None, msdb_path_for_fullthermo = None,default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = []):
        cfba = MSCommunity(model)
        cfba.constrain(media,element_uptake_limit, kinetic_coeff, abundances, msdb_path_for_fullthermo)
        cfba.gapfill(media,target,minimize,default_gapfill_templates,default_gapfill_models,test_conditions,reaction_scores,blacklist)
        self.solution = cfba.run(target,minimize,pfba,print_lp,summary)
        cfba.compute_interactions(self.solution)
        cfba.visualize()
        return community
