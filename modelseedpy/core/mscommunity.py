import logging
import itertools
import cobra
import re
import os
from numpy import zeros
from itertools import combinations
from optlang.symbolics import Zero, add
import networkx
from cobra.core.dictlist import DictList
from cobra.core import Reaction
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.fbapkg.gapfillingpkg import default_blacklist

logger = logging.getLogger(__name__)

#Adding a few exception classes to handle different types of errors
class CommunityError(Exception):
    """Error in community FBA"""
    pass

class NoSolutionExists(Exception):
    """Error execution in gapfilling."""   
    def __str__(self):
        return 'The simulation lacks a solution, and cannot be gapfilled.'

class CommunityModelSpecies:
    def __init__(self,community,biocpd,names=[]):
        self.community = community
        self.biomass_cpd = biocpd
        self.index = int(self.biomass_cpd.compartment[1:])
        self.abundance = 0
        if biocpd in community.primary_biomass.metabolites:
            self.abundance = abs(community.primary_biomass.metabolites[biocpd])
        self.biomass_drain = None
        if self.index <= len(names) and names[self.index-1] != None:
            self.id = names[self.index-1]
        else:
            if "species_name" in self.biomass_cpd.annotation:
                self.id = self.biomass_cpd.annotation["species_name"]
            else:
                self.id = "Species"+str(self.index)
        logger.info("Making atp hydrolysis reaction for species: "+self.id)
        output = FBAHelper.add_atp_hydrolysis(community.model,"c"+str(self.index))
        self.atp_hydrolysis = output["reaction"]
        self.biomasses = []
        for reaction in self.community.model.reactions:
            if biocpd in reaction.metabolites:
                if reaction.metabolites[biocpd] == 1 and len(reaction.metabolites) > 1:
                    self.biomasses.append(reaction)
                elif len(reaction.metabolites) == 1 and reaction.metabolites[biocpd] < 0:
                    self.biomass_drain = reaction
        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species "+self.id)
        if self.biomass_drain == None:
            logger.info("Making biomass drain reaction for species: "+self.id)
            self.biomass_drain = Reaction(id="DM_"+biocpd.id,
                                      name="DM_" + biocpd.name,
                                      lower_bound=0, 
                                      upper_bound=100)
            community.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({biocpd : -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'
    
    def disable_species(self):
        for reaction in self.community.model.reactions:
            if int(FBAHelper.rxn_compartment(reaction)[1:]) == self.index:
                reaction.upper_bound = 0
                reaction.lower_bound = 0
    
    def compute_max_biomass(self):
        if len(self.biomasses) == 0:
            logger.critical("No biomass found for species:"+self.id)
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.biomasses[0].forward_variable:1})
        return self.community.model.optimize()
    
    def compute_max_atp(self):
        if self.atp_hydrolysis == None:
            logger.critical("No ATP hydrolysis found for species:"+self.id)
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
        return self.community.model.optimize()

#class MSModelAdapter
    
#    def __init__(self,model):
 #       self.model = model
    
        
        

class MSCommunity:

    def __init__(self,model,names=[],abundances=None):
        self.model = MSModelAdapter(model)
        self.biomass_cpd = None
        self.primary_biomass = None
        self.biomass_drain = None
        self.gapfilling  = None
        self.constrained = False
        self.compartments = None
        self.production = None
        self.consumption = None
        self.solution = None
        self.suffix = ""
        self.verbose = True
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.species = DictList()
        id_hash = FBAHelper.msid_hash(model)
        if "cpd11416" not in id_hash:
            logger.critical("Could not find biomass compound")
        other_biocpds = []
        for biocpd in id_hash["cpd11416"]:
            if biocpd.compartment == "c0":
                self.biomass_cpd = biocpd
                for reaction in model.reactions:
                    if self.biomass_cpd in reaction.metabolites:
                        if reaction.metabolites[self.biomass_cpd] == 1 and len(reaction.metabolites) > 1:
                            self.primary_biomass = reaction
                        elif len(reaction.metabolites) == 1 and reaction.metabolites[self.biomass_cpd] < 0:
                            self.biomass_drain = reaction
            else:
                other_biocpds.append(biocpd)
        for biocpd in other_biocpds:        
            species_obj = CommunityModelSpecies(self,biocpd,names)
            self.species.append(species_obj)
        if abundances != None:
            self.set_abundance(abundances)
    
    def set_abundance(self,abundances):
        #First ensure normalization
        totalabundance = 0
        for species in abundances:
            totalabundance += abundances[species]
        #Next map abundances to all species
        for species in abundances:
            abundances[species] = abundances[species]/totalabundance
            if species in self.species:
                self.species.get_by_id(species).abundance = abundances[species]
        #Finally, remake the primary biomass reaction based on abundances
        if self.primary_biomass == None:
            logger.critical("Primary biomass reaction not found in community model")
        all_metabolites = {self.biomass_cpd:1}
        for species in self.species:
            all_metabolites[species.biomass_cpd] = -1*abundances[species.id]
        self.primary_biomass.add_metabolites(all_metabolites,combine=False)
    
    def test_individual_species(self,media = None,allow_interaction=True,run_atp=True,run_biomass=True):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        #Iterating over species and running tests
        output = {}
        for individual in self.species:
            output[individual.id] = {}
            #Running with model so changes are discarded after each simulation
            with self.model:
                #If no interaction allowed, iterate over all other species and disable them
                if not allow_interaction:
                    for indtwo in self.species:
                        if indtwo != individual:
                            indtwo.disable_species()
                #If testing biomass, setting objective to individual species biomass and optimizing
                if run_biomass:
                    output[individual.id]["biomass"] = individual.compute_max_biomass()
                #If testing atp, setting objective to individual species atp and optimizing
                if run_atp:
                    output[individual.id]["atp"] = individual.compute_max_atp()
        if self.verbose == True:
            line = "Species"
            if run_biomass:
                line += "\tBiomass"
            if run_atp:
                line += "\tATP"   
            print("\n"+line) 
            for individual in self.species:
                line  = individual.id
                if run_biomass:
                    line += "\t"+str(output[individual.id]["biomass"].objective_value)
                if run_atp:
                    line += "\t"+str(output[individual.id]["atp"].objective_value)    
                print(line)
            print("")
        return output
        #Convert table to pandas?
        
    def set_objective(self,target = None,minimize = False):
        #Setting objective function
        if target == None:
            target = self.primary_biomass.id
        sense = "max"
        if minimize:
            sense = "min"
        self.model.objective = self.model.problem.Objective(
            self.model.reactions.get_by_id(target).flux_expression, 
            direction=sense
        )
    
    def atp_correction(self,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
        atpcorrect = MSATPCorrection(self.model,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0)
        
    
    def gapfill(self, media = None, target = None, minimize = False,default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = []):
        if target == None:
            target = self.primary_biomass.id
        self.set_objective(target,minimize)
        gfname = FBAHelper.medianame(media)+"-"+target+"-"+self.suffix
        self.gapfillings[gfname] = MSGapfill(self.model, default_gapfill_templates, default_gapfill_models, test_conditions, reaction_scores, blacklist) 
        gfresults = self.gapfillings[gfname].run_gapfilling(media,target)
        if gfresults is None:
            logger.critical("Gapfilling failed with the specified model, media, and target reaction.")
        return self.gapfillings[gfname].integrate_gapfill_solution(gfresults)
    
    def constrain(self,element_uptake_limit = None, kinetic_coeff = None,msdb_path_for_fullthermo = None, verbose = True):
        # applying uptake constraints
        element_contraint_name = ''
        if element_uptake_limit is not None:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
            element_contraint_name = 'eup'
        # applying kinetic constraints
        kinetic_contraint_name = ''
        if kinetic_coeff is not None:
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
            kinetic_contraint_name = 'ckp'
        # applying FullThermo constraints
        thermo_contraint_name = ''
        if msdb_path_for_fullthermo is not None:
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_path':msdb_path_for_fullthermo}, verbose)
            thermo_contraint_name = 'ftp'
        self.suffix = '_'.join([self.model.id, thermo_contraint_name, kinetic_contraint_name, element_contraint_name])+".lp"
    
    def predict_abundances(self,media=None,kinetic_coeff = 2000,print_lp=False):
        #Using with to ensure the changes to the model are not retained
        output = {"solution":None,"abundances":{}}
        with self.model:
            self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
            objcoef = {}
            for species in self.species:
                objcoef[species.biomasses[0].forward_variable] = 1
            new_objective = self.model.problem.Objective(Zero,direction="max")
            self.model.objective = new_objective
            new_objective.set_linear_coefficients(objcoef)
            if not print_lp:
                self.print_lp(self.suffix+"-"+FBAHelper.medianame(media)+"-Coef"+str(kinetic_coeff)+"-PredAbund")
            output["solution"] = self.model.optimize()
            totalgrowth = 0
            for species in self.species:
                output["abundances"][species.id] = output["solution"].fluxes[species.biomasses[0].id]
                totalgrowth += output["abundances"][species.id]
            if totalgrowth > 0:
                for species in self.species:
                    output["abundances"][species.id] = output["abundances"][species.id]/totalgrowth
            if self.verbose == True:
                print("\nSpecies\tAbundance")
                for species in self.species:
                    print(species.id+"\t"+str(output["abundances"][species.id]))
                print("")
        return output
    
    def run(self,media,target = None,minimize = False,pfba = True,print_lp = False):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        # conditionally print the LP file of the model
        if print_lp:
            self.print_lp()
        #Solving model
        self.set_objective(target,minimize)
        solution = None
        if pfba:
            solution = cobra.flux_analysis.pfba(self.model)
        else:
            solution = self.model.optimize()
        #Checking that an optimal solution exists
        if solution.status != 'optimal':
            logger.warning("No solution found for the simulation.")
            return None
        if self.verbose == True:
            print('\nModel objective value:', solution.objective_value,'\n')
            print(self.model.summary())
            print("")
        return solution
    
    def print_lp(self,suffix=None):
        if suffix == None:
            suffix = self.suffix
        count_iteration = 0
        filename = "Community-"+suffix+"-"+str(count_iteration)+".lp"
        while os.path.exists(filename):
            count_iteration += 1
            filename = "Community-"+suffix+"-"+str(count_iteration)+".lp"
        with open(filename, 'w') as out:
            out.write(str(self.model.solver))
            out.close()
        
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
        ''' VISUALIZE FLUXES ''' 
        # graph the community network
        if graph:
            graph = networkx.Graph()
            for num in self.compartments:
                graph.add_node(num)
            for com in combinations(self.compartments, 2):
                species_1 = int(com[0])-1
                species_2 = int(com[1])-1

                interaction_net_flux = round(self.production[species_1][species_2] - self.consumption[species_1][species_2])
                if species_1 < species_2:
                    graph.add_edge(com[0],com[1],flux = interaction_net_flux)
                elif species_1 > species_2:
                    graph.add_edge(com[0],com[1],flux = -interaction_net_flux)

            pos = networkx.circular_layout(graph)
            networkx.draw_networkx(graph,pos)
            labels = networkx.get_edge_attributes(graph,'flux')
            networkx.draw_networkx_edge_labels(graph,pos,edge_labels=labels)
        
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

    @staticmethod
    def community_fba(model,media = None,target = None,minimize = False,pfba = True,print_lp = False,summary = False,element_uptake_limit = None, kinetic_coeff = None,names=[],abundances = None, msdb_path_for_fullthermo = None,default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = []):
        cfba = MSCommunity(model,names,abundances)
        cfba.constrain(media,element_uptake_limit, kinetic_coeff, abundances, msdb_path_for_fullthermo)
        cfba.gapfill(media,target,minimize,default_gapfill_templates,default_gapfill_models,test_conditions,reaction_scores,blacklist)
        self.solution = cfba.run(target,minimize,pfba,print_lp,summary)
        cfba.compute_interactions(self.solution)
        cfba.visualize()
        return community
