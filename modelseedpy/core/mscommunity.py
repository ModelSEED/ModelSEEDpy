import logging
import itertools
import cobra
import re
import os
from numpy import zeros
from itertools import combinations
from optlang.symbolics import Zero, add
import networkx
import pandas as pd
from cobra.core.dictlist import DictList
from cobra.core import Reaction
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.fbapkg.gapfillingpkg import default_blacklist

logger = logging.getLogger(__name__)

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
        FBAHelper.add_autodrain_reactions_to_community_model(community.model)
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
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_Biomass")
        return self.community.model.optimize()
    
    def compute_max_atp(self):
        if self.atp_hydrolysis == None:
            logger.critical("No ATP hydrolysis found for species:"+self.id)
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_ATP")
        return self.community.model.optimize()

class MSCommunity:

    def __init__(self,model,names=[],abundances=None,pfba = True,lp_filename = None):
        #Setting model and package manager
        self.model = model
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        #Set this attribute to a filename to cause lp files to be created
        self.lp_filename = lp_filename
        #Set this to false to turn off pfba
        self.pfba = pfba
        #Solution from the last run FBA analysis is always stashed here
        self.solution = None
        #Data attributes
        self.biomass_cpd = None
        self.primary_biomass = None
        self.biomass_drain = None
        self.msgapfill  = None
        self.element_uptake_limit = None
        self.kinetic_coeff = None
        self.msdb_path_for_fullthermo = None
        self.species = DictList()
        #Computing data from model
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
    
    #Manipulation functions
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
            
    def constrain(self,element_uptake_limit = None, kinetic_coeff = None,msdb_path_for_fullthermo = None):
        # applying uptake constraints
        self.element_uptake_limit = element_uptake_limit
        if element_uptake_limit is not None:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        # applying kinetic constraints
        self.kinetic_coeff = kinetic_coeff
        if kinetic_coeff is not None:
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
        # applying FullThermo constraints
        self.msdb_path_for_fullthermo = msdb_path_for_fullthermo
        if msdb_path_for_fullthermo is not None:
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_path':msdb_path_for_fullthermo}, verbose)
    
    #Utility functions
    def print_lp(self,filename = None):
        if filename == None:
            filename = self.lp_filename
        if filename != None:
            with open(filename+".lp", 'w') as out:
                out.write(str(self.model.solver))
                out.close()
    
    def compute_interactions(self,solution = None,threshold=1):
        #Checking for existance of solution
        if solution == None:
            solution = self.solution
        if solution == None:
            logger.warning("No feasible solution!")
            return None
        #Initializing datafram
        metabolite_data = {}
        species_data = {"Environment":{}}
        species_list = {"Environment":{}}
        species_array = [None for i in range(1000)]
        met_array = []
        data = {"IDs":[],"Metabolites/Species":[]}
        #Cycling through metabolites and flagging external metabolites
        for met in self.model.metabolites:
            if met.compartment == "e0":
                met_array.append(met)
                data["IDs"].append(met.id)
                data["Metabolites/Species"].append(met.name)
                metabolite_data[met] = {} 
                for individual in self.species:
                    metabolite_data[met][individual.id] = 0
                metabolite_data[met]["Environment"] = 0
        for individual in self.species:
            species_data[individual.id] = {}
            species_list[individual.id] = {}
            data[individual.id] = []
            data["IDs"].append(individual.index)
            data["Metabolites/Species"].append(individual.id)
            for other in self.species:
                species_data[individual.id][other.id] = 0
                species_list[individual.id][other.id] = []
            species_data["Environment"][individual.id] = 0
            species_list["Environment"][individual.id] = []
            species_data[individual.id]["Environment"] = 0
            species_list[individual.id]["Environment"] = []
            species_array[individual.index] = individual
        data["IDs"].append("Environment")
        data["Metabolites/Species"].append("Environment")
        for individual in self.species:
            data["IDs"].append(individual.index)
            data["Metabolites/Species"].append(individual.id+" list")
        data["IDs"].append("Environment")
        data["Metabolites/Species"].append("Environment")
        data["Environment"] = []
        #Cycling through reactions and computing metabolite input and output
        for rxn in self.model.reactions:
            if rxn.id[0:3] == "EX_" and abs(solution.fluxes[rxn.id]) > Zero:
                cpd = list(rxn.metabolites.keys())[0]
                if cpd in metabolite_data:
                    metabolite_data[cpd]["Environment"] += -1*solution.fluxes[rxn.id]
            if len(rxn.id.split("_")) > 1:
                cmp = rxn.id.split("_").pop()
                comp_index = int(cmp[1:])
                for metabolite in rxn.metabolites:
                    if metabolite in metabolite_data:
                        if species_array[comp_index] != None:
                            metabolite_data[metabolite][species_array[comp_index].id] += solution.fluxes[rxn.id]*rxn.metabolites[metabolite]
        #Now translating metbaolite input and output into species interaction flux
        for met in metabolite_data:
            #Iterating through the producers of the metabolite:
            total = 0
            for individual in self.species:
                if metabolite_data[met][individual.id] > Zero:
                    total += metabolite_data[met][individual.id]
            if metabolite_data[met]["Environment"] > Zero:
                total += metabolite_data[met]["Environment"]
            for individual in self.species:
                if metabolite_data[met][individual.id] > Zero:
                    for other in self.species:
                        if metabolite_data[met][other.id] < Zero:
                            species_data[individual.id][other.id] += -1*metabolite_data[met][individual.id]*metabolite_data[met][other.id]/total
                            if -1*metabolite_data[met][individual.id]*metabolite_data[met][other.id]/total > threshold:
                                species_list[individual.id][other.id].append(met.name)
                    if metabolite_data[met]["Environment"] < Zero:
                        species_data[individual.id]["Environment"] += -1*metabolite_data[met][individual.id]*metabolite_data[met]["Environment"]/total
                        if -1*metabolite_data[met][individual.id]*metabolite_data[met]["Environment"]/total > threshold:
                            species_list[individual.id]["Environment"].append(met.name)
            if metabolite_data[met]["Environment"] > Zero:
                for other in self.species:
                    if metabolite_data[met][other.id] < Zero:
                        species_data["Environment"][other.id] += -1*metabolite_data[met]["Environment"]*metabolite_data[met][other.id]/total
                        if -1*metabolite_data[met]["Environment"]*metabolite_data[met][other.id]/total > threshold:
                            species_list["Environment"][other.id].append(met.name)
        #Now reformating all data into a dataframe
        for met in met_array:
            for individual in self.species:
                data[individual.id].append(metabolite_data[met][individual.id])
            data["Environment"].append(metabolite_data[met]["Environment"])
        for individual in self.species:
            for other in self.species:
                data[individual.id].append(species_data[individual.id][other.id])
            data[individual.id].append(species_data[individual.id]["Environment"])
        for individual in self.species:
            data["Environment"].append(species_data["Environment"][individual.id])
        data["Environment"].append(0)
        for individual in self.species:    
            for other in self.species:
                data[individual.id].append("; ".join(species_list[individual.id][other.id]))
            data[individual.id].append("; ".join(species_list[individual.id]["Environment"]))
        for individual in self.species:
            data["Environment"].append("; ".join(species_list["Environment"][individual.id]))
        data["Environment"].append("")
        df = pd.DataFrame(data)
        logger.info(df)
        return df
        
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
    
    #Analysis functions
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
    
    def test_individual_species(self,media = None,allow_interaction=True,run_atp=True,run_biomass=True):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        #Iterating over species and running tests
        data = {"Species":[],"Biomass":[],"ATP":[]}
        for individual in self.species:
            data["Species"].append(individual.id)
            #Running with model so changes are discarded after each simulation
            with self.model:
                #If no interaction allowed, iterate over all other species and disable them
                if not allow_interaction:
                    for indtwo in self.species:
                        if indtwo != individual:
                            indtwo.disable_species()
                #If testing biomass, setting objective to individual species biomass and optimizing
                if run_biomass:
                    data["Biomass"].append(individual.compute_max_biomass())
                #If testing atp, setting objective to individual species atp and optimizing
                if run_atp:
                    data["ATP"].append(individual.compute_max_atp())
        df = pd.DataFrame(data)
        logger.info(df)
        return df
    
    def atp_correction(self,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
        atpcorrect = MSATPCorrection(self.model,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0)      
    
    def predict_abundances(self,media=None,pfba=True,kinetic_coeff = None):
        with self.model:
            #Kinetic coefficients must be used for this formulation to work
            if kinetic_coeff == None:
                kinetic_coeff = self.kinetic_coeff
            if kinetic_coeff == None:
                kinetic_coeff = 2000
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
            objcoef = {}
            for species in self.species:
                objcoef[species.biomasses[0].forward_variable] = 1
            new_objective = self.model.problem.Objective(Zero,direction="max")
            self.model.objective = new_objective
            new_objective.set_linear_coefficients(objcoef)
            self.run(media,pfba)
            return self._compute_relative_abundance_from_solution()
        return None
    
    def run(self,media,pfba = None):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        self.print_lp()
        if pfba == None:
            pfba = self.pfba
        if pfba == None:
            pfba = True
        if pfba:
            self._set_solution(cobra.flux_analysis.pfba(self.model))
        else:
            self._set_solution(self.model.optimize())
        if self.solution == None:
            return None
        logger.info(self.model.summary())
        return self.solution

    #Internal functions
    def _compute_relative_abundance_from_solution(self,solution = None):
        if solution == None:
            solution = self.solution
        if sollution == None:
            logger.warning("No feasible solution!")
            return None
        data = {"Species":[],"Abundance":[]}
        totalgrowth = 0
        for species in self.species:
            totalgrowth += self.solution.fluxes[species.biomasses[0].id]
        if totalgrowth == 0:
            logger.warning("Community model did not grow!")
            return None
        for species in self.species:
            data["Species"].append(species.id)
            data["Abundance"].append(self.solution.fluxes[species.biomasses[0].id]/totalgrowth)
        df = pd.DataFrame(data)
        logger.info(df)
        return df
    
    def _set_solution(self,solution):
        self.solution = None
        if solution.status != 'optimal':
            logger.warning("No solution found for the simulation.")
            return
        self.solution = solution
