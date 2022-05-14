import logging
#import itertools
import cobra
from itertools import combinations
from optlang.symbolics import Zero
import networkx
import sigfig
import pandas as pd
from cobra.core.dictlist import DictList
from cobra.core import Reaction
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msatpcorrection import MSATPCorrection
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from matplotlib import pyplot

logger = logging.getLogger(__name__)

class CommunityModelSpecies:
    def __init__(self,
                 community,         # MSCommunity environment
                 biomass_cpd,       # metabolite in the biomass reaction
                 names=[]           # names of the community species
                 ):
        self.community, self.biomass_cpd = community, biomass_cpd
        self.species_num = int(self.biomass_cpd.compartment[1:])
        self.abundance = 0
        if self.biomass_cpd in self.community.primary_biomass.metabolites:
            self.abundance = abs(self.community.primary_biomass.metabolites[self.biomass_cpd])
        if self.species_num <= len(names) and names[self.species_num-1] != None:
            self.id = names[self.species_num-1]
        else:
            if "species_name" in self.biomass_cpd.annotation:
                self.id = self.biomass_cpd.annotation["species_name"]
            else:
                self.id = "Species"+str(self.species_num)
                
        logger.info("Making atp hydrolysis reaction for species: "+self.id)
        atp_rxn = FBAHelper.add_atp_hydrolysis(self.community.model,"c"+str(self.species_num))
        # FBAHelper.add_autodrain_reactions_to_self.community_model(self.community.model)  #!!! FIXME This FBAHelper function is not defined.
        self.atp_hydrolysis = atp_rxn["reaction"]
        self.biomass_drain = None
        self.biomasses = []
        for reaction in self.community.model.reactions:
            if self.biomass_cpd in reaction.metabolites:
                if reaction.metabolites[self.biomass_cpd] == 1 and len(reaction.metabolites) > 1:
                    self.biomasses.append(reaction)
                elif len(reaction.metabolites) == 1 and reaction.metabolites[self.biomass_cpd] < 0:
                    self.biomass_drain = reaction
                    
        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species "+self.id)
        if self.atp_hydrolysis == None:
            logger.critical("No ATP hydrolysis found for species:"+self.id)
        if self.biomass_drain == None:
            logger.info("Making biomass drain reaction for species: "+self.id)
            self.biomass_drain = Reaction(id="DM_"+self.biomass_cpd.id,
                                      name="DM_" + self.biomass_cpd.name,
                                      lower_bound=0, 
                                      upper_bound=100)
            self.community.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({self.biomass_cpd : -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'
    
    def disable_species(self):
        for reaction in self.community.model.reactions:
            if int(FBAHelper.rxn_compartment(reaction)[1:]) == self.species_num:
                reaction.upper_bound, reaction.lower_bound = 0, 0
    
    def compute_max_biomass(self):
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.biomasses[0].forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_Biomass")
        return self.community.model.optimize()
    
    def compute_max_atp(self):
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_ATP")
        return self.community.model.optimize()

class MSCommunity:

    def __init__(self, model, 
                 names=[], abundances=None,   # names and abundances of the community species
                 pfba = True,                 # specify whether parsimonious FBA will be simulated
                 lp_filename = None           # specify a filename to create an lp file 
                 ):
        #Setting model and package manager
        self.model, self.lp_filename, self.pfba = model, lp_filename, pfba
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        self.gapfillings = {}
        #Define Data attributes as None
        self.solution = self.biomass_cpd = self.primary_biomass = self.biomass_drain = self.msgapfill = self.element_uptake_limit = self.kinetic_coeff = self.modelseed_db_path = None
        self.species = DictList()
        #Computing data from model
        msid_cobraid_hash = FBAHelper.msid_hash(model)
        if "cpd11416" not in msid_cobraid_hash:
            logger.critical("Could not find biomass compound")
        other_biomass_cpds = []
        for biomass_cpd in msid_cobraid_hash["cpd11416"]:
            if biomass_cpd.compartment == "c0":
                self.biomass_cpd = biomass_cpd
                for reaction in model.reactions:
                    if self.biomass_cpd in reaction.metabolites:
                        if reaction.metabolites[self.biomass_cpd] == 1 and len(reaction.metabolites) > 1:
                            self.primary_biomass = reaction
                        elif reaction.metabolites[self.biomass_cpd] < 0 and len(reaction.metabolites) == 1:
                            self.biomass_drain = reaction
            else:
                other_biomass_cpds.append(biomass_cpd)
        for biomass_cpd in other_biomass_cpds:        
            species_obj = CommunityModelSpecies(self,biomass_cpd,names)
            self.species.append(species_obj)
        if abundances != None:
            self.set_abundance(abundances)
    
    #Manipulation functions
    def set_abundance(self,abundances):
        #ensure normalization
        total_abundance = sum([abundances[species] for species in abundances])
        #map abundances to all species
        for species in abundances:
            abundances[species] = abundances[species]/total_abundance
            if species in self.species:
                self.species.get_by_id(species).abundance = abundances[species]
        #remake the primary biomass reaction based on abundances
        if self.primary_biomass == None:
            logger.critical("Primary biomass reaction not found in community model")
        all_metabolites = {self.biomass_cpd:1}
        for species in self.species:
            all_metabolites[species.biomass_cpd] = -1*abundances[species.id]
        self.primary_biomass.add_metabolites(all_metabolites,combine=False)
        
    def set_objective(self,target = None,minimize = False):
        if target == None:
            target = self.primary_biomass.id
        sense = "max"
        if minimize:
            sense = "min"
        self.model.objective = self.model.problem.Objective(
            self.model.reactions.get_by_id(target).flux_expression, 
            direction=sense
        )
            
    def constrain(self,element_uptake_limit = None, kinetic_coeff = None,modelseed_db_path = None):
        # applying uptake constraints
        self.element_uptake_limit = element_uptake_limit
        if element_uptake_limit is not None:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        # applying kinetic constraints
        self.kinetic_coeff = kinetic_coeff
        if kinetic_coeff is not None:
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
        # applying FullThermo constraints
        self.modelseed_db_path = modelseed_db_path
        if modelseed_db_path is not None:
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_db_path':modelseed_db_path})
    
    #Utility functions
    def print_lp(self,filename = None):
        if filename is None:
            filename = self.lp_filename
        if filename is not None:
            with open(filename+".lp", 'w') as out:
                out.write(str(self.model.solver))
                out.close()
    
    def compute_interactions(self,solution = None,threshold=1, visualize = False):
        #Check for existance of solution
        if solution == None:
            solution = self.solution 
        if solution == None:
            logger.warning("No feasible solution!")
            return None
        #Initialize data 
        metabolite_data, species_data, species_collection = {}, {"Environment":{}}, {"Environment":{}}
        data = {"IDs":[],"Metabolites/Species":[], "Environment":[]}
        met_list, species_list = [], [None for i in range(1000)]

        #establish spreadsheet infrastructure for only extracellular metabolites 
        for met in self.model.metabolites:
            if met.compartment == "e0":
                met_list.append(met)
                data["IDs"].append(met.id)
                data["Metabolites/Species"].append(met.name)
                
                metabolite_data[met] = {} 
                metabolite_data[met]["Environment"] = 0
                for individual in self.species:
                    metabolite_data[met][individual.id] = 0
                    
        for individual in self.species:
            species_data[individual.id], species_collection[individual.id] = {}, {}
            species_list[individual.species_num] = individual
            data[individual.id] = []
            data["IDs"].append(individual.species_num)
            data["Metabolites/Species"].append(individual.id)
            for other in self.species:
                species_data[individual.id][other.id] = 0
                species_collection[individual.id][other.id] = []
                
            species_data["Environment"][individual.id] = species_data[individual.id]["Environment"] = 0
            species_collection["Environment"][individual.id], species_collection[individual.id]["Environment"] = [], []
            
        data["IDs"].append("Environment")
        data["Metabolites/Species"].append("Environment")
        for individual in self.species:
            data["IDs"].append(individual.species_num)
            data["Metabolites/Species"].append(individual.id+" list")
            
        # computing net metabolite flux from each reaction
        for rxn in self.model.reactions:
            if rxn.id[0:3] == "EX_" and abs(solution.fluxes[rxn.id]) > Zero:
                cpd = list(rxn.metabolites.keys())[0]   
                if cpd in metabolite_data:
                    metabolite_data[cpd]["Environment"] += -1*solution.fluxes[rxn.id]
            if len(rxn.id.split("_")) > 1:
                cmp = rxn.id.split("_").pop()
                print(cmp)
                comp_index = int(cmp[1:])
                for metabolite in rxn.metabolites:
                    if metabolite in metabolite_data:
                        if species_list[comp_index] != None:
                            metabolite_data[metabolite][species_list[comp_index].id] += solution.fluxes[rxn.id]*rxn.metabolites[metabolite]
                            
        # translating net metbaolite flux into species interaction flux
        for met in metabolite_data:
            #Iterating through the metabolite producers
            total = sum([metabolite_data[met][individual.id] for individual in self.species if metabolite_data[met][individual.id] > Zero])
            if metabolite_data[met]["Environment"] > Zero:
                total += metabolite_data[met]["Environment"]
            for individual in self.species:
                if metabolite_data[met][individual.id] > Zero:
                    for other in self.species:
                        if metabolite_data[met][other.id] < Zero:
                            normalized_flux = -1*metabolite_data[met][individual.id]*metabolite_data[met][other.id]/total
                            species_data[individual.id][other.id] += normalized_flux
                            if normalized_flux > threshold:
                                species_collection[individual.id][other.id].append(met.name)
                    if metabolite_data[met]["Environment"] < Zero:
                        normalized_flux = -1*metabolite_data[met][individual.id]*metabolite_data[met]["Environment"]/total
                        species_data[individual.id]["Environment"] += normalized_flux
                        if normalized_flux > threshold:
                            species_collection[individual.id]["Environment"].append(met.name)
            if metabolite_data[met]["Environment"] > Zero:
                for individual in self.species:
                    if metabolite_data[met][individual.id] < Zero:
                        normalized_flux = -1*metabolite_data[met]["Environment"]*metabolite_data[met][individual.id]/total
                        species_data["Environment"][individual.id] += normalized_flux
                        if normalized_flux > threshold:
                            species_collection["Environment"][individual.id].append(met.name)
                            
        #Reformat data into a dataframe
        for met in met_list:
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
                if other != individual:
                    data[individual.id].append("; ".join(species_collection[individual.id][other.id]))
            data[individual.id].append("; ".join(species_collection[individual.id]["Environment"]))
        for individual in self.species:
            data["Environment"].append("; ".join(species_collection["Environment"][individual.id]))
        self.cross_feeding_df = pd.DataFrame(data)
        logger.info(self.cross_feeding_df)
        
        # graph the network diagram
        if visualize:
            self.graph = networkx.Graph()
            for species in self.species:
                self.graph.add_node(species.num_species)
            for com in combinations(self.compartments, 2): #!!! FIXME the total fluxes between each species must be amalgamated from the metabolite-resoluton fluxes
                species_1 = int(com[0])-1
                species_2 = int(com[1])-1

                interaction_net_flux = sigfig.round((self.production[species_1][species_2] - self.consumption[species_1][species_2]), 3)
                if species_1 < species_2:
                    self.graph.add_edge(com[0],com[1],flux = interaction_net_flux)
                elif species_1 > species_2:
                    self.graph.add_edge(com[0],com[1],flux = -interaction_net_flux)

            self.pos = networkx.circular_layout(self.graph)
            networkx.draw_networkx(self.graph,self.pos)
            self.labels = networkx.get_edge_attributes(self.graph,'flux')
            networkx.draw_networkx_edge_labels(self.graph,self.pos,edge_labels=self.labels)
            pyplot.show()
        
        return self.cross_feeding_df
                   
    
    #Analysis functions
    def gapfill(self, media = None, target = None, minimize = False,default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = [], suffix = None, solver = 'glpk'):
        if target == None:
            target = self.primary_biomass.id
        self.set_objective(target,minimize)
        gfname = FBAHelper.medianame(media)+"-"+target
        if suffix is not None:
            gfname += f"-{suffix}"
        self.gapfillings[gfname] = MSGapfill(self.model, default_gapfill_templates, default_gapfill_models, test_conditions, reaction_scores, blacklist) 
        gfresults = self.gapfillings[gfname].run_gapfilling(media,target, solver = solver)
        if gfresults is None:
            logger.critical("Gapfilling failed with the specified model, media, and target reaction.")
        return self.gapfillings[gfname].integrate_gapfill_solution(gfresults)
    
    def test_individual_species(self,media = None,allow_cross_feeding=True,run_atp=True,run_biomass=True):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        #Iterating over species and running tests
        data = {"Species":[],"Biomass":[],"ATP":[]}
        for individual in self.species:
            data["Species"].append(individual.id)
            with self.model: # WITH, here, discards changes after each simulation
                #If no interaction allowed, iterate over all other species and disable them
                if not allow_cross_feeding:
                    for indtwo in self.species:
                        if indtwo != individual:
                            indtwo.disable_species()
                if run_biomass:  #If testing biomass, setting objective to individual species biomass and optimizing
                    data["Biomass"].append(individual.compute_max_biomass())
                if run_atp:      #If testing atp, setting objective to individual species atp and optimizing
                    data["ATP"].append(individual.compute_max_atp())
        df = pd.DataFrame(data)
        logger.info(df)
        return df
    
    def atp_correction(self,core_template, atp_medias, max_gapfilling=None, gapfilling_delta=0): 
        self.atpcorrect = MSATPCorrection(self.model,core_template, atp_medias, compartment="c0", max_gapfilling=None, gapfilling_delta=0)      
    
    def predict_abundances(self,media=None,pfba=True,kinetic_coeff = None):
        with self.model:   # WITH, here, discards changes after each simulation
            if kinetic_coeff == None:
                kinetic_coeff = self.kinetic_coeff
            if kinetic_coeff == None: #Kinetic coefficients must be used for this formulation to work
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
        if (pfba == None and self.pfba == None) or pfba == True or self.pfba == True:
            self._set_solution(cobra.flux_analysis.pfba(self.model))
        else:
            self._set_solution(self.model.optimize())
        if self.solution == None:
            return None
        logger.info(self.model.summary())
        return self.solution
    

    #Internal functions
    def _compute_relative_abundance_from_solution(self,solution = None):
        if solution == None and self.solution == None:
            logger.warning("No feasible solution!")
            return None
        data = {"Species":[],"Abundance":[]}
        totalgrowth = sum([self.solution.fluxes[species.biomasses[0].id] for species in self.species])
        if totalgrowth == 0:
            logger.warning("The community did not grow!")
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
