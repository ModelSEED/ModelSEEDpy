from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from modelseedpy.core import MSATPCorrection
from cobra import Model, Reaction, Metabolite
from cobra.core.dictlist import DictList
from cobra.io import save_matlab_model
from itertools import combinations
from optlang.symbolics import Zero
from matplotlib import pyplot
from pandas import DataFrame
import logging
#import itertools
import cobra
import networkx
import sigfig
import re, os

logger = logging.getLogger(__name__)

class CommunityModelSpecies:
    def __init__(self,
                 community,         # MSCommunity environment
                 biomass_cpd,       # metabolite in the biomass reaction
                 names=[]           # names of the community species
                 ):
        self.community, self.biomass_cpd = community, biomass_cpd
        self.index = int(self.biomass_cpd.compartment[1:])
        self.abundance = 0
        if self.biomass_cpd in self.community.primary_biomass.metabolites:
            self.abundance = abs(self.community.primary_biomass.metabolites[self.biomass_cpd])
        if self.index <= len(names) and names[self.index-1]:
            self.id = names[self.index-1]
        else:
            if "species_name" in self.biomass_cpd.annotation:
                self.id = self.biomass_cpd.annotation["species_name"]
            else:
                self.id = "Species"+str(self.index)
                
        logger.info("Making atp hydrolysis reaction for species: "+self.id)
        atp_rxn = FBAHelper.add_atp_hydrolysis(self.community.model,"c"+str(self.index))
        # FBAHelper.add_autodrain_reactions_to_self.community_model(self.community.model)  # !!! FIXME This FBAHelper function is not defined.
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
        if not self.biomass_drain:
            logger.info("Making biomass drain reaction for species: "+self.id)
            self.biomass_drain = Reaction(id="DM_"+self.biomass_cpd.id,
                                      name="DM_" + self.biomass_cpd.name,
                                      lower_bound=0, upper_bound=100)
            self.community.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({self.biomass_cpd: -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'
    
    def disable_species(self):
        for reaction in self.community.model.reactions:
            if int(FBAHelper.rxn_compartment(reaction)[1:]) == self.index:
                reaction.upper_bound = reaction.lower_bound = 0
    
    def compute_max_biomass(self):
        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species "+self.id)
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.biomasses[0].forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_Biomass")
        return self.community.model.optimize()
    
    def compute_max_atp(self):
        if not self.atp_hydrolysis:
            logger.critical("No ATP hydrolysis found for species:"+self.id)
        self.community.model.objective = self.community.model.problem.Objective(Zero,direction="max")
        self.community.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
        if self.community.lp_filename:
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
        if abundances:
            self.set_abundance(abundances)
            
    @staticmethod
    def build_from_species_models(models,mdlid=None,name=None,names=[],abundances=None):
        """Merges the input list of single species metabolic models into a community metabolic model
        
        Parameters
        ----------
        models : list<Cobra.Model>
            List of models to be merged into a community model
        mdlid : string
            String specifying community model ID
        name : string
            String specifying community model name
        names : list<string>
            List of human readable names for models being merged
        abundances : dict<string,float>
            Hash of relative abundances for input models in community model
        
        Returns
        -------
        Cobra.Model
            Community model object
            
        Raises
        ------
        """
        newmodel = Model(mdlid,name)
        newutl = MSModelUtil(newmodel)
        biomass_compounds = []
        index = 1
        biomass_index = 2
        for model in models:        
            new_metabolites = []
            new_reactions = []
            #Rename metabolites
            for met in model.metabolites:
                #Renaming compartments
                if re.search('[a-z+](\d*)$', met.compartment):
                    m = re.search('([a-z]+)(\d*)$', met.compartment)
                    if len(m[2]) == 0:
                        if m[1] == "e":
                            met.compartment += "0"
                        else:
                            met.compartment += str(index)
                    elif m[1] == "e":
                        met.compartment = m[1]+"0"
                    else:
                        met.compartment = m[1]+str(index)
                #Processing metabolite ID
                output = MSModelUtil.parse_id(met)
                if output == None:
                    if met.compartment[0] != "e":
                        met.id += str(index)
                elif output[1] != "e":
                    if len(output[2]) == 0:
                        met.id = met.id+str(index)
                    else:
                        met.id = output[0]+"_"+output[1]+str(index)
                if met.id not in newmodel.metabolites:
                    new_metabolites.append(met)
                    if met.id == "cpd11416":
                        biomass_compounds.append(met)
            #Rename reactions
            for rxn in model.reactions:
                if rxn.id[0:3] != "EX_":
                    if re.search('^(bio)(\d+)$', rxn.id) != None:
                        rxn.id = "bio"+str(biomass_index)
                        biomass_index += 1
                    else:
                        output = MSModelUtil.parse_id(rxn)
                        if output == None:
                            if rxn.compartment.id[0] != "e":
                                rxn.id += str(index)
                        elif output[1] != "e":
                            if len(output[2]) == 0:
                                rxn.id = rxn.id+str(index)
                            else:
                                rxn.id = output[0]+"_"+output[1]+str(index)
                if rxn.id not in newmodel.reactions:
                    new_reactions.append(rxn)
            #Adding new reactions and compounds to base model
            newmodel.add_reactions(new_reactions)
            newmodel.add_metabolites(new_metabolites)
            index += 1
        #Create community biomass
        comm_biomass = Metabolite("cpd11416_c0", None, "Community biomass", 0, "c0")
        metabolites = {comm_biomass : 1}
        comm_biorxn = Reaction(id="bio1", name= "bio1", lower_bound=0, upper_bound=100)
        count = len(biomass_compounds)
        for cpd in biomass_compounds:
            metabolites[cpd] = -1/count
        comm_biorxn.add_metabolites(metabolites)
        newmodel.add_reactions([comm_biorxn])
        newutl.add_exchanges_for_metabolites([comm_biomass],0,100,'SK_')
        return MSCommunity(newmodel,names,abundances)
    
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
        
    def set_objective(self,target = None,minimize = False): #!!! Mustn't a multilevel objective be set for community models?
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
        if element_uptake_limit:
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        # applying kinetic constraints
        self.kinetic_coeff = kinetic_coeff
        if kinetic_coeff:
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)
        # applying FullThermo constraints
        self.modelseed_db_path = modelseed_db_path
        if modelseed_db_path:
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_db_path':modelseed_db_path})
    
    #Utility functions
    def print_lp(self,filename = None):
        if not filename:
            filename = self.lp_filename
        if filename:
            with open(filename+".lp", 'w') as out:
                out.write(str(self.model.solver))
                out.close()
    
    def compute_interactions(self, solution = None,               # the COBRA simulation solution that will be parsed and visualized
                             threshold: int = 1,                  #!!! What is this threshold?
                             visualize: bool = True,              # specifies whether the net flux will be depicted in a network diagram
                             export_directory: str = None,        # specifies the directory to which the network diagram and associated datatable will be exported, where None does not export the content
                             node_metabolites: bool = True,       # specifies whether the metabolites of each node will be printed
                             x_offset: float = 0.15,              # specifies the x-axis buffer between each species node and its metabolite list in the network diagram
                             show_figure: bool = True             # specifies whether the figure will be printed to the console
                             ):
        #Check for solution
        if not solution:
            solution = self.solution 
        if not solution:
            logger.warning("No feasible solution!")
            return None
        
        #Initialize data 
        metabolite_data, species_data, species_collection = {}, {"Environment":{}}, {"Environment":{}}
        data = {"IDs":[],"Metabolites/Donor":[], "Environment":[]}
        met_list, species_list = [], [None for i in range(1000)]

        #establish spreadsheet infrastructure for only extracellular metabolites 
        for met in self.model.metabolites:
            if met.compartment == "e0":
                met_list.append(met)
                data["IDs"].append(met.id)
                data["Metabolites/Donor"].append(met.name)
                
                metabolite_data[met] = {} 
                metabolite_data[met]["Environment"] = 0
                for individual in self.species:
                    metabolite_data[met][individual.id] = 0
                    
        for individual in self.species:
            species_data[individual.id], species_collection[individual.id] = {}, {}
            species_list[individual.index] = individual
            data[individual.id] = []
            data["IDs"].append(individual.index)
            data["Metabolites/Donor"].append(individual.id)
            for other in self.species:
                species_data[individual.id][other.id] = 0
                species_collection[individual.id][other.id] = []
                
            species_data["Environment"][individual.id] = species_data[individual.id]["Environment"] = 0
            species_collection["Environment"][individual.id], species_collection[individual.id]["Environment"] = [], []
            
        data["IDs"].append("Environment")
        data["Metabolites/Donor"].append("Environment")
        for individual in self.species:
            data["IDs"].append(individual.index)
            data["Metabolites/Donor"].append(individual.id+" list")
            
        # computing net metabolite flux from each reaction
        for rxn in self.model.reactions:
            if rxn.id[0:3] == "EX_" and abs(solution.fluxes[rxn.id]) > Zero:
                cpd = list(rxn.metabolites.keys())[0]   
                if cpd in metabolite_data:
                    metabolite_data[cpd]["Environment"] += -1*solution.fluxes[rxn.id]
            if len(rxn.id.split("_")) > 1:
                comp_index = int(rxn.id.split("_")[-1][1:])
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
                    # calculate the total net flux between each combination of species, and track the involved metabolites
                    for other in self.species:
                        if metabolite_data[met][other.id] < Zero:
                            normalized_flux = abs(metabolite_data[met][individual.id]*metabolite_data[met][other.id])/total
                            species_data[individual.id][other.id] += normalized_flux
                            if normalized_flux > threshold:
                                species_collection[individual.id][other.id].append(met.name)
                    # calculate the total net flux between the species and the environment, and track the involved metabolites
                    if metabolite_data[met]["Environment"] < Zero:
                        normalized_flux = abs(metabolite_data[met][individual.id]*metabolite_data[met]["Environment"])/total
                        species_data[individual.id]["Environment"] += normalized_flux
                        if normalized_flux > threshold:
                            species_collection[individual.id]["Environment"].append(met.name)
            if metabolite_data[met]["Environment"] > Zero:
                for individual in self.species:
                    if metabolite_data[met][individual.id] < Zero:
                        normalized_flux = abs(metabolite_data[met]["Environment"]*metabolite_data[met][individual.id])/total
                        species_data["Environment"][individual.id] += normalized_flux
                        if normalized_flux > threshold:
                            species_collection["Environment"][individual.id].append(met.name)
                            
        # construct a dataframe
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
                data[individual.id].append("; ".join(species_collection[individual.id][other.id]))
            data[individual.id].append("; ".join(species_collection[individual.id]["Environment"]))
        for individual in self.species:
            data["Environment"].append("; ".join(species_collection["Environment"][individual.id]))
        data["Environment"].append(0), data["IDs"].append("Environment list"), data["Metabolites/Donor"].append("Environment list")
        
        self.cross_feeding_df = DataFrame(data)
        logger.info(self.cross_feeding_df)
        
        # graph the network diagram
        if visualize:
            self._visualize_cross_feeding(export_directory, node_metabolites, x_offset, show_figure)
        
        return self.cross_feeding_df
    
    def _visualize_cross_feeding(self, export_directory, node_metabolites = True, x_offset = 0.15, show_figure = True):
        # construct an efficient DataFrame of the cross-feeding interactions
        net_cross_feeding = {}
        for index, row in self.cross_feeding_df.iterrows():
            if re.search('Species\d+', row["Metabolites/Donor"]):
                net_cross_feeding[row["Metabolites/Donor"]] = row[len(self.species):]
        
        # define species and the metabolite fluxes
        net_cross_feeding = DataFrame(net_cross_feeding)
        self.graph = networkx.Graph()
        species_nums = {}
        for species in self.species:
            species_nums[species.index]= set()
            self.graph.add_node(species.index)
            for index, entry in net_cross_feeding[f'Species{species.index} list'].iteritems():
                if 'Species' in index and re.search('(\d+)', index).group() != species.index:
                    species_nums[species.index].update(entry.split('; '))
                    
        # define the net fluxes for each combination of two species
        for species_1, species_2 in combinations(list(species_nums.keys()), 2):
            species_2_to_1 = net_cross_feeding.at[f'Species{species_2}', f'Species{species_1}']
            species_1_to_2 = net_cross_feeding.at[f'Species{species_1}', f'Species{species_2}']
            interaction_net_flux = sigfig.round(species_2_to_1 - species_1_to_2, 3)
            self.graph.add_edge(species_1,species_2,flux = interaction_net_flux)  # The graph plots directionally toward the larger numbered species

        # compose the nextwork diagram of net fluxes
        self.pos = networkx.circular_layout(self.graph)
        if node_metabolites:
            for species in self.pos:
                x, y = self.pos[species]
                metabolites = '\n'.join(species_nums[species])
                pyplot.text(x+x_offset, y, metabolites)
        networkx.draw_networkx(self.graph,self.pos)
        self.labels = networkx.get_edge_attributes(self.graph,'flux')
        networkx.draw_networkx_edge_labels(self.graph,self.pos,edge_labels=self.labels)
        
        if export_directory:
            pyplot.savefig(os.path.join(export_directory, 'cross_feeding_diagram.svg'))
            self.cross_feeding_df.to_csv(os.path.join(export_directory, 'cross_feeding.csv'))
            
        if show_figure:
            pyplot.show()
    
    #Analysis functions
    def gapfill(self, media = None, target = None, minimize = False,default_gapfill_templates = [], default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = [], suffix = None, solver = 'glpk'):
        if not target:
            target = self.primary_biomass.id
        self.set_objective(target,minimize)
        gfname = FBAHelper.medianame(media)+"-"+target
        if suffix:
            gfname += f"-{suffix}"
        self.gapfillings[gfname] = MSGapfill(self.model, default_gapfill_templates, default_gapfill_models, test_conditions, reaction_scores, blacklist) 
        gfresults = self.gapfillings[gfname].run_gapfilling(media,target, solver = solver)
        if not gfresults:
            logger.critical("Gapfilling failed with the specified model, media, and target reaction.")
            return None
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
        df = DataFrame(data)
        logger.info(df)
        return df
    
    def atp_correction(self,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0): 
        self.atpcorrect = MSATPCorrection(self.model,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0)      
    
    def predict_abundances(self,media=None,pfba=True,kinetic_coeff = None):
        with self.model:   # WITH, here, discards changes after each simulation
            if not kinetic_coeff:
                kinetic_coeff = self.kinetic_coeff
            if not kinetic_coeff: #Kinetic coefficients must be used for this formulation to work
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
        save_matlab_model(self.model, self.model.name+".mat")        
        if pfba or self.pfba:
            self._set_solution(cobra.flux_analysis.pfba(self.model))
        else:
            self._set_solution(self.model.optimize())
        if not self.solution:
            return None
        logger.info(self.model.summary())
        return self.solution
    

    #Internal functions
    def _compute_relative_abundance_from_solution(self,solution = None):
        if not solution and not self.solution:
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
        df = DataFrame(data)
        logger.info(df)
        return df
    
    def _set_solution(self,solution):
        self.solution = None
        if solution.status != 'optimal':
            logger.warning("No solution found for the simulation.")
            return
        self.solution = solution

    def steady_com(self,):
        from reframed.community import SteadyCom, SteadyComVA
        
        reframed_model = FBAHelper.get_reframed_model(self.model)
        