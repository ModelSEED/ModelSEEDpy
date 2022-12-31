# -*- coding: utf-8 -*-
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.community.mssteadycom import MSSteadyCom
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.exceptions import ObjectAlreadyDefinedError, ParameterError
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
#from modelseedpy.fbapkg.gapfillingpkg import default_blacklist
from modelseedpy.core.msatpcorrection import MSATPCorrection
from cobra import Reaction
from cobra.core.dictlist import DictList
from cobra.io import save_matlab_model
from optlang.symbolics import Zero
from pandas import DataFrame
from pprint import pprint
import logging

# import itertools
import cobra
import re, os

logger = logging.getLogger(__name__)


def isNumber(string):
    try:
        float(string)
        return True
    except:
        return False


class CommunityModelSpecies:
    def __init__(self,
                 community,         # MSCommunity environment
                 biomass_cpd,       # metabolite in the biomass reaction
                 names=None,        # names of the community species
                 name=None,         # the name of a species
                 index=None         # the index of the species
                 ):
        self.community, self.biomass_cpd = community, biomass_cpd
        print(self.biomass_cpd.compartment)
        self.index = int(self.biomass_cpd.compartment[1:]) # if index is None else index
        self.abundance = 0
        if self.biomass_cpd in self.community.primary_biomass.metabolites:
            self.abundance = abs(self.community.primary_biomass.metabolites[self.biomass_cpd])
        if name:
            self.id = name
        elif names and self.index < len(names):
            self.id = names[self.index-1]
        else:
            if "species_name" in self.biomass_cpd.annotation:
                self.id = self.biomass_cpd.annotation["species_name"]
            else:
                self.id = "Species"+str(self.index)

        logger.info("Making atp hydrolysis reaction for species: "+self.id)
        atp_rxn = FBAHelper.add_atp_hydrolysis(self.community.util.model,"c"+str(self.index))
        # FBAHelper.add_autodrain_reactions_to_self.community_model(self.community.model)  # !!! FIXME This FBAHelper function is not defined.
        self.atp_hydrolysis = atp_rxn["reaction"]
        self.biomass_drain = None
        self.biomasses, self.reactions = [], []
        for rxn in self.community.util.model.reactions:
            rxnComp = FBAHelper.rxn_compartment(rxn)
            if not rxnComp:
                print(f"The reaction {rxn.id} strangely lacks a compartment.")
            elif int(rxnComp[1:]) == self.index and 'bio' not in rxn.name:
                self.reactions.append(rxn)
            if self.biomass_cpd in rxn.metabolites:
                if rxn.metabolites[self.biomass_cpd] == 1 and len(rxn.metabolites) > 1:
                    self.biomasses.append(rxn)
                elif len(rxn.metabolites) == 1 and rxn.metabolites[self.biomass_cpd] < 0:
                    self.biomass_drain = rxn

        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species " + self.id)
        if not self.biomass_drain:
            logger.info("Making biomass drain reaction for species: "+self.id)
            self.biomass_drain = Reaction(
                id="DM_"+self.biomass_cpd.id, name="DM_" + self.biomass_cpd.name, lower_bound=0, upper_bound=100)
            self.community.util.model.add_reactions([self.biomass_drain])
            self.biomass_drain.add_metabolites({self.biomass_cpd: -1})
            self.biomass_drain.annotation["sbo"] = 'SBO:0000627'

    def disable_species(self):
        for reaction in self.community.model.reactions:
            reaction_index = FBAHelper.rxn_compartment(reaction)[1:]
            if int(reaction_index) == self.index:
                reaction.upper_bound = reaction.lower_bound = 0

    def compute_max_biomass(self):
        if len(self.biomasses) == 0:
            logger.critical("No biomass reaction found for species "+self.id)
        FBAHelper.add_objective(self.community.model, Zero, coef={self.biomasses[0].forward_variable:1})
        if self.community.lp_filename != None:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_Biomass")
        return self.community.model.optimize()

    def compute_max_atp(self):
        if not self.atp_hydrolysis:
            logger.critical("No ATP hydrolysis found for species:"+self.id)
        FBAHelper.add_objective(self.community.model, Zero, coef={self.atp_hydrolysis.forward_variable:1})
        if self.community.lp_filename:
            self.community.print_lp(self.community.lp_filename+"_"+self.id+"_ATP")
        return self.community.model.optimize()

    def compute_max_atp(self):
        if not self.atp_hydrolysis:
            logger.critical("No ATP hydrolysis found for species:" + self.id)
        self.community.model.objective = self.community.model.problem.Objective(
            Zero, direction="max"
        )
        self.community.model.objective.set_linear_coefficients(
            {self.atp_hydrolysis.forward_variable: 1}
        )
        if self.community.lp_filename:
            self.community.print_lp(self.community.lp_filename + "_" + self.id + "_ATP")
        return self.community.model.optimize()


class MSCommunity:
    def __init__(self, model=None,            # the model that will be defined
                 models:list=None,            # the list of models that will be assembled into a community
                 names=None, abundances=None, # names and abundances of the community species
                 pfba = True,                 # specify whether parsimonious FBA will be simulated
                 lp_filename = None           # specify a filename to create an lp file
                 ):
        self.lp_filename, self.pfba = lp_filename, pfba
        self.gapfillings = {}

        #Define Data attributes as None
        self.solution = self.biomass_cpd = self.primary_biomass = self.biomass_drain = None
        self.msgapfill = self.element_uptake_limit = self.kinetic_coeff = self.msdb_path = None
        self.species = DictList()

        # defining the models
        model = model if not models else build_from_species_models(
            models, names=names, abundances=abundances, cobra_model=True)
        self.util = MSModelUtil(model)
        self.id = self.util.model.id
        self.pkgmgr = MSPackageManager.get_pkg_mgr(self.util.model)
        msid_cobraid_hash = FBAHelper.msid_hash(self.util.model)
        if "cpd11416" not in msid_cobraid_hash:
            raise KeyError("Could not find biomass compound for the model.")
        other_biomass_cpds = []
        for self.biomass_cpd in msid_cobraid_hash["cpd11416"]:
            print(self.biomass_cpd)
            if self.biomass_cpd.compartment == "c0":
                for reaction in self.util.model.reactions:
                    if self.biomass_cpd in reaction.metabolites:
                        print(reaction)
                        if reaction.metabolites[self.biomass_cpd] == 1 and len(reaction.metabolites) > 1:
                            if self.primary_biomass:
                                raise ObjectAlreadyDefinedError(
                                    f"The primary biomass {self.primary_biomass} is already defined,"
                                    f"hence, the {reaction} cannot be defined as the model primary biomass.")
                            print('primary biomass defined', reaction)
                            self.primary_biomass = reaction
                        elif (
                            reaction.metabolites[self.biomass_cpd] < 0
                            and len(reaction.metabolites) == 1
                        ):
                            self.biomass_drain = reaction
            elif 'c' in self.biomass_cpd.compartment:
                other_biomass_cpds.append(self.biomass_cpd)
        for biomass_cpd in other_biomass_cpds:
            print(biomass_cpd)
            self.species.append(CommunityModelSpecies(self, biomass_cpd, names))
        if abundances:
            self.set_abundance(abundances)

    #Manipulation functions
    def set_abundance(self, abundances):
        #calculate the normalized biomass
        total_abundance = sum([abundances[species] for species in abundances])
        # map abundances to all species
        for species in abundances:
            abundances[species] /= total_abundance
            if species in self.species:
                self.species.get_by_id(species).abundance = abundances[species]
        #remake the primary biomass reaction based on abundances
        if self.primary_biomass is None:
            logger.critical("Primary biomass reaction not found in community model")
        all_metabolites = {self.biomass_cpd: 1}
        for species in self.species:
            all_metabolites[species.biomass_cpd] = -abundances[species.id]
        self.primary_biomass.add_metabolites(all_metabolites,combine=False)

    def set_objective(self,target = None,minimize = False):  #!!! Mustn't a multilevel objective be set for community models?
        if target is None:
            target = self.primary_biomass.id
        sense = "max"
        if minimize:
            sense = "min"
        self.util.model.objective = self.util.model.problem.Objective(
            self.util.model.reactions.get_by_id(target).flux_expression,
            direction=sense
        )

    def constrain(self, element_uptake_limit=None, kinetic_coeff=None, msdb_path=None):
        # applying uptake constraints
        if element_uptake_limit:
            self.element_uptake_limit = element_uptake_limit
            self.pkgmgr.getpkg("ElementUptakePkg").build_package(element_uptake_limit)
        # applying kinetic constraints
        if kinetic_coeff:
            self.kinetic_coeff = kinetic_coeff
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff, self)
        # applying FullThermo constraints
        if msdb_path:
            self.msdb_path = msdb_path
            self.pkgmgr.getpkg("FullThermoPkg").build_package({'modelseed_db_path':msdb_path})

    def steadycom(self, solution=None, media=None, filename=None, export_directory=None,
                  node_metabolites=True, flux_threshold=1, visualize=True, ignore_mets=None):
        return MSSteadyCom.compute(
            self, self.util.model, None, solution or self.run(media), filename=filename, export_directory=export_directory,
            node_metabolites=node_metabolites, flux_threshold=flux_threshold, visualize=visualize,
            show_figure=True, ignore_mets=ignore_mets)

    #Utility functions
    def print_lp(self,filename = None):
        if not filename:
            filename = self.lp_filename
        if filename:
            with open(filename+".lp", 'w') as out:
                out.write(str(self.util.model.solver))
                out.close()

    #Analysis functions
    def gapfill(self, media = None, target = None, minimize = False, default_gapfill_templates=None, default_gapfill_models=None,
                test_conditions=None, reaction_scores=None, blacklist=None, suffix = None, solver = 'glpk'):
        default_gapfill_templates = default_gapfill_templates or []
        default_gapfill_models = default_gapfill_models or []
        test_conditions = test_conditions or []
        reaction_scores = reaction_scores or {}
        blacklist = blacklist or []
        if not target:
            target = self.primary_biomass.id
        self.set_objective(target, minimize)
        gfname = FBAHelper.medianame(media) + "-" + target
        if suffix:
            gfname += f"-{suffix}"
        self.gapfillings[gfname] = MSGapfill(self.util.model, default_gapfill_templates, default_gapfill_models, test_conditions, reaction_scores, blacklist)
        gfresults = self.gapfillings[gfname].run_gapfilling(media,target, solver = solver)
        if not gfresults:
            logger.critical(
                "Gapfilling failed with the specified model, media, and target reaction."
            )
            return None
        return self.gapfillings[gfname].integrate_gapfill_solution(gfresults)

    def test_individual_species(self, media=None, allow_cross_feeding=True, run_atp=True, run_biomass=True):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        # Iterating over species and running tests
        data = {"Species": [], "Biomass": [], "ATP": []}
        for individual in self.species:
            data["Species"].append(individual.id)
            with self.util.model:
                #If no interaction allowed, iterate over all other species and disable them
                if not allow_cross_feeding:
                    for indtwo in self.species:
                        if indtwo != individual:
                            indtwo.disable_species()
                if (
                    run_biomass
                ):  # If testing biomass, setting objective to individual species biomass and optimizing
                    data["Biomass"].append(individual.compute_max_biomass())
                if (
                    run_atp
                ):  # If testing atp, setting objective to individual species atp and optimizing
                    data["ATP"].append(individual.compute_max_atp())
        df = DataFrame(data)
        logger.info(df)
        return df

    def atp_correction(self,core_template, atp_medias, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
        self.atpcorrect = MSATPCorrection(self.util.model,core_template, atp_medias,
                                          atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0)

    def predict_abundances(self,media=None,pfba=True,kinetic_coeff = None):
        with self.util.model:   # WITH, here, discards changes after each simulation
            if not kinetic_coeff:
                kinetic_coeff = self.kinetic_coeff
            if not kinetic_coeff: #Kinetic coefficients must be used for this formulation to work
                kinetic_coeff = 2000
            self.pkgmgr.getpkg("CommKineticPkg").build_package(kinetic_coeff,self)

            objcoef = {}
            for species in self.species:
                objcoef[species.biomasses[0].forward_variable] = 1
            new_objective = self.util.model.problem.Objective(Zero,direction="max")
            self.util.model.objective = new_objective
            new_objective.set_linear_coefficients(objcoef)
            self.run(media, pfba)
            return self._compute_relative_abundance_from_solution()

    def run(self,media,pfba = None):  # !!! why is the media needed to execute the model?
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        self.print_lp()
        save_matlab_model(self.util.model, self.util.model.name+".mat")
        if pfba or self.pfba:
            self._set_solution(cobra.flux_analysis.pfba(self.util.model))
        else:
            self._set_solution(self.util.model.optimize())
        if not self.solution:
            return None
        logger.info(self.util.model.summary())
        return self.solution

    def _compute_relative_abundance_from_solution(self,solution = None):
        if not solution and not self.solution:
            logger.warning("No feasible solution!")
            return None
        data = {"Species": [], "Abundance": []}
        totalgrowth = sum(
            [self.solution.fluxes[species.biomasses[0].id] for species in self.species]
        )
        if totalgrowth == 0:
            logger.warning("The community did not grow!")
            return None
        for species in self.species:
            data["Species"].append(species.id)
            data["Abundance"].append(
                self.solution.fluxes[species.biomasses[0].id] / totalgrowth
            )
        df = DataFrame(data)
        logger.info(df)
        return df

    def _set_solution(self,solution):
        self.solution = None
        if solution.status != "optimal":
            logger.warning("No solution found for the simulation.")
            return None
        self.solution = solution
