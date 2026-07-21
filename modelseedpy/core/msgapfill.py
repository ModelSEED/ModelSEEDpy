#!/usr/bin/python
# -*- coding: utf-8 -*-
import logging
import cobra
import re
import json
import numpy as np
import pandas as pd
import time
from optlang.symbolics import Zero, add
from modelseedpy.core import FBAHelper  # !!! the import is never used
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.exceptions import GapfillingError
from collections import defaultdict


logger = logging.getLogger(__name__)
logger.setLevel(
    logging.WARNING
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO


class MSGapfill:
    @staticmethod
    def gapfill_count(solution):
        total = 0
        if "new" in solution:
            total += len(solution["new"])
        if "reversed" in solution:
            total += len(solution["reversed"])
        return total

    def __init__(
        self,
        model_or_mdlutl,
        default_gapfill_templates=[],
        default_gapfill_models=[],
        test_conditions=[],
        reaction_scores={},
        blacklist=[],
        atp_gapfilling=False,
        minimum_obj=0.01,
        default_excretion=100,
        default_uptake=0,
        default_target=None,
        base_media = None,
        base_media_target_element = "C"
    ):
        # Discerning input is model or mdlutl and setting internal links
        if isinstance(model_or_mdlutl, MSModelUtil):
            self.model = model_or_mdlutl.model
            self.mdlutl = model_or_mdlutl
        else:
            self.model = model_or_mdlutl
            self.mdlutl = MSModelUtil.get(model_or_mdlutl)
        # Setting gapfilling attribute in model utl so link is bidirectional
        if not atp_gapfilling:
            self.mdlutl.gfutl = self
        self.auto_sink = [
            "cpd01042",
            "cpd02701",
            "cpd11416",
            "cpd15302",
            "cpd03091",
        ]  # the cpd11416 compound is filtered during model extension with templates
        # Cloning model to create gapfilling model
        self.gfmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
        self.gfmodel.id = self.model.id + "_gf"
        self.gfmodelutl = MSModelUtil.get(self.gfmodel)
        # Getting package manager for gapfilling model
        self.gfpkgmgr = MSPackageManager.get_pkg_mgr(self.gfmodelutl)
        # Setting target from input
        if default_target:
            self.default_target = default_target
            self.gfmodel.objective = self.gfmodel.problem.Objective(
                self.gfmodel.reactions.get_by_id(default_target).flux_expression,
                direction="max",
            )
        # Setting parameters for gapfilling
        self.last_solution = None
        self.model_penalty = 1
        self.default_minimum_objective = minimum_obj
        self.default_gapfill_models = default_gapfill_models
        self.default_gapfill_templates = default_gapfill_templates
        self.gapfill_templates_by_index, self.gapfill_models_by_index = {}, {}
        self.gapfill_all_indecies_with_default_templates = True
        self.gapfill_all_indecies_with_default_models = True
        self.blacklist = list(set(blacklist))
        self.test_condition_iteration_limit = 10
        self.test_conditions = test_conditions
        self.reaction_scores = reaction_scores
        self.default_excretion = default_excretion
        self.default_uptake = default_uptake
        self.minimum_obj = minimum_obj
        self.base_media = base_media
        self.base_media_target_element = base_media_target_element
        self.cumulative_gapfilling = []
        # Building gapfilling package
        self.gfpkgmgr.getpkg("GapfillingPkg").build_package(
            {
                "auto_sink": self.auto_sink,
                "model_penalty": self.model_penalty,
                "default_gapfill_models": self.default_gapfill_models,
                "default_gapfill_templates": self.default_gapfill_templates,
                "gapfill_templates_by_index": self.gapfill_templates_by_index,
                "gapfill_models_by_index": self.gapfill_models_by_index,
                "gapfill_all_indecies_with_default_templates": self.gapfill_all_indecies_with_default_templates,
                "gapfill_all_indecies_with_default_models": self.gapfill_all_indecies_with_default_models,
                "default_excretion": default_excretion,
                "default_uptake": default_uptake,
                "minimum_obj": minimum_obj,
                "blacklist": self.blacklist,
                "reaction_scores": self.reaction_scores,
                "set_objective": 1,
                "base_media": base_media,
                "base_media_target_element":base_media_target_element
            }
        )

    def test_gapfill_database(self, media, target=None, before_filtering=True,active_reactions=[]):
        # Testing if gapfilling can work before filtering
        if target:
            self.gfpkgmgr.getpkg("GapfillingPkg").set_base_objective(target,None)
        else:
            target = str(self.gfmodel.objective)
            target = target.split(" ")[0]
            target = target[13:]
        #Setting media
        self.gfpkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if self.gfpkgmgr.getpkg("GapfillingPkg").test_gapfill_database(active_reactions):
            return True
        if self.gfpkgmgr.getpkg("GapfillingPkg").test_solution.status == 'infeasible':
            return False
        gf_sensitivity = {}
        if target != "rxn00062_c0":
            gf_sensitivity = self.mdlutl.get_attributes("gf_sensitivity", {})
        if media.id not in gf_sensitivity:
            gf_sensitivity[media.id] = {}
        if target not in gf_sensitivity[media.id]:
            gf_sensitivity[media.id][target] = {}
        filter_msg = " "
        note = "FAF"
        if before_filtering:
            filter_msg = " before filtering "
            note = "FBF"
        gf_sensitivity[media.id][target][
            note
        ] = self.mdlutl.find_unproducible_biomass_compounds(target)
        if target != "rxn00062_c0":
            self.mdlutl.save_attributes(gf_sensitivity, "gf_sensitivity")
        logger.warning(
            "No gapfilling solution found"
            + filter_msg
            + "for "
            + media.id
            + " activating "
            + target
        )
        return False

    def test_and_adjust_gapfilling_conditions(self,medias,targets,thresholds,prefilter=True):
        output = {
            "medias":[],
            "targets":[],
            "thresholds":[],
            "conditions":[],
            "active_reactions":[]
        }
        logger.debug("Testing unfiltered database")
        for i,media in enumerate(medias):
            active_reactions = []
            if self.test_gapfill_database(media,targets[i],before_filtering=True,active_reactions=active_reactions):
                output["medias"].append(media)
                output["targets"].append(targets[i])
                output["thresholds"].append(thresholds[i])
                output["active_reactions"].append(active_reactions)
                output["conditions"].append({
                    "media": media,
                    "is_max_threshold": False,
                    "threshold": thresholds[i],
                    "objective": targets[i],
                })
        # Filtering
        if prefilter:
            logger.debug("Filtering database")
            self.prefilter(growth_conditions=output["conditions"],active_reaction_sets=output["active_reactions"])
            medias = []
            targets = []
            thresholds = []
            conditions = []
            active_reaction_sets = []
            logger.debug("Testing filtered database")
            for i,media in enumerate(output["medias"]):
                active_reactions = []
                if self.test_gapfill_database(media,output["targets"][i],before_filtering=False,active_reactions=active_reactions):
                    medias.append(media)
                    targets.append(output["targets"][i])
                    thresholds.append(output["thresholds"][i])
                    conditions.append(output["conditions"][i])
                    active_reaction_sets.append(active_reactions)
            output["medias"] = medias
            output["targets"] = targets
            output["thresholds"] = thresholds
            output["conditions"] = conditions
            output["active_reactions"] = active_reaction_sets
        return output

    def prefilter(self,test_conditions=None,growth_conditions=[],use_prior_filtering=False,base_filter_only=False,active_reaction_sets=[]):
        """Prefilters the database by removing any reactions that break specified ATP tests
        Parameters
        ----------
        test_conditions : []
            List of conditions to be tested when filtering the gapfilling database. If not specified, the test_conditions attribute will be used
        """
        if not test_conditions:
            test_conditions = self.test_conditions
        if self.test_conditions:
            logger.debug(f"PREFILTERING WITH {str(len(growth_conditions))} GROWTH CONDITIONS")
            base_filter = None
            if use_prior_filtering:
                base_filter = self.mdlutl.get_attributes("gf_filter", {})
            self.gfpkgmgr.getpkg("GapfillingPkg").filter_database_based_on_tests(
                self.test_conditions,
                growth_conditions=growth_conditions,
                base_filter=base_filter,
                base_filter_only=base_filter_only,
                active_reaction_sets=active_reaction_sets
            )
            gf_filter = self.gfpkgmgr.getpkg("GapfillingPkg").modelutl.get_attributes(
                "gf_filter", {}
            )
            base_filter = self.mdlutl.get_attributes("gf_filter", {})
            for media_id in gf_filter:
                base_filter[media_id] = gf_filter[media_id]

    def run_gapfilling(
        self,
        media=None,
        target=None,
        minimum_obj=None,
        binary_check=False,
        prefilter=True,
    ):
        """Run gapfilling on a single media condition to force the model to achieve a nonzero specified objective
        Parameters
        ----------
        media : MSMedia
            Media in which the model should be gapfilled
        target : string
            Name or expression describing the reaction or combination of reactions to the optimized
        minimum_obj : double
            Value to use for the minimal objective threshold that the model must be gapfilled to achieve
        binary_check : bool
            Indicates if the solution should be checked to ensure it is minimal in the number of reactions involved
        prefilter : bool
            Indicates if the gapfilling database should be prefiltered using the tests provided in the MSGapfill constructor before running gapfilling
        """
        # Setting target and media if specified
        if not target:
            target = self.default_target
        if not minimum_obj:
            minimum_obj = self.default_minimum_objective
        self.gfpkgmgr.getpkg("GapfillingPkg").set_base_objective(target,minimum_obj)
        if media:
            self.gfpkgmgr.getpkg("GapfillingPkg").set_media(media)

        # Testing if gapfilling can work before filtering
        if not self.test_gapfill_database(media,target,before_filtering=prefilter):
            return None

        # Filtering
        if prefilter:
            self.prefilter(growth_conditions=[{
                "media": media,
                "is_max_threshold": False,
                "threshold": minimum_obj,
                "objective": target,
            }])
            if not self.test_gapfill_database(media,target,before_filtering=False):
                return None

        # Printing the gapfilling LP file
        self.mdlutl.printlp(model=self.gfmodel,filename="StandardGapfill",print=False)

        # Running gapfil/ling and checking solution
        sol = self.gfmodel.optimize()
        logger.debug(
            f"gapfill solution objective value {sol.objective_value} ({sol.status}) for media {media}"
        )
        if sol.status != "optimal":
            logger.warning("No solution found for %s", media)
            return None

        # Computing solution and ensuring all tests still pass
        self.last_solution = self.gfpkgmgr.getpkg(
            "GapfillingPkg"
        ).compute_gapfilled_solution()
        if self.test_conditions:
            self.last_solution = self.gfpkgmgr.getpkg(
                "GapfillingPkg"
            ).run_test_conditions(
                self.test_conditions,
                self.last_solution,
                self.test_condition_iteration_limit,
            )
            if self.last_solution is None:
                logger.warning(
                    "no solution could be found that satisfied all specified test conditions in specified iterations!"
                )
                return None

        # Running binary check to reduce solution to minimal reaction solution
        if binary_check:
            self.last_solution = self.gfpkgmgr.getpkg(
                "GapfillingPkg"
            ).binary_check_gapfilling_solution()

        # Setting last solution data
        self.last_solution["media"] = media
        self.last_solution["target"] = target
        self.last_solution["minobjective"] = minimum_obj
        self.last_solution["binary_check"] = binary_check
        return self.last_solution
    
    def run_global_gapfilling(
        self,
        medias,
        targets,
        thresholds,
        binary_check=False,
        prefilter=True,
    ):
        """Run gapfilling on a single media condition to force the model to achieve a nonzero specified objective
        Parameters
        ----------
        medias : [MSMedia]
            Media in which the model should be gapfilled
        targets : [string]
            Name or expression describing the reaction or combination of reactions to the optimized
        thresholds : [double]
            Value to use for the minimal objective threshold that the model must be gapfilled to achieve
        binary_check : bool
            Indicates if the solution should be checked to ensure it is minimal in the number of reactions involved
        prefilter : bool
            Indicates if the gapfilling database should be prefiltered using the tests provided in the MSGapfill constructor before running gapfilling
        check_for_growth : bool
            Indicates if the model should be checked to ensure that the resulting gapfilling solution produces a nonzero objective
        """
        start_time = time.time()
        # Testing if gapfilling can work before filtering
        test_output = self.test_and_adjust_gapfilling_conditions(medias,targets,thresholds,prefilter=prefilter)
        #If none of the media conditions can be gapfilled, then return None
        if len(test_output["medias"]) == 0:
            return None
        #Adding max flux variables
        self.gfpkgmgr.getpkg("GapfillingPkg").create_max_flux_variables()
        #Instantiating all models to be merged
        merged_model = None
        model_list = []
        pkgmgrs = {}
        for i,media in enumerate(test_output["medias"]):
            #Setting the objective
            self.gfpkgmgr.getpkg("GapfillingPkg").set_base_objective(test_output["targets"][i],test_output["thresholds"][i])
            #Setting the media
            self.gfpkgmgr.getpkg("GapfillingPkg").set_media(media)
            #Copying model and either making it the base model or adding to the model list
            model_cpy = self.gfmodel.copy()
            
            if i == 0:
                merged_model = model_cpy
            else:
                model_list.append(model_cpy)
        #Merging all models
        mergpkgmgr = MSPackageManager.get_pkg_mgr(merged_model)
        mergpkgmgr.getpkg("ProblemReplicationPkg").build_package({
            "models":model_list,
            "shared_variable_packages":{
                "GapfillingPkg" : ["rmaxf","fmaxf"]
            }
        })
        mergfpkg = mergpkgmgr.getpkg("GapfillingPkg")
        origgfpkg = self.gfpkgmgr.getpkg("GapfillingPkg")
        #Setting the objective
        reaction_objective = merged_model.problem.Objective(Zero, direction="min")
        obj_coef = dict()
        gfrxnidhash = dict()
        for rxnid in mergfpkg.variables["rmaxf"]:
            gfrxnidhash[rxnid] = {"reverse":mergfpkg.variables["rmaxf"][rxnid]}
            if rxnid in origgfpkg.gapfilling_penalties:
                if "reverse" in origgfpkg.gapfilling_penalties[rxnid]:
                    obj_coef[mergfpkg.variables["rmaxf"][rxnid]] = abs(origgfpkg.gapfilling_penalties[rxnid]["reverse"])
                else:
                    obj_coef[mergfpkg.variables["rmaxf"][rxnid]] = 1
            else:
                obj_coef[mergfpkg.variables["rmaxf"][rxnid]] = 1
        for rxnid in mergfpkg.variables["fmaxf"]:
            if rxnid not in gfrxnidhash:
                gfrxnidhash[rxnid] = {"forward":mergfpkg.variables["fmaxf"][rxnid]}
            else:
                gfrxnidhash[rxnid]["forward"] = mergfpkg.variables["fmaxf"][rxnid]
            if rxnid in origgfpkg.gapfilling_penalties:
                if "forward" in origgfpkg.gapfilling_penalties[rxnid]:
                    obj_coef[mergfpkg.variables["fmaxf"][rxnid]] = abs(origgfpkg.gapfilling_penalties[rxnid]["forward"])
                else:
                    obj_coef[mergfpkg.variables["fmaxf"][rxnid]] = 1
            else:
                obj_coef[mergfpkg.variables["fmaxf"][rxnid]] = 1
        merged_model.objective = reaction_objective
        reaction_objective.set_linear_coefficients(obj_coef)
        # Printing the gapfilling LP file
        self.mdlutl.printlp(model=merged_model,filename="GlobalGapfill",print=True)

        # Running gapfilling and checking solution
        logger.debug("Starting global optimization- %s", time.time()-start_time)
        sol = merged_model.optimize()
        logger.debug("Global optimization complete- %s", time.time()-start_time)
        logger.debug(
            f"gapfill solution objective value {sol.objective_value} ({sol.status}) for media {media}"
        )
        if sol.status != "optimal":
            logger.warning("No solution found for %s", media)
            return None

        # Computing solution and ensuring all tests still pass
        self.last_solution = {"new":{},"reversed":{},"media":test_output["medias"][0],"target":test_output["targets"][0],"minobjective":test_output["thresholds"][0],"binary_check":False}
        flux_values = {}
        for rxnid in origgfpkg.gapfilling_penalties:
            flux_values[rxnid] = {}
            flux_values[rxnid]["reverse"] = merged_model.reactions.get_by_id(rxnid).reverse_variable.primal
            flux_values[rxnid]["forward"] = merged_model.reactions.get_by_id(rxnid).forward_variable.primal
        for rxnid in gfrxnidhash:
            if rxnid not in flux_values:
                flux_values[rxnid] = {}
            penalty = 0
            if "reverse" in gfrxnidhash[rxnid]:
                if rxnid in origgfpkg.gapfilling_penalties and "reverse" in origgfpkg.gapfilling_penalties[rxnid]:
                    penalty = origgfpkg.gapfilling_penalties[rxnid]["reverse"]
                if gfrxnidhash[rxnid]["reverse"].primal > 1e-8:
                    logger.debug(f"{rxnid} reverse {gfrxnidhash[rxnid]['reverse'].primal} {penalty}")
                flux_values[rxnid]["reverse"] = gfrxnidhash[rxnid]["reverse"].primal
            penalty = 0
            if "forward" in gfrxnidhash[rxnid]:
                if rxnid in origgfpkg.gapfilling_penalties and "forward" in origgfpkg.gapfilling_penalties[rxnid]:
                    penalty = origgfpkg.gapfilling_penalties[rxnid]["forward"]
                if gfrxnidhash[rxnid]["forward"].primal > 1e-8:
                    logger.debug(f"{rxnid} forward {gfrxnidhash[rxnid]['forward'].primal} {penalty}")
                flux_values[rxnid]["forward"] = gfrxnidhash[rxnid]["forward"].primal
        global_solution = origgfpkg.compute_gapfilled_solution(flux_values)
        logger.debug(f"Global solution: {global_solution}")
        logger.debug("Global gapfilling done - %s", time.time()-start_time)
        return global_solution

    def run_multi_gapfill(
        self,
        media_list,
        target=None,
        target_hash={},
        minimum_objectives={},
        default_minimum_objective=None,
        binary_check=False,
        prefilter=True,
        check_for_growth=True,
        gapfilling_mode="Sequential",
        run_sensitivity_analysis=True,
        integrate_solutions=True,
        remove_unneeded_reactions=True
    ):
        """Run gapfilling across an array of media conditions ultimately using different integration policies: simultaneous gapfilling, independent gapfilling, cumulative gapfilling
        Parameters
        ----------
        media_list : [MSMedia]
            List of the medias in which the model should be gapfilled
        target : string
            Name or expression describing the reaction or combination of reactions to the optimized
        minimum_objectives : {string - media ID : double - minimum objective value}
            Media-specific minimal objective thresholds that the model must be gapfilled to achieve
        default_minimum_objective : double
            Default value to use for the minimal objective threshold that the model must be gapfilled to achieve
        binary_check : bool
            Indicates if the solution should be checked to ensure it is minimal in the number of reactions involved
        prefilter : bool
            Indicates if the gapfilling database should be prefiltered using the tests provided in the MSGapfill constructor before running gapfilling
        check_for_growth : bool
            Indicates if the model should be checked to ensure that the resulting gapfilling solution produces a nonzero objective
        gapfilling_mode : string
            Indicates the integration policy to be used: Global, Independent, and Cumulative
        run_sensitivity_analysis : bool
            Indicates if sensitivity analysis should be run on the gapfilling solution to determine biomass dependency
        """
        #If not integrating, backing up and replacing self.mdlutl
        oldmdlutl = self.mdlutl
        if not integrate_solutions:
            self.model = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
            self.mdlutl = MSModelUtil.get(self.model)
        #Setting the default minimum objective
        if default_minimum_objective == None:
            default_minimum_objective = self.default_minimum_objective
        self.gfpkgmgr.getpkg("GapfillingPkg").parameters["minimum_obj"] = default_minimum_objective
        # Testing if gapfilling can work before and after filtering
        targets = []
        thresholds = []
        for media in media_list:
            currtarget = target
            if media in target_hash:
                currtarget = target_hash[media]
            targets.append(currtarget)
            minimum_obj = default_minimum_objective
            if media in minimum_objectives:
                minimum_obj = minimum_objectives[media]
            thresholds.append(minimum_obj)
        test_output = self.test_and_adjust_gapfilling_conditions(media_list,targets,thresholds,prefilter=prefilter)
        #If there are no media left, don't run gapfilling
        if len(test_output["medias"]) == 0:
            return None
        #Iterating over all media and running gapfilling
        solution_dictionary = {}
        cumulative_solution = []
        for i,media in enumerate(test_output["medias"]):
            #Implementing specified gapfilling mode
            if gapfilling_mode == "Independent" or gapfilling_mode == "Sequential":
                logger.debug("Running %s gapfilling!", gapfilling_mode)
                solution = self.run_gapfilling(
                    media,
                    test_output["targets"][i],
                    test_output["thresholds"][i],
                    binary_check,
                    False,
                )
                #If there is a solution, go ahead and integrate it into the model
                if solution:
                    solution_dictionary[media] = self.integrate_gapfill_solution(
                        solution,
                        cumulative_solution=cumulative_solution,
                        remove_unneeded_reactions=remove_unneeded_reactions,
                        check_for_growth=check_for_growth,
                        gapfilling_mode=gapfilling_mode
                    )
                    #If we are doing cumulative gapfilling, then we need adjust the gapfilling objective so it no longer penalizes using the current solution reactions
                    if gapfilling_mode == "Sequential":
                        self.gfpkgmgr.getpkg("GapfillingPkg").compute_gapfilling_penalties(exclusion_solution=cumulative_solution,reaction_scores=self.reaction_scores)
                        self.gfpkgmgr.getpkg("GapfillingPkg").build_gapfilling_objective_function()
        if gapfilling_mode == "Global":
            #Now we run simultaneous gapfilling on a combination of all our various gapfilled models
            logger.debug("Running global gapfilling!")
            full_solution = self.run_global_gapfilling(
                medias=test_output["medias"],
                targets=test_output["targets"],
                thresholds=test_output["thresholds"],
                binary_check=binary_check,
                prefilter=False
            )
            #Now we integrate the full solution into the model for every media which effectively determines which reactions are needed for each media
            for i,item in enumerate(test_output["medias"]):
                copy_solution = full_solution.copy()
                copy_solution["media"] = item
                copy_solution["target"] = test_output["targets"][i]
                copy_solution["minobjective"] = test_output["thresholds"][i]
                copy_solution["binary_check"] = binary_check
                #In this case we donot remove unnneeded reactions from the model because they may be needed for other media
                solution_dictionary[item] = self.integrate_gapfill_solution(
                    copy_solution,
                    cumulative_solution=cumulative_solution,
                    remove_unneeded_reactions=False,
                    check_for_growth=check_for_growth,
                    gapfilling_mode=gapfilling_mode
                )
            #Now we remove reactions uneeded for any of the specified media conditions
            #These is a danger here that the integration step will put a reaction into a solution that subsequently gets removed at this step. This is something to look out for
            unneeded = self.mdlutl.test_solution(
                cumulative_solution,
                test_output["targets"],
                test_output["medias"],
                thresholds=test_output["thresholds"],
                remove_unneeded_reactions=True,
                do_not_remove_list=[]
            )#Returns reactions in cumulative solution that are not needed for growth
            logger.debug("Unneeded in global gapfill: %s", unneeded)
        elif gapfilling_mode == "Sequential":
            #Restoring the gapfilling objective function
            self.gfpkgmgr.getpkg("GapfillingPkg").compute_gapfilling_penalties(reaction_scores=self.reaction_scores)
            self.gfpkgmgr.getpkg("GapfillingPkg").build_gapfilling_objective_function()
        #Running sensitivity analysis once on the cumulative solution for all media
        #with open("datacache/solutions.json", 'w') as f:
            #json.dump(solution_dictionary,f,indent=4,skipkeys=True)
        if run_sensitivity_analysis:
            logger.debug(
                "Gapfilling sensitivity analysis running"
            )
            #First aggregating all unique reactions with a media for each
            reaction_media_hash = {}
            solution_rxn_types = ["new","reversed"]
            media_reaction_hash = {}
            for media in solution_dictionary:
                if solution_dictionary[media]["growth"] > 0:
                    for rxn_type in solution_rxn_types:
                        for rxn_id in solution_dictionary[media][rxn_type]:
                            if rxn_id not in reaction_media_hash:
                                reaction_media_hash[rxn_id] = {}
                            if solution_dictionary[media][rxn_type][rxn_id] not in reaction_media_hash[rxn_id]:
                                reaction_media_hash[rxn_id][solution_dictionary[media][rxn_type][rxn_id]] = media
                                if media not in media_reaction_hash:
                                    media_reaction_hash[media] = {}
                                media_reaction_hash[media][rxn_id] = solution_dictionary[media][rxn_type][rxn_id]
            #Running sensitivity analysis on minimal reactions in each media
            rxn_sensitivity_hash = {}
            for media in media_reaction_hash:
                test_solution = []
                for rxn in media_reaction_hash[media]:
                    test_solution.append([rxn, media_reaction_hash[media][rxn]])
                self.mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                sensitivity_results = self.mdlutl.find_unproducible_biomass_compounds(
                    target, test_solution
                )
                for rxn in sensitivity_results:
                    if rxn not in rxn_sensitivity_hash:
                        rxn_sensitivity_hash[rxn] = {}
                    for dir in sensitivity_results[rxn]:
                        rxn_sensitivity_hash[rxn][dir] = sensitivity_results[rxn][dir]
            #Building gapfilling sensitivity output
            gf_sensitivity = self.mdlutl.get_attributes("gf_sensitivity", {})
            for media in solution_dictionary:
                if media.id not in gf_sensitivity:
                    gf_sensitivity[media.id] = {}
                if target not in gf_sensitivity[media.id]:
                    gf_sensitivity[media.id][target] = {}
                if solution_dictionary[media]["growth"] > 0:
                    gf_sensitivity[media.id][target]["success"] = {}
                    for rxn_type in solution_rxn_types:
                        for rxn_id in solution_dictionary[media][rxn_type]:
                            if rxn_id not in gf_sensitivity[media.id][target]["success"]:
                                gf_sensitivity[media.id][target]["success"][rxn_id] = {}
                            gf_sensitivity[media.id][target]["success"][rxn_id][solution_dictionary[media][rxn_type][rxn_id]] = rxn_sensitivity_hash[rxn_id][solution_dictionary[media][rxn_type][rxn_id]]
                else:
                    gf_sensitivity[media.id][target]["failure"] = {}
            self.mdlutl.save_attributes(gf_sensitivity, "gf_sensitivity")
        #Restoring backedup model
        self.mdlutl = oldmdlutl
        self.model = oldmdlutl.model
        #Returning the solution dictionary
        return solution_dictionary

    def integrate_gapfill_solution(
        self,solution,cumulative_solution=[],remove_unneeded_reactions=False,check_for_growth=True,gapfilling_mode="Sequential"
    ):
        """Integrating gapfilling solution into model
        Parameters
        ----------
        solution : dict
            Specifies the reactions to be added to the model to implement the gapfilling solution
        cumulative_solution : list
            Optional array to cumulatively track all reactions added to the model when integrating multiple solutions
        remove_unneeded_reactions : bool
            Indicate where unneeded reactions should be removed from the model
        check_for_growth : bool
            Indicate if the model should be checked to ensure that the resulting gapfilling solution produces a nonzero objective
        gapfilling_mode : Cumulative, Independent, Simultaneous
            Specify what the gapfilling mode is because this determines how integration is performed
        """
        logger.debug(f"Initial solution: {str(solution)}")
        original_objective = self.mdlutl.model.objective
        self.mdlutl.model.objective = solution["target"]
        self.mdlutl.model.objective.direction = "max"
        #If gapfilling mode is independent, we should remove the cumulative solution from the model before integrating the current solution
        if gapfilling_mode == "Independent":
            for item in cumulative_solution:
                rxn = self.model.reactions.get_by_id(item[0])
                if item[1] == ">":
                    rxn.upper_bound = 0
                else:
                    rxn.lower_bound = 0
        new_cumulative_reactions = []
        #Converting the solution to list
        list_solution = self.mdlutl.convert_solution_to_list(solution)
        for item in list_solution:
            if item[0] not in self.model.reactions:
                logger.debug(f"adding reaction: {str(item[0])}")
                #Copying and adding the reaction to the model
                rxn = self.gfmodel.reactions.get_by_id(item[0])
                rxn = rxn.copy()
                self.model.add_reactions([rxn])
                #Clearing current bounds because we only want to add reaction in the direction it was gapfilled in
                rxn.upper_bound = 0
                rxn.lower_bound = 0
            logger.debug(f"integrating rxn: {item[0]}")
            rxn = self.model.reactions.get_by_id(item[0])
            #Setting genes if the reaction has no genes
            if len(rxn.genes) == 0:
                #Setting genes from reaction scores in we have them
                coreid = re.sub(r"_[a-z]\d+$", "", item[0])
                if coreid in self.reaction_scores:
                    logger.debug(f"Found reaction scores for coreid: {coreid}")
                    bestgene = None
                    bestscore = None
                    for gene in self.reaction_scores[coreid]:
                        score = None
                        if isinstance(self.reaction_scores[coreid][gene], dict):
                            score = self.reaction_scores[coreid][gene]["probability"]
                        else:
                            score = self.reaction_scores[coreid][gene]
                        if (
                            not bestgene
                            or score
                            > bestscore
                        ):
                            bestgene = gene
                            bestscore = score
                    rxn = self.model.reactions.get_by_id(item[0])
                    logger.debug(f"Assigning gene to reaction: {item[0]} {bestgene}")
                    rxn.gene_reaction_rule = bestgene
                    rxn.notes["new_genes"] = bestgene
            #Setting bounds according to the direction the reaction was gapfilled in
            if item[1] == ">":
                rxn.upper_bound = 100
            else:
                rxn.lower_bound = -100
            #Adding reaction to cumulative solution if it is not already there
            if not self.mdlutl.find_item_in_solution(cumulative_solution,item):
                new_cumulative_reactions.append([item[0], item[1],item[2]])
        #Testing the full cumulative solution to see which reactions are needed for current media/target
        full_solution = cumulative_solution + new_cumulative_reactions
        logger.debug(f"Full solution: {str(full_solution)}")
        #Setting up structure to store the finalized solution for this media/target
        current_media_target_solution = {"growth":0,"media":solution["media"],"target":solution["target"],"minobjective":solution["minobjective"],"binary_check":solution["binary_check"] ,"new":{},"reversed":{}}
        #If gapfilling is independent, we only check the specific solution
        if gapfilling_mode == "Independent":
            unneeded = self.mdlutl.test_solution(list_solution,[solution["target"]],[solution["media"]],[solution["minobjective"]],remove_unneeded_reactions,do_not_remove_list=cumulative_solution)#Returns reactions in input solution that are not needed for growth
            for item in list_solution:
                if not self.mdlutl.find_item_in_solution(unneeded,item):
                    current_media_target_solution[item[2]][item[0]] = item[1]
                    if not self.mdlutl.find_item_in_solution(cumulative_solution,item): 
                        cumulative_solution.append(item)
                #elif not remove_unneeded_reactions and not self.mdlutl.find_item_in_solution(cumulative_solution,item):
                #    cumulative_solution.append(item)
            logger.debug(f"Cumulative media target solution: {str(current_media_target_solution)}")
        else:
            unneeded = self.mdlutl.test_solution(full_solution,[solution["target"]],[solution["media"]],[solution["minobjective"]],remove_unneeded_reactions,do_not_remove_list=cumulative_solution)#Returns reactions in input solution that are not needed for growth
            for item in cumulative_solution:
                if not self.mdlutl.find_item_in_solution(unneeded,item):
                    current_media_target_solution[item[2]][item[0]] = item[1]
            for item in new_cumulative_reactions:
                if not self.mdlutl.find_item_in_solution(unneeded,item):
                    current_media_target_solution[item[2]][item[0]] = item[1]
                    cumulative_solution.append(item)
                #elif not remove_unneeded_reactions:
                #    cumulative_solution.append(item)
        logger.debug(f"Unneeded: {str(unneeded)}")
        logger.debug(f"Cumulative: {str(self.cumulative_gapfilling)}")
        #Checking that the final integrated model grows
        if check_for_growth:
            self.mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(solution["media"])
            current_media_target_solution["growth"] = self.mdlutl.model.slim_optimize()
            logger.debug(f"Growth: {str(current_media_target_solution['growth'])} {solution['media'].id}")
        # Adding the gapfilling solution data to the model, which is needed for saving the model in KBase
        self.mdlutl.add_gapfilling(current_media_target_solution)
        # Testing which gapfilled reactions are needed to produce each reactant in the objective function
        self.cumulative_gapfilling.extend(cumulative_solution)
        return current_media_target_solution

    def compute_reaction_weights_from_expression_data(self, omics_data, annoont):
        """Computing reaction weights based on input gene-level omics data
        Parameters
        ----------
        omics_data : pandas dataframe with genes as rows and conditions as columns
            Specifies the reactions to be added to the model to implement the gapfilling solution
        annoont : annoont object
            Contains reaction, feature id, ontologies, probabilities. Restructured into dataframe in function
        Returns :
            A dictionary with Rxns as the keys and calculated result as the value.
        """

        ### Restructure annoont into Dataframe
        rows_list = []
        for reaction, genes in annoont.get_reaction_gene_hash(feature_type="gene").items():
            for gene, gene_info in genes.items():
                # Initialize the row with 'Gene' and 'Reactions'
                row = {"Gene": gene, "Reactions": reaction}
                # Loop through each evidence in the gene's evidence list
                for evidence in gene_info["evidence"]:
                    # Construct column name from the event and ontology for uniqueness
                    column_name = f"{evidence['ontology']}"
                    if column_name in row:
                        row[column_name] = f"{row[column_name]}, {evidence['term']}"
                    else:
                        row[column_name] = evidence["term"]
                rows_list.append(row)
        restructured_anoot = pd.DataFrame(rows_list)

        ### Integrate Omics, set weights, find indexes for features
        feature_ids_set = set(omics_data.index)

        # Find indices where 'Gene' values are in 'feature_ids'
        # isin method returns a boolean series that is True where tbl_supAno['Gene'] is in feature_ids_set
        mask = restructured_anoot["Gene"].isin(feature_ids_set)
        # Get the indices of True values in the mask
        idx_measuredGene = mask[mask].index.tolist()
        # Calculate the dimensions for the measuredGeneScore array
        num_genes = len(restructured_anoot["Gene"])
        num_columns = len(restructured_anoot.columns[2:])
        # Initialize the measuredGeneScore array with zeros
        measuredGeneScore = np.zeros((num_genes, num_columns))
        measuredGeneScore[idx_measuredGene, :] = 1
        num_weights = len(restructured_anoot.columns[3:])
        w = np.repeat(1 / num_weights, num_weights)

        ### Calculate Weights and generate the reaction/weight hash
        num_cols = len(restructured_anoot.columns[2:])
        w = np.full((num_cols, 1), 1 / num_cols)
        p = np.zeros(len(restructured_anoot["Reactions"]))
        # computed_weights is the rxn_hash ({rxn: weight, ...})
        computed_weights = {}

        # Precompute gene reaction lookups
        gene_reaction_lookup = {}
        for idx, row in restructured_anoot.iterrows():
            gene = row["Gene"]
            reaction = row["Reactions"]
            if gene in gene_reaction_lookup:
                gene_reaction_lookup[gene].append(reaction)
            else:
                gene_reaction_lookup[gene] = [reaction]

        for rxn in range(0, len(restructured_anoot)):
            substr_rxns = [rxn for rxn in restructured_anoot["Reactions"][[rxn]]]
            # Get the indices of the rows where the condition is True
            mask = restructured_anoot["Reactions"] == substr_rxns[0]
            idx_gene = mask[mask].index
            nAG = 0
            nMG = 0
            nCG = 0

            if len(idx_gene) > 0:
                # number of genes that map to a reaction
                nAG = len(idx_gene)
                for iGene in range(0, nAG):
                    subset = restructured_anoot.iloc[idx_gene[iGene], 2:].to_numpy()
                    # Checking for non-empty elements in the subset
                    non_empty_check = np.vectorize(lambda x: x is not None and x == x)(
                        subset
                    )
                    # Finding the maximum value between the non-empty check and the corresponding row in measuredGeneScore
                    max_value = np.maximum(
                        non_empty_check, measuredGeneScore[idx_gene[iGene], :]
                    )
                    # Multiplying by the weight and adding to nMG
                    nMG += max(sum((max_value * w)))
                    selected_gene = restructured_anoot["Gene"].iloc[idx_gene[iGene]]

                    # Finding reactions associated with genes that contain the selected gene
                    associated_reactions = gene_reaction_lookup.get(selected_gene, [])

                    # Checking if there are more than one unique reactions
                    if len(associated_reactions) > 1:
                        nCG += 1

                p[rxn] = (nMG / nAG) * (1 / (1 + (nCG / nAG)))

            # Add item to output rxn hash dictionary
            computed_weights[restructured_anoot.iloc[rxn, 0]] = p[rxn]

        return computed_weights

    @staticmethod
    def gapfill(
        model,
        media=None,
        target_reaction="bio1",
        default_gapfill_templates=[],
        default_gapfill_models=[],
        test_conditions=[],
        reaction_scores={},
        blacklist=[],
    ):
        gapfiller = MSGapfill(
            model,
            default_gapfill_templates,
            default_gapfill_models,
            test_conditions,
            reaction_scores,
            blacklist,
        )
        gfresults = gapfiller.run_gapfilling(media, target_reaction)
        return gapfiller.integrate_gapfill_solution(gfresults)
