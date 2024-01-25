# -*- coding: utf-8 -*-
import logging
import cobra
import re
import json
from optlang.symbolics import Zero, add
from modelseedpy.core import FBAHelper  # !!! the import is never used
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.exceptions import GapfillingError

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO  # WARNING
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
            "cpd02701",
            "cpd11416",
            "cpd15302",
            "cpd03091",
        ]  # the cpd11416 compound is filtered during model extension with templates
        # Cloning model to create gapfilling model
        self.gfmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
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
        self.lp_filename = self.last_solution = None
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

    def test_gapfill_database(self, media, target=None, before_filtering=True):
        # Testing if gapfilling can work before filtering
        if target:
            self.gfpkgmgr.getpkg("GapfillingPkg").set_base_objective(target,None)
        else:
            target = str(self.gfmodel.objective)
            target = target.split(" ")[0]
            target = target[13:]
        if self.gfpkgmgr.getpkg("GapfillingPkg").test_gapfill_database():
            return True
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

    def prefilter(self,test_conditions=None):
        """Prefilters the database by removing any reactions that break specified ATP tests
        Parameters
        ----------
        test_conditions : []
            List of conditions to be tested when filtering the gapfilling database. If not specified, the test_conditions attribute will be used
        """
        if not test_conditions:
            test_conditions = self.test_conditions
        if self.test_conditions:
            self.gfpkgmgr.getpkg("GapfillingPkg").filter_database_based_on_tests(
                self.test_conditions
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
            self.prefilter()
            if not self.test_gapfill_database(media,target,before_filtering=False):
                return None

        # Printing the gapfilling LP file
        if self.lp_filename:
            with open(self.lp_filename, "w") as out:
                out.write(str(self.gfmodel.solver))

        # Running gapfilling and checking solution
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
        # Testing if gapfilling can work before filtering
        final_media = []
        final_targets = []
        final_thresholds = []
        for i,media in enumerate(medias):
            if self.test_gapfill_database(media,targets[i],before_filtering=True):
                final_media.append(media)
                final_targets.append(targets[i])
                final_thresholds.append(thresholds[i])
        # Filtering
        if prefilter:
            self.prefilter()
            medias = []
            targets = []
            thresholds = []
            for i,media in enumerate(final_media):
                if self.test_gapfill_database(media,final_targets[i],before_filtering=True):
                    medias.append(media)
                    targets.append(targets[i])
                    thresholds.append(thresholds[i])
        #If none of the media conditions can be gapfilled, then return None
        if len(medias) == 0:
            return None
        #Instantiating all models to be merged
        merged_model = None
        model_list = []
        pkgmgrs = {}
        for i,media in enumerate(medias):
            model_cpy = self.gfmodel.copy()
            pkgmgrs[model_cpy] = MSPackageManager.get_pkg_mgr(model_cpy)
            #Creating max flux variables
            pkgmgrs[model_cpy].getpkg("GapfillingPkg").create_max_flux_variables()
            #Setting the objective
            pkgmgrs[model_cpy].getpkg("GapfillingPkg").set_base_objective(targets[i],thresholds[i])
            #Setting the media
            pkgmgrs[model_cpy].getpkg("GapfillingPkg").set_media(media)
            if i == 0:
                merged_model = model_cpy
            else:
                model_list.append(model_cpy)
        #Merging all models
        gfpkg = pkgmgrs[merged_model].getpkg("GapfillingPkg")
        pkgmgrs[merged_model].getpkg("ProblemReplicationPkg").build_package({
            "models":model_list,
            "shared_variable_packages":{
                gfpkg : ["rmaxf","fmaxf"]
            }
        })
        #Setting the objective
        reaction_objective = merged_model.problem.Objective(Zero, direction="min")
        obj_coef = dict()
        for reaction in merged_model.reactions:
            if reaction.id in gfpkg.gapfilling_penalties:
                if reaction.id[0:3] != "EX_":
                    if "reverse" in gfpkg.gapfilling_penalties[reaction.id]:
                        if reaction.id in gfpkg.maxflux_variables:
                            if "reverse" in gfpkg.maxflux_variables[reaction.id]:
                                obj_coef[gfpkg.maxflux_variables[reaction.id]["reverse"]] = abs(
                                    gfpkg.gapfilling_penalties[reaction.id]["reverse"]
                                )
                    if "forward" in gfpkg.gapfilling_penalties[reaction.id]:
                        if reaction.id in gfpkg.maxflux_variables:
                            if "forward" in gfpkg.maxflux_variables[reaction.id]:
                                obj_coef[gfpkg.maxflux_variables[reaction.id]["forward"]] = abs(
                                    gfpkg.gapfilling_penalties[reaction.id]["forward"]
                                )
        merged_model.objective = reaction_objective
        reaction_objective.set_linear_coefficients(obj_coef)
        gfpkg.parameters["gfobj"] = self.model.objective
        
        # Printing the gapfilling LP file
        if self.lp_filename:
            with open(self.lp_filename, "w") as out:
                out.write(str(merged_model.solver))

        # Running gapfilling and checking solution
        sol = merged_model.optimize()
        logger.debug(
            f"gapfill solution objective value {sol.objective_value} ({sol.status}) for media {media}"
        )
        if sol.status != "optimal":
            logger.warning("No solution found for %s", media)
            return None

        # Computing solution and ensuring all tests still pass
        self.last_solution = {"new":{},"reversed":{},"media":medias[0],"target":targets[0],"minobjective":thresholds[0],"binary_check":False}
        flux_values = {}
        for rxn in self.model.reactions:
            flux_values[rxn.id] = {
                "reverse":  self.gfpkgmgr.getpkg("GapfillingPkg").maxflux_variables[reaction.id]["reverse"].primal,
                "forward":  self.gfpkgmgr.getpkg("GapfillingPkg").maxflux_variables[reaction.id]["forward"].primal
            }
        self.gfpkgmgr.getpkg("GapfillingPkg").compute_gapfilled_solution(flux_values)
        return self.last_solution

    def run_multi_gapfill(
        self,
        media_list,
        target=None,
        minimum_objectives={},
        default_minimum_objective=None,
        binary_check=False,
        prefilter=True,
        check_for_growth=True,
        gapfilling_mode="Cumulative",
        run_sensitivity_analysis=True,
        integrate_solutions=True
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
        if not default_minimum_objective:
            default_minimum_objective = self.default_minimum_objective
        #Running prefiltering once for all media if specified. Rememeber - filtering does not care about the target or media - it is just a set of tests that are run on the database
        if prefilter:
            self.prefilter()
        #Iterating over all media and running gapfilling
        solution_dictionary = {}
        cumulative_solution = []
        targets = []
        thresholds = []
        for item in media_list:
            targets.append(target)
            #Determining the minimum objective for the current media
            minimum_obj = default_minimum_objective
            if item in minimum_objectives:
                minimum_obj = minimum_objectives[item]
            thresholds.append(minimum_obj)
            #Implementing specified gapfilling mode
            if gapfilling_mode == "Independent" or gapfilling_mode == "Cumulative":           
                self.lp_filename = item.id+"-gf.lp"
                solution = self.run_gapfilling(
                    item,
                    target,
                    minimum_obj,
                    binary_check,
                    False,
                )
                #If there is a solution, go ahead and integrate it into the model
                if solution:
                    solution_dictionary[item] = self.integrate_gapfill_solution(
                        solution,
                        cumulative_solution=cumulative_solution,
                        remove_unneeded_reactions=True,
                        check_for_growth=check_for_growth,
                        gapfilling_mode=gapfilling_mode
                    )
                    #If we are doing cumulative gapfilling, then we need adjust the gapfilling objective so it no longer penalizes using the current solution reactions
                    if gapfilling_mode == "Cumulative":
                        self.gfpkgmgr.getpkg("GapfillingPkg").compute_gapfilling_penalties(exclusion_solution=cumulative_solution,reaction_scores=self.reaction_scores)
                        self.gfpkgmgr.getpkg("GapfillingPkg").build_gapfilling_objective_function()
        if gapfilling_mode == "Global":
            #Now we run simultaneous gapfilling on a combination of all our various gapfilled models
            full_solution = self.run_global_gapfilling(
                media_list,
                targets,
                thresholds,
                binary_check,
                False,
                check_for_growth,
            )
            #Now we integrate the full solution into the model for every media which effectively determines which reactions are needed for each media
            for i,item in enumerate(media_list):
                full_solution["media"] = item
                full_solution["target"] = targets[i]
                full_solution["minobjective"] = thresholds[i]
                #In this case we donot remove unnneeded reactions from the model because they may be needed for other media
                solution_dictionary[item] = self.integrate_gapfill_solution(
                    full_solution,
                    cumulative_solution=cumulative_solution,
                    remove_unneeded_reactions=False,
                    check_for_growth=check_for_growth,
                    gapfilling_mode=gapfilling_mode
                )
            #Now we remove reactions uneeded for any of the specified media conditions
            #These is a danger here that the integration step will put a reaction into a solution that subsequently gets removed at this step. This is something to look out for
            unneeded = self.mdlutl.test_solution(
                cumulative_solution,
                targets,
                media_list,
                thresholds=[0.1],
                remove_unneeded_reactions=True,
                do_not_remove_list=[]
            )#Returns reactions in cumulative solution that are not needed for growth
        elif gapfilling_mode == "Cumulative":
            #Restoring the gapfilling objective function
            self.gfpkgmgr.getpkg("GapfillingPkg").compute_gapfilling_penalties(reaction_scores=self.reaction_scores)
            self.gfpkgmgr.getpkg("GapfillingPkg").build_gapfilling_objective_function()
        #Running sensitivity analysis once on the cumulative solution for all media
        """ 
        if run_sensitivity_analysis:#TODO - need to redesign this run operate on multiple solutions at once...
            logger.info(
                "Gapfilling sensitivity analysis running on succesful run in "
                + solution["media"].id
                + " for target "
                + solution["target"]
            )
            test_solution = []
            for rxn_id in solution["reversed"]:
                test_solution.append([rxn_id, solution["reversed"][rxn_id]])
            for rxn_id in solution["new"]:
                test_solution.append([rxn_id, solution["new"][rxn_id]])
            gf_sensitivity = self.mdlutl.get_attributes("gf_sensitivity", {})
            if solution["media"].id not in gf_sensitivity:
                gf_sensitivity[solution["media"].id] = {}
            if solution["target"] not in gf_sensitivity[solution["media"].id]:
                gf_sensitivity[solution["media"].id][solution["target"]] = {}
            gf_sensitivity[solution["media"].id][solution["target"]][
                "success"
            ] = self.mdlutl.find_unproducible_biomass_compounds(
                solution["target"], test_solution
            )
            self.mdlutl.save_attributes(gf_sensitivity, "gf_sensitivity") """
        #Restoring backedup model
        self.mdlutl = oldmdlutl
        self.model = oldmdlutl.model
        #Returning the solution dictionary
        return solution_dictionary

    def integrate_gapfill_solution(
        self,solution,cumulative_solution=[],remove_unneeded_reactions=False,check_for_growth=True,gapfilling_mode="Cumulative"
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
        print("Initial solution:",str(solution))
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
                #Copying and adding the reaction to the model
                rxn = self.gfmodel.reactions.get_by_id(item[0])
                rxn = rxn.copy()
                self.model.add_reactions([rxn])
                #Clearing current bounds because we only want to add reaction in the direction it was gapfilled in
                rxn.upper_bound = 0
                rxn.lower_bound = 0
                #Setting genes from reaction scores in we have them
                coreid = re.sub(r"_[a-z]\d+$", "", item[0])
                if coreid in self.reaction_scores:
                    bestgene = None
                    for gene in self.reaction_scores[coreid]:
                        if (
                            not bestgene
                            or self.reaction_scores[coreid][gene]
                            > self.reaction_scores[coreid][bestgene]
                        ):
                            bestgene = gene
                    rxn = self.model.reactions.get_by_id(item[0])
                    rxn.gene_reaction_rule = bestgene
            rxn = self.model.reactions.get_by_id(item[0])
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
        print("Unneeded:",str(unneeded))
        #Checking that the final integrated model grows
        if check_for_growth:
            self.mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(solution["media"])
            current_media_target_solution["growth"] = self.mdlutl.model.slim_optimize()
            print("Growth:",str(current_media_target_solution["growth"]),solution["media"].id)
        # Adding the gapfilling solution data to the model, which is needed for saving the model in KBase
        self.mdlutl.add_gapfilling(current_media_target_solution)
        #If gapfilling mode is independent, restoring cumulative solution to the model
        if gapfilling_mode == "Independent":
            for item in cumulative_solution:
                rxn = self.model.reactions.get_by_id(item[0])
                if item[1] == ">":
                    rxn.upper_bound = 100
                else:
                    rxn.lower_bound = -100
        #Restoring the original objective
        self.mdlutl.model.objective = original_objective
        #Return the stripped down solution object containing only needed reactions
        return current_media_target_solution

    def compute_reaction_weights_from_expression_data(self, omics_data, conditions=[]):
        """Computing reaction weights based on input gene-level omics data
        Parameters
        ----------
        omics_data : pandas dataframe with genes as rows and conditions as columns
            Specifies the reactions to be added to the model to implement the gapfilling solution
        conditions : list
            Optional array containing the IDs of the columns in omics_data from which data should be used.
            If an empty array (or no array) is supplied, data from all columns will be used. When multiple columns are
            used, the data from those columns should be normalized first, then added together
        """
        # Validitions:
        # 1.) An conditions listed in the conditions argument should match the columns in the omics_data dataframe
        # 2.) Most (~80%) of the genes in the model should match genes in the omics_data dataframe
        # 3.) The omics_data dataframe should have at least 2 columns
        # 4.) The omics_data dataframe should have at least 2 rows
        # 5.) Logging should be used to report out which genes in the model don't match any genes in the omics_data dataframe
        pass

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
