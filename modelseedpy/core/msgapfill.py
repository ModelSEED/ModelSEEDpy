#!/usr/bin/python
# -*- coding: utf-8 -*-
import logging
import cobra
import re
import json
import numpy as np
import pandas as pd
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
        default_uptake=100,
        default_target=None,
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
            }
        )

    def test_gapfill_database(self, media, target=None, before_filtering=True):
        # Testing if gapfilling can work before filtering
        if target:
            self.gfmodel.objective = self.gfmodel.problem.Objective(
                self.gfmodel.reactions.get_by_id(target).flux_expression,
                direction="max",
            )
            self.gfpkgmgr.getpkg("GapfillingPkg").reset_original_objective()
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

    def prefilter(self, media, target):
        # Filtering breaking reactions out of the database
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

        # Testing if gapfilling can work after filtering
        if not self.test_gapfill_database(media, target, before_filtering=False):
            return False
        return True

    def run_gapfilling(
        self,
        media=None,
        target=None,
        minimum_obj=None,
        binary_check=False,
        prefilter=True,
        check_for_growth=True,
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
        check_for_growth : bool
            Indicates if the model should be checked to ensure that the resulting gapfilling solution produces a nonzero objective
        """
        # Setting target and media if specified
        if target:
            self.gfmodel.objective = self.gfmodel.problem.Objective(
                self.gfmodel.reactions.get_by_id(target).flux_expression,
                direction="max",
            )
            self.gfpkgmgr.getpkg("GapfillingPkg").reset_original_objective()
        else:
            target = self.default_target
        if media:
            self.gfpkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if not minimum_obj:
            minimum_obj = self.default_minimum_objective
        if minimum_obj:
            self.gfpkgmgr.getpkg("GapfillingPkg").set_min_objective(minimum_obj)

        # Testing if gapfilling can work before filtering
        if not self.test_gapfill_database(media, before_filtering=True):
            return None

        # Filtering
        if prefilter:
            if not self.prefilter(media, target):
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
        self.last_solution["minobjective"] = self.gfpkgmgr.getpkg(
            "GapfillingPkg"
        ).parameters["minimum_obj"]
        self.last_solution["binary_check"] = binary_check
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
        simultaneous_gapfilling=False,
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
        """

        if not default_minimum_objective:
            default_minimum_objective = self.default_minimum_objective
        solution_dictionary = {}
        if simultaneous_gapfilling:
            for item in media_list:
                pass
        else:
            first = True
            for item in media_list:
                minimum_obj = default_minimum_objective
                if item in minimum_objectives:
                    minimum_obj = minimum_objectives[item]
                if first:
                    solution_dictionary[item] = self.run_gapfilling(
                        item,
                        target,
                        minimum_obj,
                        binary_check,
                        prefilter,
                        check_for_growth,
                    )
                else:
                    solution_dictionary[item] = self.run_gapfilling(
                        item, None, minimum_obj, binary_check, False, check_for_growth
                    )
                false = False
            return solution_dictionary

    def integrate_gapfill_solution(
        self, solution, cumulative_solution=[], link_gaps_to_objective=True
    ):
        """Integrating gapfilling solution into model
        Parameters
        ----------
        solution : dict
            Specifies the reactions to be added to the model to implement the gapfilling solution
        cumulative_solution : list
            Optional array to cumulatively track all reactions added to the model when integrating multiple solutions
        """
        for rxn_id in solution["reversed"]:
            rxn = self.model.reactions.get_by_id(rxn_id)
            if solution["reversed"][rxn_id] == ">" and rxn.upper_bound <= 0:
                cumulative_solution.append([rxn_id, ">"])
                rxn.upper_bound = 100
            elif solution["reversed"][rxn_id] == "<" and rxn.lower_bound >= 0:
                cumulative_solution.append([rxn_id, "<"])
                rxn.lower_bound = -100
        for rxn_id in solution["new"]:
            if rxn_id not in self.model.reactions:
                rxn = self.gfmodel.reactions.get_by_id(rxn_id)
                rxn = rxn.copy()
                self.model.add_reactions([rxn])
                coreid = re.sub(r"_[a-z]\d+$", "", rxn_id)
                if coreid in self.reaction_scores:
                    bestgene = None
                    for gene in self.reaction_scores[coreid]:
                        if (
                            not bestgene
                            or self.reaction_scores[coreid][gene]
                            > self.reaction_scores[coreid][bestgene]
                        ):
                            bestgene = gene
                    rxn = self.model.reactions.get_by_id(rxn_id)
                    rxn.gene_reaction_rule = bestgene
                if solution["new"][rxn_id] == ">":
                    cumulative_solution.append([rxn_id, ">"])
                    rxn.upper_bound = 100
                    rxn.lower_bound = 0
                else:
                    cumulative_solution.append([rxn_id, "<"])
                    rxn.upper_bound = 0
                    rxn.lower_bound = -100

        # Sometimes for whatever reason, the solution includes useless reactions that should be stripped out before saving the final model
        unneeded = self.mdlutl.test_solution(
            solution, keep_changes=True
        )  # Strips out unneeded reactions - which undoes some of what is done above
        for item in unneeded:
            for oitem in cumulative_solution:
                if item[0] == oitem[0] and item[1] == oitem[1]:
                    cumulative_solution.remove(oitem)
                    break
        # Adding the gapfilling solution data to the model, which is needed for saving the model in KBase
        self.mdlutl.add_gapfilling(solution)
        # Testing which gapfilled reactions are needed to produce each reactant in the objective function
        if link_gaps_to_objective:
            logger.info(
                "Gapfilling sensitivity analysis running on succesful run in "
                + solution["media"].id
                + " for target "
                + solution["target"]
            )
            gf_sensitivity = self.mdlutl.get_attributes("gf_sensitivity", {})
            if solution["media"].id not in gf_sensitivity:
                gf_sensitivity[solution["media"].id] = {}
            if solution["target"] not in gf_sensitivity[solution["media"].id]:
                gf_sensitivity[solution["media"].id][solution["target"]] = {}
            gf_sensitivity[solution["media"].id][solution["target"]][
                "success"
            ] = self.mdlutl.find_unproducible_biomass_compounds(
                solution["target"], cumulative_solution
            )
            self.mdlutl.save_attributes(gf_sensitivity, "gf_sensitivity")
        self.cumulative_gapfilling.extend(cumulative_solution)

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
        Returns :
            A dictionary with Rxns as the keys and calculated result as the value.
        """
        # Validitions:
        # 1.) An conditions listed in the conditions argument should match the columns in the omics_data dataframe
        # 2.) Most (~80%) of the genes in the model should match genes in the omics_data dataframe
        # 3.) The omics_data dataframe should have at least 2 columns
        # 4.) The omics_data dataframe should have at least 2 rows
        # 5.) Logging should be used to report out which genes in the model don't match any genes in the omics_data dataframe

        # Assumptions:
        # omics_data is an input with the following columns:
        # Col 1: Gene
        # Col 2: Reactions
        # Cols 3-9: Annotations
        # Unlike with the MATLAB code, this will not check if rxn is not in draft - this is only making calculations based on the columns in omics_data.

        # Notes:
        # This has only been tested on dummy data and in this case python performs the same as MATLAB. However, this needs to be tested with an example input.
        # Two outputs are currently created.
        # A table that matches with t in the MATLAB code.
        # The requested hash (dictionary) with Rxns and floats.
        # Many more inputs were required with the MATLAB code, this has attempted to condense everything to handle the single input.
        # Conditions has not been added yet. I think as other questions are answered, the use for this will be more clear.

        # Questions
        # When might the example data be updated?
        # Are other inputs required?
        # Is this output (Dictionary with RXN: Value), correct?
        # When will the next steps for setting up the kbase jupyter notebook be ready?

        measuredGeneScore = np.zeros(
            (omics_data.shape[0], len(omics_data.columns[3:10]))
        )
        num_cols = len(omics_data.columns[3:10])
        w = np.full((num_cols, 1), 1 / num_cols)
        p = np.zeros(len(omics_data["Reactions"]))

        # t is table to match and check against MatLab code.
        t = pd.DataFrame()
        # rxn_hash is the desired output
        rxn_hash = {}

        for rxn in range(0, len(omics_data)):
            substr_rxns = [rxn for rxn in omics_data["Reactions"][[rxn]]]
            # Get the indices of the rows where the condition is True
            mask = omics_data["Reactions"].apply(
                lambda x: any(substr in x for substr in substr_rxns)
            )
            idx_gene = mask[mask].index
            nAG = 0
            nMG = 0
            nCG = 0

            if len(idx_gene) > 0:
                # number of genes that map to a reaction
                nAG = len(idx_gene)
                for iGene in range(0, nAG):
                    subset = omics_data.iloc[idx_gene[iGene], 3:9].to_numpy()
                    # Checking for non-empty elements in the subset
                    non_empty_check = np.vectorize(lambda x: x is not None and x == x)(
                        subset
                    )  # x == x checks for NaN
                    # Finding the maximum value between the non-empty check and the corresponding row in measuredGeneScore
                    max_value = np.maximum(
                        non_empty_check, measuredGeneScore[idx_gene[iGene], :]
                    )
                    # Multiplying by the weight and adding to nMG
                    nMG += max(sum((max_value * w)))
                    selected_gene = omics_data["Gene"].iloc[idx_gene[iGene]]

                    # Finding reactions associated with genes that contain the selected gene
                    associated_reactions = omics_data["Reactions"][
                        omics_data["Gene"].str.contains(selected_gene)
                    ]
                    # Checking if there are more than one unique reactions
                    if len(associated_reactions.unique()) > 1:
                        nCG += 1

                p[rxn] = (nMG / nAG) * (1 / (1 + (nCG / nAG)))

            # format table
            new_row = {
                "iRxn": rxn,
                "nMG": nMG,
                "nCG": nCG,
                "nAG": nAG,
                "Values": p[rxn],  # Values is equivalent to Var5 in the MatLab Code
            }

            # Append the new row to the table
            t = t.append(new_row, ignore_index=True)
            # Add item to output rxn dictionary
            rxn_hash[omics_data.iloc[rxn, 0]] = p[rxn]

        return rxn_hash

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
