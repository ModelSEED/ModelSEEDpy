from modelseedpy.core.exceptions import ObjectiveError, FeasibilityError
from modelseedpy.fbapkg.reactionusepkg import ReactionUsePkg
from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from itertools import combinations, permutations
from optlang import Variable, Constraint
from optlang.symbolics import Zero
from math import isclose, inf
from deepdiff import DeepDiff
from time import process_time
from pprint import pprint
import logging
import json, re

logger = logging.getLogger(__name__)


class MinimalMediaPkg:
    
    @staticmethod
    def _exchange_solution(sol_dict):
        return {rxn:flux for rxn, flux in sol_dict.items() if "EX_" in rxn.id and flux < 0}

    @staticmethod
    def _influx_objective(model):
        influxes = []
        for ex_rxn in FBAHelper.exchange_reactions(model):
            if len(ex_rxn.reactants) == 1:  # this is essentially 100% of exchanges
                influxes.append(ex_rxn.reverse_variable)
            else:  # this captures edge cases
                influxes.append(ex_rxn.forward_variable)
        return influxes

    @staticmethod
    def _varName_to_ID(var):
        rxnID = var.name
        if "_ru" in var:
            rxnID = rxnID.replace("_ru", "")
        return rxnID
    
    @staticmethod
    def minimize_flux(org_model, minimal_growth=None, printing=True):
        """minimize the total in-flux of exchange reactions in the model"""
        model = org_model.copy()
        minimal_growth = minimal_growth or model.optimize().objective_value
        FBAHelper.add_minimal_objective_cons(model, min_value=minimal_growth)
        FBAHelper.add_objective(model, sum(MinimalMediaPkg._influx_objective(model)), "min")
        sol = model.optimize()
        min_flux_media = MinimalMediaPkg._exchange_solution(FBAHelper.solution_to_rxns_dict(sol, model))
        total_flux = sum([abs(ex.flux) for ex in min_flux_media])
        # verify the medium
        model2 = org_model.copy()
        model2.medium = {rxn.id: abs(flux) for rxn, flux in min_flux_media.items()}
        simulated_sol = model2.optimize()
        if simulated_sol.status != "optimal":
            raise FeasibilityError(f"The simulation was not optimal, with a status of {simulated_sol.status}")
        if not isclose(simulated_sol.objective_value, minimal_growth):
            raise ObjectiveError(f"The assigned minimal_growth of {minimal_growth} was not maintained during the simulation,"
                                 f" where the observed growth value was {simulated_sol.objective_value}.")
        if printing:
            print(f"The minimal flux media consists of {len(min_flux_media)} compounds and a {total_flux} total influx,"
                  f" with a growth value of {simulated_sol.objective_value}")
        return min_flux_media

    @staticmethod
    def minimize_components(org_model, minimal_growth=None, printing=True):
        """minimize the quantity of metabolites that are consumed by the model"""
        model = org_model.copy()
        variables = {"ru":{}}
        FBAHelper.add_minimal_objective_cons(
            model, sum([rxn.flux_expression for rxn in model.reactions if "bio" in rxn.id]), minimal_growth)

        # define the binary variable and constraint
        for ex_rxn in FBAHelper.exchange_reactions(model):
            # define the variable
            variables["ru"][ex_rxn.id] = Variable(ex_rxn.id+"_ru", lb=0, ub=1, type="binary")
            model.add_cons_vars(variables["ru"][ex_rxn.id])
            # bin_flux: {rxn_bin}*1000 >= {rxn_rev_flux}
            FBAHelper.create_constraint(model, Constraint(Zero, lb=0, ub=None, name=ex_rxn.id+"_bin"),
                                            coef={variables["ru"][ex_rxn.id]: 1000, ex_rxn.reverse_variable: -1})
        FBAHelper.add_objective(model, sum([var for var in variables["ru"].values()]), "min")

        # determine each solution
        solution_dicts = []
        sol = model.optimize()
        while sol.status == "optimal":  # limit to 100
            print("sol_dicts length", len(solution_dicts), end="\r")
            sol_dict = FBAHelper.solution_to_variables_dict(sol, model)
            solution_dicts.append(sol_dict)
            ## omit the solution from the next search
            FBAHelper.create_constraint(model, Constraint(  # build exclusion use can be emulated
                Zero, lb=None, ub=len(sol_dict)-1,name=ex_rxn.id + f"_exclusion_sol{len(solution_dicts)}"),
                sol_dict)
            sol = model.optimize()
        if not solution_dicts:
            logger.error("No simulations were feasible.")

        # search the permutation space by omitting previously investigated solution_dicts
        interdependencies = {}
        for sol_index, sol_dict in enumerate(solution_dicts):
            sol_exchanges = [rxn for rxn in sol_dict if "EX_" in rxn.name]
            print(sol_exchanges)
            interdependencies[sol_index] = MinimalMediaPkg._examine_permutations(model, sol_exchanges, variables, sol_dict, sol_index)

    @staticmethod
    def _knockout(org_model, exVar, variables, sol_dict, sol_index):
        # knockout the specified exchange
        knocked_model = org_model.copy()
        exID = MinimalMediaPkg._varName_to_ID(exVar.name)
        coef = {variables["ru"][exID]: 0}
        coef.update({variables["ru"][MinimalMediaPkg._varName_to_ID(exVar2.name)]: 1
                     for exVar2 in sol_dict if exVar != exVar2 and "EX_" in exVar2.name})
        FBAHelper.create_constraint(knocked_model, Constraint(Zero, lb=0.1, ub=None, name=f"{exVar.name}-sol{sol_index}"), coef)
        return knocked_model.optimize()

    @staticmethod
    def _examine_permutations(model, exchange_ids_to_explore, variables, sol_dict, sol_index):
        for exID in exchange_ids_to_explore:
            sol_dict_sans_ex = sol_dict.copy()
            sol_dict_sans_ex.pop(exID)
            # interdependencies[sol_index][exID] = MinimalMediaPkg._examine_permutations(
            #     exID, sol_dict, sol_index, variables, sol_dict_sans_ex)
            interdependencies = {}
            
            ## explore permutations after removing the selected variable
            diff = DeepDiff(sol_dict_sans_ex, FBAHelper.solution_to_dict(
                MinimalMediaPkg._knockout(model, exID, variables, sol_dict, sol_index)))
            if diff:  # the addition of new exchanges or altered exchange fluxes are detected after the removed exchange
                for key, value in diff.items():
                    print(key, value)
                    # new_mets = {re.search("(?<=[\')(.+)(?=\'])", met): change for met, change in value.items()}
                    # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                    interdependencies[sol_index][exID].update(new_mets)
                    MinimalMediaPkg._examine_permutations(model, exchange_ids_to_explore, variables, sol_dict, sol_index)
                # coef = {variables["met"][exID]: 0 for cpd in new_mets.keys()}
                # coef.update({variables["met"][exID]: 1 for exID in sol_dict if exID not in new_mets.keys()})
                cpd_name = "_".join(new_mets.keys())
                BaseFBAPkg.build_constraint(self, "met", 0.1, None, coef, f"{cpd_name}-sol{sol_index}")
                new_sol = self.model.optimize()
                if new_sol.status != "optimal":
                    return interdependencies
                MinimalMediaPkg._examine_permutations(exID, new_sol, sol_index, sol_dict_sans_ex)
            return interdependencies

    @staticmethod
    def jenga_method(org_models, com_model=None, syntrophy=True, min_growth=0.1, conserved_cpds=[], export=True, printing=True,
                     compatibilize=False):
        from cobra.medium import minimal_medium

        # determine the unique combination of all species minimal media
        models = org_models
        if compatibilize:
            models = MSCompatibility.standardize(org_models, conflicts_file_name="standardization_corrections.json", printing=printing)
        media = {"community_media": {}, "members": {}}
        for org_model in models:
            model = org_model.copy()
            print(model.id)
            media["members"][model.id] = {"media": minimal_medium(org_model, min_growth, minimize_components=True).to_dict()}
            model.medium = media["members"][model.id]["media"]
            media["members"][model.id]["solution"] = FBAHelper.solution_to_dict(model.optimize())
            media["community_media"] = FBAHelper.sum_dict(model.medium, media["community_media"])

        # subtract syntrophic interactions and remove satisfied fluxes
        org_media = media["community_media"].copy()
        original_time = syntrophic_time = process_time()
        if printing:
            print(f"Initial media defined with {len(media['community_media'])} exchanges")
        changed = 0
        syntrophic_media = media["community_media"].copy()
        if syntrophy:
            for model in models:
                for rxnID, flux in media["members"][model.id]["solution"].items():
                    if rxnID in media["community_media"] and flux > 0:
                        stoich = list(model.reactions.get_by_id(rxnID).metabolites.values())[0]
                        media["community_media"][rxnID] -= flux * stoich
                        changed += 1
            media["community_media"] = {ID: flux for ID, flux in media["community_media"].items() if flux > 0}

            syntrophic_media = media["community_media"].copy()
            syntrophy_diff = DeepDiff(org_media, syntrophic_media)
            changed_quantity = 0 if not syntrophy_diff else len(list(syntrophy_diff.values())[0].values())
            syntrophic_time = process_time()
            if printing:
                print(
                    f"Syntrophic fluxes examined after {(syntrophic_time - original_time) / 60} minutes, with {changed_quantity} change(s):",
                    syntrophy_diff)

        # JENGA method of further reduction
        changed = 0
        if com_model:
            community_model = com_model
            ## identify additionally redundant compounds
            redundant_cpds = set()
            community_model.medium = media["community_media"]
            original_obj_value = com_model.optimize().objective_value
            for cpd in media["community_media"]:
                new_media = media["community_media"].copy()
                new_media.pop(cpd)
                community_model.medium = new_media
                sol = community_model.optimize()
                if isclose(sol.objective_value, original_obj_value, abs_tol=1e-4):
                    redundant_cpds.add(cpd)

            ## vet the permutations
            permuts = [p for p in permutations(redundant_cpds)]
            if printing:
                print(
                    f"The {len(permuts)} permutations of the {redundant_cpds} redundant compounds, from absolute tolerance of 1e-4, will be examined.")
            permutation_results = []
            best = 0
            failed_permutation_starts = []
            for perm_index, permut in enumerate(permuts):
                print(f"{perm_index}/{len(permuts)}", end="\r")
                permutation_segments = [permut[:index] for index in range(len(permut), 2, -1)]
                ### eliminate previously discovered failures and successes, respectively
                if not any([seg in failed_permutation_starts for seg in permutation_segments]):
                    if best < 3 or not any([set(permut[:best - 1]) == set(list(success)[:best - 1]) for success in permutation_results]):
                        successful_removal = 0
                        new_media = media["community_media"].copy()
                        for cpd in permut:
                            ### parameterize and simulate the community
                            new_media.pop(cpd)
                            community_model.medium = new_media
                            sol = community_model.optimize()
                            if isclose(sol.objective_value, original_obj_value, abs_tol=1e-7):
                                successful_removal += 1
                                continue
                            # print("objective value discrepancy:", sol.objective_value, original_obj_value)
                            failed_permutation_starts.append(permut[:successful_removal + 1])
                            break
                        if successful_removal >= best:
                            if successful_removal > best:
                                best = successful_removal
                                permutation_results = []
                            permut_removable = permut[:best]  # slice only the elements that are removable
                            if permut_removable not in permutation_results:
                                permutation_results.append(permut_removable)
                                if printing:
                                    print(permut_removable)

            ## filter to only the most minimal media
            new_combinations = []
            for permut in permutation_results:
                if permut not in permutation_results:
                    new_combinations.append(permut)

            unique_combinations, unique_paths = [], []
            for removal_path in permutation_results:
                path_permutations = permutations(removal_path)
                if all([set(path) in permutation_results for path in path_permutations]):
                    combination = combinations(removal_path, len(removal_path))
                    for com in combination:
                        if com not in unique_combinations:
                            unique_combinations.append(com)
                else:
                    unique_paths.append(removal_path)
            if unique_combinations and printing:
                print("Unique combinations:")
                print(len(unique_combinations), unique_combinations)
            if unique_paths and printing:
                print("Unique paths:")
                print(len(unique_paths), unique_paths)

            # further remove compounds from the media, while defaulting to the removal with the largest ID values
            possible_removal_tracker = {}
            possible_options = unique_combinations + unique_paths
            if conserved_cpds:
                possible_options = [opt for opt in possible_options if not any(cpd in conserved_cpds for cpd in opt)]
            best = -inf
            for possible_removal in possible_options:
                cpdID_sum = sum([int(cpd.split('_')[1].replace("cpd", "")) for cpd in possible_removal])
                if cpdID_sum > best:
                    best = cpdID_sum
                    possible_removal_tracker = {best: [possible_removal]}
                elif cpdID_sum == best:
                    possible_removal_tracker[best].append(possible_removal)
            media["community_media"] = FBAHelper.remove_media_compounds(
                media["community_media"], list(possible_removal_tracker.values())[0][0], printing)
            if printing:
                print(possible_removal_tracker)
                pprint(media["community_media"])

        jenga_media = media["community_media"].copy()
        jenga_time = process_time()
        jenga_difference = DeepDiff(syntrophic_media, jenga_media)
        changed_quantity = 0 if not jenga_difference else len(list(jenga_difference.values())[0])
        if printing:
            print(f"Jenga fluxes examined after {(jenga_time - syntrophic_time) / 60} minutes, with {changed_quantity} change(s):",
                  jenga_difference)
        if export:
            if not com_model:
                export_name = "_".join([model.id for model in models]) + "_media.json"
            else:
                export_name = com_model.id + "_media.json"
            with open(export_name, 'w') as out:
                json.dump(media, out, indent=3)
        return media
        

class MinimalMedia(BaseFBAPkg):
    """A class that determines the minimal media of a model"""  # investigate conversion to a staticmethod or single function for convenience and in-line utility
    def __init__(self, model, min_growth=0.1):
        # define the system
        BaseFBAPkg.__init__(self, model, "Minimal Media", {"met":"metabolite"}, {"met":"string", "obj":"string"})
        
        # define the exchange variables, the minimal objective, and the minimal growth value
        for cpd in FBAHelper.exchange_reactions(self.model):
            BaseFBAPkg.build_variable(self,"met",0,1,"binary",cpd)
        self.model = FBAHelper.add_objective(self.model, sum([var for var in self.variables["met"].values()]), "min")
        BaseFBAPkg.build_constraint(self, "obj", min_growth, None, {
            rxn.forward_variable:1 for rxn in FBAHelper.bio_reactions(self.model)}, "min_growth")
        
        # determine the set of media solutions
        solutions = []
        sol = self.model.optimize()
        while sol.status == "optimal":
            solutions.append(sol)
            sol_dict = FBAHelper.solution_to_variables_dict(sol, model)
            ## omit the solution from the next search
            BaseFBAPkg.build_constraint(self, "met", len(sol_dict)-1, len(sol_dict)-1, 
                                        coef=sol_dict, cobra_obj=f"sol{len(solutions)}")
            sol = self.model.optimize()
        if not solutions:
            logger.error("No simulations were feasible.")
                
        # search the permutation space by omitting previously investigated solutions
        self.interdependencies = {}
        for sol_index, sol in enumerate(solutions):
            self.interdependencies[sol_index] = {}
            for cpd in sol:
                self.interdependencies[sol_index][cpd] = {}
                coef = {self.variables["met"][cpd]:0}
                coef.update({self.variables["met"][cpd2]:1 for cpd2 in sol if cpd != cpd2})
                BaseFBAPkg.build_constraint(self, "met", sol.objective_value, 
                                            sol.objective_value, coef, f"{cpd}-sol{sol_index}")
                new_sol = self.model.optimize()
                diff = DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))
                
                ## track metabolites that fill the void from the removed compound
                while diff:
                    for key, value in diff.items():
                        new_mets = {re.search("(?<=[\')(.+)(?=\'])", met):change for met, change in value.items()}
                        # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                        self.interdependencies[sol_index][cpd].update(new_mets)
                    diff = self.test_compounds(cpd, sol, sol_index, new_mets.keys())
                        
    def test_compounds(self, cpd, sol, sol_index, zero_compounds):
        coef = {self.variables["met"][cpd]:0 for cpd in zero_compounds}
        coef.update({self.variables["met"][cpd]:1 for cpd in sol if cpd not in zero_compounds})
        cpd_name = "_".join(zero_compounds)
        BaseFBAPkg.build_constraint(self, "met", sol.objective_value, 
                                    sol.objective_value, coef, f"{cpd_name}-sol{sol_index}")
        new_sol = self.model.optimize()
        if new_sol.status != "optimal":
            return False
        return DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))