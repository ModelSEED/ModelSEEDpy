from modelseedpy.core.exceptions import ObjectiveError, FeasibilityError
from modelseedpy.fbapkg.reactionusepkg import ReactionUsePkg
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint
from cobra.medium import minimal_medium
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
        if isinstance(list(sol_dict.keys())[0], str):
            return {rxn:abs(flux) for rxn, flux in sol_dict.items() if "EX_" in rxn and flux < 0}
        return {rxn:abs(flux) for rxn, flux in sol_dict.items() if "EX_" in rxn.id and flux < 0}

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
    def _compatibilize(org_models, printing=False):
        from modelseedpy.community.mscompatibility import MSCompatibility
        return MSCompatibility.standardize(org_models, conflicts_file_name="standardization_corrections.json", printing=printing)

    @staticmethod
    def minimize_flux(org_model, minimal_growth=None, printing=True):
        """minimize the total in-flux of exchange reactions in the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
        model = org_model.copy()
        minimal_growth = minimal_growth or model.optimize().objective_value
        FBAHelper.add_minimal_objective_cons(model, min_value=minimal_growth)
        FBAHelper.add_objective(model, sum(MinimalMediaPkg._influx_objective(model)), "min")
        min_media = MinimalMediaPkg._exchange_solution(FBAHelper.solution_to_dict(model.optimize()))
        total_flux = sum([abs(flux) for flux in min_media.values()])
        # verify the medium
        model2 = org_model.copy()
        model2.medium = min_media
        simulated_sol = model2.optimize()
        if simulated_sol.status != "optimal":
            raise FeasibilityError(f"The simulation was not optimal, with a status of {simulated_sol.status}")
        if not isclose(simulated_sol.objective_value, minimal_growth):
            raise ObjectiveError(f"The assigned minimal_growth of {minimal_growth} was not maintained during the simulation,"
                                 f" where the observed growth value was {simulated_sol.objective_value}.")
        if printing:
            print(f"The minimal flux media consists of {len(min_media)} compounds and a {total_flux} total influx,"
                  f" with a growth value of {simulated_sol.objective_value}")
        return min_media

    @staticmethod
    def minimize_components(org_model, minimal_growth=None, printing=True):
        """minimize the quantity of metabolites that are consumed by the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
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
    def comm_media_est(models, min_growth=0.1, minimization_method="jenga", printing=False):
        media = {"community_media": {}, "members": {}}
        for org_model in models:
            model = org_model.copy()
            reactions = [rxn.name for rxn in model.variables]
            duplicate_reactions = DeepDiff(sorted(reactions), sorted(set(reactions)))
            if duplicate_reactions:
                logger.critical(f'CodeError: The model {model.id} contains {duplicate_reactions}'
                                f' that compromise the model.')
            if minimization_method == "minFlux":
                media["members"][model.id] = {"media": MinimalMediaPkg.minimize_flux(model, min_growth, printing)}
            elif minimization_method == "minComponent":
                media["members"][model.id] = {"media": minimal_medium(model, min_growth, minimize_components=True).to_dict()}
            elif minimization_method == "jenga":
                media["members"][model.id] = {"media": MinimalMediaPkg.jenga_method(model, printing=printing)}
            model.medium = media["members"][model.id]["media"]
            media["members"][model.id]["solution"] = FBAHelper.solution_to_dict(model.optimize())
            media["community_media"] = FBAHelper.sum_dict(model.medium, media["community_media"])
        return media

    @staticmethod
    def interacting_comm_media(models, minimization_method="jenga", media=None, printing=True):
        # define the community minimal media
        media = media or MinimalMediaPkg.comm_media_est(models, 0.1, minimization_method, printing=printing)
        org_media = media["community_media"].copy()
        original_time = process_time()
        # remove exchanges that can be satisfied by cross-feeding
        for model in models:
            for rxnID, flux in media["members"][model.id]["solution"].items():
                if rxnID in media["community_media"] and flux > 0:  ## outflux in solutions
                    stoich = list(model.reactions.get_by_id(rxnID).metabolites.values())[0]
                    media["community_media"][rxnID] += flux*stoich  ## the cytoplasmic removal is captured by negative reactant stoich
        media["community_media"] = {ID: flux for ID, flux in media["community_media"].items() if flux > 0}  # influx in media
        syntrophic_diff = DeepDiff(org_media, media["community_media"])
        changed_quantity = 0 if not syntrophic_diff else len(list(chain(*[v for v in list(dict(syntrophic_diff).values())])))
        if printing:
            print(
                f"Syntrophic fluxes examined after {(process_time() - original_time) / 60} minutes, "
                f"with {changed_quantity} change(s): {syntrophic_diff}")
        return media

    @staticmethod
    def jenga_method(org_model, member_models=None, syntrophy=True, minimal_growth=0.1,
                     conserved_cpds:list=None, export=True, printing=True, compatibilize=False):
        # copy and compatibilize the parameter objects
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
        copied_model = org_model.copy()
        models = [] if not member_models else member_models[:]
        if compatibilize:
            copied_model = MinimalMediaPkg._compatibilize(copied_model)
        original_media = MinimalMediaPkg.minimize_flux(copied_model, printing=False)
        if models:
            media = MinimalMediaPkg.comm_media_est(models, minimal_growth)
            if compatibilize:
                models = MinimalMediaPkg._compatibilize(member_models, printing)
            # subtract syntrophic interactions from media requirements
            if syntrophy:
                media["community_media"] = MinimalMediaPkg.interacting_comm_media(
                    models, "jenga", media, printing)["community_media"]
            original_media = media["community_media"].copy()

        # identify removal=ble compounds
        original_time = process_time()
        copied_model.medium = original_media
        original_obj_value = org_model.optimize().objective_value
        redundant_cpds = set()
        print(original_media)
        for cpd in original_media:
            new_media = original_media.copy()
            new_media.pop(cpd)
            copied_model.medium = new_media
            sol_obj_val = copied_model.slim_optimize()
            if isclose(sol_obj_val, original_obj_value, abs_tol=1e-4):
                redundant_cpds.add(cpd)
            else:
                logger.debug(f"The {sol_obj_val} objective value after the removal of {cpd} "
                             f"does not match the original objective value of {original_obj_value}.")
        if not redundant_cpds:
            logger.debug("None of the media components were determined to be removable.")
            return original_media

        # vet all permutation removals of the redundant compounds
        permuts = [p for p in permutations(redundant_cpds)]
        if printing:
            print(f"The {len(permuts)} permutations of the {redundant_cpds} redundant compounds, "
                  "from absolute tolerance of 1e-4, will be examined.")
        permut_results, failed_permut_starts = [], []
        best = 0
        print(len(permuts))
        for perm_index, permut in enumerate(permuts):
            print(f"{perm_index+1}/{len(permuts)}", end="\r")
            successful_removal = 0
            permut_segments = [permut[:index] for index in range(len(permut), 2, -1)]
            ## eliminate previously discovered failures and successes, respectively
            if any([seg in failed_permut_starts for seg in permut_segments]):
                continue
            if best >= len(permut)/2 and any([set(permut[:best-1]) == set(
                    list(success)[:best-1]) for success in permut_results]):
                continue
            new_media = original_media.copy()
            for cpd in permut:
                ### parameterize and simulate the community
                new_media.pop(cpd)
                copied_model.medium = new_media
                sol = copied_model.optimize()
                if not isclose(sol.objective_value, original_obj_value, abs_tol=1e-7):
                    failed_permut_starts.append(permut[:successful_removal + 1])
                    break
                successful_removal += 1

            if successful_removal >= best:
                if successful_removal > best:
                    best = successful_removal
                    permut_results = []
                permut_removable = permut[:best]  # slice only the elements that are removable
                if permut_removable not in permut_results:
                    permut_results.append(permut_removable)
                    if printing:
                        print(permut_removable)
                        print("best:", best)

        # filter to only the most minimal media
        unique_combinations, unique_paths = [], []
        for removal_path in permut_results:
            path_permutations = permutations(removal_path)
            if all([path in permut_results for path in path_permutations]):
                for com in combinations(removal_path, len(removal_path)):
                    com = set(com)
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
        best_removals = {}
        possible_removals = unique_combinations + unique_paths
        if conserved_cpds:
            possible_removals = [opt for opt in possible_removals if not any(cpd in conserved_cpds for cpd in opt)]
        best = -inf
        for removal in possible_removals:
            cpdID_sum = sum([int(cpd.split('_')[1].replace("cpd", "") if "cpd" in cpd else 500) for cpd in removal])
            if cpdID_sum > best:
                best = cpdID_sum
                best_removals = {best: [removal]}
            elif cpdID_sum == best:
                best_removals[best].append(removal)
        ## arbitrarily select the first removal from those that both maximize the summed cpdID and avoid conserved compounds
        comm_media = FBAHelper.remove_media_compounds(
            original_media, list(best_removals.values())[0][0], printing)
        if printing:
            print(best_removals)
            pprint(comm_media)
        if models:
            media["community_media"] = comm_media

        # communicate results
        jenga_media = comm_media.copy()
        jenga_difference = DeepDiff(original_media, jenga_media)
        changed_quantity = 0 if not jenga_difference else len(list(jenga_difference.values())[0])
        if printing:
            print(f"Jenga fluxes examined after {(process_time()-original_time)/60} minutes, "
                  f"with {changed_quantity} change(s): {jenga_difference}")
        final_media = comm_media if "media" not in locals() else media
        if export:
            if member_models:
                export_name = "_".join([model.id for model in models]) + "_media.json"
            else:
                export_name = copied_model.id + "_media.json"
            with open(export_name, 'w') as out:
                json.dump(final_media, out, indent=3)
        return final_media
        

# class MinimalMedia(BaseFBAPkg):
#     """A class that determines the minimal media of a model"""  # investigate conversion to a staticmethod or single function for convenience and in-line utility
#     def __init__(self, model, min_growth=0.1):
#         # define the system
#         BaseFBAPkg.__init__(self, model, "Minimal Media", {"met":"metabolite"}, {"met":"string", "obj":"string"})
#
#         # define the exchange variables, the minimal objective, and the minimal growth value
#         for cpd in FBAHelper.exchange_reactions(self.model):
#             BaseFBAPkg.build_variable(self,"met",0,1,"binary",cpd)
#         self.model = FBAHelper.add_objective(self.model, sum([var for var in self.variables["met"].values()]), "min")
#         BaseFBAPkg.build_constraint(self, "obj", min_growth, None, {
#             rxn.forward_variable:1 for rxn in FBAHelper.bio_reactions(self.model)}, "min_growth")
#
#         # determine the set of media solutions
#         solutions = []
#         sol = self.model.optimize()
#         while sol.status == "optimal":
#             solutions.append(sol)
#             sol_dict = FBAHelper.solution_to_variables_dict(sol, model)
#             ## omit the solution from the next search
#             BaseFBAPkg.build_constraint(self, "met", len(sol_dict)-1, len(sol_dict)-1,
#                                         coef=sol_dict, cobra_obj=f"sol{len(solutions)}")
#             sol = self.model.optimize()
#         if not solutions:
#             logger.error("No simulations were feasible.")
#
#         # search the permutation space by omitting previously investigated solutions
#         self.interdependencies = {}
#         for sol_index, sol in enumerate(solutions):
#             self.interdependencies[sol_index] = {}
#             for cpd in sol:
#                 self.interdependencies[sol_index][cpd] = {}
#                 coef = {self.variables["met"][cpd]:0}
#                 coef.update({self.variables["met"][cpd2]:1 for cpd2 in sol if cpd != cpd2})
#                 BaseFBAPkg.build_constraint(self, "met", sol.objective_value,
#                                             sol.objective_value, coef, f"{cpd}-sol{sol_index}")
#                 new_sol = self.model.optimize()
#                 diff = DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))
#
#                 ## track metabolites that fill the void from the removed compound
#                 while diff:
#                     for key, value in diff.items():
#                         new_mets = {re.search("(?<=[\')(.+)(?=\'])", met):change for met, change in value.items()}
#                         # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
#                         self.interdependencies[sol_index][cpd].update(new_mets)
#                     diff = self.test_compounds(cpd, sol, sol_index, new_mets.keys())
#
#     def test_compounds(self, cpd, sol, sol_index, zero_compounds):
#         coef = {self.variables["met"][cpd]:0 for cpd in zero_compounds}
#         coef.update({self.variables["met"][cpd]:1 for cpd in sol if cpd not in zero_compounds})
#         cpd_name = "_".join(zero_compounds)
#         BaseFBAPkg.build_constraint(self, "met", sol.objective_value,
#                                     sol.objective_value, coef, f"{cpd_name}-sol{sol_index}")
#         new_sol = self.model.optimize()
#         if new_sol.status != "optimal":
#             return False
#         return DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))