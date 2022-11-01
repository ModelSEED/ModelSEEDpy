from modelseedpy.core.exceptions import ObjectiveError, FeasibilityError
from modelseedpy.fbapkg.reactionusepkg import ReactionUsePkg
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint
from cobra.medium import minimal_medium
from optlang.symbolics import Zero
from math import isclose, inf, factorial
from deepdiff import DeepDiff
from time import process_time
from pprint import pprint
import logging
import json, re

logger = logging.getLogger(__name__)


class MSMinimalMedia:

    @staticmethod
    def _exchange_solution(sol_dict):
        if isinstance(list(sol_dict.keys())[0], str):
            return {rxn:abs(flux) for rxn, flux in sol_dict.items() if "EX_" in rxn and flux < 0}
        elif hasattr(list(sol_dict.keys())[0], "id"):
            return {rxn:abs(flux) for rxn, flux in sol_dict.items() if "EX_" in rxn.id and flux < 0}
        return {rxn:abs(flux) for rxn, flux in sol_dict.items() if "EX_" in rxn.name and flux < 0}

    @staticmethod
    def _influx_objective(model):
        influxes = []
        for ex_rxn in MSMinimalMedia.media_exchange_variables(model):
            if len(ex_rxn.reactants) == 1:  # this is essentially every exchange
                influxes.append(ex_rxn.reverse_variable)
            else:  # this captures edge cases
                influxes.append(ex_rxn.forward_variable)
        return influxes

    @staticmethod
    def media_exchange_variables(model):
        return [exRXN for exRXN in FBAHelper.exchange_reactions(model) if exRXN.id in model.medium]

    @staticmethod
    def _var_to_ID(var):
        rxnID = var.name
        if "_ru" in rxnID:
            rxnID = rxnID.replace("_ru", "")
        return rxnID

    @staticmethod
    def _compatibilize(org_models, printing=False):
        from modelseedpy.community.mscompatibility import MSCompatibility
        return MSCompatibility.standardize(org_models, conflicts_file_name="standardization_corrections.json", printing=printing)

    @staticmethod
    def minimize_flux(org_model, minimal_growth=None, environment=None, printing=True):
        """minimize the total in-flux of exchange reactions in the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
        model = org_model.copy()
        model.medium = environment or model.medium
        # define the MILP
        minimal_growth = minimal_growth or model.optimize().objective_value
        FBAHelper.add_minimal_objective_cons(model, min_value=minimal_growth)
        media_exchanges = MSMinimalMedia._influx_objective(model)
        FBAHelper.add_objective(model, sum(media_exchanges), "min")
        for exRXN in MSMinimalMedia.media_exchange_variables(model):
            FBAHelper.create_constraint(model, Constraint(
                exRXN.reverse_variable, lb=None,ub=environment[exRXN.id], name=exRXN.id + "_maxFlux"))
        # parse the minimal media
        min_media = MSMinimalMedia._exchange_solution(FBAHelper.solution_to_dict(model.optimize()))
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
    def _model_growth(sol_dict):
        return sum([flux for var, flux in sol_dict.items() if "bio" in var.name])

    @staticmethod
    def _define_min_objective(model, interacting):
        rxns = FBAHelper.exchange_reactions(model)
        if not interacting:
            rxns = FBAHelper.transport_reactions(model)
        vars = {}
        for rxn in rxns:
            # define the variable
            vars[rxn.id] = Variable(rxn.id+"_ru", lb=0, ub=1, type="binary")
            model.add_cons_vars(vars[rxn.id])
            model.solver.update()
            # bin_flux: {rxn_bin}*1000 >= {rxn_rev_flux}
            FBAHelper.create_constraint(model, Constraint(Zero, lb=0, ub=None, name=rxn.id+"_bin"),
                                        coef={vars[rxn.id]: 1000, rxn.reverse_variable: -1})
        FBAHelper.add_objective(model, sum(list(vars.values())), "min")
        return model, vars

    @staticmethod
    def minimize_components(org_model, minimal_growth=0.1, solution_limit=100, interacting=True, printing=True):
        """minimize the quantity of metabolites that are consumed by the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
        model = org_model.copy()
        time1 = process_time()
        variables = {"ru":{}}
        FBAHelper.add_minimal_objective_cons(
            model, sum([rxn.flux_expression for rxn in model.reactions if "bio" in rxn.id]), minimal_growth)
        time2 = process_time()
        # define the binary variable and constraint
        model, variables["ru"] = MSMinimalMedia._define_min_objective(model, interacting)
        # print(model.objective)
        print(f"objective definition: {(time2-time1)/60} mins")

        # determine each solution
        interdependencies = {}
        solution_dicts = []
        sol = model.optimize()
        minMedia = [0]*100
        time3 = process_time()
        print(f"initial optimize complete: {(time3-time2)/60} mins")
        count = 0
        while sol.status == "optimal" and count < solution_limit:
            print("sol_dicts length", len(solution_dicts), end="\r")
            sol_dict = FBAHelper.solution_to_variables_dict(sol, model)
            solution_dicts.append(sol_dict)
            time4 = process_time()
            print(f"solution parsed to variables: {(time4-time3)/60} mins")
            if not isclose(minimal_growth, MSMinimalMedia._model_growth(sol_dict),rel_tol=0.00001*minimal_growth):
                raise ObjectiveError(f"The parameterized minimal growth of {minimal_growth} was not maintained, and "
                                     f"the simulation achieved {MSMinimalMedia._model_growth(sol_dict)} growth.")
            sol_media = MSMinimalMedia._exchange_solution(sol_dict)
            minMedia = sol_media if len(sol_media) < len(minMedia) else minMedia
            ## omit the solution from the next search
            FBAHelper.create_constraint(model, Constraint(  # build exclusion use can be emulated
                Zero, lb=None, ub=len(sol_dict)-1,name=f"exclude_sol{len(solution_dicts)}"), sol_dict)

            # search the permutation space by omitting previously investigated solution_dicts
            # sol_exchanges = [rxn for rxn in sol_dict if "EX_" in rxn.name]
            # interdependencies[count] = MSMinimalMedia._examine_permutations(
            #     model, sol_exchanges, variables, sol_dict, count, interacting)
            time5 = process_time()
            print(f"minimal media is calculated and prevented from future iterations: {(time5-time4)/60} mins")

            sol = model.optimize()
            count += 1
        if not solution_dicts:
            logger.error("No simulations were feasible.")
        time6 = process_time()
        print(f"final optimization: {(time6-time5)/60} mins")

        return {rxn.id:flux for rxn, flux in minMedia.items()}

    @staticmethod
    def _knockout(org_model, rxnVar, variables, sol_dict, sol_index, interacting):
        # knockout the specified exchange
        knocked_model = org_model.copy()
        knocked_model, vars = MSMinimalMedia._define_min_objective(knocked_model, {"ru":{}}, interacting)
        coef = {rxnVar: 0}
        if interacting:
            coef.update({variables["ru"][MSMinimalMedia._var_to_ID(rxnVar2)]: 1
                         for rxnVar2 in sol_dict if rxnVar != rxnVar2 and "EX_" in rxnVar2.name})
        else:
            coef.update({variables["ru"][MSMinimalMedia._var_to_ID(rxnVar2)]: 1
                         for rxnVar2 in sol_dict if (
                                 rxnVar != rxnVar2 and any(["_e0" in met.id for met in rxnVar2.metabolites]))
                         })
        FBAHelper.create_constraint(knocked_model, Constraint(Zero, lb=0.1, ub=None, name=f"{rxnVar.name}-sol{sol_index}"), coef)
        return knocked_model.optimize()

    @staticmethod
    def _examine_permutations(model, exchange_ids_to_explore, variables, sol_dict, sol_index, interacting):
        for index, ex in enumerate(exchange_ids_to_explore):
            print(f"{ex.name}: {index}/{len(exchange_ids_to_explore)-1} exchanges to explore")
            sol_dict_sans_ex = sol_dict.copy()
            sol_dict_sans_ex.pop(ex)
            # interdependencies[sol_index][exID] = MSMinimalMedia._examine_permutations(
            #     exID, sol_dict, sol_index, variables, sol_dict_sans_ex)
            interdependencies = {}

            ## explore permutations after removing the selected variable
            diff = DeepDiff(sol_dict_sans_ex, FBAHelper.solution_to_dict(MSMinimalMedia._knockout(
                model, ex, variables, sol_dict, sol_index, interacting)))
            if diff:  # the addition of new exchanges or altered exchange fluxes are detected after the removed exchange
                print(diff)
                for key, changes in diff.items():
                    # for change in changes:
                    #     print(change)
                    changed_reactions = [re.search("(?<=\[\')(.+)(?=\'\])", change).group() for change in changes]
                    # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                    for exchange in [rxn for rxn in changed_reactions if "EX_" in rxn]:
                        interdependencies[sol_index][exchange] = MSMinimalMedia._examine_permutations(
                            model, exchange_ids_to_explore, variables, sol_dict,sol_index+1, interacting)
                # coef = {variables["met"][exID]: 0 for cpd in new_mets.keys()}
                # coef.update({variables["met"][exID]: 1 for exID in sol_dict if exID not in new_mets.keys()})
                cpd_name = "_".join(new_mets.keys())
                new_sol = model.optimize()
                new_sol_dict = FBAHelper.solution_to_variables_dict(new_sol, model)
                new_sol_exchanges = [rxn for rxn in sol_dict if "EX_" in rxn.name]
                if new_sol.status != "optimal":
                    return interdependencies
                MSMinimalMedia._examine_permutations(model, new_sol_exchanges, variables, new_sol_dict, sol_index+1, interacting)
            return interdependencies

    @staticmethod
    def comm_media_est(models, comm_model, min_growth=0.1, minimization_method="jenga", environment=None, printing=False):
        media = {"community_media": {}, "members": {}}
        for org_model in models:
            model = org_model.copy()
            model.medium = environment or model.medium
            reactions = [rxn.name for rxn in model.variables]
            duplicate_reactions = DeepDiff(sorted(reactions), sorted(set(reactions)))
            if duplicate_reactions:
                logger.critical(f'CodeError: The model {model.id} contains {duplicate_reactions}'
                                f' that compromise the model.')
            if minimization_method == "minFlux":
                media["members"][model.id] = {"media": MSMinimalMedia.minimize_flux(model, min_growth, printing)}
            elif minimization_method == "minComponents":
                media["members"][model.id] = {"media": minimal_medium(model, min_growth, minimize_components=True).to_dict()}
            elif minimization_method == "jenga":
                media["members"][model.id] = {"media": MSMinimalMedia.jenga_method(model, printing=printing)}
                media["community_media"] = FBAHelper.sum_dict(media["members"][model.id]["media"], media["community_media"])
            model.medium = media["members"][model.id]["media"]
            media["members"][model.id]["solution"] = FBAHelper.solution_to_dict(model.optimize())
        if comm_model:
            if environment:
                comm_model.medium = environment
            if minimization_method == "jenga":
                print("Community models are too excessive for direct assessment via the JENGA method; "
                      "thus, the community minimal media is estimated as the combination of member media.")
            elif minimization_method == "minFlux":
                media["community_media"] = MSMinimalMedia.minimize_flux(comm_model, min_growth, printing)
            elif minimization_method == "minComponents":
                media["community_media"] = minimal_medium(comm_model, min_growth, minimize_components=True).to_dict()
        return media

    @staticmethod
    def interacting_comm_media(models, comm_model, minimization_method="jenga", min_growth=0.1,
                               media=None, environment=None, printing=True):
        # define the community minimal media
        media = media or MSMinimalMedia.comm_media_est(
            models, comm_model, min_growth, minimization_method, environment, printing=printing)
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
    def jenga_method(org_model, org_media=None, conserved_cpds:list=None,
                     export=True, printing=True, compatibilize=False, environment=None):
        # copy and compatibilize the parameter objects
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                                 "which is incompatible with minimal media computations.")
        copied_model = org_model.copy()
        copied_model.medium = environment or copied_model.medium
        if compatibilize:
            copied_model = MSMinimalMedia._compatibilize(copied_model)
        # TODO - the COBRA method must eventually be replaced with MSMinimalMedia.minimize_components(copied_model, printing=False)
        original_media = org_media or minimal_medium(copied_model, minimize_components=True).to_dict()
        # {cpd.replace("EX_", ""): flux for cpd, flux in .items()}

        # identify removal=ble compounds
        original_time = process_time()
        copied_model.medium = original_media
        original_obj_value = org_model.optimize().objective_value
        redundant_cpds = set()
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
        if len(redundant_cpds) > 9:
            import sigfig
            num_permuts = sigfig.round(factorial(len(redundant_cpds)), sigfigs=2, format='sci')
            raise FeasibilityError(f"The model {copied_model.id} contains {len(redundant_cpds)} removable"
                                   f" compounds, which yields {num_permuts} permutations and is untenable for computation."
                                   " Select a different minimal media method such as 'minFlux' or 'minComponents'.")

        # vet all permutation removals of the redundant compounds
        permuts = [p for p in permutations(redundant_cpds)]
        if printing:
            print(f"The {len(permuts)} permutations of the {redundant_cpds} redundant compounds, "
                  "from absolute tolerance of 1e-4, will be examined.")
        permut_results, failed_permut_starts = [], []
        best = 0
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
        media = FBAHelper.remove_media_compounds(
            original_media, list(best_removals.values())[0][0], printing)
        if printing:
            print(best_removals)
            pprint(media)

        # communicate results
        jenga_media = media.copy()
        jenga_difference = DeepDiff(original_media, jenga_media)
        changed_quantity = 0 if not jenga_difference else len(list(jenga_difference.values())[0])
        if printing:
            print(f"Jenga fluxes examined after {(process_time()-original_time)/60} minutes, "
                  f"with {changed_quantity} change(s): {jenga_difference}")
        if export:
            export_name = copied_model.id + "_media.json"
            with open(export_name, 'w') as out:
                json.dump(media, out, indent=3)
        return media

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