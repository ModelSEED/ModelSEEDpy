from modelseedpy.core.exceptions import ObjectiveError, FeasibilityError
from modelseedpy.fbapkg.reactionusepkg import ReactionUsePkg
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint
from cobra.medium import minimal_medium
from optlang.symbolics import Zero
from math import isclose, inf, factorial
#from deepdiff import DeepDiff
from time import process_time
from pprint import pprint
import logging
import json, re

logger = logging.getLogger(__name__)


def _exchange_solution(sol_dict):
    if isinstance(list(sol_dict.keys())[0], str):
        return {
            rxn: abs(flux)
            for rxn, flux in sol_dict.items()
            if "EX_" in rxn and flux < 0
        }
    elif hasattr(list(sol_dict.keys())[0], "id"):
        return {
            rxn.id: abs(flux)
            for rxn, flux in sol_dict.items()
            if "EX_" in rxn.id and flux < 0
        }
    return {
        rxn.name: abs(flux)
        for rxn, flux in sol_dict.items()
        if "EX_" in rxn.name and flux < 0
    }


def _model_growth(sol_dict):
    return sum(
        [flux for var, flux in sol_dict.items() if re.search(r"(^bio\d+$)", var.name)]
    )


def _var_to_ID(var):
    rxnID = var.name
    if "_ru" in rxnID:
        rxnID = rxnID.replace("_ru", "")
    return rxnID


def _compatibilize(org_models, printing=False):
    from commscores import GEMCompatibility

    return GEMCompatibility.standardize(
        org_models,
        conflicts_file_name="standardization_corrections.json",
        printing=printing,
    )


def verify(org_model, min_media):
    model2 = org_model.copy()
    model2.medium = min_media
    return model2.optimize()


def bioFlux_check(model, sol=None, sol_dict=None, min_growth=0.1):
    sol_dict = sol_dict or FBAHelper.solution_to_variables_dict(sol, model)
    # print({k:v for k,v in sol_dict.items() if v > 1E-8})
    simulated_growth = max(
        sum(
            [
                flux
                for var, flux in sol_dict.items()
                if re.search(r"(^bio\d+$)", var.name)
            ]
        ),
        sol.objective_value,
    )
    if simulated_growth < min_growth * 0.9999:
        raise ObjectiveError(
            f"The assigned minimal_growth of {min_growth} was not maintained during the simulation,"
            f" where the observed growth value was {simulated_growth}."
        )
    if sol.status != "optimal":
        display(sol)
    return sol_dict


def minimizeFlux_withGrowth(model_util, min_growth, obj):
    model_util.add_minimal_objective_cons(min_growth)
    model_util.add_objective(obj, "min")
    # print(model_util.model.objective)
    # print([(cons.lb, cons.expression) for cons in model_util.model.constraints if "min" in cons.name])
    sol = model_util.model.optimize()
    # print(sol.objective_value)
    sol_dict = bioFlux_check(model_util.model, sol, min_growth=min_growth)
    return sol, sol_dict


class MSMinimalMedia:

    @staticmethod
    def _influx_objective(model_util, interacting):
        rxns = (
            model_util.exchange_list() if interacting else model_util.transport_list()
        )
        influxes = []
        for rxn in rxns:
            if any(
                ["e0" in met.id for met in rxn.reactants]
            ):  # this is essentially every exchange
                influxes.append(rxn.reverse_variable)
            elif any(
                ["e0" in met.id for met in rxn.products]
            ):  # this captures edge cases or transporters
                influxes.append(rxn.forward_variable)
            else:
                logger.critical(
                    f"The reaction {rxn} lacks exchange metabolites, which indicates an error."
                )
        return influxes

    @staticmethod
    def minimize_flux(
        org_model, min_growth=None, environment=None, interacting=True, printing=True
    ):
        """minimize the total in-flux of exchange reactions in the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(
                f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                "which is incompatible with minimal media computations."
            )
        model_util = MSModelUtil(org_model, True)
        model_util.add_medium(environment or model_util.model.medium)
        # define the MILP
        min_growth = (
            model_util.model.slim_optimize()
            if min_growth is None
            else min(min_growth, model_util.model.slim_optimize())
        )
        # min_flux = MSMinimalMedia._min_consumption_objective(model_util, interacting)
        media_exchanges = MSMinimalMedia._influx_objective(model_util, interacting)
        # parse the minimal media
        sol, sol_dict = minimizeFlux_withGrowth(
            model_util, min_growth, sum(media_exchanges)
        )
        min_media = _exchange_solution(sol_dict)
        total_flux = sum([abs(flux) for flux in min_media.values()])
        simulated_sol = verify(org_model, min_media)
        if simulated_sol.status != "optimal":
            raise FeasibilityError(
                f"The simulation was not optimal, with a status of {simulated_sol.status}"
            )
        if printing:
            print(
                f"The minimal flux media for {org_model.id} consists of {len(min_media)} compounds and a {total_flux} total influx,"
                f" with a growth value of {simulated_sol.objective_value}"
            )
        return min_media, sol

    @staticmethod
    def _min_consumption_objective(model_util, interacting):
        rxns = (
            model_util.exchange_list() if interacting else model_util.transport_list()
        )
        vars = {}
        for rxn in rxns:
            cons_name = rxn.id + "_bin"
            if cons_name in model_util.model.constraints:
                print(
                    f"The {cons_name} constraint already exists in "
                    f"{model_util.model.id} and thus is skipped.\n"
                )
                continue

            # define the variable
            var_name = rxn.id + "_ru"
            if var_name in model_util.model.variables:
                print(
                    f"The {var_name} variable already exists in "
                    f"{model_util.model.id} and thus is skipped.\n"
                )
                continue
            vars[rxn.id] = Variable(var_name, lb=0, ub=1, type="binary")
            model_util.add_cons_vars([vars[rxn.id]])
            # bin_flux: {rxn_bin}*1000 >= {rxn_rev_flux}
            model_util.create_constraint(
                Constraint(Zero, lb=0, ub=None, name=cons_name),
                coef={vars[rxn.id]: 1000, rxn.reverse_variable: -1},
            )
        return vars

    @staticmethod
    def conserved_exchanges():
        pass

    @staticmethod
    def relative_media():
        pass

    @staticmethod
    def minimize_components(
        org_model,
        min_growth=0.1,
        environment=None,
        interacting=True,
        solution_limit=5,
        printing=True,
    ):
        """minimize the quantity of metabolites that are consumed by the model"""
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(
                f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                "which is incompatible with minimal media computations."
            )
        model_util = MSModelUtil(org_model, True)
        model_util.add_timeout(10)
        print("Minimal Components media")
        if environment:
            model_util.add_medium(environment)
        # ic(org_model, min_growth, solution_limit)
        model_util.add_minimal_objective_cons(
            min_growth
        )  # , model_util.model.reactions.bio1.flux_expression)
        # print(model_util.model.constraints[-1])
        # define the binary variable and constraint
        time1 = process_time()
        variables = {
            "ru": MSMinimalMedia._min_consumption_objective(model_util, interacting)
        }
        model_util.add_objective(sum(variables["ru"].values()), "min")
        time2 = process_time()
        print(f"\nDefinition of minimum objective time: {(time2 - time1)/60} mins")

        # determine each solution
        # interdependencies = {}
        solution_dicts, min_media = [], [0] * 1000
        sol = (
            model_util.model.optimize()
        )  # TODO This is the troublesome line that occasionally refuses to solve
        if "optimal" not in sol.status:
            raise FeasibilityError(
                f"The simulation for minimal uptake in {model_util.model.id} was {sol.status}."
            )
        time3 = process_time()
        broken = False
        while (
            not broken
            and sol.status == "optimal"
            and len(solution_dicts) < solution_limit
        ):
            print(f"Iteration {len(solution_dicts)}", end="\r")
            sol_dict = FBAHelper.solution_to_variables_dict(sol, model_util.model)
            ## ensure that the minimal growth is respected
            simulated_growth = _model_growth(sol_dict)
            if simulated_growth < min_growth * 0.9999:
                raise ObjectiveError(
                    f"The minimal growth of {min_growth} was not maintained; "
                    f"the simulation achieved {simulated_growth} growth."
                )
            sol_rxns_dict = FBAHelper.solution_to_rxns_dict(sol, model_util.model)
            solution_dicts.append(sol_dict)
            sol_media = _exchange_solution(sol_rxns_dict)
            min_media = sol_media if len(sol_media) < len(min_media) else min_media
            ## omit the solution from future searches
            model_util.create_constraint(
                Constraint(
                    Zero,
                    lb=None,
                    ub=len(sol_dict) - 1,
                    name=f"exclude_sol{len(solution_dicts)}",
                ),
                sol_dict,
            )

            # search the permutation space by omitting previously investigated solution_dicts
            # sol_exchanges = [rxn for rxn in sol_dict if "EX_" in rxn.name]
            # interdependencies[count] = MSMinimalMedia._examine_permutations(
            #     model, sol_exchanges, variables, sol_dict, count, interacting)
            try:
                sol = model_util.model.optimize()
            except:
                broken = True
            if broken:
                break
        if not solution_dicts:
            raise FeasibilityError("The model was not feasibly simulated.")
        min_media = {rxn: flux for rxn, flux in min_media.items()}
        simulated_sol = verify(org_model, min_media)
        if simulated_sol.status != "optimal":
            raise FeasibilityError(
                f"The predicted media {min_media} is not compatible with its model {org_model.id}, "
                f"and possesses a(n) {simulated_sol.status} status."
            )
        time6 = process_time()
        print(f"Optimization time: {(time6-time3)/60} mins")
        return min_media

    @staticmethod
    def _knockout(org_model, rxnVar, variables, sol_dict, sol_index, interacting):
        # knockout the specified exchange
        knocked_model_utl = MSModelUtil(org_model, True)
        knocked_model_utl, vars = MSMinimalMedia._min_consumption_objective(
            knocked_model_utl, interacting
        )
        coef = {rxnVar: 0}
        if interacting:
            coef.update(
                {
                    variables["ru"][_var_to_ID(rxnVar2)]: 1
                    for rxnVar2 in sol_dict
                    if rxnVar != rxnVar2 and "EX_" in rxnVar2.name
                }
            )
        else:
            coef.update(
                {
                    variables["ru"][_var_to_ID(rxnVar2)]: 1
                    for rxnVar2 in sol_dict
                    if (
                        rxnVar != rxnVar2
                        and any(["_e0" in met.id for met in rxnVar2.metabolites])
                    )
                }
            )
        knocked_model_utl.create_constraint(
            Constraint(Zero, lb=0.1, ub=None, name=f"{rxnVar.name}-sol{sol_index}"),
            coef,
        )
        return knocked_model_utl.optimize()

    @staticmethod
    def _examine_permutations(
        model, exchange_ids_to_explore, variables, sol_dict, sol_index, interacting
    ):
        for index, ex in enumerate(exchange_ids_to_explore):
            print(
                f"{ex.name}: {index}/{len(exchange_ids_to_explore)-1} exchanges to explore"
            )
            sol_dict_sans_ex = sol_dict.copy()
            sol_dict_sans_ex.pop(ex)
            # interdependencies[sol_index][exID] = MSMinimalMedia._examine_permutations(
            #     exID, sol_dict, sol_index, variables, sol_dict_sans_ex)
            interdependencies = {}

            ## explore permutations after removing the selected variable
            diff = DeepDiff(
                sol_dict_sans_ex,
                FBAHelper.solution_to_dict(
                    MSMinimalMedia._knockout(
                        model, ex, variables, sol_dict, sol_index, interacting
                    )
                ),
            )
            if (
                diff
            ):  # the addition of new exchanges or altered exchange fluxes are detected after the removed exchange
                print(diff)
                for key, changes in diff.items():
                    # for change in changes:
                    #     print(change)
                    changed_reactions = [
                        re.search("(?<=\[')(.+)(?='\])", change).group()
                        for change in changes
                    ]
                    # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                    for exchange in [rxn for rxn in changed_reactions if "EX_" in rxn]:
                        interdependencies[sol_index][exchange] = (
                            MSMinimalMedia._examine_permutations(
                                model,
                                exchange_ids_to_explore,
                                variables,
                                sol_dict,
                                sol_index + 1,
                                interacting,
                            )
                        )
                # coef = {variables["met"][exID]: 0 for cpd in new_mets.keys()}
                # coef.update({variables["met"][exID]: 1 for exID in sol_dict if exID not in new_mets.keys()})
                # cpd_name = "_".join(new_mets.keys())
                new_sol = model.optimize()
                new_sol_dict = FBAHelper.solution_to_variables_dict(new_sol, model)
                new_sol_exchanges = [rxn for rxn in sol_dict if "EX_" in rxn.name]
                if new_sol.status != "optimal":
                    return interdependencies
                MSMinimalMedia._examine_permutations(
                    model,
                    new_sol_exchanges,
                    variables,
                    new_sol_dict,
                    sol_index + 1,
                    interacting,
                )
            return interdependencies

    @staticmethod
    def determine_min_media(
        model,
        minimization_method="minFlux",
        min_growth=None,
        environment=None,
        interacting=True,
        solution_limit=5,
        printing=True,
    ):
        if minimization_method == "minFlux":
            return MSMinimalMedia.minimize_flux(
                model, min_growth, environment, interacting, printing
            )
        if minimization_method == "minComponents":
            return minimal_medium(model, min_growth, minimize_components=True)
            # return MSMinimalMedia.minimize_components(
            #     model, min_growth, environment, interacting, solution_limit, printing)
        if minimization_method == "jenga":
            return MSMinimalMedia.jenga_method(model, printing=printing)

    @staticmethod
    def comm_media_est(
        models,
        comm_model,
        minimization_method="minComponents",
        min_growth=0.1,
        environment=None,
        interacting=True,
        n_solutions=5,
        printing=False,
    ):
        media = {"community_media": {}, "members": {}}
        # print("com_media_est")
        for org_model in models:
            model_util = MSModelUtil(org_model, True)
            # print(model_util.model.optimize())
            if environment:
                # print(environment)
                model_util.add_medium(environment)
            # reactions = [rxn.name for rxn in model.variables]
            # duplicate_reactions = DeepDiff(sorted(set(reactions)), sorted(reactions))
            # if duplicate_reactions:
            #     logger.critical(f'CodeError: The model {model.id} contains {duplicate_reactions}'
            #                     f' that compromise the model.')
            media["members"][model_util.model.id] = {
                "media": MSMinimalMedia.determine_min_media(
                    model_util.model,
                    minimization_method,
                    min_growth,
                    environment,
                    interacting,
                    n_solutions,
                    printing,
                ),
                "solution": FBAHelper.solution_to_dict(model_util.model.optimize()),
            }
            if minimization_method == "jenga":
                media["community_media"] = FBAHelper.sum_dict(
                    media["members"][model_util.model.id]["media"],
                    media["community_media"],
                )
        if comm_model:
            comm_util = MSModelUtil(comm_model)
            if environment:
                comm_util.add_medium(environment)
            # if minimization_method == "jenga":
            #     print("Community models are too excessive for direct assessment via the JENGA method; "
            #           "thus, the community minimal media is estimated as the combination of member media.")
            #     return media
            media["community_media"] = MSMinimalMedia.determine_min_media(
                comm_model,
                minimization_method,
                min_growth,
                environment,
                interacting,
                n_solutions,
                printing,
            )
        return media

    @staticmethod
    def interacting_comm_media(
        models,
        comm_model,
        minimization_method="jenga",
        min_growth=0.1,
        media=None,
        environment=None,
        printing=True,
    ):
        # define the community minimal media
        media = media or MSMinimalMedia.comm_media_est(
            models,
            comm_model,
            min_growth,
            minimization_method,
            environment,
            printing=printing,
        )
        org_media = media["community_media"].copy()
        original_time = process_time()
        # remove exchanges that can be satisfied by cross-feeding
        for model in models:
            for rxnID, flux in media["members"][model.id]["solution"].items():
                if (
                    rxnID in media["community_media"] and flux > 0
                ):  ## outflux in solutions
                    stoich = list(
                        model.reactions.get_by_id(rxnID).metabolites.values()
                    )[0]
                    media["community_media"][rxnID] += (
                        flux * stoich
                    )  ## the cytoplasmic removal is captured by negative reactant stoich
        media["community_media"] = {
            ID: flux for ID, flux in media["community_media"].items() if flux > 0
        }  # influx in media
        syntrophic_diff = DeepDiff(org_media, media["community_media"])
        changed_quantity = (
            0
            if not syntrophic_diff
            else len(list(chain(*[v for v in list(dict(syntrophic_diff).values())])))
        )
        if printing:
            print(
                f"Syntrophic fluxes examined after {(process_time() - original_time) / 60} minutes, "
                f"with {changed_quantity} change(s): {syntrophic_diff}"
            )
        return media

    @staticmethod
    def jenga_method(
        org_model,
        org_media=None,
        conserved_cpds: list = None,
        export=True,
        printing=True,
        compatibilize=False,
        environment=None,
    ):
        # copy and compatibilize the parameter objects
        if org_model.slim_optimize() == 0:
            raise ObjectiveError(
                f"The model {org_model.id} possesses an objective value of 0 in complete media, "
                "which is incompatible with minimal media computations."
            )
        copied_model = org_model.copy()
        copied_model.medium = environment or copied_model.medium
        if compatibilize:
            copied_model = _compatibilize(copied_model)
        original_media = org_media or MSMinimalMedia.minimize_components(copied_model)
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
                logger.debug(
                    f"The {sol_obj_val} objective value after the removal of {cpd} "
                    f"does not match the original objective value of {original_obj_value}."
                )
        if not redundant_cpds:
            logger.debug(
                "None of the media components were determined to be removable."
            )
            return original_media
        if len(redundant_cpds) > 9:
            import sigfig

            num_permuts = sigfig.round(
                factorial(len(redundant_cpds)), sigfigs=2, format="sci"
            )
            raise FeasibilityError(
                f"The model {copied_model.id} contains {len(redundant_cpds)} removable"
                f" compounds, which yields {num_permuts} permutations and is untenable for computation."
                " Select a different minimal media method such as 'minFlux' or 'minComponents'."
            )

        # vet all permutation removals of the redundant compounds
        permuts = [p for p in permutations(redundant_cpds)]
        if printing:
            print(
                f"The {len(permuts)} permutations of the {redundant_cpds} redundant compounds, "
                "from absolute tolerance of 1e-4, will be examined."
            )
        permut_results, failed_permut_starts = [], []
        best = 0
        for perm_index, permut in enumerate(permuts):
            print(f"{perm_index+1}/{len(permuts)}", end="\r")
            successful_removal = 0
            permut_segments = [permut[:index] for index in range(len(permut), 2, -1)]
            ## eliminate previously discovered failures and successes, respectively
            if any([seg in failed_permut_starts for seg in permut_segments]):
                continue
            if best >= len(permut) / 2 and any(
                [
                    set(permut[: best - 1]) == set(list(success)[: best - 1])
                    for success in permut_results
                ]
            ):
                continue
            new_media = original_media.copy()
            for cpd in permut:
                ### parameterize and simulate the community
                new_media.pop(cpd)
                copied_model.medium = new_media
                sol = copied_model.optimize()
                if not isclose(sol.objective_value, original_obj_value, abs_tol=1e-7):
                    failed_permut_starts.append(permut[: successful_removal + 1])
                    break
                successful_removal += 1

            if successful_removal >= best:
                if successful_removal > best:
                    best = successful_removal
                    permut_results = []
                permut_removable = permut[
                    :best
                ]  # slice only the elements that are removable
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
            possible_removals = [
                opt
                for opt in possible_removals
                if not any(cpd in conserved_cpds for cpd in opt)
            ]
        best = -inf
        for removal in possible_removals:
            cpdID_sum = sum(
                [
                    int(cpd.split("_")[1].replace("cpd", "") if "cpd" in cpd else 500)
                    for cpd in removal
                ]
            )
            if cpdID_sum > best:
                best = cpdID_sum
                best_removals = {best: [removal]}
            elif cpdID_sum == best:
                best_removals[best].append(removal)
        ## arbitrarily select the first removal from those that both maximize the summed cpdID and avoid conserved compounds
        media = FBAHelper.remove_media_compounds(
            original_media, list(best_removals.values())[0][0], printing
        )
        if printing:
            print(best_removals)
            pprint(media)

        # communicate results
        jenga_media = media.copy()
        jenga_difference = DeepDiff(original_media, jenga_media)
        changed_quantity = (
            0 if not jenga_difference else len(list(jenga_difference.values())[0])
        )
        if printing:
            print(
                f"Jenga fluxes examined after {(process_time()-original_time)/60} minutes, "
                f"with {changed_quantity} change(s): {jenga_difference}"
            )
        if export:
            export_name = copied_model.id + "_media.json"
            with open(export_name, "w") as out:
                json.dump(media, out, indent=3)
        return media
