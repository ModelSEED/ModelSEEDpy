from modelseedpy.core.msminimalmedia import minimizeFlux_withGrowth, bioFlux_check
from modelseedpy.core.exceptions import NoFluxError, ObjectiveError
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper
from cobra import Model, Reaction, Metabolite
from cobra.medium import minimal_medium

# from commscores import GEMCompatibility
from cobra.flux_analysis import pfba
from collections import OrderedDict
from optlang.symbolics import Zero
from optlang import Constraint
from math import inf, isclose
from pandas import DataFrame
from pprint import pprint
from numpy import mean
import re


def strip_comp(ID):
    ID = ID.replace("-", "~")
    return re.sub("(\_\w\d)", "", ID)


def export_lp(model, name):
    with open(f"{name}.lp", "w") as out:
        out.write(model.solver.to_lp())


def correct_nonMSID(nonMSobject, output, model_index):
    name, compartment = output
    index = 0 if compartment == "e" else model_index
    nonMSobject.compartment = compartment + str(index)
    comp = re.search(r"(_[a-z]\d+$)", nonMSobject.id)
    if comp is None and rf"[{compartment}]" in nonMSobject.id:
        return nonMSobject.id.replace(
            rf"[{compartment}]", f"_{nonMSobject.compartment}"
        )
    elif comp is None:
        return nonMSobject.id + f"_{nonMSobject.compartment}"
    return "_".join([nonMSobject.id.replace(comp.group(), ""), nonMSobject.compartment])


def build_from_species_models(
    org_models,
    model_id=None,
    name=None,
    abundances=None,
    standardize=False,
    MSmodel=True,
    commkinetics=True,
    copy_models=True,
    printing=False,
):
    """Merges the input list of single species metabolic models into a community metabolic model

    Parameters
    ----------
    org_models : list<Cobra.Model> to be merged into a community model
    model_id : string specifying community model ID
    name : string specifying community model name
    names : list<string>  human-readable names for models being merged
    abundances : dict<string,float> relative abundances for input models in community model
    cobra_model : bool for whether the raw COBRA model is returned
    standardize: bool for whether the exchanges of each member model will be standardized (True) or just aligned.

    Returns
    -------
    Cobra.Model for the desired Community

    Raises
    ------
    """
    # construct the new model
    models = org_models  # if not standardize else GEMCompatibility.standardize(
    # org_models, exchanges=True, conflicts_file_name='exchanges_conflicts.json')
    biomass_indices = []
    biomass_index = minimal_biomass_index = 2
    new_metabolites, new_reactions = set(), set()
    member_biomasses = {}
    for model_index, org_model in enumerate(models):
        model_util = MSModelUtil(org_model, copy=copy_models)
        model_reaction_ids = [rxn.id for rxn in model_util.model.reactions]
        model_index += 1
        # if MSmodel:
        # Rename metabolites
        for met in model_util.model.metabolites:
            # Renaming compartments
            output = MSModelUtil.parse_id(met)
            if printing:
                print(met, output)
            if output is None:
                if printing:
                    print(
                        f"The {met.id} ({output}; {hasattr(met, 'compartment')}) is unpredictable."
                    )
                met.id = correct_nonMSID(met, (met.id, "c"), model_index)
            elif len(output) == 2:
                met.id = correct_nonMSID(met, output, model_index)
            elif len(output) == 3:
                name, compartment, out_index = output
                index = 0 if compartment == "e" else model_index
                if out_index == "":
                    met.id += str(index)
                    met.compartment += str(index)
                elif compartment == "e":
                    met.compartment = "e0"
                else:
                    met.compartment = compartment + str(index)
                    met.id = name + "_" + met.compartment
            new_metabolites.add(met)
            if "cpd11416_c" in met.id or "biomass" in met.id:
                member_biomasses[org_model.id] = met
        # Rename reactions
        for (
            rxn
        ) in (
            model_util.model.reactions
        ):  # !!! all reactions should have a non-zero compartment index
            if rxn.id[0:3] != "EX_":
                ## biomass reactions
                if re.search("^(bio)(\d+)$", rxn.id):
                    index = int(re.sub(r"(^bio)", "", rxn.id))
                    if biomass_index == 2:
                        while f"bio{biomass_index}" in model_reaction_ids:
                            biomass_index += 1
                    if index not in biomass_indices and index >= minimal_biomass_index:
                        biomass_indices.append(index)
                    else:  # biomass indices can be decoupled from the respective reaction indices of the same model
                        rxn.id = "bio" + str(biomass_index)
                        if rxn.id not in model_reaction_ids:
                            biomass_indices.append(biomass_index)
                        else:
                            index = minimal_biomass_index
                            rxn.id = "bio" + str(index)
                            while (
                                rxn.id not in model_reaction_ids
                                and index not in biomass_indices
                            ):
                                index += 1
                                rxn.id = "bio" + str(index)
                            biomass_indices.append(index)
                    biomass_index += 1
                ## non-biomass reactions
                else:
                    initialID = str(rxn.id)
                    output = MSModelUtil.parse_id(rxn)
                    if output is None:
                        if printing:
                            print(
                                f"The {rxn.id} ({output}; {hasattr(rxn, 'compartment')}) is unpredictable."
                            )
                        try:
                            rxn.id = correct_nonMSID(rxn, (rxn.id, "c"), model_index)
                        except ValueError:
                            pass
                    elif len(output) == 2:
                        rxn.id = correct_nonMSID(rxn, output, model_index)
                    elif len(output) == 3:
                        name, compartment, index = output
                        if compartment != "e":
                            rxn.name = f"{name}_{compartment}{model_index}"
                            rxn_id = re.search(r"(.+\_\w)(?=\d+)", rxn.id).group()
                            if index == "":
                                rxn.id += str(model_index)
                            else:
                                rxn.id = rxn_id + str(model_index)
                    finalID = str(rxn.id)
                    string_diff = ""
                    for index, let in enumerate(finalID):
                        if (
                            index >= len(initialID)
                            or index < len(initialID)
                            and let != initialID[index]
                        ):
                            string_diff += let
                    if string_diff != f"_{compartment}{model_index}" and printing:
                        print(
                            f"The ID {initialID} is changed with {string_diff} to create the final ID {finalID}"
                        )
            new_reactions.add(rxn)
        # else:
        #     # TODO develop a method for compartmentalizing models without editing all reaction IDs or assuming their syntax
        #     pass
    # adds only unique reactions and metabolites to the community model
    newmodel = Model(
        model_id or "+".join([model.id for model in models]),
        name or "+".join([model.name for model in models]),
    )
    newmodel.add_reactions(FBAHelper.filter_cobra_set(new_reactions))
    newmodel.add_metabolites(FBAHelper.filter_cobra_set(new_metabolites))

    # Create community biomass
    comm_biomass = Metabolite("cpd11416_c0", None, "Community biomass", 0, "c0")
    metabolites = {comm_biomass: 1}
    ## constrain the community abundances
    if abundances:
        abundances = {
            met: abundances[memberID] for memberID, met in member_biomasses.items()
        }
    else:
        abundances = {
            cpd: -1 / len(member_biomasses) for cpd in member_biomasses.values()
        }
    ## define community biomass components
    metabolites.update(abundances)
    comm_biorxn = Reaction(id="bio1", name="bio1", lower_bound=0, upper_bound=1000)
    comm_biorxn.add_metabolites(metabolites)
    newmodel.add_reactions([comm_biorxn])
    # update model components
    newutl = MSModelUtil(newmodel)
    newutl.add_objective(comm_biorxn.flux_expression)
    newutl.model.add_boundary(
        comm_biomass, "sink"
    )  # Is a sink reaction for reversible cpd11416_c0 consumption necessary?
    ## proportionally limit the fluxes to their abundances
    if commkinetics:
        add_commkinetics(newutl, models, member_biomasses, abundances)
    # add the metadata of community composition
    if hasattr(newutl.model, "_context"):
        newutl.model._contents.append(member_biomasses)
    elif hasattr(newutl.model, "notes"):
        newutl.model.notes.update(member_biomasses)
    # print([cons.name for cons in newutl.model.constraints])
    return newutl.model


def add_commkinetics(util, models, member_biomasses, abundances):
    # TODO this creates an error with the member biomass reactions not being identified in the model
    coef = {}
    for model in models:
        coef[member_biomasses[model.id]] = -abundances[member_biomasses[model.id]]
        for rxn in model.reactions:
            if rxn.id[:3] == "rxn":
                coef[rxn.forward_variable] = coef[rxn.reverse_variable] = 1
    util.create_constraint(
        Constraint(Zero, name="member_flux_limit"), coef=coef, printing=True
    )


def phenotypes(community_members, phenotype_flux_threshold=0.1, solver: str = "glpk"):
    # log information of each respective model
    models = OrderedDict()
    solutions = []
    media_conc = set()
    # calculate all phenotype profiles for all members
    comm_members = community_members.copy()
    # print(community_members)
    for (
        org_model,
        content,
    ) in (
        community_members.items()
    ):  # community_members excludes the stationary phenotype
        print("\n", org_model.id)
        org_model.solver = solver
        all_phenotypes = "phenotypes" not in content
        model_util = MSModelUtil(org_model, True)
        if "org_coef" not in locals():
            org_coef = {
                model_util.model.reactions.get_by_id(
                    "EX_cpd00007_e0"
                ).reverse_variable: -1
            }
        model_util.standard_exchanges()
        models[org_model.id] = {
            "exchanges": model_util.exchange_list(),
            "solutions": {},
            "name": content["name"],
        }
        phenotypes = (
            {
                met.name: {"consumed": met.id.replace("EX_", "").replace("_e0", "")}
                for met in model_util.carbon_exchange_mets_list(include_unknown=False)
            }
            if all_phenotypes
            else content["phenotypes"]
        )
        # print(phenotypes)
        models[org_model.id]["phenotypes"] = ["stationary"] + [
            content["phenotypes"].keys() for member, content in comm_members.items()
        ]
        phenoRXNs = [
            pheno_cpd
            for pheno, pheno_cpds in content["phenotypes"].items()
            for pheno_cpd in pheno_cpds["consumed"]
        ]
        media = {cpd: 100 for cpd, flux in model_util.model.medium.items()}
        # TODO correct or remove the media, since it seems to be overwritten by the optimization of all carbon exchanges
        ### eliminate hydrogen absorption
        media.update({"EX_cpd11640_e0": 0})
        past_phenoRXNs = []
        for name, phenoCPDs in phenotypes.items():
            pheno_util = MSModelUtil(model_util.model, True)
            metID = phenoCPDs["consumed"][0]
            try:
                phenoRXN = pheno_util.model.reactions.get_by_id(f"EX_{metID}_e0")
                if past_phenoRXNs:
                    del media[past_phenoRXNs[-1]]
            except Exception as e:
                print(e, f"\nEX_{metID}_e0 is not in the model {org_model.id}")
                continue
            media.update({phenoRXN.id: 100})
            pheno_util.add_medium(media)
            print(phenoRXN.id)
            pheno_util.model.solver = solver
            ### define an oxygen absorption relative to the phenotype carbon source
            # O2_consumption: EX_cpd00007_e0 <= phenotype carbon source    # formerly <= 2 * sum(primary carbon fluxes)
            coef = org_coef.copy()
            coef.update({phenoRXN.reverse_variable: 1})
            pheno_util.create_constraint(
                Constraint(Zero, lb=0, ub=None, name="EX_cpd00007_e0_limitation"),
                coef=coef,
            )

            ## minimize the influx of all carbonaceous exchanges, mostly non-phenotype compounds, at a fixed biomass growth
            min_growth = float(1)  # arbitrarily assigned minimal growth
            pheno_util.add_minimal_objective_cons(min_growth)
            phenoRXN.upper_bound = 0
            for ex in pheno_util.carbon_exchange_list():
                exMet = ex.id.replace("EX_", "").replace("_e0", "")
                if exMet in phenoRXNs and exMet != metID:
                    ex.lower_bound = 0
                # print(f"The new bounds of {exMet} exchange are: {ex.bounds}")
            pheno_util.add_objective(
                Zero,
                "min",
                coef={
                    ex.reverse_variable: 1000 if ex.id != phenoRXN.id else 1
                    for ex in pheno_util.carbon_exchange_list()
                },
            )
            # export_lp(pheno_util.model, f"minimize_cInFlux_{phenoRXN.id}")
            sol = pheno_util.model.optimize()
            if sol.status != "optimal":
                pheno_util.model.remove_cons_vars(["EX_cpd00007_e0_limitation"])
                coef.update({phenoRXN.reverse_variable: 5})
                pheno_util.create_constraint(
                    Constraint(Zero, lb=0, ub=None, name="EX_cpd00007_e0_limitation"),
                    coef=coef,
                )
                sol = pheno_util.model.optimize()
            bioFlux_check(pheno_util.model, sol)
            ### limit maximum consumption to the values from the previous minimization
            for ex in pheno_util.carbon_exchange_list():
                #### (limiting the reverse_variable is more restrictive than the net flux variable)
                if ex.id != phenoRXN.id:
                    ex.reverse_variable.ub = abs(min(0, sol.fluxes[ex.id]))

            ## maximize the phenotype yield with the previously defined growth and constraints
            pheno_util.add_objective(phenoRXN.reverse_variable, "min")
            # export_lp(pheno_util.model, f"maximize_phenoYield_{phenoRXN.id}")
            pheno_sol = pheno_util.model.optimize()
            bioFlux_check(pheno_util.model, pheno_sol)
            pheno_influx = pheno_sol.fluxes[phenoRXN.id]
            if pheno_influx >= 0:
                if not all_phenotypes:
                    print(
                        f"The phenotype carbon source has a flux of {pheno_sol.fluxes[phenoRXN.id]}."
                    )
                    pprint(
                        {
                            rxn: flux
                            for rxn, flux in pheno_sol.fluxes.items()
                            if flux != 0
                        }
                    )
                    # TODO gapfill the model in media the non-functioning carbon source
                    raise NoFluxError(
                        f"The (+) net flux of {pheno_influx} for the {phenoRXN.id} phenotype"
                        f" indicates that it is an implausible phenotype."
                    )
                print(
                    f"NoFluxError: The (+) net flux of {pheno_influx} for the {phenoRXN.id}"
                    " phenotype indicates that it is an implausible phenotype."
                )
                continue
            phenoRXN.lower_bound = phenoRXN.upper_bound = pheno_influx

            ## maximize excretion of all potential carbon byproducts whose #C's < phenotype source #C's
            phenotype_source_carbons = FBAHelper.rxn_mets_list(phenoRXN)[0].elements[
                "C"
            ]
            minimum_fluxes = {}
            for carbon_source in pheno_util.carbon_exchange_list(include_unknown=False):
                if (
                    0
                    < FBAHelper.rxn_mets_list(carbon_source)[0].elements["C"]
                    < phenotype_source_carbons
                ):
                    pheno_util.add_objective(carbon_source.flux_expression, "max")
                    minObj = pheno_util.model.slim_optimize()
                    # print(carbon_source.reaction, "\t", carbon_source.flux_expression, "\t", minObj)
                    if minObj > phenotype_flux_threshold:
                        minimum_fluxes[carbon_source.id] = minObj
            # TODO limit the possible excreted compounds to only those that are defined in the media
            excreted_compounds = list(
                [exID for exID in minimum_fluxes.keys() if exID != "EX_cpd00011_e0"]
            )
            # minimum_fluxes_df = DataFrame(data=list(minimum_fluxes.values()), index=excreted_compounds, columns=["min_flux"])
            # max_excretion_cpd = minimum_fluxes_df["minimum"].idxmin()
            ### optimize the excretion of the discovered phenotype excreta
            if "excreted" in phenoCPDs:
                phenoCPDs["excreted"] = [
                    f"EX_{cpd}_e0" for cpd in phenoCPDs["excreted"]
                ]
                phenoCPDs["excreted"].extend(excreted_compounds)
            else:
                phenoCPDs["excreted"] = excreted_compounds
            pheno_excreta = [
                pheno_util.model.reactions.get_by_id(excreta)
                for excreta in phenoCPDs["excreted"]
            ]
            pheno_util.add_objective(
                sum([ex.flux_expression for ex in pheno_excreta]), "max"
            )
            # export_lp(pheno_util.model, "maximize_excreta")
            sol = pheno_util.model.optimize()
            bioFlux_check(pheno_util.model, sol)
            for ex in pheno_excreta:
                ex.lower_bound = ex.upper_bound = sol.fluxes[ex.id]

            ## minimize flux of the total simulation flux through pFBA
            # TODO discover why some phenotypes are infeasible with pFBA
            try:
                pheno_sol = pfba(pheno_util.model)
            # pheno_util.add_objective(sum([rxn.flux_expression for rxn in pheno_util.e]), "min")
            # pheno_sol = pheno_util.model.optimize()
            except Exception as e:
                print(
                    f"The {phenoRXN.id} phenotype of the {pheno_util.model} model is "
                    f"unable to be simulated with pFBA and yields a < {e} > error."
                )
            sol_dict = FBAHelper.solution_to_variables_dict(pheno_sol, pheno_util.model)
            simulated_growth = sum(
                [
                    flux
                    for var, flux in sol_dict.items()
                    if re.search(r"(^bio\d+$)", var.name)
                ]
            )
            if not isclose(simulated_growth, min_growth):
                display(
                    [
                        (rxn, flux)
                        for rxn, flux in pheno_sol.fluxes.items()
                        if "EX_" in rxn and flux != 0
                    ]
                )
                raise ObjectiveError(
                    f"The assigned minimal_growth of {min_growth} was not optimized"
                    f" during the simulation, where the observed growth was {simulated_growth}."
                )

            ## store solution fluxes and update the community_members phenotypes
            met_name = strip_comp(name).replace(" ", "-")
            col = content["name"] + "_" + met_name
            models[pheno_util.model.id]["solutions"][col] = pheno_sol
            solutions.append(
                models[pheno_util.model.id]["solutions"][col].objective_value
            )
            met_name = met_name.replace("_", "-").replace("~", "-")
            if all_phenotypes:
                if "phenotypes" not in comm_members[org_model]:
                    comm_members[org_model]["phenotypes"] = {
                        met_name: {"consumed": [strip_comp(metID)]}
                    }
                if met_name not in comm_members[org_model]["phenotypes"]:
                    comm_members[org_model]["phenotypes"].update(
                        {met_name: {"consumed": [strip_comp(metID)]}}
                    )
                else:
                    comm_members[org_model]["phenotypes"][met_name]["consumed"] = [
                        strip_comp(metID)
                    ]
                met_pheno = content["phenotypes"][met_name]
                if (
                    "excreted" in met_pheno
                    and strip_comp(metID) in met_pheno["excreted"]
                ):
                    comm_members[org_model]["phenotypes"][met_name].update(
                        {"excreted": met_pheno}
                    )
            past_phenoRXNs.append(phenoRXN.id)

    # construct the parsed table of all exchange fluxes for each phenotype
    cols = {}
    ## biomass row
    cols["rxn"] = ["bio"]
    for content in models.values():
        for col in content["solutions"]:
            cols[col] = [0]
            if col not in content["solutions"]:
                continue
            bio_rxns = [x for x in content["solutions"][col].fluxes.index if "bio" in x]
            flux = mean(
                [
                    content["solutions"][col].fluxes[rxn]
                    for rxn in bio_rxns
                    if content["solutions"][col].fluxes[rxn] != 0
                ]
            )
            cols[col] = [flux]
    ## exchange reactions rows
    looped_cols = cols.copy()
    looped_cols.pop("rxn")
    for content in models.values():
        for ex_rxn in content["exchanges"]:
            cols["rxn"].append(ex_rxn.id)
            for col in looped_cols:
                ### reactions that are not present in the columns are ignored
                flux = (
                    0
                    if (
                        col not in content["solutions"]
                        or ex_rxn.id not in list(content["solutions"][col].fluxes.index)
                    )
                    else content["solutions"][col].fluxes[ex_rxn.id]
                )
                cols[col].append(flux)
    ## construct the DataFrame
    fluxes_df = DataFrame(data=cols)
    fluxes_df.index = fluxes_df["rxn"]
    fluxes_df.drop("rxn", axis=1, inplace=True)
    fluxes_df = fluxes_df.groupby(fluxes_df.index).sum()
    fluxes_df = fluxes_df.loc[(fluxes_df != 0).any(axis=1)]
    fluxes_df.astype(str)
    # fluxes_df.to_csv("fluxes.csv")
    return fluxes_df, comm_members
