# -*- coding: utf-8 -*-
from __future__ import absolute_import

import logging
from chemicals import periodic_table
import re
from optlang import Objective
from cobra.core import Gene, Metabolite, Model, Reaction   # !!! Gene and Model are never used
from optlang import Model
from cobra.util import solver as sutil  # !!! sutil is never used
import time
from modelseedpy.biochem import from_local
from optlang import Constraint
from scipy.odr import Output  # !!! Output is never used
from typing import Iterable
from chemw import ChemMW
from numpy import nan
from warnings import warn

# from Carbon.Aliases import false

logger = logging.getLogger(__name__)


class FBAHelper:
    @staticmethod
    def add_autodrain_reactions_to_community_model(
        model, auto_sink=["cpd02701", "cpd15302"]
    ):
        # Adding missing drains in the base model
        drain_reactions = []
        for metabolite in model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid in auto_sink:
                if metabolite.compartment == "c0":
                    met_id = metabolite.id
                    if all(
                        [
                            rxn not in model.reactions
                            for rxn in [f"EX_{met_id}", f"DM_{met_id}", f"SK_{met_id}"]
                        ]
                    ):
                        drain_reaction = FBAHelper.add_drain_from_metabolite_id(
                            model, metabolite.id, 0, 100, "DM_"
                        )
                        if not drain_reaction:
                            logger.info("Adding " + met_id + " DM")
                            drain_reactions.append(drain_reaction)
        model.add_reactions(drain_reactions)

    @staticmethod
    def add_drain_from_metabolite_id(
        model, cpd_id, uptake, excretion, prefix="EX_", prefix_name="Exchange for "
    ):
        """
        :param model:
        :param cpd_id:
        :param uptake:
        :param excretion:
        :param prefix:
        :param prefix_name:
        :return:
        """
        if cpd_id in model.metabolites:
            cobra_metabolite = model.metabolites.get_by_id(cpd_id)
            drain_reaction = Reaction(
                id=f"{prefix}{cpd_id}",
                name=prefix_name + cobra_metabolite.name,
                lower_bound=-uptake,
                upper_bound=excretion,
            )
            drain_reaction.add_metabolites({cobra_metabolite: -1})
            drain_reaction.annotation["sbo"] = "SBO:0000627"
            # model.add_reactions([drain_reaction])
            return drain_reaction
        return None
    
    @staticmethod
    def test_condition_list(model, condition_list, pkgmgr):
        for condition in condition_list:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(condition["media"])
            model.objective = condition["objective"]
            if condition["is_max_threshold"]:
                model.objective.direction = "max"
            else:
                model.objective.direction = "min"
            objective = model.slim_optimize()
            if model.solver.status != 'optimal':
                with open("debug.lp", 'w') as out:
                    out.write(str(model.solver))
                    out.close()
                logger.critical("Infeasible problem - LP file printed to debug!")
                return False
            if objective >= condition["threshold"] and condition["is_max_threshold"]:
                logger.info("FAILED")
                return False
            elif objective <= condition["threshold"] and not condition["is_max_threshold"]:
                logger.info("FAILED")
                return False
        return True
        
    @staticmethod
    def reaction_expansion_test(model, reaction_list, condition_list, pkgmgr):
        # First knockout all reactions in the input list and save original bounds
        original_bound = []
        for item in reaction_list:
            if item[1] == ">":
                original_bound.append(item[0].upper_bound)
                item[0].upper_bound = 0
            else:
                original_bound.append(item[0].lower_bound)
                item[0].lower_bound = 0
        # Now restore reactions one at a time
        filtered_list = []
        for index, item in enumerate(reaction_list):
            logger.info("Testing "+item[0].id)
            if item[1] == ">":
                item[0].upper_bound = original_bound[index]
                if not FBAHelper.test_condition_list(model, condition_list, pkgmgr):
                    item[0].upper_bound = 0
                    filtered_list.append(item)
            else:
                item[0].lower_bound = original_bound[index]
                if not FBAHelper.test_condition_list(model, condition_list, pkgmgr):
                    item[0].lower_bound = 0
                    filtered_list.append(item)
        return filtered_list

    @staticmethod
    def set_reaction_bounds_from_direction(reaction, direction, add=False):
        if direction == "<":
            reaction.lower_bound = -100
            if not add:
                reaction.upper_bound = 0
        if direction == ">":
            reaction.upper_bound = 100
            if not add:
                reaction.lower_bound = 0
        reaction.update_variable_bounds()

    @staticmethod
    def set_objective_from_target_reaction(model, target_reaction, minimize=False):
        target_reaction = model.reactions.get_by_id(target_reaction)
        sense = "max"
        if minimize:
            sense = "min"
        model.objective = model.problem.Objective(
            target_reaction.flux_expression, direction=sense
        )
        return target_reaction

    @staticmethod
    def modelseed_id_from_cobra_metabolite(metabolite):
        if re.search("^(cpd\d+)", metabolite.id):
            m = re.search("^(cpd\d+)", metabolite.id)
            return m[1]
        # TODO: should check to see if ModelSEED ID is in the annotations for the compound
        return None

    @staticmethod
    def modelseed_id_from_cobra_reaction(reaction):
        if re.search("^(rxn\d+)", reaction.id):
            m = re.search("^(rxn\d+)", reaction.id)
            return m[1]
        # TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None

    @staticmethod
    def metabolite_mw(metabolite):
        try:
            chem_mw = ChemMW()
            chem_mw.mass(metabolite.formula)
            return chem_mw.raw_mw
        except:
            warn(
                "The compound "
                + metabolite.id
                + " possesses an unconventional formula {metabolite.formula}; hence, the MW cannot be computed."
            )

    @staticmethod
    def elemental_mass():
        return {element.symbol: element.MW for element in periodic_table}

    @staticmethod
    def get_modelseed_db_api(modelseed_path):
        return from_local(modelseed_path)

    @staticmethod
    def is_ex(reaction):
        # TODO: check for SBO
        if len(reaction.id) > 3 and reaction.id[0:3] in ["EX_", "DM_", "SK_"]:
            return True
        return False

    @staticmethod
    def is_biomass(reaction):
        # TODO: check for SBO
        return reaction.id[0:3] == "bio"

    @staticmethod
    def exchange_hash(model):  #!!! This function is pointless?
        exchange_hash = {}  # !!! this variable is never used
        for reaction in model.reactions:
            if len(reaction.metabolites) == 1:
                for metabolite in reaction.metabolites:
                    (base, comp, index) = FBAHelper.parse_id(metabolite)
                    # exchange_hash[base][comp]

    @staticmethod
    def rxn_mets_list(rxn):
        return [met for met in rxn.reactants+rxn.products]

    @staticmethod
    def find_reaction(model, stoichiometry):
        reaction_strings = FBAHelper.stoichiometry_to_string(stoichiometry)
        atpstring = reaction_strings[0]
        rxn_hash = FBAHelper.rxn_hash(model)
        if atpstring in rxn_hash:
            return rxn_hash[atpstring]
        return None

    @staticmethod
    def msid_hash(model):
        output = {}
        for met in model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(met)
            if msid is not None:
                if msid not in output:
                    output[msid] = []
                output[msid].append(met)
        return output

    @staticmethod
    def sum_dict(d1,d2):
        for key, value in d1.items():
            if key in d2:
                d2[key] += value
            else:
                d2[key] = value
        return d2
    
    @staticmethod
    def rxn_hash(model): 
        output = {}
        for rxn in model.reactions:
            reaction_strings = FBAHelper.stoichiometry_to_string(rxn.metabolites)
            output[reaction_strings[0]] = [rxn, 1]
            output[reaction_strings[1]] = [rxn, -1]
        return output

    @staticmethod
    def rxn_compartment(reaction):
        compartments = list(reaction.compartments)
        if len(compartments) == 1:
            return compartments[0]
        cytosol = None
        for comp in compartments:
            if comp[0:1] == "c":
                cytosol = comp
            elif comp[0:1] != "e":
                return comp
        return cytosol
    
    @staticmethod
    def stoichiometry_to_string(stoichiometry):
        reactants, products = [], []
        for met in stoichiometry:
            stoich = stoichiometry[met]
            if not isinstance(met, str):
                met = (
                    None
                    if FBAHelper.modelseed_id_from_cobra_metabolite(met) == "cpd00067"
                    else met.id
                )
            if met:
                if stoich < 0:
                    reactants.append(met)
                else:
                    products.append(met)
        return [
            "+".join(sorted(reactants)) + "=" + "+".join(sorted(products)),
            "+".join(sorted(products)) + "=" + "+".join(sorted(reactants)),
        ]

    @staticmethod
    def add_atp_hydrolysis(model, compartment):
        # Searching for ATP hydrolysis compounds
        coefs = {
            "cpd00002": [-1, compartment],
            "cpd00001": [-1, compartment],
            "cpd00008": [1, compartment],
            "cpd00009": [1, compartment],
            "cpd00067": [1, compartment],
        }
        stoichiometry = {}
        id_hash = FBAHelper.msid_hash(model)
        for msid, content in coefs.items():
            if msid not in id_hash:
                logger.warning("Compound " + msid + " not found in model!")
                return None
            else:
                for cpd in id_hash[msid]:
                    if cpd.compartment == content[1]:
                        stoichiometry[cpd] = content[0]
        output = FBAHelper.find_reaction(model, stoichiometry)
        if (
            output and output[1] == 1
        ):  # !!! the second element of the output is 1/0 and not a direction string
            return {"reaction": output[0], "direction": ">", "new": False}
        cobra_reaction = Reaction(
            "rxn00062_" + compartment,
            name="ATP hydrolysis",
            lower_bound=0,
            upper_bound=1000,
        )
        cobra_reaction.annotation.update(
            {"sbo": "SBO:0000176", "seed.reaction": "rxn00062"}
        )  # biochemical reaction
        cobra_reaction.add_metabolites(stoichiometry)
        model.add_reactions([cobra_reaction])
        return {"reaction": cobra_reaction, "direction": ">", "new": True}

    @staticmethod
    def parse_id(cobra_obj):
        if re.search("(.+)_([a-z])(\d+)$", cobra_obj.id):
            m = re.search("(.+)_([a-z])(\d+)$", cobra_obj.id)
            return (m[1], m[2], int(m[3]))
        return None

    @staticmethod
    def id_from_ref(ref):
        array = ref.split("/")
        return array[-1]

    @staticmethod
    def medianame(media):
        if media == None:
            return "Complete"
        return media.id

    @staticmethod
    def validate_dictionary(dictionary, required_keys, optional_keys={}):
        for item in required_keys:
            if item not in dictionary:
                raise ValueError("Required key " + item + " is missing!")
        for key in optional_keys:
            if key not in dictionary:
                dictionary[key] = optional_keys[key]
        return dictionary

    @staticmethod
    def parse_media(media):
        return [cpd.id for cpd in media.data["mediacompounds"]]

    def get_reframed_model(
        kbase_model,
    ):
        from reframed import from_cobrapy

        reframed_model = from_cobrapy(kbase_model)
        if hasattr(kbase_model, "id"):
            reframed_model.id = kbase_model.id
        for comp in reframed_model.compartments:
            if 'e' in comp:
                reframed_model.compartments[comp].external = True

        return reframed_model

    @staticmethod
    def add_vars_cons(model, vars_cons):
        model.add_cons_vars(vars_cons)
        model.solver.update()
        return model

    @staticmethod
    def update_model_media(model, media):
        medium = {}
        model_reactions = [rxn.id for rxn in model.reactions]
        for cpd in media.data["mediacompounds"]:
            ex_rxn = f"EX_{cpd.id}"
            if ex_rxn not in model_reactions:
                model.add_boundary(
                    metabolite=Metabolite(id=cpd.id, name=cpd.name, compartment="e0"),
                    type="exchange",
                    lb=cpd.minFlux,
                    ub=cpd.maxFlux,
                )
            medium[ex_rxn] = cpd.maxFlux
        model.medium = medium
        return model

    @staticmethod
    def filter_cobra_set(cobra_set):
        unique_ids = set(obj.id for obj in cobra_set)
        unique_objs = set()
        for obj in cobra_set:
            if obj.id in unique_ids:
                unique_objs.add(obj)
                unique_ids.remove(obj.id)
        return unique_objs

    @staticmethod
    def get_reframed_model(
        kbase_model,
    ):
        from reframed import from_cobrapy

        reframed_model = from_cobrapy(kbase_model)
        if hasattr(kbase_model, "id"):
            reframed_model.id = kbase_model.id
        reframed_model.compartments.e0.external = True
        return reframed_model

    @staticmethod
    def parse_media(media):
        return [cpd.id for cpd in media.data['mediacompounds']]
    
    @staticmethod
    def parse_df(df, float_values=True):
        if isinstance(df, tuple):
            return df
        from collections import namedtuple
        dataframe = namedtuple("DataFrame", ("index", "columns", "values"))
        df.dropna(inplace=True)
        values = df.to_numpy()
        if float_values:
            values = values.astype("float64")
        return dataframe(list(df.index), list(df.columns), values)
    
    @staticmethod
    def add_cons_vars(model, vars_cons, sloppy=False):
        model.add_cons_vars(vars_cons, sloppy=sloppy)
        model.solver.update()
    
    @staticmethod
    def remove_cons_vars(model, vars_cons):
        model.remove_cons_vars(vars_cons)
        model.solver.update()
        
    @staticmethod
    def create_constraint(model, constraint, coef=None):
        model.add_cons_vars(constraint)
        model.solver.update()
        if coef:
            constraint.set_linear_coefficients(coef)
            model.solver.update()
    
    @staticmethod
    def add_objective(model, objective, direction="max", coef=None):
        model.objective = Objective(objective, direction=direction)
        model.solver.update()
        if coef:
            model.objective.set_linear_coefficients(coef)
            model.solver.update()

    @staticmethod
    def add_minimal_objective_cons(model, min_value=0.1, objective_expr=None):
        objective_expr = objective_expr or model.objective.expression
        FBAHelper.create_constraint(model, Constraint(objective_expr, lb=min_value, ub=None, name="min_value"))
    
    @staticmethod
    def add_exchange_to_model(model, cpd, rxnID):
        model.add_boundary(metabolite=Metabolite(id=cpd.id, name=cpd.name, compartment="e0"), 
            reaction_id=rxnID, type="exchange", lb=cpd.minFlux, ub=cpd.maxFlux)
    
    @staticmethod
    def update_model_media(model, media):
        medium = model.medium
        model_reactions = [rxn.id for rxn in model.reactions]
        for cpd in media.data["mediacompounds"]:
            ex_rxn = f"EX_{cpd.id}_e0"
            if ex_rxn not in model_reactions:
                model = FBAHelper.add_exchange_to_model(model, cpd, ex_rxn)
            medium[ex_rxn] = cpd.maxFlux
        model.medium = medium
        
        return model
    
    @staticmethod
    def filter_cobra_set(cobra_set):
        unique_ids = set(obj.id for obj in cobra_set)
        unique_objs = set()
        for obj in cobra_set:
            if obj.id in unique_ids:
                unique_objs.add(obj)
                unique_ids.remove(obj.id)
        return unique_objs
    
    @staticmethod
    def solution_to_dict(solution):
        return {key:flux for key, flux in solution.fluxes.items()}
    
    @staticmethod
    def solution_to_rxns_dict(solution, model):
        return {model.reactions.get_by_id(key):flux for key, flux in solution.fluxes.items()}
        
    @staticmethod
    def solution_to_variables_dict(solution, model):
        return {model.variables.get(key):flux for key, flux in solution.fluxes.items()}
    
    @staticmethod
    def remove_media_compounds(media_dict, compounds, printing=True):
        edited_dic = media_dict.copy()
        for cpd in compounds:
            if cpd in edited_dic:
                edited_dic.pop(cpd)
                if printing:
                    print(f"{cpd} removed")
            else:
                print(f"ERROR: The {cpd} is not located in the media.")
        return edited_dic

    @staticmethod
    def IDRxnMets(rxn):
        if not isinstance(rxn, dict):
            return {met.id: stoich for met, stoich in rxn.metabolites.items()}
        else:
            return {met.id: stoich for met, stoich in rxn.items()}

    @staticmethod
    def convert_kbase_media(kbase_media):
        return {"EX_"+exID: -bound[0] for exID, bound in kbase_media.get_media_constraints().items()}
