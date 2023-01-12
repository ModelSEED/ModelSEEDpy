# -*- coding: utf-8 -*-
from __future__ import absolute_import
import logging
from chemicals import periodic_table
import re
from cobra.core import (
    Gene,
    Metabolite,
    Model,
    Reaction,
)  # !!! Gene, Metabolite, and Model are never used
from cobra.util import solver as sutil  # !!! sutil is never used
import time
from scipy.odr.odrpack import Output  # !!! Output is never used
from chemw import ChemMW
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
            if not metabolite.formula:
                return 0
            formula = re.sub("R\d*", "", metabolite.formula)
            chem_mw = ChemMW(printing=False)
            chem_mw.mass(formula)
            return chem_mw.raw_mw
        except:
            warn(
                "The compound "
                + metabolite.id
                + " possesses an unconventional formula {metabolite.formula}; hence, the MW cannot be computed."
            )
            return 0

    @staticmethod
    def elemental_mass():
        return {element.symbol: element.MW for element in periodic_table}

    @staticmethod
    def get_modelseed_db_api(modelseed_path):
        from modelseedpy.biochem import from_local

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
        cytosol = othercomp = None
        for comp in compartments:
            if comp[0:1] == "c":
                cytosol = comp
            elif comp[0:1] != "e":
                othercomp = comp
        return othercomp or cytosol

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
        reframed_model.compartments.e0.external = True
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
    def parse_df(df):
        from numpy import array

        return array(
            dtype=object, object=[array(df.index), array(df.columns), df.to_numpy()]
        )
