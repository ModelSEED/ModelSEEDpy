# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
import json
from optlang.symbolics import Zero, add
from cobra import Model, Reaction, Metabolite
from cobra.io import (
    load_json_model,
    save_json_model,
    load_matlab_model,
    save_matlab_model,
    read_sbml_model,
    write_sbml_model,
)
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

logger = logging.getLogger(__name__)

base_blacklist = {
    "rxn10157": "<",
    "rxn09295": "<",
    "rxn05938": "<",
    "rxn08628": ">",
    "rxn10155": "<",
    "rxn01353": "<",
    "rxn05683": "<",
    "rxn09193": "<",
    "rxn09003": "<",
    "rxn01128": ">",
    "rxn08655": "<",
    "rxn09272": "<",
    "rxn05313": "<",
    "rxn01510": ">",
    "rxn05297": ">",
    "rxn00507": "<",
    "rxn05596": "<",
    "rxn01674": "<",
    "rxn01679": "<",
    "rxn00778": ">",
    "rxn05206": ">",
    "rxn00239": "<",
    "rxn05937": "<",
    "rxn00715": "<",
    "rxn05638": ">",
    "rxn05289": ">",
    "rxn00839": "<",
    "rxn08866": "<",
    "rxn10901": "<",
    "rxn09331": "<",
    "rxn05242": "<",
    "rxn12549": "<",
    "rxn13143": "<",
    "rxn12498": "<",
    "rxn08373": "<",
    "rxn05208": "<",
    "rxn09372": "<",
    "rxn00571": ">",
    "rxn08104": "<",
    "rxn08704": "<",
    "rxn07191": "<",
    "rxn09672": "<",
    "rxn01048": ">",
    "rxn11267": ">",
    "rxn08290": "<",
    "rxn09307": "<",
    "rxn05676": ">",
    "rxn09653": "<",
    "rxn11277": "<",
    "rxn00976": "<",
    "rxn02520": "<",
    "rxn08275": "<",
    "rxn09121": "<",
    "rxn08999": "<",
    "rxn08633": "<",
    "rxn08610": "<",
    "rxn09218": "<",
    "rxn05626": "<",
    "rxn11320": "<",
    "rxn10058": ">",
    "rxn08544": "<",
    "rxn12539": "<",
    "rxn08990": "<",
    "rxn09348": "<",
    "rxn00378": "<",
    "rxn05243": "<",
    "rxn02154": "<",
    "rxn12587": "<",
    "rxn00125": "<",
    "rxn05648": "<",
    "rxn13722": "<",
    "rxn10910": ">",
    "rxn05308": ">",
    "rxn08585": "<",
    "rxn14207": "<",
    "rxn08682": "<",
    "rxn10895": "<",
    "rxn09655": "<",
    "rxn11934": "<",
    "rxn01742": ">",
    "rxn05222": ">",
    "rxn09942": "<",
    "rxn13753": ">",
    "rxn10857": "<",
    "rxn03468": "<",
    "rxn04942": "<",
    "rxn10990": ">",
    "rxn08639": "<",
    "rxn09248": "<",
    "rxn11935": ">",
    "rxn00870": ">",
    "rxn08314": "<",
    "rxn09378": "<",
    "rxn09269": "<",
    "rxn10057": ">",
    "rxn13702": ">",
    "rxn00517": "<",
    "rxn09221": ">",
    "rxn01505": ">",
    "rxn13692": ">",
    "rxn05573": "<",
    "rxn10123": ">",
    "rxn09005": "<",
    "rxn05244": "<",
    "rxn05940": "<",
    "rxn10124": ">",
    "rxn06202": ">",
    "rxn09660": "<",
    "rxn02260": ">",
    "rxn08912": "<",
    "rxn05760": ">",
    "rxn05580": ">",
    "rxn02181": ">",
    "rxn09339": "<",
    "rxn00767": "<",
    "rxn09118": "<",
    "rxn05303": "<",
    "rxn06110": "<",
    "rxn12800": "<",
    "rxn10966": "<",
    "rxn12561": "<",
    "rxn04678": ">",
    "rxn10818": "<",
    "rxn08166": "<",
    "rxn02044": ">",
    "rxn12623": "<",
    "rxn13392": ">",
    "rxn02283": "<",
    "rxn13647": ">",
    "rxn08653": "<",
    "rxn05218": ">",
    "rxn11676": ">",
    "rxn00197": "<",
    "rxn00697": "<",
    "rxn12575": ">",
    "rxn08188": "<",
    "rxn01215": "<",
    "rxn08730": ">",
    "rxn08519": ">",
    "rxn08642": "<",
    "rxn05245": "<",
    "rxn04042": "<",
    "rxn01443": ">",
    "rxn08535": "<",
    "rxn03983": "<",
    "rxn08317": "<",
    "rxn14173": ">",
    "rxn08868": "<",
    "rxn05893": ">",
    "rxn00435": ">",
    "rxn13724": "<",
    "rxn09681": "<",
    "rxn00572": ">",
    "rxn05942": "<",
    "rxn11158": "<",
    "rxn05562": "<",
    "rxn10868": "<",
    "rxn10426": "<",
    "rxn00941": ">",
    "rxn08240": "<",
    "rxn05220": ">",
    "rxn01228": ">",
    "rxn12540": "<",
    "rxn10618": ">",
    "rxn09659": "<",
    "rxn08985": ">",
    "rxn05523": "<",
    "rxn00421": "<",
    "rxn09385": "<",
    "rxn08542": "<",
    "rxn09658": "<",
    "rxn01173": "<",
    "rxn10977": "<",
    "rxn05216": "<",
    "rxn13748": ">",
    "rxn10769": ">",
    "rxn00451": "<",
    "rxn01639": "<",
    "rxn08661": "<",
    "rxn09308": "<",
    "rxn09260": "<",
    "rxn00253": "<",
    "rxn05207": "<",
    "rxn01667": "<",
    "rxn08063": "<",
    "rxn01508": ">",
    "rxn09657": "<",
    "rxn01209": ">",
    "rxn00548": ">",
    "rxn12617": "<",
    "rxn08747": ">",
    "rxn08096": "<",
    "rxn11951": "<",
    "rxn09061": "<",
    "rxn10978": "<",
    "rxn02748": ">",
    "rxn09663": "<",
    "rxn08737": "<",
    "rxn13127": "<",
    "rxn09366": "<",
    "rxn05634": "<",
    "rxn05554": "<",
    "rxn09266": ">",
    "rxn04676": ">",
    "rxn11078": ">",
    "rxn04932": "<",
    "rxn00607": ">",
    "rxn08856": "<",
    "rxn12624": "<",
    "rxn05215": "<",
    "rxn13686": "<",
    "rxn12529": "<",
    "rxn00234": "<",
    "rxn13689": ">",
    "rxn08117": "<",
    "rxn05315": ">",
    "rxn08865": "<",
    "rxn11678": ">",
    "rxn00518": "<",
    "rxn00195": "<",
    "rxn10054": "<",
    "rxn12532": "<",
    "rxn05902": ">",
    "rxn12777": "<",
    "rxn12822": ">",
    "rxn13735": ">",
    "rxn00427": "<",
    "rxn13196": "<",
    "rxn08284": "<",
    "rxn10576": ">",
    "rxn00891": "<",
    "rxn08293": "<",
    "rxn00374": ">",
    "rxn08795": "<",
    "rxn12583": "<",
    "rxn00918": ">",
    "rxn08525": "<",
    "rxn10427": ">",
    "rxn09271": "<",
    "rxn10860": "<",
    "rxn10600": ">",
    "rxn13729": ">",
    "rxn01375": "<",
    "rxn13726": ">",
    "rxn10587": "<",
    "rxn08672": "<",
    "rxn10588": ">",
    "rxn08152": ">",
    "rxn09306": "<",
    "rxn00635": "<",
    "rxn08427": "<",
    "rxn05225": ">",
    "rxn00680": ">",
    "rxn08786": ">",
    "rxn08721": "<",
    "rxn11339": "<",
    "rxn05749": "<",
    "rxn01187": ">",
    "rxn08625": "<",
    "rxn06677": "<",
    "rxn12302": ">",
    "rxn02770": "<",
    "rxn05628": "<",
    "rxn13706": ">",
    "rxn12739": "<",
    "rxn00177": "<",
    "rxn09896": ">",
    "rxn12574": "<",
    "rxn12533": ">",
    "rxn08537": ">",
    "rxn05651": ">",
    "rxn08170": "<",
    "rxn05240": "<",
    "rxn00663": ">",
    "rxn12589": "<",
    "rxn09299": "<",
    "rxn02059": "<",
    "rxn12217": ">",
    "rxn06592": "<",
    "rxn05939": ">",
    "rxn08581": "<",
    "rxn00430": "<",
    "rxn09283": ">",
    "rxn08919": "<",
    "rxn13660": "<",
    "rxn08065": "<",
    "rxn08428": ">",
    "rxn10936": ">",
    "rxn05238": ">",
    "rxn05685": "<",
    "rxn08920": ">",
    "rxn07193": "<",
    "rxn08265": "<",
    "rxn12554": "<",
    "rxn08094": "<",
    "rxn13727": ">",
    "rxn04158": "<",
    "rxn09839": "<",
    "rxn10820": "<",
    "rxn00869": ">",
    "rxn00331": ">",
    "rxn09034": "<",
    "rxn01136": "<",
    "rxn09247": "<",
    "rxn08302": "<",
    "rxn10594": "<",
    "rxn08670": ">",
    "rxn11334": "<",
    "rxn09941": "<",
    "rxn02919": "<",
    "rxn09670": "<",
    "rxn10892": "<",
    "rxn09794": "<",
    "rxn02332": ">",
    "rxn00244": ">",
    "rxn08030": "<",
    "rxn12526": "<",
    "rxn13150": ">",
    "rxn05486": "<",
    "rxn10852": ">",
    "rxn13790": ">",
    "rxn06348": ">",
    "rxn09172": ">",
    "rxn03653": ">",
    "rxn05213": "<",
    "rxn01869": "<",
    "rxn08142": "<",
    "rxn12606": "<",
    "rxn11916": ">",
    "rxn05748": "<",
    "rxn08543": "<",
    "rxn01107": ">",
    "rxn05708": "<",
    "rxn08169": "<",
    "rxn06641": ">",
    "rxn12578": "<",
    "rxn01172": "<",
    "rxn02120": ">",
    "rxn05669": "<",
    "rxn11322": "<",
    "rxn12630": "<",
    "rxn00698": "<",
    "rxn05507": ">",
    "rxn12530": "<",
    "rxn09304": "<",
    "rxn05532": ">",
    "rxn03644": ">",
    "rxn08733": "<",
    "rxn13733": "<",
    "rxn10044": ">",
    "rxn00176": ">",
    "rxn01364": ">",
    "rxn02198": ">",
    "rxn06990": "<",
    "rxn08424": "<",
    "rxn08069": "<",
    "rxn05611": "<",
    "rxn11973": "<",
    "rxn12665": ">",
    "rxn05241": "<",
    "rxn08982": ">",
    "rxn00542": ">",
    "rxn12588": "<",
    "rxn03517": ">",
    "rxn01805": "<",
    "rxn13203": ">",
    "rxn08614": "<",
    "rxn12200": ">",
    "rxn13811": "<",
    "rxn08377": "<",
    "rxn11342": ">",
    "rxn02976": "<",
    "rxn08217": "<",
    "rxn07921": ">",
    "rxn09944": ">",
    "rxn02401": "<",
    "rxn08429": ">",
    "rxn00905": "<",
    "rxn08196": "<",
    "rxn03054": "<",
    "rxn08643": "<",
    "rxn01874": "<",
    "rxn08028": "<",
    "rxn01641": ">",
    "rxn03442": "<",
    "rxn02172": "<",
    "rxn10692": ">",
    "rxn10613": ">",
    "rxn12928": ">",
    "rxn12994": ">",
    "rxn13843": ">",
    "rxn12942": ">",
    "rxn12934": ">",
    "rxn16827": ">",
    "rxn12941": ">",
    "rxn01736": ">",
    "rxn14109": ">",
    "rxn15060": ">",
    "rxn15064": ">",
    "rxn30685": ">",
    "rxn10095": ">",
    "rxn16143": ">",
    "rxn25271": ">",
    "rxn25160": ">",
    "rxn30917": ">",
    "rxn16843": ">",
    "rxn08921": ">",
    "rxn09390": ">",
    "rxn27362": ">",
    "rxn02664": ">",
    "rxn24638": ">",
    "rxn24613": ">",
    "rxn24611": ">",
    "rxn14428": ">",
    "rxn03079": ">",
    "rxn03020": ">",
    "rxn10471": "<",
}


class GapfillingPkg(BaseFBAPkg):
    """ """

    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "gapfilling", {}, {})
        self.gapfilling_penalties = None

    def build(self, template, minimum_objective=0.01):
        parameters = {
            "default_gapfill_templates": [template],
            "gapfill_all_indecies_with_default_templates": 1,
            "minimum_obj": minimum_objective,
            "set_objective": 1,
        }
        self.build_package(parameters)

    def get_model_index_hash(self):
        """
        Determine all indices that should be gap filled
        :return:
        """
        index_hash = {"none": 0}
        for metabolite in self.model.metabolites:
            if re.search("_[a-z]\d+$", metabolite.id) is not None:
                m = re.search("_([a-z])(\d+)$", metabolite.id)
                if m[1] != "e":
                    if m[2] not in index_hash:
                        index_hash[m[2]] = 0
                    index_hash[m[2]] += 1
            else:
                index_hash["none":0]
                # Iterating over all indecies with more than 10 intracellular compounds:
        return index_hash

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {
                "auto_sink": ["cpd02701", "cpd11416", "cpd15302", "cpd03091"],
                "extend_with_template": 1,
                "model_penalty": 1,
                "default_gapfill_models": [],
                "default_gapfill_templates": [],
                "gapfill_templates_by_index": {},
                "gapfill_models_by_index": {},
                "reaction_scores": {},
                "gapfill_all_indecies_with_default_templates": 1,
                "gapfill_all_indecies_with_default_models": 1,
                "default_excretion": 100,
                "default_uptake": 100,
                "minimum_obj": 0.01,
                "set_objective": 1,
                "minimize_exchanges": False,
                "blacklist": [],
            },
        )
        # Adding model reactions to original reaction list
        self.parameters["original_reactions"] = []
        for rxn in self.model.reactions:
            if FBAHelper.is_ex(rxn):
                continue
            if FBAHelper.is_biomass(rxn):
                continue
            if rxn.lower_bound < 0:
                self.parameters["original_reactions"].append([rxn, "<"])
            if rxn.upper_bound > 0:
                self.parameters["original_reactions"].append([rxn, ">"])
        # Adding constraint for target reaction
        self.parameters["origobj"] = self.model.objective
        self.pkgmgr.getpkg("ObjConstPkg").build_package(
            self.parameters["minimum_obj"], None
        )

        # Determine all indecies that should be gapfilled
        indexhash = self.get_model_index_hash()

        # Iterating over all indecies with more than 10 intracellular compounds:
        self.gapfilling_penalties = dict()
        for index in indexhash:
            if indexhash[index] > 10:
                if index == "none":
                    for template in self.parameters["default_gapfill_templates"]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(
                            template, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(
                            gfmdl, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
                if index in self.parameters["gapfill_templates_by_index"]:
                    for template in self.parameters["gapfill_templates_by_index"][
                        index
                    ]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(
                            template, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
                if index in self.parameters["gapfill_models_by_index"]:
                    for gfmdl in self.parameters["gapfill_models_by_index"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(
                            gfmdl, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
                if self.parameters["gapfill_all_indecies_with_default_templates"]:
                    for template in self.parameters["default_gapfill_templates"]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(
                            template, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
                if self.parameters["gapfill_all_indecies_with_default_models"]:
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(
                            gfmdl, index
                        )
                        self.gapfilling_penalties.update(new_penalties)
        # Rescaling penalties by reaction scores and saving genes
        for reaction in self.gapfilling_penalties:
            rxnid = reaction.split("_")[0]
            if rxnid in self.parameters["reaction_scores"]:
                highest_score = 0
                for gene in self.parameters["reaction_scores"][rxnid]:
                    if highest_score < self.parameters["reaction_scores"][rxnid][gene]:
                        highest_score = self.parameters["reaction_scores"][rxnid][gene]
                factor = 0.1
                if "reverse" in self.gapfilling_penalties[reaction]:
                    self.gapfilling_penalties[reaction]["reverse"] = (
                        factor * self.gapfilling_penalties[reaction]["reverse"]
                    )
                if "forward" in self.gapfilling_penalties[reaction]:
                    self.gapfilling_penalties[reaction]["forward"] = (
                        factor * self.gapfilling_penalties[reaction]["forward"]
                    )

        self.model.solver.update()
        if self.parameters["set_objective"] == 1:
            reaction_objective = self.model.problem.Objective(Zero, direction="min")
            obj_coef = dict()
            for reaction in self.model.reactions:
                if reaction.id in self.gapfilling_penalties:
                    if (
                        self.parameters["minimize_exchanges"]
                        or reaction.id[0:3] != "EX_"
                    ):
                        # Minimizing gapfilled reactions
                        if "reverse" in self.gapfilling_penalties[reaction.id]:
                            obj_coef[reaction.reverse_variable] = abs(
                                self.gapfilling_penalties[reaction.id]["reverse"]
                            )
                        if "forward" in self.gapfilling_penalties[reaction.id]:
                            obj_coef[reaction.forward_variable] = abs(
                                self.gapfilling_penalties[reaction.id]["forward"]
                            )
                else:
                    obj_coef[reaction.forward_variable] = 0
                    obj_coef[reaction.reverse_variable] = 0
            self.model.objective = reaction_objective
            reaction_objective.set_linear_coefficients(obj_coef)

    def extend_model_with_model_for_gapfilling(self, source_model, index):
        new_metabolites = {}
        new_reactions = {}
        new_exchange = []
        new_demand = []
        new_penalties = dict()
        local_remap = {}
        # Adding metabolites from source model to gapfill model
        for cobra_metabolite in source_model.metabolites:
            original_id = cobra_metabolite.id
            if re.search("(.+)_([a-z])\d+$", cobra_metabolite.id) != None:
                m = re.search("(.+)_([a-z])\d+$", cobra_metabolite.id)
                if m[2] == "e":
                    cobra_metabolite.compartment = "e0"
                    cobra_metabolite.id = m[1] + "_e0"
                else:
                    cobra_metabolite.compartment = m[2] + index
                    cobra_metabolite.id = m[1] + "_" + m[2] + index
                if (
                    cobra_metabolite.id not in self.model.metabolites
                    and cobra_metabolite.id not in new_metabolites
                ):
                    new_metabolites[cobra_metabolite.id] = cobra_metabolite
                    local_remap[original_id] = cobra_metabolite
                    if m[1] + "_" + m[2] in self.parameters["auto_sink"]:
                        new_demand.append(cobra_metabolite)
                    if m[2] == "e":
                        new_exchange.append(cobra_metabolite)
        # Adding all metabolites to model prior to adding reactions
        self.model.add_metabolites(new_metabolites.values())
        # Adding reactions from source model to gapfill model
        for modelreaction in source_model.reactions:
            if re.search("(.+)_([a-z])\d+$", modelreaction.id) != None:
                m = re.search("(.+)_([a-z])\d+$", modelreaction.id)
                if m[1] not in self.parameters["blacklist"]:
                    if m[1] in base_blacklist:
                        if base_blacklist[m[1]] == ">" or base_blacklist[m[1]] == "=":
                            cobra_reaction.upper_bound = 0
                        if base_blacklist[m[1]] == "<" or base_blacklist[m[1]] == "=":
                            cobra_reaction.lower_bound = 0
                    cobra_reaction = modelreaction.copy()
                    cobra_reaction.id = groups[1] + "_" + groups[2] + index
                    if (
                        cobra_reaction.id not in self.model.reactions
                        and cobra_reaction.id not in new_reactions
                    ):
                        new_reactions[cobra_reaction.id] = cobra_reaction
                        new_penalties[cobra_reaction.id] = dict()
                        new_penalties[cobra_reaction.id]["added"] = 1
                        if cobra_reaction.lower_bound < 0:
                            new_penalties[cobra_reaction.id][
                                "reverse"
                            ] = self.parameters["model_penalty"]
                        if cobra_reaction.upper_bound > 0:
                            new_penalties[cobra_reaction.id][
                                "forward"
                            ] = self.parameters["model_penalty"]
                        # Updating metabolites in reaction to new model
                        metabolites = cobra_reaction.metabolites
                        new_stoichiometry = {}
                        for metabolite in metabolites:
                            # Adding new coefficient:
                            new_stoichiometry[local_remap[metabolite.id]] = metabolites[
                                metabolite
                            ]
                            # Zeroing out current coefficients
                            if local_remap[metabolite.id] != metabolite:
                                new_stoichiometry[metabolite] = 0
                        cobra_reaction.add_metabolites(new_stoichiometry, combine=False)
                    elif (
                        cobra_reaction.lower_bound < 0
                        and self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).lower_bound
                        == 0
                    ):
                        self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).lower_bound = cobra_reaction.lower_bound
                        self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["reverse"] = self.parameters[
                            "model_penalty"
                        ]
                        new_penalties[cobra_reaction.id]["reversed"] = 1
                    elif (
                        cobra_reaction.upper_bound > 0
                        and self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).upper_bound
                        == 0
                    ):
                        self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).upper_bound = cobra_reaction.upper_bound
                        self.model.reactions.get_by_id(
                            cobra_reaction.id
                        ).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["forward"] = model_penalty
                        new_penalties[cobra_reaction.id]["reversed"] = 1

                        # Only run this on new exchanges so we don't readd for all exchanges
        self.modelutl.add_exchanges_for_metabolites(
            new_exchange,
            self.parameters["default_uptake"],
            self.parameters["default_excretion"],
        )
        # Only run this on new demands so we don't readd for all exchanges
        self.modelutl.add_exchanges_for_metabolites(
            new_demand,
            self.parameters["default_uptake"],
            self.parameters["default_excretion"],
            "DM_",
        )
        # Adding all new reactions to the model at once (much faster than one at a time)
        self.model.add_reactions(new_reactions.values())
        return new_penalties

    def extend_model_with_template_metabolites(self, template, index="0"):
        """
        Add all missing template compartment compounds to the model
        :param template:
        :param index:
        :return:
        """
        new_metabolites = {}
        new_exchange = []
        new_demand = []
        for template_compound in template.compcompounds:
            compartment = template_compound.compartment
            compartment_index = "0" if compartment == "e" else index
            cobra_metabolite = template_compound.to_metabolite(compartment_index)
            # cobra_metabolite = self.convert_template_compound(template_compound, compartment_index, template)  # TODO: move function out
            if (
                cobra_metabolite.id not in self.model.metabolites
                and cobra_metabolite.id not in new_metabolites
            ):
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                # self.model.add_metabolites([cobra_metabolite])
                msid = FBAHelper.modelseed_id_from_cobra_metabolite(cobra_metabolite)
                if msid in self.parameters["auto_sink"]:
                    if msid != "cpd11416" or cobra_metabolite.compartment == "c0":
                        new_demand.append(cobra_metabolite)
                if compartment == "e":
                    new_exchange.append(cobra_metabolite)
        # Adding all metabolites to model prior to adding reactions
        self.model.add_metabolites(new_metabolites.values())

        return new_exchange, new_demand

    # Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def extend_model_with_template_for_gapfilling(self, template, index):
        logger.debug(f"extend model with template: {template}, index: {index}")

        new_reactions = {}
        new_penalties = dict()

        # Adding all metabolites to model prior to adding reactions
        new_exchange, new_demand = self.extend_model_with_template_metabolites(
            template, index
        )

        for template_reaction in template.reactions:
            if template_reaction.reference_id in self.parameters["blacklist"]:
                continue
            cobra_reaction = self.convert_template_reaction(
                template_reaction, index, template, 1
            )  # TODO: move function out
            if template_reaction.reference_id in base_blacklist:
                if (
                    base_blacklist[template_reaction.reference_id] == ">"
                    or base_blacklist[template_reaction.reference_id] == "="
                ):
                    cobra_reaction.upper_bound = 0
                if (
                    base_blacklist[template_reaction.reference_id] == "<"
                    or base_blacklist[template_reaction.reference_id] == "="
                ):
                    cobra_reaction.lower_bound = 0
            new_penalties[cobra_reaction.id] = dict()
            if (
                cobra_reaction.id not in self.model.reactions
                and cobra_reaction.id not in new_reactions
            ):
                # Adding any template reactions missing from the present model
                new_reactions[cobra_reaction.id] = cobra_reaction
                if cobra_reaction.lower_bound < 0:
                    new_penalties[cobra_reaction.id]["reverse"] = (
                        template_reaction.base_cost + template_reaction.reverse_penalty
                    )
                if cobra_reaction.upper_bound > 0:
                    new_penalties[cobra_reaction.id]["forward"] = (
                        template_reaction.base_cost + template_reaction.forward_penalty
                    )
                new_penalties[cobra_reaction.id]["added"] = 1
            elif template_reaction.GapfillDirection == "=":
                # Adjusting directionality as needed for existing reactions
                model_reaction = self.model.reactions.get_by_id(cobra_reaction.id)
                new_penalties[cobra_reaction.id]["reversed"] = 1
                if model_reaction.lower_bound == 0:
                    model_reaction.lower_bound = template_reaction.lower_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_reaction.id]["reverse"] = (
                        template_reaction.base_cost + template_reaction.reverse_penalty
                    )
                if model_reaction.upper_bound == 0:
                    model_reaction.upper_bound = template_reaction.upper_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_reaction.id]["forward"] = (
                        template_reaction.base_cost + template_reaction.forward_penalty
                    )
        # Only run this on new exchanges so we don't read for all exchanges
        exchanges = self.modelutl.add_exchanges_for_metabolites(
            new_exchange,
            self.parameters["default_uptake"],
            self.parameters["default_excretion"],
        )
        for ex in exchanges:
            new_penalties[ex.id] = {"added": 1, "reverse": 1, "forward": 1}

        # Only run this on new demands so we don't readd for all exchanges
        exchanges = self.modelutl.add_exchanges_for_metabolites(
            new_demand,
            self.parameters["default_uptake"],
            self.parameters["default_excretion"],
            "DM_",
        )
        for ex in exchanges:
            new_penalties[ex.id] = {"added": 1, "reverse": 1, "forward": 1}

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.model.add_reactions(new_reactions.values())
        return new_penalties

    def convert_template_compound(self, template_compound, index, template):
        base_id = template_compound.id.split("_")[0]
        base_compound = template.compounds.get_by_id(base_id)
        new_id = template_compound.id
        new_id += str(index)
        compartment = template_compound.compartment
        compartment += str(index)

        met = Metabolite(
            new_id,
            formula=base_compound.formula,
            name=base_compound.name,
            charge=template_compound.charge,
            compartment=compartment,
        )

        met.annotation[
            "sbo"
        ] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        met.annotation["seed.compound"] = base_id
        return met

    def convert_template_reaction(
        self, template_reaction, index, template, for_gapfilling=1
    ):
        array = template_reaction.id.split("_")
        base_id = array[0]
        new_id = template_reaction.id
        new_id += str(index)

        lower_bound = template_reaction.lower_bound
        upper_bound = template_reaction.upper_bound

        direction = template_reaction.GapfillDirection
        if for_gapfilling == 0:
            direction = template_reaction.direction

        if direction == ">":
            lower_bound = 0
        elif direction == "<":
            upper_bound = 0

        cobra_reaction = Reaction(
            new_id,
            name=template_reaction.name,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )

        object_stoichiometry = {}
        for m, value in template_reaction.metabolites.items():
            metabolite_id = m.id
            template_compound = template.compcompounds.get_by_id(m.id)
            compartment = template_compound.compartment
            if compartment == "e":
                metabolite_id = m.id + "0"
            else:
                metabolite_id = m.id + str(index)
            metabolite = self.model.metabolites.get_by_id(metabolite_id)
            object_stoichiometry[metabolite] = value

        cobra_reaction.add_metabolites(object_stoichiometry)

        cobra_reaction.annotation["sbo"] = "SBO:0000176"  # biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = template_reaction.reference_id

        return cobra_reaction

    def binary_check_gapfilling_solution(self, solution=None, flux_values=None):
        if solution is None:
            solution = self.compute_gapfilled_solution()
        if flux_values is None:
            flux_values = self.modelutl.compute_flux_values_from_variables()
        filter = {}
        for rxn_id in solution["reversed"]:
            filter[rxn_id] = solution["reversed"][rxn_id]
        for rxn_id in solution["new"]:
            filter[rxn_id] = solution["new"][rxn_id]
        self.pkgmgr.getpkg("ReactionUsePkg").build_package(filter)
        objcoef = {}
        for rxnid in filter:
            if filter[rxnid] == ">":
                objcoef[self.pkgmgr.getpkg("ReactionUsePkg").variables["fu"][rxnid]] = 1
            if filter[rxnid] == "<":
                objcoef[self.pkgmgr.getpkg("ReactionUsePkg").variables["ru"][rxnid]] = 1
        new_solution = {}
        with self.model:
            # Setting all gapfilled reactions not in the solution to zero
            self.knockout_gf_reactions_outside_solution(solution, flux_values)
            # Setting the objective to be minimization of sum of binary variables
            min_reaction_objective = self.model.problem.Objective(Zero, direction="min")
            self.model.objective = min_reaction_objective
            min_reaction_objective.set_linear_coefficients(objcoef)
            self.model.optimize()
            new_solution = self.compute_gapfilled_solution()
        return new_solution

    def knockout_gf_reactions_outside_solution(self, solution=None, flux_values=None):
        """
        This function is designed to KO all gap filled reactions not included in the solution
        """
        if solution == None:
            solution = self.compute_gapfilled_solution()
        if flux_values == None:
            flux_values = self.modelutl.compute_flux_values_from_variables()
        for rxnobj in self.model.reactions:
            if rxnobj.id in self.gapfilling_penalties:
                if (
                    "reverse" in self.gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["reverse"] <= Zero
                ):
                    rxnobj.lower_bound = 0
                if (
                    "forward" in self.gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["forward"] <= Zero
                ):
                    rxnobj.upper_bound = 0
                rxnobj.update_variable_bounds()

    def run_test_conditions(self, condition_list, solution=None, max_iterations=10):
        if solution == None:
            solution = self.compute_gapfilled_solution()
        reaction_list = []
        for rxnid in solution["reversed"]:
            reaction_list.append(
                [self.model.reactions.get_by_id(rxnid), solution["reversed"][rxnid]]
            )
        for rxnid in solution["new"]:
            reaction_list.append(
                [self.model.reactions.get_by_id(rxnid), solution["new"][rxnid]]
            )
        filtered_list = []
        with self.model:
            # Setting all gapfilled reactions not in the solution to zero
            self.knockout_gf_reactions_outside_solution(solution)
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            for condition in condition_list:
                condition["change"] = True
            filtered_list = self.modelutl.reaction_expansion_test(
                reaction_list, condition_list
            )
            for condition in condition_list:
                condition["change"] = False
        if len(filtered_list) > 0:
            if max_iterations > 0:
                print("Gapfilling test failed " + str(11 - max_iterations))
                # Forcing filtered reactions to zero
                for item in filtered_list:
                    if item[1] == ">":
                        self.model.reactions.get_by_id(item[0].id).upper_bound = 0
                    else:
                        self.model.reactions.get_by_id(item[0].id).lower_bound = 0
                # Restoring lower bound on biomass constraint
                self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"][
                    "1"
                ].lb = self.parameters["minimum_obj"]
                # Reoptimizing
                self.model.optimize()
                return self.run_test_conditions(
                    condition_list, None, max_iterations - 1
                )
            return None
        return solution

    def filter_database_based_on_tests(self, test_conditions):
        # Preserving the gapfilling objective function
        gfobj = self.model.objective
        # Setting the minimal growth constraint to zero
        self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
        # Setting the objective to the original default objective for the model
        self.model.objective = self.parameters["origobj"]
        # Testing if the minimal objective can be achieved before filtering
        solution = self.model.optimize()
        print(
            "Objective before filtering:",
            solution.objective_value,
            "; min objective:",
            self.parameters["minimum_obj"],
        )
        with open("debuggf.lp", "w") as out:
            out.write(str(self.model.solver))
        if solution.objective_value < self.parameters["minimum_obj"]:
            save_json_model(self.model, "gfdebugmdl.json")
            logger.critical(
                "Model cannot achieve the minimum objective even before filtering!"
            )
        # Filtering the database of any reactions that violate the specified tests
        filetered_list = []
        with self.model:
            rxnlist = []
            for reaction in self.model.reactions:
                if reaction.id in self.gapfilling_penalties:
                    if "reverse" in self.gapfilling_penalties[reaction.id]:
                        rxnlist.append([reaction, "<"])
                    if "forward" in self.gapfilling_penalties[reaction.id]:
                        rxnlist.append([reaction, ">"])

            filtered_list = self.modelutl.reaction_expansion_test(
                rxnlist, test_conditions
            )
        # Now constraining filtered reactions to zero
        for item in filtered_list:
            logger.debug("Filtering:", item[0].id, item[1])
            if item[1] == ">":
                self.model.reactions.get_by_id(item[0].id).upper_bound = 0
            else:
                self.model.reactions.get_by_id(item[0].id).lower_bound = 0
        # Now testing if the gapfilling minimum objective can still be achieved
        solution = self.model.optimize()
        print(
            "Objective after filtering:",
            solution.objective_value,
            "; min objective:",
            self.parameters["minimum_obj"],
        )
        # Now we need to restore a minimal set of filtered reactions such that we permit the minimum objective to be reached
        if solution.objective_value < self.parameters["minimum_obj"]:
            # Restoring the minimum objective constraint
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"][
                "1"
            ].lb = self.parameters["minimum_obj"]
            new_objective = self.model.problem.Objective(Zero, direction="min")
            filterobjcoef = dict()
            for item in filtered_list:
                rxn = self.model.reactions.get_by_id(item[0].id)
                if item[1] == ">":
                    filterobjcoef[rxn.forward_variable] = item[3]
                    rxn.upper_bound = item[2]
                else:
                    filterobjcoef[rxn.reverse_variable] = item[3]
                    rxn.lower_bound = item[2]
            self.model.objective = new_objective
            new_objective.set_linear_coefficients(filterobjcoef)
            solution = self.model.optimize()
            count = len(filtered_list)
            for item in filtered_list:
                rxn = self.model.reactions.get_by_id(item[0].id)
                if solution.fluxes[rxn.id] > 0.0000001:
                    if item[1] == "<":
                        count += -1
                        rxn.lower_bound = 0
                elif solution.fluxes[rxn.id] < -0.0000001:
                    if item[1] == ">":
                        count += -1
                        rxn.upper_bound = 0
                else:
                    if item[1] == ">":
                        count += -1
                        rxn.upper_bound = 0
                    else:
                        count += -1
                        rxn.lower_bound = 0
            print("Reactions unfiltered:", count)
            # Checking for model reactions that can be removed to enable all tests to pass
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            filtered_list = self.modelutl.reaction_expansion_test(
                self.parameters["original_reactions"], test_conditions
            )
            for item in filtered_list:
                logger.debug("Filtering:", item[0].id, item[1])
                if item[1] == ">":
                    self.model.reactions.get_by_id(item[0].id).upper_bound = 0
                else:
                    self.model.reactions.get_by_id(item[0].id).lower_bound = 0
        # Restoring gapfilling objective function and minimal objective constraint
        self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = self.parameters[
            "minimum_obj"
        ]
        self.model.objective = gfobj

    def compute_gapfilled_solution(self, flux_values=None):
        if flux_values is None:
            flux_values = self.modelutl.compute_flux_values_from_variables()
        output = {"reversed": {}, "new": {}}
        for reaction in self.model.reactions:
            if reaction.id in self.gapfilling_penalties:
                if (
                    flux_values[reaction.id]["forward"] > Zero
                    and "forward" in self.gapfilling_penalties[reaction.id]
                ):
                    if "added" in self.gapfilling_penalties[reaction.id]:
                        output["new"][reaction.id] = ">"
                    else:
                        output["reversed"][reaction.id] = ">"
                elif (
                    flux_values[reaction.id]["reverse"] > Zero
                    and "reverse" in self.gapfilling_penalties[reaction.id]
                ):
                    if "added" in self.gapfilling_penalties[reaction.id]:
                        output["new"][reaction.id] = "<"
                    else:
                        output["reversed"][reaction.id] = "<"
        return output
