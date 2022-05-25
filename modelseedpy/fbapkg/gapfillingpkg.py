# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero
from cobra import Reaction, Metabolite
from modelseedpy.fbapkg import BaseFBAPkg
from modelseedpy.core import FBAHelper

logger = logging.getLogger(__name__)

default_blacklist = ["rxn12985", "rxn00238", "rxn07058", "rxn05305", "rxn00154", "rxn09037", "rxn10643",
                     "rxn11317", "rxn05254", "rxn05257", "rxn05258", "rxn05259", "rxn05264", "rxn05268",
                     "rxn05269", "rxn05270", "rxn05271", "rxn05272", "rxn05273", "rxn05274", "rxn05275",
                     "rxn05276", "rxn05277", "rxn05278", "rxn05279", "rxn05280", "rxn05281", "rxn05282",
                     "rxn05283", "rxn05284", "rxn05285", "rxn05286", "rxn05963", "rxn05964", "rxn05971",
                     "rxn05989", "rxn05990", "rxn06041", "rxn06042", "rxn06043", "rxn06044", "rxn06045",
                     "rxn06046", "rxn06079", "rxn06080", "rxn06081", "rxn06086", "rxn06087", "rxn06088",
                     "rxn06089", "rxn06090", "rxn06091", "rxn06092", "rxn06138", "rxn06139", "rxn06140",
                     "rxn06141", "rxn06145", "rxn06217", "rxn06218", "rxn06219", "rxn06220", "rxn06221",
                     "rxn06222", "rxn06223", "rxn06235", "rxn06362", "rxn06368", "rxn06378", "rxn06474",
                     "rxn06475", "rxn06502", "rxn06562", "rxn06569", "rxn06604", "rxn06702", "rxn06706",
                     "rxn06715", "rxn06803", "rxn06811", "rxn06812", "rxn06850", "rxn06901", "rxn06971",
                     "rxn06999", "rxn07123", "rxn07172", "rxn07254", "rxn07255", "rxn07269", "rxn07451",
                     "rxn09037", "rxn10018", "rxn10077", "rxn10096", "rxn10097", "rxn10098", "rxn10099",
                     "rxn10101", "rxn10102", "rxn10103", "rxn10104", "rxn10105", "rxn10106", "rxn10107",
                     "rxn10109", "rxn10111", "rxn10403", "rxn10410", "rxn10416", "rxn11313", "rxn11316",
                     "rxn11318", "rxn11353", "rxn05224", "rxn05795", "rxn05796", "rxn05797", "rxn05798",
                     "rxn05799", "rxn05801", "rxn05802", "rxn05803", "rxn05804", "rxn05805", "rxn05806",
                     "rxn05808", "rxn05812", "rxn05815", "rxn05832", "rxn05836", "rxn05851", "rxn05857",
                     "rxn05869", "rxn05870", "rxn05884", "rxn05888", "rxn05896", "rxn05898", "rxn05900",
                     "rxn05903", "rxn05904", "rxn05905", "rxn05911", "rxn05921", "rxn05925", "rxn05936",
                     "rxn05947", "rxn05956", "rxn05959", "rxn05960", "rxn05980", "rxn05991", "rxn05992",
                     "rxn05999", "rxn06001", "rxn06014", "rxn06017", "rxn06021", "rxn06026", "rxn06027",
                     "rxn06034", "rxn06048", "rxn06052", "rxn06053", "rxn06054", "rxn06057", "rxn06059",
                     "rxn06061", "rxn06102", "rxn06103", "rxn06127", "rxn06128", "rxn06129", "rxn06130",
                     "rxn06131", "rxn06132", "rxn06137", "rxn06146", "rxn06161", "rxn06167", "rxn06172",
                     "rxn06174", "rxn06175", "rxn06187", "rxn06189", "rxn06203", "rxn06204", "rxn06246",
                     "rxn06261", "rxn06265", "rxn06266", "rxn06286", "rxn06291", "rxn06294", "rxn06310",
                     "rxn06320", "rxn06327", "rxn06334", "rxn06337", "rxn06339", "rxn06342", "rxn06343",
                     "rxn06350", "rxn06352", "rxn06358", "rxn06361", "rxn06369", "rxn06380", "rxn06395",
                     "rxn06415", "rxn06419", "rxn06420", "rxn06421", "rxn06423", "rxn06450", "rxn06457",
                     "rxn06463", "rxn06464", "rxn06466", "rxn06471", "rxn06482", "rxn06483", "rxn06486",
                     "rxn06492", "rxn06497", "rxn06498", "rxn06501", "rxn06505", "rxn06506", "rxn06521",
                     "rxn06534", "rxn06580", "rxn06585", "rxn06593", "rxn06609", "rxn06613", "rxn06654",
                     "rxn06667", "rxn06676", "rxn06693", "rxn06730", "rxn06746", "rxn06762", "rxn06779",
                     "rxn06790", "rxn06791", "rxn06792", "rxn06793", "rxn06794", "rxn06795", "rxn06796",
                     "rxn06797", "rxn06821", "rxn06826", "rxn06827", "rxn06829", "rxn06839", "rxn06841",
                     "rxn06842", "rxn06851", "rxn06866", "rxn06867", "rxn06873", "rxn06885", "rxn06891",
                     "rxn06892", "rxn06896", "rxn06938", "rxn06939", "rxn06944", "rxn06951", "rxn06952",
                     "rxn06955", "rxn06957", "rxn06960", "rxn06964", "rxn06965", "rxn07086", "rxn07097",
                     "rxn07103", "rxn07104", "rxn07105", "rxn07106", "rxn07107", "rxn07109", "rxn07119",
                     "rxn07179", "rxn07186", "rxn07187", "rxn07188", "rxn07195", "rxn07196", "rxn07197",
                     "rxn07198", "rxn07201", "rxn07205", "rxn07206", "rxn07210", "rxn07244", "rxn07245",
                     "rxn07253", "rxn07275", "rxn07299", "rxn07302", "rxn07651", "rxn07723", "rxn07736",
                     "rxn07878", "rxn11417", "rxn11582", "rxn11593", "rxn11597", "rxn11615", "rxn11617",
                     "rxn11619", "rxn11620", "rxn11624", "rxn11626", "rxn11638", "rxn11648", "rxn11651",
                     "rxn11665", "rxn11666", "rxn11667", "rxn11698", "rxn11983", "rxn11986", "rxn11994",
                     "rxn12006", "rxn12007", "rxn12014", "rxn12017", "rxn12022", "rxn12160", "rxn12161",
                     "rxn01267", "rxn05294", "rxn04656"]


class GapfillingPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "gapfilling", {}, {})
        self.gapfilling_penalties = None

    def build(self, template, minimum_objective=0.01):
        parameters = {
            "default_gapfill_templates": [template],
            "gapfill_all_indecies_with_default_templates": 1,
            "minimum_obj": minimum_objective,
            "set_objective": 1
        }
        self.build_package(parameters)

    def build_package(self, parameters):
        self.validate_parameters(parameters, [], {
            "auto_sink": ["cpd02701", "cpd11416", "cpd15302"],
            "extend_with_template":1,
            "model_penalty":1,
            "default_gapfill_models":[],
            "default_gapfill_templates":[],
            "gapfill_templates_by_index":{},
            "gapfill_models_by_index":{},
            "reaction_scores":{},
            "gapfill_all_indecies_with_default_templates":1,
            "gapfill_all_indecies_with_default_models":1,
            "default_excretion":100,
            "default_uptake":-100,
            "minimum_obj":0.01,
            "set_objective":1,
            "blacklist":default_blacklist
        })
        # Adding constraint for target reaction
        self.pkgmgr.getpkg("ObjConstPkg").build_package(self.parameters["minimum_obj"],None)
                
        # Determine all indecies that should be gapfilled
        indexhash = self._get_model_index_hash()

        # Iterating over all indecies with more than 10 intracellular compounds:
        self.gapfilling_penalties = dict()
        for index, val in indexhash.items():
            if val > 10:
                if index == "none":
                    for template in self.parameters["default_gapfill_templates"]:
                        self.gapfilling_penalties.update(self.extend_model_with_template_for_gapfilling(template, index))
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        self.gapfilling_penalties.update(self.extend_model_with_model_for_gapfilling(gfmdl, index))
                if index in self.parameters["gapfill_templates_by_index"]:
                    for template in self.parameters["gapfill_templates_by_index"][index]:
                        self.gapfilling_penalties.update(self.extend_model_with_template_for_gapfilling(template, index))
                if index in self.parameters["gapfill_models_by_index"]:
                    for gfmdl in self.parameters["gapfill_models_by_index"]:
                        self.gapfilling_penalties.update(self.extend_model_with_model_for_gapfilling(gfmdl, index))
                if self.parameters["gapfill_all_indecies_with_default_templates"]:
                    for template in self.parameters["default_gapfill_templates"]:
                        self.gapfilling_penalties.update(self.extend_model_with_template_for_gapfilling(template, index))
                if self.parameters["gapfill_all_indecies_with_default_models"]:
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        self.gapfilling_penalties.update(self.extend_model_with_model_for_gapfilling(gfmdl, index))
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
                    self.gapfilling_penalties[reaction]["reverse"] = factor * self.gapfilling_penalties[reaction]["reverse"]
                if "forward" in self.gapfilling_penalties[reaction]:
                    self.gapfilling_penalties[reaction]["forward"] = factor * self.gapfilling_penalties[reaction]["forward"]

        self.model.solver.update()
        if self.parameters["set_objective"] == 1:
            reaction_objective = self.model.problem.Objective(Zero, direction="min")
            obj_coef = dict()
            for reaction in self.model.reactions:
                if reaction.id in self.gapfilling_penalties:
                    # Minimizing gapfilled reactions
                    if "reverse" in self.gapfilling_penalties[reaction.id]:
                        obj_coef[reaction.reverse_variable] = abs(self.gapfilling_penalties[reaction.id]["reverse"])
                    # elif default_penalty != 0:
                    #    obj_coef[reaction.reverse_variable] = 0
                    if "forward" in self.gapfilling_penalties[reaction.id]:
                        obj_coef[reaction.forward_variable] = abs(self.gapfilling_penalties[reaction.id]["forward"])
                    # elif default_penalty != 0:
                    #    obj_coef[reaction.forward_variable] = 0
                else:
                    obj_coef[reaction.forward_variable] = obj_coef[reaction.reverse_variable] = 0
            self.model.objective = reaction_objective
            reaction_objective.set_linear_coefficients(obj_coef)

    def _get_model_index_hash(self):
        """
        Determine all indices that should be gap filled
        :return:
        """
        index_hash = {"none": 0}
        for metabolite in self.model.metabolites:
            if re.search('_[a-z]\d+$', metabolite.id) is not None:
                m = re.search('_([a-z])(\d+)$', metabolite.id)
                if m[1] != "e":
                    if m[2] not in index_hash:
                        index_hash[m[2]] = 0
                    index_hash[m[2]] += 1
            else:
                index_hash["none":0]
                # Iterating over all indecies with more than 10 intracellular compounds:
        return index_hash

    def extend_model_with_model_for_gapfilling(self, source_model, index):
        self.new_metabolites, self.new_reactions, local_remap, new_penalties = {}, {}, {}, {}
        new_exchange, new_demand = [], []
        # Adding metabolites from source model to gapfill model
        for cobra_met in source_model.metabolites:
            original_id = cobra_met.id
            if re.search('(.+)_([a-z])\d+$', cobra_met.id):
                m = re.search('(.+)_([a-z])\d+$', cobra_met.id)
                if m[2] == "e":
                    cobra_met.compartment = "e0"
                    cobra_met.id = m[1] + "_e0"
                else:
                    cobra_met.compartment = m[2] + index
                    cobra_met.id = m[1] + "_" + m[2] + index
                if cobra_met.id not in self.model.metabolites and cobra_met.id not in self.new_metabolites:
                    self.new_metabolites[cobra_met.id] = cobra_met
                    local_remap[original_id] = cobra_met
                    if m[1] + "_" + m[2] in self.parameters["auto_sink"]:
                        new_demand.append(cobra_met)
                    if m[2] == "e":
                        new_exchange.append(cobra_met)
        # Adding all metabolites to model prior to adding reactions
        self.model.add_metabolites(self.new_metabolites.values())
        # Adding reactions from source model to gapfill model
        for modelreaction in source_model.reactions:
            if re.search('(.+)_([a-z])\d+$', modelreaction.id) != None:
                m = re.search('(.+)_([a-z])\d+$', modelreaction.id)
                if m[1] not in self.parameters["blacklist"]:
                    cobra_rxn = modelreaction.copy()
                    cobra_rxn.id = m[1] + "_" + m[2] + index    
                    if cobra_rxn.id not in self.model.reactions and cobra_rxn.id not in self.new_reactions:
                        self.new_reactions[cobra_rxn.id] = cobra_rxn
                        new_penalties[cobra_rxn.id] = {}
                        new_penalties[cobra_rxn.id]["added"] = 1
                        if cobra_rxn.lower_bound < 0:
                            new_penalties[cobra_rxn.id]["reverse"] = self.parameters["model_penalty"]
                        if cobra_rxn.upper_bound > 0:
                            new_penalties[cobra_rxn.id]["forward"] = self.parameters["model_penalty"]
                        # Updating metabolites in reaction to new model
                        new_stoichiometry = {}
                        for met in cobra_rxn.metabolites:
                            # Adding new coefficient:
                            new_stoichiometry[local_remap[met.id]] = cobra_rxn.metabolites[met]
                            # Zeroing out current coefficients
                            if local_remap[met.id] != met:
                                new_stoichiometry[met] = 0
                        cobra_rxn.add_metabolites(new_stoichiometry, combine=False)
                    elif cobra_rxn.lower_bound < 0 and self.model.reactions.get_by_id(cobra_rxn.id).lower_bound == 0:
                        self.model.reactions.get_by_id(cobra_rxn.id).lower_bound = cobra_rxn.lower_bound
                        self.model.reactions.get_by_id(cobra_rxn.id).update_variable_bounds()
                        new_penalties[cobra_rxn.id]["reverse"] = self.parameters["model_penalty"]
                        new_penalties[cobra_rxn.id]["reversed"] = 1
                    elif cobra_rxn.upper_bound > 0 and self.model.reactions.get_by_id(cobra_rxn.id).upper_bound == 0:
                        self.model.reactions.get_by_id(cobra_rxn.id).upper_bound = cobra_rxn.upper_bound
                        self.model.reactions.get_by_id(cobra_rxn.id).update_variable_bounds()
                        new_penalties[cobra_rxn.id]["forward"] = self.parameters["model_penalty"]
                        new_penalties[cobra_rxn.id]["reversed"] = 1

        # Only run this on new exchanges so we don't readd for all exchanges
        for cpd in new_exchange:
            drain_reaction = FBAHelper.add_drain_from_metabolite_id(
                self.model, cpd.id, self.parameters["default_uptake"], self.parameters["default_excretion"])
            if drain_reaction.id not in self.cobramodel.reactions and drain_reaction.id not in self.new_reactions:
                self.new_reactions[drain_reaction.id] = drain_reaction

        # Only run this on new demands so we don't read for all exchanges
        for cpd_id in new_demand:
            drain_reaction = FBAHelper.add_drain_from_metabolite_id(
                self.model, cpd.id, self.parameters["default_uptake"], self.parameters["default_excretion"])
            if drain_reaction.id not in self.cobramodel.reactions and drain_reaction.id not in self.new_reactions:
                self.new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.model.add_reactions(self.new_reactions.values())
        return new_penalties

    def _extend_model_with_template_metabolites(self, template, index='0'):
        self.new_metabolites = {}
        new_exchange, new_demand = [], []
        for template_compound in template.compcompounds:
            compartment_index = "0" if template_compound.compartment == 'e' else index
            cobra_met = self._convert_template_compound(template_compound, compartment_index, template)  # TODO: move function out
            if cobra_met.id not in self.model.metabolites and cobra_met.id not in self.new_metabolites:
                self.new_metabolites[cobra_met.id] = cobra_met
                #self.model.add_metabolites([cobra_met])
                msid = FBAHelper.modelseed_id_from_cobra_metabolite(cobra_met)
                if msid in self.parameters["auto_sink"]:
                    if msid != "cpd11416" or cobra_met.compartment == "c0":
                        new_demand.append(cobra_met.id)
                if template_compound.compartment == "e":
                    new_exchange.append(cobra_met.id)
        # Adding all metabolites to model prior to adding reactions
        self.model.add_metabolites(self.new_metabolites.values())

        return new_exchange, new_demand

    # Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def extend_model_with_template_for_gapfilling(self, template, index):
        # Adding all metabolites to model prior to adding reactions
        self.new_reactions, new_penalties = {}, {}
        new_exchange, new_demand = self._extend_model_with_template_metabolites(template, index)

        for template_reaction in template.reactions:
            if template_reaction.reference_id in self.parameters["blacklist"]:
                continue
            cobra_rxn = self._convert_template_reaction(template_reaction, index, template, 1)  # TODO: move function out
            new_penalties[cobra_rxn.id] = dict()
            if cobra_rxn.id not in self.model.reactions and cobra_rxn.id not in self.new_reactions:
                # Adding any template reactions missing from the present model
                self.new_reactions[cobra_rxn.id] = cobra_rxn
                if cobra_rxn.lower_bound < 0:
                    new_penalties[cobra_rxn.id]["reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                if cobra_rxn.upper_bound > 0:
                    new_penalties[cobra_rxn.id]["forward"] = template_reaction.base_cost + template_reaction.forward_penalty
                new_penalties[cobra_rxn.id]["added"] = 1
            elif template_reaction.GapfillDirection == "=":
                # Adjusting directionality as needed for existing reactions
                model_reaction = self.model.reactions.get_by_id(cobra_rxn.id)
                new_penalties[cobra_rxn.id]["reversed"] = 1
                if model_reaction.lower_bound == 0:
                    model_reaction.lower_bound = template_reaction.lower_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_rxn.id]["reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                if model_reaction.upper_bound == 0:
                    model_reaction.upper_bound = template_reaction.upper_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_rxn.id]["forward"] = template_reaction.base_cost + template_reaction.forward_penalty
        # Only run this on new exchanges so we don't read for all exchanges
        for cpd_id in new_exchange:
            drain_reaction = FBAHelper.add_drain_from_metabolite_id(
                self.model, cpd_id, self.parameters["default_uptake"], self.parameters["default_excretion"])
            if drain_reaction is not None and drain_reaction.id not in self.new_reactions:
                new_penalties[drain_reaction.id] = {
                    'added': 1,
                    'reverse': 1,
                    'forward': 1
                }
                self.new_reactions[drain_reaction.id] = drain_reaction
        # Only run this on new demands so we don't read for all exchanges
        for cpd_id in new_demand:
            drain_reaction = FBAHelper.add_drain_from_metabolite_id(
                self.model, cpd_id, self.parameters["default_uptake"], self.parameters["default_excretion"], "DM_")
            if drain_reaction is not None and drain_reaction.id not in self.new_reactions:
                new_penalties[drain_reaction.id] = {
                    'added': 1,
                    'reverse': 1,
                    'forward': 1
                }
                self.new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.model.add_reactions(self.new_reactions.values())
        return new_penalties

    def _convert_template_compound(self, template_compound, index, template):
        base_id = template_compound.id.split("_")[0]
        base_compound = template.compounds.get_by_id(base_id)
        new_id = template_compound.id
        new_id += str(index)
        compartment = template_compound.compartment
        compartment += str(index)

        met = Metabolite(new_id, formula=base_compound.formula, name=base_compound.name,
                         charge=template_compound.charge, compartment=compartment)
        met.annotation["sbo"] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        met.annotation["seed.compound"] = base_id
        return met

    def _convert_template_reaction(self, template_reaction, index, template, for_gapfilling=1):
        new_id = template_reaction.id+str(index)
        lower_bound = template_reaction.lower_bound
        upper_bound = template_reaction.upper_bound
        direction = template_reaction.GapfillDirection
        if for_gapfilling == 0:
            direction = template_reaction.direction
        if direction == ">":
            lower_bound = 0
        elif direction == "<":
            upper_bound = 0

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

        cobra_rxn = Reaction(new_id, name=template_reaction.name, lower_bound=lower_bound, upper_bound=upper_bound)
        cobra_rxn.add_metabolites(object_stoichiometry)
        cobra_rxn.annotation["sbo"] = "SBO:0000176"  # biochemical reaction
        cobra_rxn.annotation["seed.reaction"] = template_reaction.reference_id
        return cobra_rxn
    
    def binary_check_gapfilling_solution(self, solution=None, flux_values=None):
        if solution is None:
            solution = self._compute_gapfilled_solution(flux_values)
        rxn_filter = {rxn_id:solution["reversed"][rxn_id] for rxn_id in solution["reversed"]}
        rxn_filter.update({rxn_id:solution["new"][rxn_id] for rxn_id in solution["new"]})
        self.pkgmgr.getpkg("ReactionUsePkg").build_package(rxn_filter)
        objcoef = {}
        for rxnid in rxn_filter:
            if rxn_filter[rxnid] == ">":
                objcoef[self.pkgmgr.getpkg("ReactionUsePkg").variables["fu"][rxnid]] = 1
            if rxn_filter[rxnid] == "<":
                objcoef[self.pkgmgr.getpkg("ReactionUsePkg").variables["ru"][rxnid]] = 1
        new_solution = {}
        with self.model: # to prevent the model for permanently assuming the zeroed reactions
            # Setting all gapfilled reactions not in the solution to zero
            self._knockout_gf_reactions_outside_solution(solution,flux_values)
            # Setting the objective to the minimum sum of binary variables
            self.model.objective = self.model.problem.Objective(Zero, direction="min")
            self.model.objective.set_linear_coefficients(objcoef)
            self.model.optimize()
            new_solution = self._compute_gapfilled_solution(flux_values)
        return new_solution
    
    #This function is designed to KO all gapfilled reactions not included in the solution
    def _knockout_gf_reactions_outside_solution(self,solution = None,flux_values = None):
        if solution == None:
            solution = self._compute_gapfilled_solution(flux_values)
        if flux_values == None:
            flux_values = FBAHelper.compute_flux_values_from_variables(self.model)
        for rxnobj in self.model.reactions:
            if rxnobj.id in self.gapfilling_penalties:
                if "reverse" in self.gapfilling_penalties[rxnobj.id] and flux_values[rxnobj.id]["reverse"] <= Zero:
                    rxnobj.lower_bound = 0
                if "forward" in self.gapfilling_penalties[rxnobj.id] and flux_values[rxnobj.id]["forward"] <= Zero:
                    rxnobj.upper_bound = 0
                rxnobj.update_variable_bounds()
        
    def run_test_conditions(self, condition_list, solution = None, max_iterations = 10):
        reaction_list, filtered_list = [], []
        if solution == None:
            solution = self._compute_gapfilled_solution(flux_values)
        for rxnid in solution["reversed"]:
            reaction_list.append([self.model.reactions.get_by_id(rxnid),solution["reversed"][rxnid]])
        for rxnid in solution["new"]:
            reaction_list.append([self.model.reactions.get_by_id(rxnid),solution["new"][rxnid]])
        with self.model: # to prevent the model for permanently assuming the zeroed reactions
            #Setting all gapfilled reactions not in the solution to zero
            self._knockout_gf_reactions_outside_solution(solution)
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            filtered_list = FBAHelper.reaction_expansion_test(self.model,reaction_list,condition_list,self.pkgmgr)
        if len(filtered_list) > 0:
            if max_iterations > 0:
                print("Gapfilling test failed "+str(11-max_iterations))
                #Forcing filtered reactions to zero
                for item in filtered_list:
                    if item[1] == ">":
                        self.model.reactions.get_by_id(item[0].id).upper_bound = 0
                    else:
                        self.model.reactions.get_by_id(item[0].id).lower_bound = 0
                #Restoring lower bound on biomass constraint
                self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = self.parameters["minimum_obj"]
                #Reoptimizing
                self.model.optimize()
                return self.run_test_conditions(condition_list,None,max_iterations-1)
            return None
        return solution
    
    def filter_database_based_on_tests(self,test_conditions):
        with self.model:
            rxnlist = []
            for reaction in self.model.reactions:
                if reaction.id in self.gapfilling_penalties:
                    if "reverse" in self.gapfilling_penalties[reaction.id]:
                        reaction.lower_bound = 0
                        rxnlist.append([reaction,"<"])
                    if "forward" in self.gapfilling_penalties[reaction.id]:
                        reaction.upper_bound = 0
                        rxnlist.append([reaction,">"])
                    reaction.update_variable_bounds()
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            filtered_list = FBAHelper.reaction_expansion_test(self.model,rxnlist,test_conditions,self.pkgmgr)
        #Now constraining filtered reactions to zero
        for item in filtered_list:
            print("Filtering:",item[0].id,item[1])
            if item[1] == ">":
                self.model.reactions.get_by_id(item[0].id).upper_bound = 0
            else:
                self.model.reactions.get_by_id(item[0].id).lower_bound = 0
    
    def _compute_gapfilled_solution(self, flux_values=None):
        if flux_values is None:
            flux_values = FBAHelper.compute_flux_values_from_variables(self.model)
        output = {"reversed": {}, "new": {}}
        for reaction in self.model.reactions:
            if reaction.id in self.gapfilling_penalties:
                if flux_values[reaction.id]["forward"] > Zero and "forward" in self.gapfilling_penalties[reaction.id]:
                    if "added" in self.gapfilling_penalties[reaction.id]:
                        output["new"][reaction.id] = ">"
                    else:
                        output["reversed"][reaction.id] = ">"
                elif flux_values[reaction.id]["reverse"] > Zero and "reverse" in self.gapfilling_penalties[reaction.id]:
                    if "added" in self.gapfilling_penalties[reaction.id]:
                        output["new"][reaction.id] = "<"
                    else:
                        output["reversed"][reaction.id] = "<"
        return output