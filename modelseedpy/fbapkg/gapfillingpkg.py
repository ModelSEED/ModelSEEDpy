# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
import json
from optlang.symbolics import Zero, add
from cobra import Model, Reaction, Metabolite
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

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
    """

    """

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

    def get_model_index_hash(self):
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
        self.pkgmgr.getpkg("ObjConstPkg").build_package(self.parameters["minimum_obj"],None)
                
        # Determine all indecies that should be gapfilled
        indexhash = self.get_model_index_hash()

        # Iterating over all indecies with more than 10 intracellular compounds:
        self.gapfilling_penalties = dict()
        for index in indexhash:
            if indexhash[index] > 10:
                if index == "none":
                    for template in self.parameters["default_gapfill_templates"]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(template, index)
                        self.gapfilling_penalties.update(new_penalties)
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(gfmdl, index)
                        self.gapfilling_penalties.update(new_penalties)
                if index in self.parameters["gapfill_templates_by_index"]:
                    for template in self.parameters["gapfill_templates_by_index"][index]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(template, index)
                        self.gapfilling_penalties.update(new_penalties)
                if index in self.parameters["gapfill_models_by_index"]:
                    for gfmdl in self.parameters["gapfill_models_by_index"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(gfmdl, index)
                        self.gapfilling_penalties.update(new_penalties)
                if self.parameters["gapfill_all_indecies_with_default_templates"]:
                    for template in self.parameters["default_gapfill_templates"]:
                        new_penalties = self.extend_model_with_template_for_gapfilling(template, index)
                        self.gapfilling_penalties.update(new_penalties)
                if self.parameters["gapfill_all_indecies_with_default_models"]:
                    for gfmdl in self.parameters["default_gapfill_models"]:
                        new_penalties = self.extend_model_with_model_for_gapfilling(gfmdl, index)
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
                    self.gapfilling_penalties[reaction]["reverse"] = factor * self.gapfilling_penalties[reaction]["reverse"]
                if "forward" in self.gapfilling_penalties[reaction]:
                    self.gapfilling_penalties[reaction]["forward"] = factor * self.gapfilling_penalties[reaction]["forward"]

        self.model.solver.update()
        if self.parameters["set_objective"] == 1:
            reaction_objective = self.model.problem.Objective(
                Zero,
                direction="min")
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
            if re.search('(.+)_([a-z])\d+$', cobra_metabolite.id) != None:
                m = re.search('(.+)_([a-z])\d+$', cobra_metabolite.id)
                if m[2] == "e":
                    cobra_metabolite.compartment = "e0"
                    cobra_metabolite.id = m[1] + "_e0"
                else:
                    cobra_metabolite.compartment = m[2] + index
                    cobra_metabolite.id = m[1] + "_" + m[2] + index
                if cobra_metabolite.id not in self.model.metabolites and cobra_metabolite.id not in new_metabolites:
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
            if re.search('(.+)_([a-z])\d+$', modelreaction.id) != None:
                m = re.search('(.+)_([a-z])\d+$', modelreaction.id)
                if m[1] not in self.parameters["blacklist"]:
                    cobra_reaction = modelreaction.copy()
                    cobra_reaction.id = groups[1] + "_" + groups[2] + index
                    if cobra_reaction.id not in self.model.reactions and cobra_reaction.id not in new_reactions:
                        new_reactions[cobra_reaction.id] = cobra_reaction
                        new_penalties[cobra_reaction.id] = dict();
                        new_penalties[cobra_reaction.id]["added"] = 1
                        if cobra_reaction.lower_bound < 0:
                            new_penalties[cobra_reaction.id]["reverse"] = self.parameters["model_penalty"]
                        if cobra_reaction.upper_bound > 0:
                            new_penalties[cobra_reaction.id]["forward"] = self.parameters["model_penalty"]
                        # Updating metabolites in reaction to new model
                        metabolites = cobra_reaction.metabolites;
                        new_stoichiometry = {}
                        for metabolite in metabolites:
                            # Adding new coefficient:
                            new_stoichiometry[local_remap[metabolite.id]] = metabolites[metabolite]
                            # Zeroing out current coefficients
                            if local_remap[metabolite.id] != metabolite:
                                new_stoichiometry[metabolite] = 0
                        cobra_reaction.add_metabolites(new_stoichiometry, combine=False)
                    elif cobra_reaction.lower_bound < 0 and self.model.reactions.get_by_id(
                            cobra_reaction.id).lower_bound == 0:
                        self.model.reactions.get_by_id(cobra_reaction.id).lower_bound = cobra_reaction.lower_bound
                        self.model.reactions.get_by_id(cobra_reaction.id).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["reverse"] = self.parameters["model_penalty"]
                        new_penalties[cobra_reaction.id]["reversed"] = 1
                    elif cobra_reaction.upper_bound > 0 and self.model.reactions.get_by_id(
                            cobra_reaction.id).upper_bound == 0:
                        self.model.reactions.get_by_id(cobra_reaction.id).upper_bound = cobra_reaction.upper_bound
                        self.model.reactions.get_by_id(cobra_reaction.id).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["forward"] = model_penalty
                        new_penalties[cobra_reaction.id]["reversed"] = 1

                        # Only run this on new exchanges so we don't readd for all exchanges
        self.modelutl.add_exchanges_for_metabolites(new_exchange,self.parameters["default_uptake"],self.parameters["default_excretion"])
        # Only run this on new demands so we don't readd for all exchanges
        self.modelutl.add_exchanges_for_metabolites(new_demand,self.parameters["default_uptake"],self.parameters["default_excretion"],"DM_")
        # Adding all new reactions to the model at once (much faster than one at a time)
        self.model.add_reactions(new_reactions.values())
        return new_penalties

    def extend_model_with_template_metabolites(self, template, index='0'):
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
            compartment_index = "0" if compartment == 'e' else index
            cobra_metabolite = self.convert_template_compound(template_compound, compartment_index, template)  # TODO: move function out
            if cobra_metabolite.id not in self.model.metabolites and cobra_metabolite.id not in new_metabolites:
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                #self.model.add_metabolites([cobra_metabolite])
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
        new_reactions = {}
        new_penalties = dict()

        # Adding all metabolites to model prior to adding reactions
        new_exchange, new_demand = self.extend_model_with_template_metabolites(template, index)

        for template_reaction in template.reactions:
            if template_reaction.reference_id in self.parameters["blacklist"]:
                continue
            cobra_reaction = self.convert_template_reaction(template_reaction, index, template, 1)  # TODO: move function out
            new_penalties[cobra_reaction.id] = dict()
            if cobra_reaction.id not in self.model.reactions and cobra_reaction.id not in new_reactions:
                # Adding any template reactions missing from the present model
                new_reactions[cobra_reaction.id] = cobra_reaction
                if cobra_reaction.lower_bound < 0:
                    new_penalties[cobra_reaction.id][
                        "reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                if cobra_reaction.upper_bound > 0:
                    new_penalties[cobra_reaction.id][
                        "forward"] = template_reaction.base_cost + template_reaction.forward_penalty
                new_penalties[cobra_reaction.id]["added"] = 1
            elif template_reaction.GapfillDirection == "=":
                # Adjusting directionality as needed for existing reactions
                model_reaction = self.model.reactions.get_by_id(cobra_reaction.id)
                new_penalties[cobra_reaction.id]["reversed"] = 1
                if model_reaction.lower_bound == 0:
                    model_reaction.lower_bound = template_reaction.lower_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_reaction.id][
                        "reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                if model_reaction.upper_bound == 0:
                    model_reaction.upper_bound = template_reaction.upper_bound
                    model_reaction.update_variable_bounds()
                    new_penalties[cobra_reaction.id][
                        "forward"] = template_reaction.base_cost + template_reaction.forward_penalty
        # Only run this on new exchanges so we don't read for all exchanges
        exchanges = self.modelutl.add_exchanges_for_metabolites(new_exchange,self.parameters["default_uptake"],self.parameters["default_excretion"])
        for ex in exchanges:
            new_penalties[ex.id] = {
                'added': 1,
                'reverse': 1,
                'forward': 1
            }
            
        # Only run this on new demands so we don't readd for all exchanges
        exchanges = self.modelutl.add_exchanges_for_metabolites(new_demand,self.parameters["default_uptake"],self.parameters["default_excretion"],"DM_")
        for ex in exchanges:
            new_penalties[ex.id] = {
                'added': 1,
                'reverse': 1,
                'forward': 1
            }

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

        met = Metabolite(new_id,
                         formula=base_compound.formula,
                         name=base_compound.name,
                         charge=template_compound.charge,
                         compartment=compartment)

        met.annotation["sbo"] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        met.annotation["seed.compound"] = base_id
        return met

    def convert_template_reaction(self, template_reaction, index, template, for_gapfilling=1):
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

        cobra_reaction = Reaction(new_id,
                                  name=template_reaction.name,
                                  lower_bound=lower_bound,
                                  upper_bound=upper_bound)

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
            self.knockout_gf_reactions_outside_solution(solution,flux_values)
            # Setting the objective to be minimization of sum of binary variables
            min_reaction_objective = self.model.problem.Objective(Zero, direction="min")
            self.model.objective = min_reaction_objective
            min_reaction_objective.set_linear_coefficients(objcoef)
            self.model.optimize()
            new_solution = self.compute_gapfilled_solution()
        return new_solution
    
    #This function is designed to KO all gapfilled reactions not included in the solution
    def knockout_gf_reactions_outside_solution(self,solution = None,flux_values = None):
        if solution == None:
            solution = self.compute_gapfilled_solution()
        if flux_values == None:
            flux_values = self.modelutl.compute_flux_values_from_variables()
        for rxnobj in self.model.reactions:
            if rxnobj.id in self.gapfilling_penalties:
                if "reverse" in self.gapfilling_penalties[rxnobj.id] and flux_values[rxnobj.id]["reverse"] <= Zero:
                    rxnobj.lower_bound = 0
                if "forward" in self.gapfilling_penalties[rxnobj.id] and flux_values[rxnobj.id]["forward"] <= Zero:
                    rxnobj.upper_bound = 0
                rxnobj.update_variable_bounds()
        
    def run_test_conditions(self,condition_list,solution = None,max_iterations = 10):
        if solution == None:
            solution = self.compute_gapfilled_solution()
        reaction_list = []
        for rxnid in solution["reversed"]:
            reaction_list.append([self.model.reactions.get_by_id(rxnid),solution["reversed"][rxnid]])
        for rxnid in solution["new"]:
            reaction_list.append([self.model.reactions.get_by_id(rxnid),solution["new"][rxnid]])
        filtered_list = []
        with self.model:
            #Setting all gapfilled reactions not in the solution to zero
            self.knockout_gf_reactions_outside_solution(solution)
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            for condition in condition_list:
                condition["change"] = True
            filtered_list = self.modelutl.reaction_expansion_test(reaction_list,condition_list)
            for condition in condition_list:
                condition["change"] = False
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
        filetered_list = []
        with self.model:
            rxnlist = []
            for reaction in self.model.reactions:
                if reaction.id in self.gapfilling_penalties:
                    if "reverse" in self.gapfilling_penalties[reaction.id]:
                        rxnlist.append([reaction,"<"])
                    if "forward" in self.gapfilling_penalties[reaction.id]:
                        rxnlist.append([reaction,">"])
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            filtered_list = self.modelutl.reaction_expansion_test(rxnlist,test_conditions)
        #Now constraining filtered reactions to zero
        for item in filtered_list:
            logger.debug("Filtering:",item[0].id,item[1])
            if item[1] == ">":
                self.model.reactions.get_by_id(item[0].id).upper_bound = 0
            else:
                self.model.reactions.get_by_id(item[0].id).lower_bound = 0
        #Now testing if the gapfilling minimum objective can still be achieved
        gfobj = self.model.objective
        self.model.objective = self.parameters["origobj"]
        solution = self.model.optimize()
        #Restoring the minimum objective constraint
        self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = self.parameters["minimum_obj"]
        print("Objective after filtering:",solution.objective_value,"; min objective:",self.parameters["minimum_obj"])
        if solution.objective_value < self.parameters["minimum_obj"]:
            #Now we need to restore a minimal set of filtered reactions such that we permit the minimum objective to be reached
            new_objective = self.model.problem.Objective(
                Zero,
                direction="min")
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
            print("Reactions unfiltered:",count)
            #Checking for model reactions that can be removed to enable all tests to pass
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = 0
            filtered_list = self.modelutl.reaction_expansion_test(self.parameters["original_reactions"],test_conditions)
            for item in filtered_list:
                logger.debug("Filtering:",item[0].id,item[1])
                if item[1] == ">":
                    self.model.reactions.get_by_id(item[0].id).upper_bound = 0
                else:
                    self.model.reactions.get_by_id(item[0].id).lower_bound = 0
            self.pkgmgr.getpkg("ObjConstPkg").constraints["objc"]["1"].lb = self.parameters["minimum_obj"]
        self.model.objective = gfobj
    
    def compute_gapfilled_solution(self, flux_values=None):
        if flux_values is None:
            flux_values = self.modelutl.compute_flux_values_from_variables()
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
