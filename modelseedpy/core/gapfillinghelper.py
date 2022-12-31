# -*- coding: utf-8 -*-
import logging
logger = logging.getLogger(__name__)

import re
import copy

from cobra.core.dictlist import DictList
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from cobrakbase.core.kbaseobject import AttrDict
from cobrakbase.annotation_ontology_api.annotation_ontology_apiServiceClient import (
    annotation_ontology_api,
)
from numpy.f2py.cfuncs import f90modhooks

def build_cpd_id(string):
    if string.startswith("M_") or string.startswith("M-"):
        string = string[2:]
    if not string == string.replace('-', '__DASH__'):
        logger.debug('[Species] rename: [%s] -> [%s]', string, string.replace('-', '__DASH__'))
    return string

def build_rxn_id(string):
    if string.startswith("R_") or string.startswith("R-"):
        string = string[2:]
    if not string == string.replace('-', '__DASH__'):
        logger.debug('[Reaction] rename: [%s] -> [%s]', string, string.replace('-', '__DASH__'))
    return string


# Adding a few exception classes to handle different types of errors
class ObjectError(Exception):
    """Error in the construction of a base KBase object"""

    pass


class FeasibilityError(Exception):
    """Error in FBA formulation"""

    pass

# This class caries functions designed to more easily make standard modifications to an input model - no model state is retained in this class
class GapfillingHelper:
    def __init__(
        self, blacklist=[], auto_sink=["cpd02701_c", "cpd11416_c0", "cpd15302_c"]
    ):
        self.blacklist = ["rxn12985","rxn00238","rxn07058","rxn05305","rxn00154","rxn09037","rxn10643",
        "rxn11317","rxn05254","rxn05257","rxn05258","rxn05259","rxn05264","rxn05268",
        "rxn05269","rxn05270","rxn05271","rxn05272","rxn05273","rxn05274","rxn05275",
        "rxn05276","rxn05277","rxn05278","rxn05279","rxn05280","rxn05281","rxn05282",
        "rxn05283","rxn05284","rxn05285","rxn05286","rxn05963","rxn05964","rxn05971",
        "rxn05989","rxn05990","rxn06041","rxn06042","rxn06043","rxn06044","rxn06045",
        "rxn06046","rxn06079","rxn06080","rxn06081","rxn06086","rxn06087","rxn06088",
        "rxn06089","rxn06090","rxn06091","rxn06092","rxn06138","rxn06139","rxn06140",
        "rxn06141","rxn06145","rxn06217","rxn06218","rxn06219","rxn06220","rxn06221",
        "rxn06222","rxn06223","rxn06235","rxn06362","rxn06368","rxn06378","rxn06474",
        "rxn06475","rxn06502","rxn06562","rxn06569","rxn06604","rxn06702","rxn06706",
        "rxn06715","rxn06803","rxn06811","rxn06812","rxn06850","rxn06901","rxn06971",
        "rxn06999","rxn07123","rxn07172","rxn07254","rxn07255","rxn07269","rxn07451",
        "rxn09037","rxn10018","rxn10077","rxn10096","rxn10097","rxn10098","rxn10099",
        "rxn10101","rxn10102","rxn10103","rxn10104","rxn10105","rxn10106","rxn10107",
        "rxn10109","rxn10111","rxn10403","rxn10410","rxn10416","rxn11313","rxn11316",
        "rxn11318","rxn11353","rxn05224","rxn05795","rxn05796","rxn05797","rxn05798",
        "rxn05799","rxn05801","rxn05802","rxn05803","rxn05804","rxn05805","rxn05806",
        "rxn05808","rxn05812","rxn05815","rxn05832","rxn05836","rxn05851","rxn05857",
        "rxn05869","rxn05870","rxn05884","rxn05888","rxn05896","rxn05898","rxn05900",
        "rxn05903","rxn05904","rxn05905","rxn05911","rxn05921","rxn05925","rxn05936",
        "rxn05947","rxn05956","rxn05959","rxn05960","rxn05980","rxn05991","rxn05992",
        "rxn05999","rxn06001","rxn06014","rxn06017","rxn06021","rxn06026","rxn06027",
        "rxn06034","rxn06048","rxn06052","rxn06053","rxn06054","rxn06057","rxn06059",
        "rxn06061","rxn06102","rxn06103","rxn06127","rxn06128","rxn06129","rxn06130",
        "rxn06131","rxn06132","rxn06137","rxn06146","rxn06161","rxn06167","rxn06172",
        "rxn06174","rxn06175","rxn06187","rxn06189","rxn06203","rxn06204","rxn06246",
        "rxn06261","rxn06265","rxn06266","rxn06286","rxn06291","rxn06294","rxn06310",
        "rxn06320","rxn06327","rxn06334","rxn06337","rxn06339","rxn06342","rxn06343",
        "rxn06350","rxn06352","rxn06358","rxn06361","rxn06369","rxn06380","rxn06395",
        "rxn06415","rxn06419","rxn06420","rxn06421","rxn06423","rxn06450","rxn06457",
        "rxn06463","rxn06464","rxn06466","rxn06471","rxn06482","rxn06483","rxn06486",
        "rxn06492","rxn06497","rxn06498","rxn06501","rxn06505","rxn06506","rxn06521",
        "rxn06534","rxn06580","rxn06585","rxn06593","rxn06609","rxn06613","rxn06654",
        "rxn06667","rxn06676","rxn06693","rxn06730","rxn06746","rxn06762","rxn06779",
        "rxn06790","rxn06791","rxn06792","rxn06793","rxn06794","rxn06795","rxn06796",
        "rxn06797","rxn06821","rxn06826","rxn06827","rxn06829","rxn06839","rxn06841",
        "rxn06842","rxn06851","rxn06866","rxn06867","rxn06873","rxn06885","rxn06891",
        "rxn06892","rxn06896","rxn06938","rxn06939","rxn06944","rxn06951","rxn06952",
        "rxn06955","rxn06957","rxn06960","rxn06964","rxn06965","rxn07086","rxn07097",
        "rxn07103","rxn07104","rxn07105","rxn07106","rxn07107","rxn07109","rxn07119",
        "rxn07179","rxn07186","rxn07187","rxn07188","rxn07195","rxn07196","rxn07197",
        "rxn07198","rxn07201","rxn07205","rxn07206","rxn07210","rxn07244","rxn07245",
        "rxn07253","rxn07275","rxn07299","rxn07302","rxn07651","rxn07723","rxn07736",
        "rxn07878","rxn11417","rxn11582","rxn11593","rxn11597","rxn11615","rxn11617",
        "rxn11619","rxn11620","rxn11624","rxn11626","rxn11638","rxn11648","rxn11651",
        "rxn11665","rxn11666","rxn11667","rxn11698","rxn11983","rxn11986","rxn11994",
        "rxn12006","rxn12007","rxn12014","rxn12017","rxn12022","rxn12160","rxn12161",
        "rxn01267","rxn05294","rxn04656"]
        for item in blacklist:
            if item not in self.blacklist:
                self.blacklist.append(item)
        self.auto_sink = []
        full_id = re.compile("\d+$")
        for id in auto_sink:
            if full_id.search(id):
                self.auto_sink.append(id)
            else:
                for i in range(0, 100):
                    newid = id + str(i)
                    self.auto_sink.append(newid)

        self.auto_exchange = "e0"
        self.COBRA_0_BOUND = 0
        self.COBRA_DEFAULT_LB = -1000
        self.COBRA_DEFAULT_UB = 1000

    # FBA macro analyses
    def test_reaction_additions_againt_limits(self, model, reactions, tests):
        filtered = DictList()
        with model:
            for rxn in reactions:
                if rxn.id() in model.reactions:
                    rxn_obj = model.reactions.get_by_id(rxn.id())
                else:
                    rxn_obj = model.add_reactions([rxn])
                self.set_reaction_bounds_from_direction(rxn_obj,reactions[rxn]) 
                for test in tests:
                    testmodel = model
                    with testmodel:
                        self.apply_media_to_model(testmodel,test["media"],test["default_uptake"],test["default_excretion"])  #!!! where is this function defined?
                        self.set_objective_from_target_reaction(testmodel,test["target"],test["maximize"])  #!!! where is this function defined?
                        solution = testmodel.optimize()
                        if test.maximize == 1:
                            if testmodel.objective.value() > test.limit:
                                filtered.append(test)
        return filtered

    def build_model_extended_for_gapfilling(self,extend_with_template = 1, source_models = [], input_templates = [], model_penalty = 1, reaction_scores = {}):
        #Determine all indecies that should be gapfilled
        indices = [0]*1000
        compounds = self.fbamodel["modelcompounds"]
        for compound in compounds:
            compartment = compound["modelcompartment_ref"].split("/").pop()
            basecomp = compartment[0:1]
            if not basecomp == "e":
                indices[int(compartment[1:])] += 1

        # Iterating over all indecies with more than 10 intracellular compounds:
        gapfilling_penalties = dict()
        for i, val in enumerate(indices):
            if val > 10:
                if extend_with_template == 1:
                    new_penalties = self.temp_extend_model_index_for_gapfilling(
                        i, input_templates
                    )
                    gapfilling_penalties.update(new_penalties)
                if i < len(source_models) and source_models[i] != None:
                    new_penalties = self.mdl_extend_model_index_for_gapfilling(
                        i, source_models[i], model_penalty
                    )
                    gapfilling_penalties.update(new_penalties)
        # Rescaling penalties by reaction scores and saving genes
        for reaction in gapfilling_penalties:
            rxnid = reaction.split("_")[0]
            if rxnid in reaction_scores:
                highest_score = 0
                for gene in reaction_scores[rxnid]:
                    if highest_score < reaction_scores[rxnid][gene]:
                        highest_score = reaction_scores[rxnid][gene]
                factor = 1 - 0.9 * highest_score
                if "reverse" in gapfilling_penalties[reaction]:
                    gapfilling_penalties[reaction.id]["reverse"] = factor*gapfilling_penalties[reaction.id]["reverse"]
                if "forward" in gapfilling_penalties[reaction]:
                    gapfilling_penalties[reaction.id]["forward"] = factor*gapfilling_penalties[reaction.id]["forward"]
        self.cobramodel.solver.update()
        return gapfilling_penalties

    #Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def mdl_extend_model_index_for_gapfilling(self, model, index, source_model, model_penalty):
        new_metabolites, new_reactions, new_penalties, local_remap = {}, {}, {}, {}
        new_exchange, new_demand = [], []
        comp = re.compile('(.*_*)(.)\d+$')
        for modelcompound in source_model.metabolites:
            cobra_metabolite = self.convert_modelcompound(modelcompound)
            original_id = cobra_metabolite.id
            groups = comp.match(cobra_metabolite.compartment)
            if groups[2] == "e":
                cobra_metabolite.compartment = groups[1] + groups[2] + "0"
                groups = comp.match(cobra_metabolite.id)
                cobra_metabolite.id = groups[1] + groups[2] + "0"
            else:
                cobra_metabolite.compartment = groups[1] + groups[2] + str(index)
                groups = comp.match(cobra_metabolite.id)
                cobra_metabolite.id = groups[1] + groups[2] + str(index)
            if (
                cobra_metabolite.id not in self.cobramodel.metabolites
                and cobra_metabolite.id not in new_metabolites
            ):
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                if original_id in self.auto_sink:
                    self.demand_compounds.add(cobra_metabolite.id)  #!!! where is demand_compounds defined?
                    new_demand.append(cobra_metabolite)
                if cobra_metabolite.compartment == self.auto_exchange:
                    self.exchange_compounds.add(cobra_metabolite.id)
                    new_exchange.append(cobra_metabolite)
            if cobra_metabolite.id in self.cobramodel.metabolites: #!!! where is cobramodel defined?
                cobra_metabolite = self.cobramodel.metabolites.get_by_id(
                    cobra_metabolite.id
                )
            else:  # Just in case the same compound is added twice - we want to switch the metabolite to the first new version
                cobra_metabolite = new_metabolites[cobra_metabolite.id]
            local_remap[original_id] = cobra_metabolite
        # Adding all metabolites to model prior to adding reactions
        self.cobramodel.add_metabolites(new_metabolites.values())

        for modelreaction in source_model.reactions:
            if modelreaction.id.split("_")[0] in self.blacklist:
                continue
            # cobra_reaction = self.convert_modelreaction(modelreaction)
            cobra_reaction = modelreaction.copy()
            groups = comp.match(cobra_reaction.id)
            cobra_reaction.id = groups[1] + groups[2] + str(index)
            new_penalties[cobra_reaction.id] = dict()
            # Updating metabolites in reaction to new model
            metabolites = cobra_reaction.metabolites
            new_stoichiometry = {}
            for metabolite in metabolites:
                # Adding new coefficient:
                new_stoichiometry[local_remap[metabolite.id]] = metabolites[metabolite]
                # Zeroing out current coefficients
                if local_remap[metabolite.id] != metabolite:
                    new_stoichiometry[metabolite] = 0
            cobra_reaction.add_metabolites(new_stoichiometry, combine=False)
            if (
                cobra_reaction.id not in self.cobramodel.reactions
                and cobra_reaction.id not in new_reactions
            ):
                new_reactions[cobra_reaction.id] = cobra_reaction
                new_penalties[cobra_reaction.id]["added"] = 1
                if cobra_reaction.lower_bound < 0:
                    new_penalties[cobra_reaction.id]["reverse"] = model_penalty
                    new_penalties[cobra_reaction.id]["reversed"] = True
                elif cobra_reaction.upper_bound > 0 and rxn.upper_bound == 0:
                    rxn.upper_bound = cobra_reaction.upper_bound
                    rxn.update_variable_bounds()
                    new_penalties[cobra_reaction.id]["forward"] = model_penalty
            elif (
                cobra_reaction.lower_bound < 0
                and self.cobramodel.reactions.get_by_id(cobra_reaction.id).lower_bound
                == 0
            ):
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).lower_bound = cobra_reaction.lower_bound
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).update_variable_bounds()
                new_penalties[cobra_reaction.id]["reverse"] = model_penalty
                new_penalties[cobra_reaction.id]["reversed"] = 1
            elif (
                cobra_reaction.upper_bound > 0
                and self.cobramodel.reactions.get_by_id(cobra_reaction.id).upper_bound
                == 0
            ):
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).upper_bound = cobra_reaction.upper_bound
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).update_variable_bounds()
                new_penalties[cobra_reaction.id]["forward"] = model_penalty
                new_penalties[cobra_reaction.id]["reversed"] = 1

        # Only run this on new exchanges so we don't readd for all exchanges
        for cpd in new_exchange:
            drain_reaction = self.add_drain_from_metabolite_id(cpd.id)
            if (
                drain_reaction.id not in self.cobramodel.reactions
                and drain_reaction.id not in new_reactions
            ):
                new_reactions[drain_reaction.id] = drain_reaction

        # Only run this on new demands so we don't readd for all exchanges
        for cpd_id in new_demand:
            drain_reaction = self.add_drain_from_metabolite_id(
                cpd_id,
                lower_bound=self.COBRA_0_BOUND,
                upper_bound=self.COBRA_DEFAULT_UB,
                prefix="DM_",
                prefix_name="Demand for ",
                sbo="SBO:0000627",
            )
            if (
                drain_reaction.id not in self.cobramodel.reactions
                and drain_reaction.id not in new_reactions
            ):
                new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.cobramodel.add_reactions(new_reactions.values())
        return new_penalties

    #Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def temp_extend_model_index_for_gapfilling(self,index,input_templates = []):
        new_metabolites, new_reactions, new_penalties = {}, {}, {}
        new_exchange, new_demand = [], []
        template = None
        if index < len(input_templates):
            template = input_templates[index]
        elif index in self.fbamodel['template_refs']:  #!!! where is fbamodel defined?
            template = self.kbapi.get_from_ws(self.fbamodel['template_refs'][index])
        else:
            template = self.kbapi.get_from_ws(self.fbamodel['template_ref'])
        if template.info.type != "KBaseFBA.NewModelTemplate":
            raise ObjectError(
                template.info.type + " loaded when KBaseFBA.NewModelTemplate expected"
            )

        for template_compound in template.compcompounds:
            tempindex = index
            if template_compound.templatecompartment_ref.split("/").pop() == "e":
                tempindex = 0

            cobra_metabolite = self.convert_template_compound(
                template_compound, tempindex, template
            )
            if (
                cobra_metabolite.id not in self.cobramodel.metabolites
                and cobra_metabolite.id not in new_metabolites
            ):
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                self.cobramodel.add_metabolites([cobra_metabolite])
                if cobra_metabolite.id in self.auto_sink:
                    self.demand_compounds.add(cobra_metabolite.id)
                    new_demand.append(cobra_metabolite.id)
                if cobra_metabolite.compartment == self.auto_exchange:
                    new_exchange.append(cobra_metabolite.id)
                    self.exchange_compounds.add(cobra_metabolite.id)
        # Adding all metabolites to model prior to adding reactions
        self.cobramodel.add_metabolites(new_metabolites.values())

        for template_reaction in template.reactions:
            if template_reaction.id.split("_")[0] not in self.blacklist:      
                cobra_reaction = self.convert_template_reaction(template_reaction,index,template,1)
                new_penalties[cobra_reaction.id] = dict()
                if cobra_reaction.id not in (self.cobramodel.reactions and new_reactions):
                    #Adding any template reactions missing from the present model
                    new_reactions[cobra_reaction.id] = cobra_reaction
                    if cobra_reaction.lower_bound < 0:
                        new_penalties[cobra_reaction.id]["reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                    if cobra_reaction.upper_bound > 0:
                        new_penalties[cobra_reaction.id]["forward"] = template_reaction.base_cost + template_reaction.forward_penalty
                    new_penalties[cobra_reaction.id]["added"] = True
                elif template_reaction.GapfillDirection == "=":
                    #Adjusting directionality as needed for existing reactions
                    new_penalties[cobra_reaction.id]["reversed"] = True
                    if self.cobramodel.reactions.get_by_id(cobra_reaction.id).lower_bound == 0:
                        self.cobramodel.reactions.get_by_id(cobra_reaction.id).lower_bound = template_reaction.maxrevflux
                        self.cobramodel.reactions.get_by_id(cobra_reaction.id).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["reverse"] = template_reaction.base_cost + template_reaction.reverse_penalty
                    if self.cobramodel.reactions.get_by_id(cobra_reaction.id).upper_bound == 0:
                        self.cobramodel.reactions.get_by_id(cobra_reaction.id).upper_bound = template_reaction.maxforflux
                        self.cobramodel.reactions.get_by_id(cobra_reaction.id).update_variable_bounds()
                        new_penalties[cobra_reaction.id]["forward"] = template_reaction.base_cost + template_reaction.forward_penalty

        # Only run this on new exchanges so we don't readd for all exchanges
        for cpd_id in new_exchange:
            drain_reaction = self.add_drain_from_metabolite_id(cpd_id)
            if drain_reaction and drain_reaction.id not in new_reactions:
                new_reactions[drain_reaction.id] = drain_reaction

        # Only run this on new demands so we don't readd for all exchanges
        for cpd_id in new_demand:
            drain_reaction = self.add_drain_from_metabolite_id(
                cpd_id, self.COBRA_0_BOUND, self.COBRA_DEFAULT_UB, "DM_", "Demand for "
            )
            if drain_reaction is not None and drain_reaction.id not in new_reactions:
                new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.cobramodel.add_reactions(new_reactions.values())
        return new_penalties

    def convert_modelreaction(self, reaction, bigg=False):
        lower_bound, upper_bound = reaction.get_reaction_constraints()
        rxn_id = build_rxn_id(mr_id)
        if bigg and "bigg.reaction" in annotation:
            rxn_id = annotation["bigg.reaction"]

        cobra_reaction = Reaction(
            rxn_id, name=name, lower_bound=lower_bound, upper_bound=upper_bound
        )
        cobra_reaction.annotation[
            self.SBO_ANNOTATION
        ] = "SBO:0000176"  # biochemical reaction
        cobra_reaction.annotation.update(annotation)

        if rxn_id.startswith("rxn"):
            cobra_reaction.annotation["seed.reaction"] = rxn_id.split("_")[0]

        cobra_reaction.add_metabolites(
            self.convert_modelreaction_stoichiometry(reaction)
        )

        cobra_reaction.gene_reaction_rule = reaction.gene_reaction_rule

        for genes in reaction.get_gpr():
            for gene in genes:
                if not gene in self.genes:
                    self.genes[gene] = gene

        return cobra_reaction

    def convert_modelcompound(self, metabolite, bigg=False):
        formula = metabolite.formula
        name = metabolite.name
        charge = metabolite.charge
        mc_id = metabolite.id
        compartment = metabolite.compartment
        annotation = metabolite.annotation

        id = build_cpd_id(mc_id)

        if bigg and "bigg.metabolite" in annotation:
            id = annotation["bigg.metabolite"] + "_" + compartment
            # print(id)

        met = Metabolite(
            id, formula=formula, name=name, charge=charge, compartment=compartment
        )

        met.annotation[
            self.SBO_ANNOTATION
        ] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        if id.startswith("cpd"):
            met.annotation["seed.compound"] = id.split("_")[0]
        met.annotation.update(annotation)
        return met

    def convert_modelreaction_stoichiometry(self, reaction):
        object_stoichiometry = {}
        s = reaction.stoichiometry
        for metabolite_id in s:
            if metabolite_id in self.metabolites_remap:
                object_stoichiometry[
                    self.cobramodel.metabolites.get_by_id(
                        self.metabolites_remap[metabolite_id]
                    )
                ] = s[metabolite_id]
        return object_stoichiometry

    def create_binary_variables(self, rxnobj, forward=1, reverse=1):
        if rxnobj.id not in self.binary_flux_variables:
            self.binary_flux_variables[rxnobj.id] = dict()
            self.binary_flux_constraints[rxnobj.id] = dict()
        if (
            forward == 1
            and rxnobj.upper_bound > 0
            and "forward" not in self.binary_flux_variables[rxnobj.id]
        ):
            self.binary_flux_variables[rxnobj.id][
                "forward"
            ] = self.cobramodel.problem.Variable(
                rxnobj.id + "_fb", lb=0, ub=1, type="binary"
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_variables[rxnobj.id]["forward"]
            )
            self.binary_flux_constraints[rxnobj.id][
                "forward"
            ] = self.cobramodel.problem.Constraint(
                1000 * self.binary_flux_variables[rxnobj.id]["forward"]
                - rxnobj.forward_variable,
                lb=0,
                ub=None,
                name=rxnobj.id + "_fb",
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_constraints[rxnobj.id]["forward"]
            )
        if (
            reverse == 1
            and rxnobj.lower_bound < 0
            and "reverse" not in self.binary_flux_variables[rxnobj.id]
        ):
            self.binary_flux_variables[rxnobj.id][
                "reverse"
            ] = self.cobramodel.problem.Variable(
                rxnobj.id + "_bb", lb=0, ub=1, type="binary"
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_variables[rxnobj.id]["reverse"]
            )
            self.binary_flux_constraints[rxnobj.id][
                "reverse"
            ] = self.cobramodel.problem.Constraint(
                1000 * self.binary_flux_variables[rxnobj.id]["reverse"]
                - rxnobj.forward_variable,
                lb=0,
                ub=None,
                name=rxnobj.id + "_bb",
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_constraints[rxnobj.id]["reverse"]
            )

    def binary_check_gapfilling_solution(
        self, gapfilling_penalties, add_solution_exclusion_constraint
    ):
        objcoef = {}
        flux_values = self.compute_flux_values_from_variables()
        for rxnobj in self.cobramodel.reactions:
            if rxnobj.id in gapfilling_penalties:
                if (
                    "reverse" in gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["reverse"] > Zero
                ):
                    self.create_binary_variables(rxnobj, 0, 1)
                    objcoef[self.binary_flux_variables[rxnobj.id]["reverse"]] = 1
                if (
                    "forward" in gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["forward"] > Zero
                ):
                    self.create_binary_variables(rxnobj, 1, 0)
                    objcoef[self.binary_flux_variables[rxnobj.id]["forward"]] = 1
        with self.cobramodel:
            # Setting all gapfilled reactions not in the solution to zero
            min_reaction_objective = self.cobramodel.problem.Objective(
                Zero, direction="min"
            )
            for rxnobj in self.cobramodel.reactions:
                if rxnobj.id in gapfilling_penalties:
                    if (
                        "reverse" in gapfilling_penalties[rxnobj.id]
                        and flux_values[rxnobj.id]["reverse"] <= Zero
                    ):
                        rxnobj.lower_bound = 0
                    if (
                        "forward" in gapfilling_penalties[rxnobj.id]
                        and flux_values[rxnobj.id]["forward"] <= Zero
                    ):
                        rxnobj.upper_bound = 0
                    rxnobj.update_variable_bounds()
            # Setting the objective to be minimization of sum of binary variables
            self.cobramodel.objective = min_reaction_objective
            self.cobramodel.objective.set_linear_coefficients(objcoef)
            with open("GapfillBinary.lp", "w") as out:
                out.write(str(self.cobramodel.solver))
            self.cobramodel.optimize()
            flux_values = self.compute_flux_values_from_variables()
        if add_solution_exclusion_constraint:
            self.add_binary_solution_exclusion_constraint(flux_values)
        return flux_values

    # Adds a constraint that eliminates a gapfilled solution from feasibility so a new solution can be obtained
    def add_binary_solution_exclusion_constraint(self, flux_values):
        count = len(self.solution_exclusion_constraints)
        solution_coef = {}
        solution_size = 0
        for reaction, direction in self.binary_flux_variables.items():
            if flux_values[reaction][direction] > Zero:
                solution_size += 1
                solution_coef[self.binary_flux_variables[reaction][direction]] = 1
        if len(solution_coef) > 0:
            new_exclusion_constraint = self.cobramodel.problem.Constraint(
                Zero,
                lb=None,
                ub=(solution_size - 1),
                name="exclusion." + str(count + 1),
            )
            self.cobramodel.add_cons_vars(new_exclusion_constraint)
            self.cobramodel.solver.update()
            new_exclusion_constraint.set_linear_coefficients(solution_coef)
            self.solution_exclusion_constraints.append(new_exclusion_constraint)
            return new_exclusion_constraint
        return None

    # Takes gapfilled penalties and creates and objective function minimizing gapfilled reactions
    def create_minimal_reaction_objective(self, penalty_hash, default_penalty=0):
        reaction_objective = self.cobramodel.problem.Objective(Zero, direction="min")
        obj_coef = dict()
        for reaction in self.cobramodel.reactions:
            obj_coef[reaction.forward_variable] = default_penalty
            obj_coef[reaction.reverse_variable] = default_penalty
            if reaction.id in penalty_hash:
                # Minimizing gapfilled reactions
                if "reverse" in penalty_hash[reaction.id]:
                    obj_coef[reaction.reverse_variable] = abs(
                        penalty_hash[reaction.id]["reverse"]
                    )
                elif default_penalty != 0:
                    obj_coef[reaction.reverse_variable] = default_penalty
                if "forward" in penalty_hash[reaction.id]:
                    obj_coef[reaction.forward_variable] = abs(
                        penalty_hash[reaction.id]["forward"]
                    )
                elif default_penalty != 0:
                    obj_coef[reaction.forward_variable] = default_penalty
            else:
                obj_coef[reaction.forward_variable] = default_penalty
                obj_coef[reaction.reverse_variable] = default_penalty

        self.cobramodel.objective = reaction_objective
        self.cobramodel.objective.set_linear_coefficients(obj_coef)

    # Required this function to add gapfilled compounds to a KBase model for saving gapfilled model
    def convert_cobra_compound_to_kbcompound(self, cpd, kbmodel, add_to_model=1):
        refid = "cpd00000"
        if re.search("cpd\d+_[a-z]+", cpd.id):
            refid = cpd.id
            refid = re.sub("_[a-z]\d+$", "", refid)
        cpd_data = {
            "aliases": [],
            "charge": cpd.charge,
            "compound_ref": "~/template/compounds/id/" + refid,
            "dblinks": {},
            "formula": cpd.formula,
            "id": cpd.id,
            "inchikey": "ALYNCZNDIQEVRV-UHFFFAOYSA-M",
            "modelcompartment_ref": "~/modelcompartments/id/" + cpd.id.split("_").pop(),
            "name": cpd.name(),
            "numerical_attributes": {},
            "string_attributes": {},
        }
        cpd_data = AttrDict(cpd_data)
        if kbmodel:
            kbmodel.modelcompounds.append(cpd_data)
        return cpd_data

    # Required this function to add gapfilled reactions to a KBase model for saving gapfilled model
    def convert_cobra_reaction_to_kbreaction(
        self, rxn, kbmodel, direction="=", add_to_model=1
    ):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search("rxn\d+_[a-z]+", rxn.id):
            rxnref = "~/template/reactions/id/" + rxn.id
            rxnref = re.sub("\d+$", "", rxnref)
        rxn_data = {
            "id": rxn.id,
            "aliases": [],
            "dblinks": {},
            "direction": direction,
            "edits": {},
            "gapfill_data": {},
            "maxforflux": 1000000,
            "maxrevflux": 1000000,
            "modelReactionProteins": [],
            "modelReactionReagents": [],
            "modelcompartment_ref": "~/modelcompartments/id/" + rxn.id.split("_").pop(),
            "name": rxn.name,
            "numerical_attributes": {},
            "probability": 0,
            "protons": 0,
            "reaction_ref": rxnref,
            "string_attributes": {},
        }
        rxn_data = AttrDict(rxn_data)
        for cpd in rxn.metabolites:
            if cpd.id not in kbmodel.modelcompounds:
                convert_cobra_compound_to_kbcompound(cpd, kbmodel, 1)
            rxn_data.modelReactionReagents.append(
                {
                    "coefficient": rxn.metabolites[cpd],
                    "modelcompound_ref": "~/modelcompounds/id/" + cpd.id,
                }
            )
        if add_to_model == 1:
            kbmodel.modelreactions.append(rxn_data)
        return rxn_data

    def convert_objective_to_constraint(self, lower_bound, upper_bound):
        old_obj_variable = self.cobramodel.problem.Variable(
            name="old_objective_variable", lb=lower_bound, ub=upper_bound
        )
        old_obj_constraint = self.cobramodel.problem.Constraint(
            self.cobramodel.solver.objective.expression - old_obj_variable,
            lb=0, ub=0, name="old_objective_constraint")
        self.cobramodel.add_cons_vars([old_obj_variable, old_obj_constraint])

    def compute_flux_values_from_variables(self):
        flux_values = {}
        for rxnobj in self.cobramodel.reactions:
            flux_values[rxnobj.id] = {}
            flux_values[rxnobj.id]["reverse"] = rxnobj.reverse_variable.primal
            flux_values[rxnobj.id]["forward"] = rxnobj.forward_variable.primal
        return flux_values

    def compute_gapfilled_solution(self, penalties, flux_values=None):
        if flux_values is None:
            flux_values = self.compute_flux_values_from_variables()
        output = {"reversed": {}, "new": {}}
        for reaction in self.cobramodel.reactions:
            if reaction.id in penalties:
                if (
                    flux_values[reaction.id]["forward"] > Zero
                    and "forward" in penalties[reaction.id]
                ):
                    if "added" in penalties[reaction.id]:
                        output["new"][reaction.id] = ">"
                    else:
                        output["reversed"][reaction.id] = ">"
                elif (
                    flux_values[reaction.id]["reverse"] > Zero
                    and "reverse" in penalties[reaction.id]
                ):
                    if "added" in penalties[reaction.id]:
                        output["new"][reaction.id] = "<"
                    else:
                        output["reversed"][reaction.id] = "<"
        return output

    def add_gapfilling_solution_to_kbase_model(self, newmodel, penalties, media_ref):
        largest_index = max([gapfilling.id.split(".").pop() for gapfilling in newmodel.gapfillings])
        gfid = "gf." + str(largest_index+1)
        newmodel.gapfillings.append(
            {
                "gapfill_id": newmodel.id + "." + gfid,
                "id": gfid,
                "integrated": 1,
                "integrated_solution": "0",
                "media_ref": media_ref,
            }
        )
        for reaction in self.cobramodel.reactions:
            if reaction.id in penalties:
                if (
                    reaction.forward_variable.primal > Zero
                    and "forward" in penalties[reaction.id]
                ):
                    if reaction.id not in newmodel.modelreactions:
                        self.convert_cobra_reaction_to_kbreaction(
                            reaction, newmodel, ">", 1
                        )
                    gfrxn = newmodel.modelreactions.get_by_id(reaction.id)
                    gfrxn.gapfill_data[gfid] = dict()
                    gfrxn.gapfill_data[gfid]["0"] = [">", 1, []]
                elif (
                    reaction.forward_variable.primal > Zero
                    and "reverse" in penalties[reaction.id]
                ):
                    if reaction.id not in newmodel.modelreactions:
                        self.convert_cobra_reaction_to_kbreaction(
                            reaction, newmodel, "<", 1
                        )
                    gfrxn = newmodel.modelreactions.get_by_id(reaction.id)
                    gfrxn.gapfill_data[gfid] = dict()
                    gfrxn.gapfill_data[gfid]["0"] = ["<", 1, []]

    def compute_reaction_scores(self, weigh_all_events_equally=1, weights=None):
        reaction_genes = {}
        if "genome_ref" in self.fbamodel:
            anno_api = annotation_ontology_api()
            events = anno_api.get_annotation_ontology_events(
                {
                    "input_ref": self.fbamodel["genome_ref"],
                }
            )
            for event in events:
                for gene in event["ontology_terms"]:
                    if "modelseed_ids" in event["ontology_terms"][gene]:
                        for rxn in event["ontology_terms"][gene]["modelseed_ids"]:
                            newrxn = re.sub("^MSRXN:", "", rxn)
                            if newrxn not in reaction_genes:
                                reaction_genes[newrxn] = {}
                            if gene not in reaction_genes[newrxn]:
                                reaction_genes[newrxn][gene] = 0
                            if weigh_all_events_equally == 1 or weights is None:
                                reaction_genes[newrxn][gene] += 1
                            elif event["description"] in weights:
                                reaction_genes[newrxn][gene] += weights[
                                    event["description"]
                                ]
                            elif event["event_id"] in weights:
                                reaction_genes[newrxn][gene] += weights[
                                    event["event_id"]
                                ]
                            elif event["id"] in weights:
                                reaction_genes[newrxn][gene] += weights[event["id"]]
        return reaction_genes

    def replicate_model(self, count):
        newmodel = Model(self.cobramodel.id + "_rep" + str(count))
        utilities = KBaseFBAUtilities(
            newmodel,
            newmodel,
            self.kbapi,
            self.media,
            default_uptake=self.default_uptake,
            default_excretion=self.default_excretion,
            blacklist=self.blacklist,
        )
        metabolites = []
        reactions = []
        metabolite_hash = {}
        for i in range(0, count):
            for metabolite in self.cobramodel.metabolites:
                metabolite = metabolite.copy()
                metabolite.id = metabolite.id + "__" + str(i)
                metabolite_hash[metabolite.id] = metabolite
                metabolites.append(metabolite)
            for reaction in self.cobramodel.reactions:
                reaction = reaction.copy()
                reaction.id = reaction.id + "__" + str(i)
                input_metabolites = {}
                for metabolite in reaction.metabolites:
                    newid = metabolite.id + "__" + str(i)
                    input_metabolites[metabolite_hash[newid]] = reaction.metabolites[
                        metabolite
                    ]
                reaction.add_metabolites(input_metabolites, combine=False)
                reactions.append(reaction)
        newmodel.add_metabolites(metabolites)
        newmodel.add_reactions(reactions)
        return utilities

    def test_reaction_additions_againt_limits(self, reactions, directions, tests):
        filtered_rxn = []
        filtered_direction = []
        # Using "with" to ensure we don't alter the model with these tests
        model = self.cobramodel
        with model: # conserve the original model through WITH
            for rxn in reactions:
                if rxn.id in self.cobramodel.reactions:
                    rxn_obj = self.cobramodel.reactions.get_by_id(rxn.id)
                else:
                    rxn_obj = self.cobramodel.add_reactions([rxn])
                self.set_reaction_bounds_from_direction(rxn_obj, reactions[rxn])
                for test in tests:
                    testmodel = model
                    with testmodel:
                        self.apply_media_to_model(
                            test["media"],
                            test["default_uptake"],
                            test["default_excretion"],
                            testmodel,
                        )
                        self.set_objective_from_target_reaction(
                            test["target"], test["maximize"], testmodel
                        )
                        solution = self.cobramodel.optimize()
                        if test.maximize == 1:
                            if testmodel.objective.value() > test.limit:
                                filtered_tests.update({rxn_obj:reactions[rxn]})
        return filtered_tests
    
    def set_reaction_bounds_from_direction(self, reaction, direction, add=0):
        if direction == "<":
            reaction.lower_bound = -100
            if add:
                reaction.upper_bound = 0
        if direction == ">":
            reaction.upper_bound = 100
            if add:
                reaction.lower_bound = 0
        reaction.update_variable_bounds()
