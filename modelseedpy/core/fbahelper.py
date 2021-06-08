import logging

import re
from cobra.core import Gene, Metabolite, Model, Reaction

logger = logging.getLogger(__name__)

class FBAHelper:
    @staticmethod
    def add_drain_from_metabolite_id(model,cpd_id,uptake,excretion,prefix = 'EX_', prefix_name = 'Exchange for '):
        if cpd_id in model.metabolites:
            id = prefix + cpd_id
            cobra_metabolite = model.metabolites.get_by_id(cpd_id)
            object_stoichiometry = {cobra_metabolite : -1}
            drain_reaction = Reaction(id=id, 
                                      name= prefix_name + cobra_metabolite.name, 
                                      lower_bound=-1*uptake, 
                                      upper_bound=excretion)
            drain_reaction.add_metabolites(object_stoichiometry)
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            return drain_reaction
        return None

    @staticmethod
    def set_reaction_bounds_from_direction(reaction,direction,add=0):
        if direction == "<":
            reaction.lower_bound = -100
            if add == 0:
                reaction.upper_bound = 0
        if direction == ">":
            reaction.upper_bound = 100
            if add == 0:
                reaction.lower_bound = 0
        reaction.update_variable_bounds()

    @staticmethod
    def set_objective_from_target_reaction(model,target_reaction,maximize = 1):
        target_reaction = model.reactions.get_by_id(target_reaction)
        sense = "max"
        if maximize == 0:
            sense = "min"
        target_objective = model.problem.Objective(
            1 * target_reaction.flux_expression,
            direction=sense)
        model.objective = target_objective
        return target_reaction

    @staticmethod
    def compute_flux_values_from_variables(model):
        flux_values = {}
        for rxnobj in model.reactions:
            flux_values[rxnobj.id] = {}
            flux_values[rxnobj.id]["reverse"] = rxnobj.reverse_variable.primal
            flux_values[rxnobj.id]["forward"] = rxnobj.forward_variable.primal
        return flux_values