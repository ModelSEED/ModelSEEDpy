import logging

import re
from cobra.core import Gene, Metabolite, Model, Reaction

logger = logging.getLogger(__name__)

#Parsing functions
def interpret_id(id):
    output = {
        "base_id" : None,
        "exchange" : 0,
        "model_instance" : 0,
        "compartment" : None,
        "compartment_index" : None
    }
    if id[0:3].lower() == "ex_":
        output["exchange"] = 1
        id = id[:-3]
    array = id.split("__")
    if len(array) > 1 and re.search('^\d+$',array[-1]):
        output["model_instance"] = int(array[-1])
        id = array[0]
        for i in range(1,len(array)-1):
            id += "_"+array[i]
    array = id.split("_")
    if len(array) > 1 and re.search('^[a-z]\d+$',array[-1]):
        output["compartment"] = array[-1][0:1]
        output["compartment_index"] = int(array[-1][1:])
        id = array[0]
        for i in range(1,len(array)-1):
            id += "_"+array[i]
    output["base_id"] = id
    return output

#Model manipulation
def add_drain_from_metabolite_id(model,cpd_id,lower_bound,upper_bound,media_hash = {},prefix = 'EX_', prefix_name = 'Exchange for ', sbo = 'SBO:0000627'):
    if lower_bound == None:
        lower_bound = -1*self.default_uptake
        if cpd_id in media_hash:
            lower_bound = media_hash[cpd_id]["lb"] 
    if upper_bound == None:
        upper_bound = self.default_excretion
        if cpd_id in media_hash:
            upper_bound = media_hash[cpd_id]["ub"]                            
    if cpd_id in model.metabolites:
        id = prefix + cpd_id
        cobra_metabolite = model.metabolites.get_by_id(cpd_id)
        object_stoichiometry = {cobra_metabolite : -1}
        drain_reaction = Reaction(id=id, 
                                  name= prefix_name + cobra_metabolite.name, 
                                  lower_bound=lower_bound, 
                                  upper_bound=upper_bound)
        drain_reaction.add_metabolites(object_stoichiometry)
        drain_reaction.annotation[self.SBO_ANNOTATION] = sbo    
        return drain_reaction
    return None

def replicate_model(model,count):
    newmodel = Model(model.id+"_rep"+str(count))
    metabolites = []
    reactions = []
    metabolite_hash = {}
    for i in range(0,count):
        for metabolite in model.metabolites:
            metabolite = metabolite.copy()
            metabolite.id = metabolite.id + "__" + str(i)
            metabolite_hash[metabolite.id] = metabolite
            metabolites.append(metabolite)
        for reaction in model.reactions:
            reaction = reaction.copy()
            reaction.id = reaction.id + "__" + str(i)
            input_metabolites = {}
            for metabolite in reaction.metabolites:
                newid = metabolite.id + "__" + str(i)
                input_metabolites[metabolite_hash[newid]] = reaction.metabolites[metabolite]
            reaction.add_metabolites(input_metabolites,combine=False)
            reactions.append(reaction)
    newmodel.add_metabolites(metabolites)
    newmodel.add_reactions(reactions)
    return newmodel

#Manipulating reaction bounds
def apply_media_to_model(model,media_hash = {},default_uptake,default_excretion):
    for reaction in model.reactions:
        output = interpret_id(reaction.id)
        if output["exchange"] == 1:
            if output["base_id"] in media_hash:
                reaction.lower_bound = media_hash[output["base_id"]]["lb"]
                reaction.upper_bound = media_hash[output["base_id"]]["ub"]
            else:
                reaction.lower_bound = -1*default_uptake
                reaction.upper_bound = default_excretion
            reaction.update_variable_bounds()

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

#Manipulating model objective
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

def convert_objective_to_constraint(model,lower_bound,upper_bound):
    old_obj_variable = model.problem.Variable(
        name="old_objective_variable",
        lb=lower_bound,ub=upper_bound
    )
    old_obj_constraint = model.problem.Constraint(
        self.cobramodel.solver.objective.expression - old_obj_variable,
        lb=0,
        ub=0,
        name="old_objective_constraint",
    )
    self.cobramodel.add_cons_vars([old_obj_variable, old_obj_constraint])

#Data retrieval functions
def compute_flux_values_from_variables(model):
    flux_values = {}
    for rxnobj in model.reactions:
        flux_values[rxnobj.id] = {}
        flux_values[rxnobj.id]["reverse"] = rxnobj.reverse_variable.primal
        flux_values[rxnobj.id]["forward"] = rxnobj.forward_variable.primal
    return flux_values