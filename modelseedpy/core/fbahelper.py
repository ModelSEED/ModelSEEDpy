from __future__ import absolute_import

import logging
from chemicals import periodic_table
import re
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.biochem import from_local
#from Carbon.Aliases import false

logger = logging.getLogger(__name__)

elementmass = {}
for element in periodic_table:
    elementmass[element.symbol] = element.MW


class FBAHelper:

    @staticmethod
    def add_autodrain_reactions_to_community_model(model,auto_sink = ["cpd02701", "cpd11416", "cpd15302"]):
        #Adding missing drains in the base model
        for metabolite in model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid in auto_sink:
                if msid != "cpd11416" or metabolite.compartment == "c0":
                    if "EX_"+metabolite.id not in self.model.reactions and "DM_"+metabolite.id not in self.model.reactions and "SK_"+metabolite.id not in self.model.reactions:
                        drain_reaction = FBAHelper.add_drain_from_metabolite_id(self.model,metabolite.id,0,100,"DM_")

    @staticmethod
    def add_drain_from_metabolite_id(model, cpd_id, uptake, excretion, prefix='EX_', prefix_name='Exchange for '):
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
            drain_reaction = Reaction(id=f'{prefix}{cpd_id}',
                                      name=prefix_name + cobra_metabolite.name,
                                      lower_bound=-1*uptake, 
                                      upper_bound=excretion)
            drain_reaction.add_metabolites({cobra_metabolite : -1})
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            return drain_reaction
        return None
    
    @staticmethod
    def test_condition_list(model, condition_list, pkgmgr):
        for condition in condition_list:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(condition["media"])
            model.objective = condition["objective"]
            sol = model.optimize()
            if sol.objective_value >= condition["threshold"] and condition["is_max_threshold"]:
                logger.info("FAILED")
                return False
            elif sol.objective_value <= condition["threshold"] and not condition["is_max_threshold"]:
                logger.info("FAILED")
                return False
        return True
        
    @staticmethod
    def reaction_expansion_test(model, reaction_list, condition_list):
        # First knockout all reactions in the input list and save original bounds
        original_bound = []
        for item in reaction_list:
            if item[1] == ">":
                original_bound.append(item[0].upper_bound)
                item[0].upper_bound = 0
            else:
                original_bound.append(-1*item[0].lower_bound)
                item[0].lower_bound = 0
        # Now restore reactions one at a time
        count = 0
        filtered_list = []
        for item in reaction_list:
            if item[1] == ">":
                item[0].upper_bound = original_bound[count]
                if not FBAHelper.test_condition_list(model, condition_list):
                    item[0].upper_bound = 0
                    filtered_list.append(item)
            else:
                item[0].lower_bound = original_bound[count]
                if not FBAHelper.test_condition_list(model, condition_list):
                    item[0].lower_bound = 0
                    filtered_list.append(item)
            count += 1
        return filtered_list

    @staticmethod
    def set_reaction_bounds_from_direction(reaction, direction, add=0):
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
        for rxn in model.reactions:
            flux_values[rxn.id] = {
                'reverse': rxn.reverse_variable.primal,
                'forward': rxn.forward_variable.primal
            }

        return flux_values
    
    @staticmethod
    def modelseed_id_from_cobra_metabolite(metabolite):
        if re.search('^(cpd\d+)', metabolite.id) != None:
            m = re.search('^(cpd\d+)', metabolite.id)
            return m[1]
        #TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None
        
    def modelseed_id_from_cobra_reaction(reaction):
        if re.search('^(rxn\d+)', reaction.id) != None:
            m = re.search('^(rxn\d+)', reaction.id)
            return m[1]
        #TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None
    
    @staticmethod
    def metabolite_mw(metabolite):
        mw = 0
        elements = metabolite.elements
        for element in elements:
            if element not in elementmass:
                print("Missing mass for element "+element+" in compound "+metabolite.id+". Element will be ignored when computing MW")
            else:
                mw += elements[element]*elementmass[element]
        return mw
    
    @staticmethod    
    def elemental_mass():
        # Source:https://stackoverflow.com/questions/16699180/how-to-get-molecular-weight-of-a-compound-in-python/45557858
        return elementmass
    
    @staticmethod
    def get_modelseed_db_api(modelseed_path):
        return from_local(modelseed_path)
    
    @staticmethod
    def is_ex(reaction):
        # TODO: check for SBO
        if len(reaction.id) > 3 and (reaction.id[0:3] == "EX_" or reaction.id[0:3] == "DM_" or reaction.id[0:3] == "SK_"):
            return True
        return False

    @staticmethod
    def is_biomass(reaction):
        # TODO: check for SBO
        return reaction.id[0:3] == "bio"
