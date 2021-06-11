import logging

import re
import copy
import cobra
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg import GapfillingPkg, KBaseMediaPkg

#from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)


class MSEditorAPI:

    @staticmethod
    def remove_reactions(model,rxn_id_list = []):
        # throws a key error if not in model reactions
        #for rxnid in rxn_id_list:
            #model.remove_reactions(model.reactions.get_by_id(rxnid))
        model.remove_reactions(rxn_id_list)
        
    @staticmethod    
    def edit_reaction(model,rxn_id,direction = None,gpr = None,genome = None):
        # Direction: =>, <=, or <=>
        if direction is not None:
            lower_bound = model.reactions.get_by_id(rxn_id).lower_bound
            upper_bound = model.reactions.get_by_id(rxn_id).upper_bound

            if lower_bound < 0 and upper_bound > 0: # rxn_id is reversible
                if direction == "=>":
                    model.reactions.get_by_id(rxn_id).lower_bound = 0
                elif direction == "<=":
                    model.reactions.get_by_id(rxn_id).upper_bound = 0
            elif lower_bound == 0 and upper_bound > 0: # rxn_id is forward only
                if direction == "<=":
                    model.reactions.get_by_id(rxn_id).lower_bound = -1*upper_bound
                    model.reactions.get_by_id(rxn_id).upper_bound = 0
                elif direction == "<=>":
                    model.reactions.get_by_id(rxn_id)._lower_bound = -1*upper_bound
            elif lower_bound < 0 and upper_bound == 0: # rxn_id is reverse only
                if direction == "=>":
                    model.reactions.get_by_id(rxn_id).lower_bound = 0
                    model.reactions.get_by_id(rxn_id).upper_bound = -1*lower_bound
                elif direction == "<=>":
                    model.reactions.get_by_id(rxn_id)._upper_bound = -1*lower_bound
        
        # Specify GPR as a string with boolean conditions (e.g. "(b0001 and b0002) or b1010").
        if gpr is not None:
            model.reactions.get_by_id(rxn_id).gene_reaction_rule = gpr   
    
    @staticmethod
    def edit_biomass_compound(model,biomass_id,cpd_id,new_coef,rescale = 1):
        model.reactions.get_by_id(biomass_id).add_metabolites({model.metabolties.get_by_id(cpd_id): new_coef}, combine=False)

    @staticmethod

    def add_custom_reaction(model,rxn_id,MSEquation,gpr = None,genome = None):
        pass
        #mseqn = MSEquation.build_from_palsson_string("cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004")
    
    @staticmethod  
    def add_ms_reaction(model,rxn_id,compartments,modelseed): #Andrew
        pass
        
    @staticmethod  
    def copy_model_reactions(model,source_model,rxn_id_list = []):
        for rxnid in rxn_id_list:
            model.add_reactions(source_model.reactions.get_by_id(rxnid))


class MSEquation:

    def __init__(self, stoichiometry):
        pass

    @staticmethod
    def build_from_palsson_string(equation_string):
        #cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004
        # get substrings for either side of the equation
        reactants_substring_list = equation_string[0:equation_string.find('>') - 1].split('+')
        products_substring_list = equation_string[equation_string.find('>') + 1:len(equation_string)].split('+')
        # clean up our substrings:
        for i in range(len(reactants_substring_list)):
            # remove whitespace from the front
            while (reactants_substring_list[i][0] == ' '):
                reactants_substring_list[i] = reactants_substring_list[i][1:]
            # remove whitespace from the back
            while (reactants_substring_list[i][-1] == ' '):
                reactants_substring_list[i] = reactants_substring_list[i][:-1]
        for i in range(len(products_substring_list)):
            # remove whitespace from the front
            while (products_substring_list[i][0] == ' '):
                products_substring_list[i] = products_substring_list[i][1:]
            # remove whitespace from the back
            while (products_substring_list[i][-1] == ' '):
                products_substring_list[i] = products_substring_list[i][:-1]
        variables = {}
        # add reactants to the dictionary
        for reactant in reactants_substring_list:
            coeficient = -1
            if reactant[0] == '(':
                number = ''
                position = 1
                while reactant[position] != ')':
                    number += reactant[position]
                    position += 1
                coeficient = -1 * int(number)
                reactant = reactant[position + 1:]
            identifier = reactant[0]
            if reactant[-1] == ']':
                s = ''
                position = -2
                while reactant[position] != '[':
                    s = reactant[position] + s
                    position -= 1
                identifier = s
                reactant = reactant[:position]
            variables[(reactant, identifier)] = coeficient
        # add the products to the dictionary
        for reactant in products_substring_list:
            coeficient = 1
            if reactant[0] == '(':
                number = ''
                position = 1
                while reactant[position] != ')':
                    number += reactant[position]
                    position += 1
                coeficient = int(number)
                reactant = reactant[position + 1:]
            identifier = reactant[0]
            if reactant[-1] == ']':
                s = ''
                position = -2
                while reactant[position] != '[':
                    s = reactant[position] + s
                    position -= 1
                identifier = s
                reactant = reactant[:position]
            variables[(reactant, identifier)] = coeficient
        return variables
