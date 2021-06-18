import logging

from modelseedpy.fbapkg import GapfillingPkg, KBaseMediaPkg
from cobra.core import Gene, Metabolite, Model, Reaction
from optlang.symbolics import Zero, add
from modelseedpy.core import FBAHelper
import cobra
import copy
import re


#from modelseedpy.core.msgenome import MSGenome
# FIXME: fix unused imports
import re
import copy
import cobra
from cobra.core import Gene, Metabolite, Model, Reaction

logger = logging.getLogger(__name__)


class MSEditorAPI:

    @staticmethod
    def remove_reactions(model, rxn_id_list = []):
        for rxn_id in rxn_id_list:
            if not model.reactions.has_id(rxn_id):
                raise Exception('Reaction', rxn_id, 'is not in the model.')
            model.remove_reactions([rxn_id])

    # edit reaction progam
    # ASSUMPTIONS:
    # an arrow will exist in the program, either =>, <=>, or <=
    @staticmethod
    def edit_reaction(model, rxn_id, direction=None, gpr=None, genome=None):
        # Direction: =>, <=, or <=>
        if direction is not None:
            lower_bound = model.reactions.get_by_id(rxn_id).lower_bound
            upper_bound = model.reactions.get_by_id(rxn_id).upper_bound

            if lower_bound < 0 and upper_bound > 0:  # rxn_id is reversible
                if direction == "=>":
                    model.reactions.get_by_id(rxn_id).lower_bound = 0
                elif direction == "<=":
                    model.reactions.get_by_id(rxn_id).upper_bound = 0
            elif lower_bound == 0 and upper_bound > 0:  # rxn_id is forward only
                if direction == "<=":
                    model.reactions.get_by_id(rxn_id).lower_bound = -1 * upper_bound
                    model.reactions.get_by_id(rxn_id).upper_bound = 0
                elif direction == "<=>":
                    model.reactions.get_by_id(rxn_id).lower_bound = -1 * upper_bound
            elif lower_bound < 0 and upper_bound == 0:  # rxn_id is reverse only
                if direction == "=>":
                    model.reactions.get_by_id(rxn_id).lower_bound = 0
                    model.reactions.get_by_id(rxn_id).upper_bound = -1 * lower_bound
                elif direction == "<=>":
                    model.reactions.get_by_id(rxn_id).upper_bound = -1 * lower_bound

        # Specify GPR as a string with boolean conditions (e.g. "(b0001 and b0002) or b1010").
        try:
            if gpr is not None:
                model.reactions.get_by_id(rxn_id).gene_reaction_rule = gpr
        except:
            raise Exception('invalid gpr statement, check parentheses')  # not working, unsure exactly why
    
    @staticmethod
    def edit_biomass_compound(model,biomass_id,cpd_id,new_coef,rescale = 1):
        if biomass_id in model.reactions:
            if cpd_id in model.metabolites:
                model.reactions.get_by_id(biomass_id).add_metabolites({model.metabolites.get_by_id(cpd_id): new_coef},
                                                                      combine=False)
            else:
                raise Exception('Metabolite', cpd_id, ' is not in the model.')
        else:  # if there is no biomass reaction
            biomass_rxn = Reaction(biomass_id)
            model.add_reaction(biomass_rxn)
            if cpd_id in model.metabolites:
                biomass_rxn.add_metabolites({model.metabolites.get_by_id(cpd_id): new_coef})
            else:
                raise Exception('Metabolite ', cpd_id, ' is not in the model.')

    @staticmethod
    def compute_molecular_weight(model, metabolite_id):
        if metabolite_id not in model.metabolites:
            raise Exception('Error, metabolite', metabolite_id, 'not found in the model')
        return model.metabolites.get_by_id(metabolite_id).formula_weight

    @staticmethod

    def add_custom_reaction(model,rxn_id,MSEquation,gpr = None,genome = None):
        new_rxn = Reaction(id=rxn_id)
        # going on the assumption that all metabolites are present in the model
        metabolites = {}
        for key in MSEquation.equation:
            met_id = key[0] + '_' + key[1]
            if met_id in model.metabolites:
                metabolites[met_id] = MSEquation.equation[key]
            else:
                raise Exception("Error,", met_id, "not in model metabolites list")
        model.add_reaction(new_rxn)
        # new_rxn.gene_reaction_rule = gpr
        new_rxn.add_metabolites(metabolites)

        # adjust the bounds based on the arrow direction  -1000, 1000, 0
        if MSEquation.direction == 'left':
            new_rxn.lower_bound = -1000
            new_rxn.upper_bound = 0
        if MSEquation.direction == 'reversable':
            new_rxn.lower_bound = -1000
            new_rxn.upper_bound = 1000
    
    @staticmethod  
    def add_ms_reaction(model, rxn_id, modelseed, compartment_equivalents = {'0':'c0', '1':'e0'}, direction = 'forward'):#Andrew
        ''' Add a reaction with ModelSEED parameters to the FBA simulation
        "model" (COBRA object): The metabolic model that is defined by COBRApy
        "rxn_id" (Python object, string): The ModelSEED reaction id that will be added to the model
        "Compartment_equivalents" (Python object, dictionary): The compartment codes that are used in formating the passed reactions 
        "direction" (Python object, string): The direction of the defined reaction
        "modelseed" (ModelSEED object): The ModelSEED database that describes the reaction and metabolites of the argument    
        '''
        modelseed_reaction = modelseed.get_seed_reaction(rxn_id)
        reaction_stoich = modelseed_reaction.cstoichiometry
        cobra_reaction = cobra.Reaction(rxn_id)
        cobra_reaction.name = modelseed_reaction.data['name']

        metabolites_to_add = {}
        for metabolite, stoich in reaction_stoich.items():
            id = metabolite[0]
            compound = modelseed.get_seed_compound(id).data
            compartment_number = metabolite[1]
            compartment_string = compartment_equivalents[compartment_number]        

            metabolites_to_add[cobra.Metabolite(id, name = compound['name'], compartment = compartment_string)] = stoich

        cobra_reaction.add_metabolites(metabolites_to_add)
        cobra_reaction.reaction

        if direction == 'reversible':
            cobra_reaction.lower_bound = -1000
        elif direction == 'backward':
            cobra_reaction.lower_bound = -1000
            cobra_reaction.upper_bound = 0
        elif direction == 'forward':
            pass
        else:
            directions = ['forward', 'backward', 'reversible']
            print('ERROR: The \'direction\' argument is not among the accepted {}, {}, and {} values.'.format(directions[0], directions[1], directions[2]))
            direction = input('What is the direction of the reaction?')
            while direction not in directions:
                print('ERROR: The \'direction\' argument is not among the accepted {}, {}, and {} values.'.format(directions[0], directions[1], directions[2]))
                direction = input('What is the direction of the reaction?')
            MSEditorAPI.add_ms_reaction(model, rxn_id, modelseed, direction = direction)


        model.add_reactions([cobra_reaction])
        
    @staticmethod  
    def copy_model_reactions(model,source_model,rxn_id_list = []):
        for rxnid in rxn_id_list:
            if rxnid in source_model.reactions:
                model.add_reaction(source_model.reactions.get_by_id(rxnid))
            else:
                raise Exception('The reaction', rxnid, 'in the reaction list is not in the model.')

    @staticmethod
    def copy_all_model_reactions(model,source_model):  #new method that copies all reactions, may not be necessary
        for rxnid in source_model.reactions:
            model.add_reaction(source_model.reactions.get_by_id(rxnid))


class MSEquation:

    def __init__(self, stoichiometry, d):
        self.equation = stoichiometry
        self.direction = d

    @staticmethod
    def build_from_palsson_string(equation_string, default_group='c'):  # add default group

        def clean_ends(lst):
            for i in range(len(lst)):
                # remove whitespace from the front
                while (lst[i][0] == ' '):
                    lst[i] = lst[i][1:]
                # remove whitespace from the back
                while (lst[i][-1] == ' '):
                    lst[i] = lst[i][:-1]
            return lst

        def get_coef_and_group(lst, return_dict,
                               side):  # for side variable, -1 is left side, 1 is right side, for coeficients
            for reactant in lst:
                coeficient = side
                if reactant.find('(') != -1 or reactant.find(
                        ')') != -1:  # if one is present, check to make sure both there
                    if equation_string.find(')') == -1:
                        raise Exception("Error, ')' character missing in string", reactant)
                    if equation_string.find('(') == -1:
                        raise Exception("Error, '(' character missing in string", reactant)
                    number = ''
                    position = 1
                    while reactant[position] != ')':
                        number += reactant[position]
                        position += 1
                    coeficient = side * int(number)
                    reactant = reactant[position + 1:]

                identifier = default_group
                if reactant.find('[') != -1 or reactant.find(
                        ']') != -1:  # if one is present, check to make sure both there
                    # check to see both are present
                    if equation_string.find(']') == -1:
                        raise Exception("Error, ']' character missing in string", reactant)
                    if equation_string.find('[') == -1:
                        raise Exception("Error, '[' character missing in string", reactant)
                    s = ''
                    position = -2
                    while reactant[position] != '[':
                        s = reactant[position] + s
                        position -= 1
                    identifier = s
                    reactant = reactant[:position]

                return_dict[(reactant, identifier)] = coeficient
            return return_dict

        # check for the '=' character, throw exception otherwise
        if equation_string.find('=') == -1:
            raise Exception("Error, '=' character missing, unable to split string", equation_string)

        # check direction
        reversible = False
        right = False
        left = False
        ret_str = ''
        reversible = equation_string.find('<=>') != -1
        if reversible:
            ret_str = '='
        else:  # not two ways, so check right
            right = equation_string.find('=>') != -1
            if right:
                ret_str = '>'
            else:  # if not right, check left
                left = equation_string.find('<=') != -1
                if left:
                    ret_str = '<'
                else:  # if not left, error
                    ret_str = "?"

        # cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004
        # get substrings for either side of the euqation
        reactants_substring_list = equation_string[0:equation_string.find('=') - 1].split('+')
        products_substring_list = equation_string[equation_string.find('=') + 2:len(equation_string)].split('+')

        # clean up our substrings:
        clean_ends(reactants_substring_list)
        clean_ends(products_substring_list)

        variables = {}
        # add reactants to the dictionary

        get_coef_and_group(reactants_substring_list, variables, -1)
        get_coef_and_group(products_substring_list, variables, 1)

        ret_mse = MSEquation(variables, ret_str)

        return ret_mse
