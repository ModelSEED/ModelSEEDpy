# -*- coding: utf-8 -*-
import logging

logger = logging.getLogger(__name__)
from cobra.core import Reaction
import cobra
import re


class MSEditorAPI:
    @staticmethod
    def remove_reactions(model, rxn_id_list=[]):
        model_reactions = " ".join(
            [rxn.id for rxn in model.reactions]
        )  # removed from loop for greater efficiency
        for rxn_id in rxn_id_list:
            if not model.reactions.has_id(rxn_id):
                compartment = re.search(f"(?<={rxn_id})(\_\w\d)", model_reactions)
                if not compartment:
                    raise Exception("Reaction", rxn_id, "is not in the model.")
                else:
                    rxn_id += compartment.group()
            model.remove_reactions([rxn_id])

    # edit reaction progam
    # ASSUMPTIONS:
    # an arrow will exist in the program, either =>, <=>, or <=
    @staticmethod
    def edit_reaction(model, rxn_id, direction=None, gpr=None):
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
            if gpr:
                model.reactions.get_by_id(rxn_id).gene_reaction_rule = gpr
        except:
            raise Exception(
                f'The gpr {gpr} is invalid. Perhaps check parentheses.'
            )  # not working, unsure exactly why

    @staticmethod
    def edit_biomass_compound(model, biomass_id, cpd_id, new_coef, rescale=1):
        if biomass_id in model.reactions:
            if cpd_id in model.metabolites:
                model.reactions.get_by_id(biomass_id).add_metabolites(
                    {model.metabolites.get_by_id(cpd_id): new_coef}, combine=False
                )
            else:
                raise Exception("Metabolite", cpd_id, " is not in the model.")
        else:  # if there is no biomass reaction
            biomass_rxn = Reaction(biomass_id)
            model.add_reaction(biomass_rxn)
            if cpd_id in model.metabolites:
                biomass_rxn.add_metabolites(
                    {model.metabolites.get_by_id(cpd_id): new_coef}
                )
            else:
                raise Exception("Metabolite ", cpd_id, " is not in the model.")

    @staticmethod
    def compute_molecular_weight(model, metabolite_id):
        if metabolite_id not in model.metabolites:
            raise Exception(
                "Error, metabolite", metabolite_id, "not found in the model"
            )
        return model.metabolites.get_by_id(metabolite_id).formula_weight

    @staticmethod
    def add_custom_reaction(model, rxn_id, MSEquation, gpr=None, genome=None):
        new_rxn = Reaction(id=rxn_id)
        # going on the assumption that all metabolites are present in the model
        metabolites = {}
        for key in MSEquation.equation:
            met_id = key[0] + "_" + key[1]
            if met_id in model.metabolites:
                metabolites[met_id] = stoich
            else:
                raise ValueError(f"The {met_id} is not in model.")
        model.add_reaction(new_rxn)
        new_rxn.add_metabolites(metabolites)
        if gpr:
            new_rxn.gene_reaction_rule = gpr

        # adjust the bounds based on the arrow direction  -1000, 1000, 0
        if MSEquation.direction == "left":
            new_rxn.lower_bound = -1000
            new_rxn.upper_bound = 0
        if MSEquation.direction == "reversable":
            new_rxn.lower_bound = -1000

    @staticmethod  
    def add_ms_reaction(model, rxn_id, modelseed, compartment_equivalents = {'0':'c0', '1':'e0'}, direction = '>'):#Andrew
        cobra_reaction = cobra.Reaction(rxn_id)
        modelseed_reaction = modelseed.get_seed_reaction(rxn_id)
        cobra_reaction.name = modelseed_reaction.data['name']
        reaction_stoich = modelseed_reaction.cstoichiometry

        metabolites_to_add = {}
        for metabolite, stoich in reaction_stoich.items():
            metabolites_to_add[cobra.Metabolite(
                    metabolite[0], name = modelseed.get_seed_compound(metabolite[0]).data['name'], 
                    compartment = compartment_equivalents[metabolite[1]])
                ] = stoich

        cobra_reaction.add_metabolites(metabolites_to_add)
        cobra_reaction.lower_bound = 0
        cobra_reaction.upper_bound = 1000
        if direction == '=':
            cobra_reaction.lower_bound = -1000
        elif direction == '<':
            cobra_reaction.lower_bound = -1000
            cobra_reaction.upper_bound = 0
        model.add_reactions([cobra_reaction])

    @staticmethod
    def copy_model_reactions(model, source_model, rxn_id_list=[]):
        for rxnid in rxn_id_list:
            if rxnid in source_model.reactions:
                model.add_reactions([source_model.reactions.get_by_id(rxnid)])
            else:
                raise ValueError(f'The {rxnid} reaction ID is not in the source model, and thus cannot be added to the model.')

    @staticmethod
    def copy_all_model_reactions(model,source_model):  #new method that copies all reactions, may not be necessary
        model.add_reactions([source_model.reactions.get_by_id(rxn.id) for rxn in source_model.reactions if rxn not in model.reactions])

class MSEquation:
    def __init__(self, stoichiometry, direction = None):
        self.equation = stoichiometry; self.direction = direction

    @staticmethod
    def build_from_palsson_string(equation_string, default_group='c'):  # add default group
        def clean_ends(lst):
            return [i.strip() for i in lst]

        def get_coef_and_group(lst, return_dict, side):
            # for side variable, -1 is left side, 1 is right side, for coeficients
            for reagent in lst:
                coeficient = side
                identifier = default_group
                if '(' in reagent and ')' in reagent:  
                    number = ''
                    position = 1
                    while reagent[position] != ')':
                        number += reagent[position]
                        position += 1
                    coeficient = side * float(number)
                    reagent = reagent[position+1: ]
                elif '[' in reagent and ']' in reagent: 
                    s = ''
                    position = -2
                    while reagent[position] != '[':
                        s = reagent[position] + s
                        position -= 1
                    identifier = s
                    reagent = reagent[:position]
                elif any([x in reagent for x in ['(', ')', '[', ']']]):
                    raise ValueError("A closing or opening parentheses or bracket is missing in the reaction string", reagent)
                return_dict[(reagent.strip(), identifier)] = coeficient
            return return_dict

        # check for the '=' character, throw exception otherwise
        if '=' not in equation_string:
            raise ValueError(f"Error: The '=' character is missing; hence, the reaction string {equation_string} cannot be split.")
        if '<=>' in equation_string:
            direction = '='
        elif '=>' in equation_string: 
            direction = '>'
        elif '<=' in equation_string:
            direction = '<'
        else:
            direction = '?'

        # get substrings for either side of the equation
        reactants_substring_list = equation_string[0:equation_string.find('=') - 1].split('+')
        products_substring_list = equation_string[equation_string.find('=') + 2:len(equation_string)].split('+')
        reactant_dict = get_coef_and_group([x.strip() for x in reactants_substring_list], {}, -1)
        products_dict = get_coef_and_group([x.strip() for x in products_substring_list], reactant_dict, 1)
        return MSEquation(reactant_dict.update(products_dict), direction)
