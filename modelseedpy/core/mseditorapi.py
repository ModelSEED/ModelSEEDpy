import logging

from modelseedpy.fbapkg import GapfillingPkg, KBaseMediaPkg
from cobra.core import Gene, Metabolite, Model, Reaction
from optlang.symbolics import Zero, add
from modelseedpy.core import FBAHelper
import cobra
import copy
import re


#from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)


class MSEditorAPI:

    @staticmethod
    def remove_reactions(model,rxn_id_list = []):
        for rxnid in rxn_id_list:
            pass
        pass
        
    @staticmethod    
    def edit_reaction(model,rxn_id,direction = None,gpr = None,genome = None):
        pass
    
    @staticmethod
    def edit_biomass_compound(model,biomass_id,cpd_id,new_coef,rescale = 1):
        pass

    @staticmethod

    def add_custom_reaction(model,rxn_id,MSEquation,gpr = None,genome = None):
        pass
    
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
        pass


class MSEquation:

    def __init__(self, stoichiometry):
        pass

    @staticmethod
    def build_from_palsson_string(equation_string):
        #cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004
        pass
