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
        for rxnid in rxn_id_list:
        
    @staticmethod    
    def edit_reaction(model,rxn_id,direction = None,gpr = None,genome = None]):    
    
    @staticmethod
    def edit_biomass_compound(model,biomass_id,cpd_id,new_coef,rescale = 1):
        
    @staticmethod
    def add_custom_reaction(model,rxn_id,MSEquation,gpr = None,genome = None):
    
    @staticmethod  
    def add_ms_reaction(model,rxn_id,compartments,modelseed):#Andrew
        
    @staticmethod  
    def copy_model_reactions(model,source_model,rxn_id_list = []):

class MSEquation:
    @staticmethod
    def build_from_palsson_string(equation_string):
        #cpd00001 + cpd00002[e] => (2)cpd00003 + cpd00004
        
