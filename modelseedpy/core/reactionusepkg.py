# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
import json as _json
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core.basefbapkg import BaseFBAPkg

#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation"""
    pass

#Base class for FBA packages
class ReactionUsePkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,"reaction use",{"fu":["reaction",0,1],"ru":["reaction",0,1]},{"fu":"reaction","ru":"reaction","exclusion":"other"})

    def build_package(self,filter = None):
        for reaction in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if filter == None:
                self.build_variable(reaction,"=")
                self.build_constraint(reaction)
            elif reaction.id in filter:
                self.build_variable(reaction,filter[reaction.id])
                self.build_constraint(reaction)
    
    def build_variable(self,object,direction):
        variable = None
        if (direction == ">" or direction == "=") and object.upper_bound > 0 and object.id not in self.variables["fu"]:
            variable = BaseFBAPkg.build_variable(self,"fu",0,1,"binary",object)
        if (direction == "<" or direction == "=") and object.lower_bound < 0 and object.id not in self.variables["ru"]:
            variable = BaseFBAPkg.build_variable(self,"ru",0,1,"binary",object)
        return variable
        
    def build_constraint(self,object = None):
        constraint = None
        if object.id not in self.variables["fu"]:
            constraint = BaseFBAPkg.build_constraint(self,"fu",0,None,{self.variables["fu"][object.id]:1000,object.forward_variable:-1},object)
        if object.id not in self.variables["ru"]:
            constraint = BaseFBAPkg.build_constraint(self,"ru",0,None,{self.variables["ru"][object.id]:1000,object.reverse_variable:-1},object)
        return constraint
    
    def build_exclusion_constraint(self,flux_values):
        count = len(self.constraints["exclusion"])
        solution_coef = {}
        solution_size = 0
        for rxnid in flux_values:
            if flux_values[rxnid] > Zero:
                solution_size += 1
                solution_coef[self.variables["fu"][rxnid]] = 1
            elif flux_values[rxnid] < -1*Zero:
                solution_size += 1
                solution_coef[self.variables["ru"][rxnid]] = 1            
        if len(solution_coef) > 0:
            const_name = "exclusion."+str(count+1)
            self.constraints["exclusion"][const_name] = self.model.problem.Constraint(
                Zero,lb=None,ub=(solution_size-1),name=const_name
            )
            self.cobramodel.add_cons_vars(self.constraints["exclusion"][const_name])
            self.cobramodel.solver.update()
            self.constraints["exclusion"][const_name].set_linear_coefficients(solution_coef)
            return self.constraints["exclusion"][const_name]
        return None