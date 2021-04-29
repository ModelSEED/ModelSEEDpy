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
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,"reaction use",{"revbin":["reaction",0,1],"potential":["metabolite",0,1000]},{"revbinF":"reaction","revbinR":"reaction","thermo":"reaction"})

    def build_package(self,filter = None):
        for metabolite in self.model.metabolites:
            self.build_variable("potential",metabolite)
        for reaction in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if filter == None or reaction.id in filter:
                self.build_variable("revbin",reaction)
                self.build_constraint("revbin",reaction)
                self.build_constraint("thermo",reaction)
    
    def build_variable(self,type,object = None):
        if type == "potential":
            return BaseFBAPkg.build_variable(self,type,0,1000,None,object)
        elif type == "revbin":
            return BaseFBAPkg.build_variable(self,type,0,1,"binary",object)
    
    def build_constraint(self,type,object = None):
        if type == "thermo":
            coef = {self.variables["revbin"][reaction.id] : 1000}
            for metabolite in object.metabolites:
                coef[self.variables["potential"][metabolite.id]] = object.metabolites[metabolite]
            return BaseFBAPkg.build_constraint(self,type,0,1000,coef,object)
        elif type == "revbin":
            BaseFBAPkg.build_constraint(self,type+"F",-1000,None,{self.variables["revbin"][reaction.id]:1000,object.forward_variable:-1})
            return BaseFBAPkg.build_constraint(self,type+"R",0,None,{self.variables["revbin"][reaction.id]:-1000,object.reverse_variable:-1})