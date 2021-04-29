# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
import json as _json
from cobra.core import Gene, Metabolite, Model, Reaction

#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation"""
    pass

#Base class for FBA packages
class BaseFBAPkg():
    def __init__(self,model,name,variable_types={},constraint_types={},parent = None):
        self.model = model
        self.name = name
        self.childpkgs = dict()
        self.parentpkg = parent
        self.constraints = dict()
        self.variables = dict()
        self.variable_types = vartypes
        self.constraint_types = consttypes
        for type in variable_types:
            self.variables[type] = dict()
        for type in constraint_types:
            self.constraints[type] = dict()
    
    def build_variable(self,type,lower_bound,upper_bound,vartype,object = None):
        name = None
        if self.variable_types[type][0] == "none":
            count = len(self.variables[type])
            name = str(count+1)
        else:
            name = object.id
        if name not in self.variables[type]:
            self.variables[type][name] = self.cobramodel.problem.Variable(name+"_"+type, lb=lower_bound,ub=upper_bound,type=vartype)
            self.model.add_cons_vars(self.variables[type][name])
        return self.variables[type][name]
        
    def build_constraint(self,type,lower_bound,upper_bound,coef = {},object = None):
        name = None
        if self.variable_types[type][0] == "none":
            count = len(self.variables[type])
            name = str(count+1)
        else:
            name = object.id
        if name not in self.constraints[type][name]:
            self.constraints[type][name] = self.model.problem.Constraint(
                Zero,lb=lower_bound,ub=upper_bound,name=name+"_"+type
            )
            self.model.add_cons_vars(self.constraints[type][name])
            if len(coef) > 0:
                self.constraints[type][name].set_linear_coefficients(coef)
            self.model.solver.update()
        return self.constraints[type][name]