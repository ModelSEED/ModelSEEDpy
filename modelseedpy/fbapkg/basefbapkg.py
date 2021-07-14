# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
import json as _json
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation"""
    pass

#Base class for FBA packages
class BaseFBAPkg:
    def __init__(self, model, name, variable_types={}, constraint_types={}, parent=None):
        self.model = model
        self.name = name
        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        if self.pkgmgr == None:
            self.pkgmgr = MSPackageManager.get_pkg_mgr(model,1)
        self.pkgmgr.addpkgobj(self)
        self.constraints = dict()
        self.variables = dict()
        self.parameters = dict()
        self.new_reactions = dict()
        self.variable_types = variable_types
        self.constraint_types = constraint_types
        
        for type in variable_types:
            self.variables[type] = dict()
        for type in constraint_types:
            self.constraints[type] = dict()
    
    def validate_parameters(self,params,required,defaults):
        for item in required:
            if item not in params:
                raise ValueError('Required argument '+item+' is missing!')
        for key in params:
            self.parameters[key] = params[key]
        for key in defaults:
            if key not in params:
                self.parameters[key] = defaults[key]
    
    def clear(self):
        objects = []
        for type in self.variables:
            for object in self.variables[type]:
                objects.append(self.variables[type][object])
        for type in self.constraints:
            for object in self.constraints[type]:
                objects.append(self.constraints[type][object])
        self.model.remove_cons_vars(objects)
        self.variables = {}
        self.constraints = {}
        
    def build_variable(self,type,lower_bound,upper_bound,vartype,object = None):
        name = None
        if self.variable_types[type] == "none":
            count = len(self.variables[type])
            name = str(count+1)
        elif self.variable_types[type] == "string":
            name = object
        else:
            name = object.id
        if name not in self.variables[type]:
            self.variables[type][name] = self.model.problem.Variable(name+"_"+type, lb=lower_bound,ub=upper_bound,type=vartype)
            self.model.add_cons_vars(self.variables[type][name])
        return self.variables[type][name]
        
    def build_constraint(self,type,lower_bound,upper_bound,coef = {},object = None):
        name = None
        if self.constraint_types[type] == "none":
            count = len(self.constraints[type])
            name = str(count+1)
        elif self.constraint_types[type] == "string":
            name = object
        else:
            name = object.id
        if name not in self.constraints[type]:
            self.constraints[type][name] = self.model.problem.Constraint(
                Zero,lb=lower_bound,ub=upper_bound,name=name+"_"+type
            )
            self.model.add_cons_vars(self.constraints[type][name])
            self.model.solver.update()
            if len(coef) > 0:
                self.constraints[type][name].set_linear_coefficients(coef)
            self.model.solver.update()
        return self.constraints[type][name]
    
    def all_variables(self):
        return self.pkgmgr.all_variables()
    
    def all_constraints(self):
        return self.pkgmgr.all_constraints()