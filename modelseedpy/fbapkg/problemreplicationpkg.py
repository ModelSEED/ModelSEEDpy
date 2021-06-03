# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.fbapkg.revbinpkg import RevBinPkg
from modelseedpy.fbapkg.totalfluxpkg import TotalFluxPkg

#Base class for FBA packages
class ProblemReplicationPkg(BaseFBAPkg):
    def __init__(self,model):
        clonemodel = model.copy()
        BaseFBAPkg.__init__(self,clonemodel,"problem replication",{},{})
        
    def build_package(self,parameters):
        self.validate_parameters(parameters,["models"],{
            "shared_variable_packages":{}
        })
        #First loading shared variables into a hash
        shared_var_hash = {}
        for pkg in self.parameters["shared_variable_packages"]:
            for type in self.parameters["shared_variable_packages"][pkg]:
                if type in pkg.variables:
                    for objid in pkg.variables[type]:
                        shared_var_hash[pkg.variables[type][objid].name] = pkg.variables[type][objid]
        #Now copying over variables and constraints from other models and replacing shared variables
        count = 0
        for othermdl in self.parameters["models"]:
            newobj = []
            new_var_hash = {}
            for var in othermdl.variables:
                if var.name not in shared_var_hash:
                    newvar = var.copy()
                    newvar.name = var.name+"."+str(count)
                    new_var_hash[var.name] = newvar
                    newobj.append(newvar)       
            self.model.add_cons_vars(newobj)
            newobj = []
            for const in othermdl.constraints:
                substitutions = {}
                for var in const.variables:
                    if var.name not in shared_var_hash:
                        substitutions[var] = shared_var_hash[var.name]
                    else:
                        substitutions[var] = new_var_hash[var.name]
                expression = const.expression.xreplace(substitutions)
                newconst = self.model.problem.Constraint(
                    expression,lb=const.lb,ub=const.ub,name=const.name+"."+str(count)
                )
                newobj.append(newconst)       
            self.model.add_cons_vars(newobj)
            count += 1