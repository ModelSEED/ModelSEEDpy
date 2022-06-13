# -*- coding: utf-8 -*-

from __future__ import absolute_import

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class ObjConstPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"objective constraint",{},{"objc":"none"})
        
    def build_package(self,lower_bound,upper_bound):
        coef = self.model.solver.objective.get_linear_coefficients(self.model.solver.objective.variables)
        return BaseFBAPkg.build_constraint(self,"objc",lower_bound,upper_bound,coef,None)
