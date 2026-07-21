# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class ObjConstPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "objective constraint", {}, {"objc": "none"})

    def build_package(self, lower_bound, upper_bound):
        return self.build_constraint(lower_bound, upper_bound)

    def build_constraint(self, lower_bound, upper_bound):
        coef = self.model.solver.objective.get_linear_coefficients(
            self.model.solver.objective.variables
        )
        #Check if the constraint already exists and if so, just updating bounds in place
        for name in self.constraints["objc"]:
            constraint = self.constraints["objc"][name]
            existing_coef = constraint.get_linear_coefficients(
                constraint.variables
            )
            if coef == existing_coef:
                constraint.lb = lower_bound
                constraint.ub = upper_bound
                return constraint
        return BaseFBAPkg.build_constraint(
            self, "objc", lower_bound, upper_bound, coef, None
        )
