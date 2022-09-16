# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class ReactionUsePkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "reaction use",
            {"fu": "reaction", "ru": "reaction"},
            {
                "fu": "reaction",
                "ru": "reaction",
                "exclusion": "none",
                "urev": "reaction",
            },
        )

    def build_package(self, filter=None, reversibility=0):
        for reaction in self.model.reactions:
            # Checking that reaction passes input filter if one is provided
            if filter == None:
                self.build_variable(reaction, "=")
                self.build_constraint(reaction, reversibility)
            elif reaction.id in filter:
                self.build_variable(reaction, filter[reaction.id])
                self.build_constraint(reaction, reversibility)

    def build_variable(self, object, direction):
        variable = None
        if (
            (direction == ">" or direction == "=")
            and object.upper_bound > 0
            and object.id not in self.variables["fu"]
        ):
            variable = BaseFBAPkg.build_variable(self, "fu", 0, 1, "binary", object)
        if (
            (direction == "<" or direction == "=")
            and object.lower_bound < 0
            and object.id not in self.variables["ru"]
        ):
            variable = BaseFBAPkg.build_variable(self, "ru", 0, 1, "binary", object)
        return variable

    def build_constraint(self, object, reversibility):
        constraint = None
        if (
            object.id not in self.constraints["fu"]
            and object.id in self.variables["fu"]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "fu",
                0,
                None,
                {self.variables["fu"][object.id]: 1000, object.forward_variable: -1},
                object,
            )
        if (
            object.id not in self.constraints["ru"]
            and object.id in self.variables["ru"]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "ru",
                0,
                None,
                {self.variables["ru"][object.id]: 1000, object.reverse_variable: -1},
                object,
            )
        if (
            reversibility == 1
            and object.id in self.variables["ru"]
            and object.id in self.variables["fu"]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "urev",
                None,
                1,
                {
                    self.variables["ru"][object.id]: 1,
                    self.variables["fu"][object.id]: 1,
                },
                object,
            )
        return constraint

    def build_exclusion_constraint(self, flux_values=None):
        if flux_values == None:
            flux_values = FBAHelper.compute_flux_values_from_variables(self.model)
        count = len(self.constraints["exclusion"])
        solution_coef = {}
        solution_size = 0
        for rxnid in flux_values:
            if flux_values[rxnid] > Zero:
                solution_size += 1
                solution_coef[self.variables["fu"][rxnid]] = 1
            elif flux_values[rxnid] < -1 * Zero:
                solution_size += 1
                solution_coef[self.variables["ru"][rxnid]] = 1
        if len(solution_coef) > 0:
            const_name = "exclusion." + str(count + 1)
            self.constraints["exclusion"][const_name] = self.model.problem.Constraint(
                Zero, lb=None, ub=(solution_size - 1), name=const_name
            )
            self.model.add_cons_vars(self.constraints["exclusion"][const_name])
            self.model.solver.update()
            self.constraints["exclusion"][const_name].set_linear_coefficients(
                solution_coef
            )
            return self.constraints["exclusion"][const_name]
        return None
