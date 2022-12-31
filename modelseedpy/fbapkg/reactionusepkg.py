# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging

logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! add is never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

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

    def build_package(self, rxn_filter=None, reversibility=False):
        for rxn in self.model.reactions:
            # Checking that reaction passes input filter if one is provided
            if rxn_filter == None:
                self.build_variable(rxn, "=")
                self.build_constraint(rxn, reversibility)
            elif rxn.id in rxn_filter:
                self.build_variable(rxn, rxn_filter[rxn.id])
                self.build_constraint(rxn, reversibility)

    def build_variable(self, cobra_obj, direction):
        variable = None
        if (
            (direction == ">" or direction == "=")
            and cobra_obj.upper_bound > 0
            and cobra_obj.id not in self.variables["fu"]
        ):
            variable = BaseFBAPkg.build_variable(self, "fu", 0, 1, "binary", cobra_obj)
        if (
            (direction == "<" or direction == "=")
            and cobra_obj.lower_bound < 0
            and cobra_obj.id not in self.variables["ru"]
        ):
            variable = BaseFBAPkg.build_variable(self, "ru", 0, 1, "binary", cobra_obj)
        return variable

    def build_constraint(self, cobra_obj, reversibility):
        constraint = None
        if (
            cobra_obj.id not in self.constraints["fu"]
            and cobra_obj.id in self.variables["fu"]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "fu",
                0,
                None,
                {
                    self.variables["fu"][cobra_obj.id]: 1000,
                    cobra_obj.forward_variable: -1,
                },
                cobra_obj,
            )
        if (
            cobra_obj.id not in self.constraints["ru"]
            and cobra_obj.id in self.variables["ru"]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "ru",
                0,
                None,
                {
                    self.variables["ru"][cobra_obj.id]: 1000,
                    cobra_obj.reverse_variable: -1,
                },
                cobra_obj,
            )
        if all(
            [
                reversibility,
                cobra_obj.id in self.variables["ru"],
                cobra_obj.id in self.variables["fu"],
            ]
        ):
            constraint = BaseFBAPkg.build_constraint(
                self,
                "urev",
                None,
                1,
                {
                    self.variables["ru"][cobra_obj.id]: 1,
                    self.variables["fu"][cobra_obj.id]: 1,
                },
                cobra_obj,
            )
        return constraint

    def build_exclusion_constraint(self, flux_values=None):
        flux_values = flux_values or FBAHelper.compute_flux_values_from_variables(
            self.model
        )
        count = len(self.constraints["exclusion"])
        solution_coef = {}
        solution_size = 0
        for rxnid, flux in flux_values.items():
            if flux > Zero:
                solution_size += 1
                solution_coef[self.variables["fu"][rxnid]] = 1
            elif flux < -1 * Zero:
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