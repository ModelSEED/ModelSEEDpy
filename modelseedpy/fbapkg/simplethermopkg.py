# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from optlang.symbolics import Zero

# Base class for FBA packages
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "simple thermo",
            {"potential": "metabolite", "dgbinF": "reaction", "dgbinR": "reaction"},
            {"thermo": "reaction"},
        )
        self.pkgmgr.addpkgs(["RevBinPkg"])

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {
                "filter": None,
                "min_potential": 0,
                "max_potential": 1000,
                "dgbin": False,
                "reduced_constraints": False,
            },
        )
        self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
        for metabolite in self.model.metabolites:
            self.build_variable(metabolite)
        for reaction in self.model.reactions:
            if reaction.id[:3] not in ["EX_", "SK_", "DM_"]:
                # determine the range of Delta_rG values
                objective_coefficient = {}
                for metabolite in reaction.metabolites:
                    objective_coefficient[
                        self.variables["potential"][metabolite.id]
                    ] = reaction.metabolites[metabolite]

                # define the maximum progression
                self.model.objective = self.model.problem.Objective(
                    Zero, direction="max"
                )
                self.model.objective.set_linear_coefficients(objective_coefficient)
                solution = self.model.optimize()
                max_value = solution.objective_value

                # define the minimum progression
                self.model.objective = self.model.problem.Objective(
                    Zero, direction="min"
                )
                self.model.objective.set_linear_coefficients(objective_coefficient)
                solution = self.model.optimize()
                min_value = solution.objective_value

                # determine the maximum Delta_rG magnitude
                reaction_energy_range = [min_value, max_value]
                max_energy_magnitude = max(
                    abs(energy) for energy in reaction_energy_range
                )

                # build constraints for the filtered reactions
                if (
                    self.parameters["filter"] == None
                    or reaction.id in self.parameters["filter"]
                ):
                    self.build_constraint(reaction, max_energy_magnitude)

        if self.parameters["dgbin"]:
            # define the model objective as the sum of the dgbin variables
            self.optimize_dgbin()

    def build_variable(self, object):
        return BaseFBAPkg.build_variable(
            self,
            "potential",
            self.parameters["min_potential"],
            self.parameters["max_potential"],
            "continuous",
            object,
        )

    def build_constraint(self, object, max_energy_magnitude):
        # Gibbs: dg = Sum(st(i,j)*p(j))
        # 0 <= max_energy_magnitude*revbin(i) - max_energy_magnitude*dgbinR + max_energy_magnitude*dgbinF + Sum(st(i,j)*p(j)) <= max_energy_magnitude

        coef = {}
        for metabolite in object.metabolites:
            coef[self.variables["potential"][metabolite.id]] = object.metabolites[
                metabolite
            ]

        if not self.parameters["reduced_constraints"]:
            coef[
                self.pkgmgr.getpkg("RevBinPkg").variables["revbin"][object.id]
            ] = max_energy_magnitude
            if self.parameters["dgbin"]:
                # build the dgbin variables
                BaseFBAPkg.build_variable(self, "dgbinF", 0, 1, "binary", object)
                BaseFBAPkg.build_variable(self, "dgbinR", 0, 1, "binary", object)

                # define the dgbin coefficients
                coef[self.variables["dgbinF"][object.id]] = max_energy_magnitude
                coef[self.variables["dgbinR"][object.id]] = -max_energy_magnitude

            # build the constraint
            built_constraint = BaseFBAPkg.build_constraint(
                self, "thermo", 0, max_energy_magnitude, coef, object
            )
        else:
            built_constraint = None

        return built_constraint

    def optimize_dgbin(self):
        # create the sum of dgbin variables
        dgbin_sum_coef = {}
        for reaction in self.variables["dgbinF"]:
            print(f"{self.model.solver.status} status for {reaction}")
            try:
                dgbin_sum_coef[self.variables["dgbinF"][reaction].primal] = 1
            except:
                print("--> ERROR: The simulation lack a solution.")
        for reaction in self.variables["dgbinR"]:
            print(f"{self.model.solver.status} status for {reaction}")
            try:
                dgbin_sum_coef[self.variables["dgbinR"][reaction].primal] = 1
            except:
                print("--> ERROR: The simulation lack a solution.")

        # set the dgbin sum as the model objective
        self.model.objective = self.model.problem.Objective(Zero, direction="max")
        self.model.objective.set_linear_coefficients(dgbin_sum_coef)
