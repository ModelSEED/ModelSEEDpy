# -*- coding: utf-8 -*-

from __future__ import absolute_import

from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class ChangeOptPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "change optimization",
            {"bgoal": "string"},
            {"bgoalc": "reaction"},
        )

    def build_package(self, target_values={}, build_objective=True):
        self.validate_parameters({}, [], {"target_values": target_values})
        # Creating objective but won't set objective unless build_objective is true
        goal_objective = self.model.problem.Objective(Zero, direction="max")
        obj_coef = dict()
        for rxnid in target_values:
            print(rxnid)
            if rxnid in self.model.reactions:
                print("FOUND!")
                rxn = self.model.reactions.get_by_id(rxnid)
                var = self.build_variable("bgoal", rxn)
                obj_coef[var] = target_values[rxnid]["objcoef"]
                const = self.build_constraint("bgoalc", rxn)
        # Now setting the goal objective if build_objective is true
        if build_objective:
            self.model.objective = goal_objective
            goal_objective.set_linear_coefficients(obj_coef)

    def build_variable(self, type, object):
        if type == "bgoal":
            goal = self.parameters["target_values"][object.id]
            return BaseFBAPkg.build_variable(
                self, type, 0, 1, "binary", object.id + goal["direction"]
            )

    def build_constraint(self, type, object):
        if type == "bgoalc":
            # For lower: goal - vi + 1000 - 1000 gi >= 0
            # For higher: vi - goal  + 1000 - 1000 gi >= 0
            lb = -1000
            goal = self.parameters["target_values"][object.id]
            var = self.variables["bgoal"][object.id + goal["direction"]]
            fluxvar = None
            if goal["flux"] < 0:
                fluxvar = object.reverse_variable
            else:
                fluxvar = object.forward_variable
            coef = {var: -1000}
            if goal["direction"] == "low":
                coef[fluxvar] = -1
                lb += -1 * abs(goal["flux"])
            else:
                coef[fluxvar] = 1
                lb += abs(goal["flux"])
            return BaseFBAPkg.build_constraint(self, type, lb, None, coef, object)
