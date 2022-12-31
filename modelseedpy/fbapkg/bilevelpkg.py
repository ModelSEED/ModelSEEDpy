# -*- coding: utf-8 -*-

from __future__ import absolute_import

import re
from optlang.symbolics import Zero, add  # !!! Neither import is used
from cobra.core import (
    Gene,
    Metabolite,
    Model,
    Reaction,
)  # !!! None of these imports are used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
import re

# Base class for FBA packages
class BilevelPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "reaction use",
            {"dualconst": "string", "dualub": "string", "duallb": "string"},
            {"dualvar": "string", "objective": "string", "dualbin": "string"},
        )

    def build_package(self, filter=None, binary_variable_count=0):
        self.validate_parameters(
            {}, [], {"binary_variable_count": binary_variable_count}
        )
        print("binary_variable_count:", binary_variable_count)
        varhash, coefficients, obj_coef = {}, {}, {}
        objective = self.model.solver.objective

        # Creating new objective coefficient and bound variables
        if self.parameters["binary_variable_count"] > 0:
            for reaction in self.model.reactions:
                var = self.build_variable("flxcmp", reaction, None)
        # Retrieving model data with componenent flux variables
        # Using the JSON calls because get_linear_coefficients is REALLY slow  #!!! get_linear_coefficients is still used?
        mdldata = self.model.solver.to_json()
        consthash = {}
        for const in mdldata["constraints"]:
            consthash[const["name"]] = const
        variables = list(self.model.solver.variables)
        objterms = objective.get_linear_coefficients(variables)

        # Adding binary variables and constraints which should not be included in dual formulation
        if self.parameters["binary_variable_count"] > 0:
            for reaction in self.model.reactions:
                self.build_variable("bflxcmp", reaction, None)

        # Now implementing dual variables and constraints
        for var in variables:
            varhash[var.name] = var
        for const in list(self.model.solver.constraints):
            var = self.build_variable("dualconst", const, obj_coef)
            if all(
                [
                    var,
                    const.name in consthash,
                    "expression" in consthash[const.name],
                    "args" in consthash[const.name]["expression"],
                ]
            ):
                for item in consthash[const.name]["expression"]["args"]:
                    if all(
                        [
                            "args" in item,
                            len(item["args"]) >= 2,
                            item["args"][1]["name"] in varhash,
                        ]
                    ):
                        var_name = varhash[item["args"][1]["name"]]
                        if var_name not in coefficients:
                            coefficients[var_name] = {}
                        coefficients[var_name][var] = item["args"][0]["value"]
        for var in variables:
            if var.type == "continuous":
                dvar = self.build_variable("duallb", var, obj_coef)   #!!! why is this repeated twice?
                if dvar != None:
                    if var not in coefficients:
                        coefficients[var] = {}
                    coefficients[var][dvar] = 1
                dvar = self.build_variable("dualub", var, obj_coef)
                if dvar != None:
                    if var not in coefficients:
                        coefficients[var] = {}
                    coefficients[var][dvar] = 1
                self.build_constraint("dualvar", var, objective, objterms, coefficients)
        self.build_constraint("objective", None, objective, objterms, obj_coef)

    def build_variable(self, obj_type, cobra_obj, obj_coef):
        if obj_type == "dualconst":
            lb = -1000000
            ub = -lb
            coef = 0
            if cobra_obj.lb == None:
                lb = 0
                coef = cobra_obj.ub
            if cobra_obj.ub == None:
                ub = 0
                coef = cobra_obj.lb
            var = BaseFBAPkg.build_variable(
                self, obj_type, lb, ub, "continuous", cobra_obj.name
            )
            obj_coef[var] = coef
            return var
        if (
            obj_type == "dualub"
        ):  # constrain this variable to zero when the binary variable is zero
            var = BaseFBAPkg.build_variable(
                self, obj_type, 0, 1000000, "continuous", cobra_obj.name
            )
            if re.search("(.+)_(fflxcmp\d+)$", cobra_obj.name) is not None:
                match = re.search("(.+)_(fflxcmp\d+)$", cobra_obj.name)
                bvar = self.variables[match[2]][match[1]]
                BaseFBAPkg.build_constraint(
                    self, "dualbin", None, 0, {var: 1, bvar: -1000000}, cobra_obj.name
                )
            obj_coef[var] = cobra_obj.ub
            return var
        if obj_type == "duallb":
            var = BaseFBAPkg.build_variable(
                self, obj_type, -1000000, 0, "continuous", cobra_obj.name
            )
            # if re.search('(.+)_(fflxcmp\d+)$', cobra_obj.name) is not None:
            # m = re.search('(.+)_(fflxcmp\d+)$', metabolite.id)
            # bvar = self.variables[m[2]][m[1]]
            # BaseFBAPkg.build_constraint(self,cobra_obj.name+"_lbdualbin",None,0,{var:-1,bvar:-1000000},cobra_obj)
            obj_coef[var] = cobra_obj.lb
            return var
        if obj_type == "flxcmp" and self.parameters["binary_variable_count"] > 0:
            denominator = 2 ** self.parameters["binary_variable_count"] - 1
            coefs = [{}, {}]
            for i in range(0, self.parameters["binary_variable_count"]):
                value = 2**i
                if cobra_obj.lower_bound < 0:
                    self.add_variable_type("rflxcmp" + str(i), "reaction")
                    var = BaseFBAPkg.build_variable(
                        self,
                        "rflxcmp" + str(i),
                        0,
                        -1 * value * cobra_obj.lower_bound / denominator,
                        "continuous",
                        cobra_obj,
                    )
                    coefs[0][var] = -1
                if cobra_obj.upper_bound > 0:
                    self.add_variable_type("fflxcmp" + str(i), "reaction")
                    var = BaseFBAPkg.build_variable(
                        self,
                        "fflxcmp" + str(i),
                        0,
                        value * cobra_obj.upper_bound / denominator,
                        "continuous",
                        cobra_obj,
                    )
                    coefs[1][var] = -1
            if cobra_obj.lower_bound < 0:
                # flux - flux_comp_0 - flux_comp_n = 0 - restriction of reverse fluxes by component fluxes
                self.add_constraint_type("rflxcmpc", "reaction")
                coefs[0][cobra_obj.reverse_variable] = 1
                BaseFBAPkg.build_constraint(self, "rflxcmpc", 0, 0, coefs[0], cobra_obj)
            if cobra_obj.upper_bound > 0:
                # flux - flux_comp_0 - flux_comp_n = 0 - restriction of forward fluxes by component fluxes
                self.add_constraint_type("fflxcmpc", "reaction")
                coefs[1][cobra_obj.forward_variable] = 1
                BaseFBAPkg.build_constraint(self, "fflxcmpc", 0, 0, coefs[1], cobra_obj)
            return None
        if obj_type == "bflxcmp" and self.parameters["binary_variable_count"] > 0:
            for i in range(0, self.parameters["binary_variable_count"]):
                if cobra_obj.lower_bound < 0:
                    self.add_variable_type("brflxcmp" + str(i), "reaction")
                    var = BaseFBAPkg.build_variable(
                        self, "brflxcmp" + str(i), 0, 1, "binary", cobra_obj
                    )
                    othervar = self.variables["rflxcmp" + str(i)][cobra_obj.id]
                    self.add_constraint_type("brflxcmpc" + str(i), "reaction")
                    BaseFBAPkg.build_constraint(
                        self,
                        "brflxcmpc" + str(i),
                        None,
                        0,
                        {othervar: 1, var: -1000},
                        cobra_obj,
                    )
                if cobra_obj.upper_bound > 0:
                    self.add_variable_type("bfflxcmp" + str(i), "reaction")
                    var = BaseFBAPkg.build_variable(
                        self, "bfflxcmp" + str(i), 0, 1, "binary", cobra_obj
                    )
                    othervar = self.variables["fflxcmp" + str(i)][cobra_obj.id]
                    self.add_constraint_type("bfflxcmpc" + str(i), "reaction")
                    BaseFBAPkg.build_constraint(
                        self,
                        "bfflxcmpc" + str(i),
                        None,
                        0,
                        {othervar: 1, var: -1000},
                        cobra_obj,
                    )
            return None

    def build_constraint(self, obj_type, cobra_obj, objective, objterms, coefficients):
        if obj_type == "dualvar":
            coef = {}
            lb = ub = 0
            objsign = 1
            if objective.direction == "min":
                objsign = -1
            if cobra_obj in objterms:
                lb = ub = objterms[cobra_obj]
            if cobra_obj in coefficients:
                for var in coefficients[cobra_obj]:
                    coef[var] = coefficients[cobra_obj][var]
            if cobra_obj.lb == 0:
                ub = None
            elif cobra_obj.ub == 0:
                lb = None
            return BaseFBAPkg.build_constraint(
                self, obj_type, lb, ub, coef, cobra_obj.name
            )
        elif obj_type == "objective":
            coef = {}
            objsign = 1
            if objective.direction == "min":
                objsign = -1
            for var in objterms:
                coef[var] = objsign * objterms[var]
            for dvar in coefficients:
                coef[dvar] = -1 * coefficients[dvar]
            return BaseFBAPkg.build_constraint(
                self, obj_type, 0, 0, coef, "dualobjconst"
            )
