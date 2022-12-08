# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re  # !!! import is never used
from optlang.symbolics import Zero, add  # !!! add is never used
import json as _json  # !!! import is never used
from cobra.core import (
    Gene,
    Metabolite,
    Model,
    Reaction,
)  # !!! none of these imports are used
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.exceptions import FeasibilityError

logger = logging.getLogger(__name__)

# Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation"""

    pass


class BaseFBAPkg:
    """
    Base class for FBA packages
    """

    def __init__(
        self, model, name, variable_types={}, constraint_types={}, reaction_types={}
    ):
        self.model = model
        self.modelutl = MSModelUtil.get(model)
        self.name = name

        self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
        if self.pkgmgr is None:
            self.pkgmgr = MSPackageManager.get_pkg_mgr(model, 1)
        self.pkgmgr.addpkgobj(self)

        self.constraints, self.variables, self.parameters, self.new_reactions = (
            {},
            {},
            {},
            {},
        )
        self.variable_types = variable_types
        self.constraint_types = constraint_types

        for type in variable_types:
            self.variables[type] = dict()
        for type in constraint_types:
            self.constraints[type] = dict()

    def validate_parameters(self, params, required, defaults):
        for item in required:
            if item not in params:
                raise ValueError(f"Required argument {item} is missing!")
        self.parameters.update(defaults)  # we assign all defaults
        self.parameters.update(params)  # replace defaults with params

    def clear(self):
        cobra_objs = []
        for obj_type in self.variables:
            for cobra_obj in self.variables[obj_type]:
                cobra_objs.append(self.variables[obj_type][cobra_obj])
        for obj_type in self.constraints:
            for cobra_obj in self.constraints[obj_type]:
                cobra_objs.append(self.constraints[obj_type][cobra_obj])
        self.model.remove_cons_vars(cobra_objs)
        self.variables = {}
        self.constraints = {}

    def build_variable(
        self, obj_type, lower_bound, upper_bound, vartype, cobra_obj=None
    ):
        name = None
        if self.variable_types[obj_type] == "none":
            count = len(self.variables[obj_type])
            name = str(count + 1)
        elif self.variable_types[obj_type] == "string":
            name = cobra_obj
        else:
            name = cobra_obj.id
        if name not in self.variables[obj_type]:
            self.variables[obj_type][name] = self.model.problem.Variable(
                name + "_" + obj_type, lb=lower_bound, ub=upper_bound, type=vartype
            )
            self.model.add_cons_vars(self.variables[obj_type][name])
        return self.variables[obj_type][name]

    def build_constraint(
        self, obj_type, lower_bound, upper_bound, coef={}, cobra_obj=None
    ):
        name = None
        if self.constraint_types[obj_type] == "none":
            count = len(self.constraints[obj_type])
            name = str(count + 1)
        elif self.constraint_types[obj_type] == "string":
            name = cobra_obj
        else:
            name = cobra_obj.id
        if name in self.constraints[obj_type]:
            self.model.remove_cons_vars(self.constraints[obj_type][name])
        self.constraints[obj_type][name] = self.model.problem.Constraint(
            Zero, lb=lower_bound, ub=upper_bound, name=name + "_" + obj_type
        )
        self.model.add_cons_vars(self.constraints[obj_type][name])
        self.model.solver.update()
        if len(coef) > 0:
            self.constraints[obj_type][name].set_linear_coefficients(coef)
        self.model.solver.update()
        return self.constraints[obj_type][name]

    # Utility functions
    def print_lp(self, filename=None):
        if filename is None:
            filename = self.lp_filename
        if filename is not None:
            with open(filename + ".lp", "w") as out:
                complete_line = ""
                for line in str(self.model.solver).splitlines():
                    if ":" in line:
                        if complete_line != "":
                            out.write(complete_line)
                        complete_line = ""
                    else:
                        complete_line += line

    def all_variables(self):
        return self.pkgmgr.all_variables()

    def all_constraints(self):
        return self.pkgmgr.all_constraints()

    def add_variable_type(self, name, type):
        if name not in self.variables:
            self.variables[name] = dict()
        if name not in self.variable_types:
            self.variable_types[name] = type

    def add_constraint_type(self, name, type):
        if name not in self.constraints:
            self.constraints[name] = dict()
        if name not in self.constraint_types:
            self.constraint_types[name] = type
