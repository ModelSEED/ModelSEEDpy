# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging

logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! Zero is never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class ReferenceFluxPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self, model, "reference flux", {"refv": "reaction","preferr":"reaction","nrefer":"reaction"}, {"refvc": "reaction"}
        )

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {
                "coef": {},
                "fixed": False,
                "set_objective": True
                "quadratic": True,
                "objective_coef":{},
                "default_objective_coef": 1.0
            },
        )
        for rxnid in self.parameters["coef"]:
            if rxnid in self.model.reactions:
                rxnobj = self.model.reactions.get_by_id(rxnid)
                self.build_variables(rxnobj)
                self.build_constraints(rxnobj)
        objvars = {} 
        if self.parameters["set_objective"]:
            if self.parameters["fixed"]:
                logger.warning("In fixed mode, no objective will be created")
            else:
                for rxnid in self.variables["preferr"]:
                    objcoef = self.parameters["default_objective_coef"]
                    if rxnid in self.parameters["objective_coef"]:
                        objcoef = self.parameters["objective_coef"][rxnid]
                    if self.parameters["quadratic"]:
                        objvars[self.variables["preferr"][rxnid]] = objcoef
                    else:
                        objvars[self.variables["preferr"][rxnid]] = objcoef
                        objvars[self.variables["nrefer"][rxnid]] = objcoef
                self.model.objective = self.model.problem.Objective(
                    add(objvars), direction="min", sloppy=True
                )

    def build_variables(self, object,coef=None):
        BaseFBAPkg.build_variable(
            self, "refv", object.lower_bound, object.upper_bound, "continuous", object
        )
        lower_bound = -1000
        if not self.parameters["quadratic"]:
            BaseFBAPkg.build_variable(
                self, "nrefer", 0, 1000, "continuous", object
            )
            lower_bound = 0
        BaseFBAPkg.build_variable(
            self, "preferr", lower_bound, 1000, "continuous", object
        )

    def build_constraints(self, object):
        # Variable: preferr(i) - nrefer(i) = refv(i) - forward_flux(i) + reverse_flux(i)
        # Fixed: refv(i) = forward_flux(i) - reverse_flux(i)
        coef = {
            self.variables["vfit"][object.id]: 1,
            object.forward_variable: -1,
            object.reverse_variable: 1
        }
        if not self.parameters["fixed"]:
            coef[self.variables["preferr"][object.id]] = -1
            if not self.parameters["quadratic"]:
                coef[self.variables["nrefer"][object.id]] = 1
        return BaseFBAPkg.build_constraint(
            self, "refvc", 0, 0, coef, object
        )