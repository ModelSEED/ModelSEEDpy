# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging

logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! Zero is never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class FoldChangeFittingPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self, model, "fold change fitting", {}, {}
        )

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            ["msexp"],
            {
                "quadratic": True,
                "condition_models": {}
            },
        )
        ref_flux_pkg_args = { 
            "coef": None,
            "fixed": True,
            "set_objective": False,
            "quadratic": self.parameters["quadratic"],
            "objective_coef":{},
            "default_objective_coef": 1.0
        }
        rxn_expression = self
        if isinstance(self.object, MSGenome):
            rxn_expression = self.build_reaction_expression(self.model.model)

        for feature in self.parameters["msexp"].features:
            if feature.id in self.model.reactions:
                ref_flux_pkg_args["objective_coef"][feature.id] = 1
                ref_flux_pkg_args["coef"][feature.id] = 1
        #Adding the reference flux variables to the main model
        self.pkgmgr.getpkg("ReferenceFluxPkg").build_package(ref_flux_pkg_args)
        ref_flux_pkg_args["fixed"] = False
        ref_flux_pkg_args["set_objective"] = True
        other_models = []
        for condition in self.parameters["condition_models"]:
            if condition not in self.parameters["msexp"].conditions:
                raise ValueError("Condition "+condition+" not found in expression data")
            mdl = self.parameters["condition_models"][condition]
            #Adding the reference flux variables to each model
            mdl.pkgmgr.getpkg("ReferenceFluxPkg").build_package(ref_flux_pkg_args)
            other_models.append(mdl)
        shared_var_pkgs = {"ReferenceFluxPkg": ["refv"]}
        #Replicating the other models into the main model
        self.pkgmgr.getpkg("ProblemReplicationPkg").build_package(
            {
                "models": other_models,
                "shared_variable_packages": shared_var_pkgs
            }
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