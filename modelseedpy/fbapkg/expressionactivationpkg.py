# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging

logger = logging.getLogger(__name__)
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

# Base class for FBA packages
class ExpressionActivationPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "ExpressionActivation",
            {},
            {}
        )
        self.pkgmgr.addpkgs(["ReactionActivationPkg"])

    def build_package(self,on_hash,off_hash,on_coeff=None,off_coeff=None,other_coef=0.1,max_value=0.001,use_activation_constraints=False):
        activation_filter = {}
        for rxn in on_hash:
            activation_filter[rxn] = 1
        if use_activation_constraints:
            self.pkgmgr.getpkg("ReactionActivationPkg").build_package(rxn_filter=activation_filter,max_value=max_value)
        expression_objective = self.model.problem.Objective(0, direction="min")
        obj_coef = dict()
        for rxn in self.model.reactions:
            if rxn.id in on_hash and use_activation_constraints:
                coef = on_coeff
                if coef == None:
                    coef = on_hash[rxn.id]
                obj_coef[self.pkgmgr.getpkg("ReactionActivationPkg").variables["fra"][rxn.id]] = -1*coef
                obj_coef[self.pkgmgr.getpkg("ReactionActivationPkg").variables["rra"][rxn.id]] = -1*coef
            elif rxn.id in off_hash:
                coef = off_coeff
                if coef == None:
                    coef = off_hash[rxn.id]
                obj_coef[rxn.forward_variable] = coef
                obj_coef[rxn.reverse_variable] = coef
            elif rxn.id[0:3] == "bio" or rxn.id[0:3] == "EX_" or rxn.id[0:3] == "SK_" or rxn.id[0:3] == "DM_":
                pass
            else:
                obj_coef[rxn.forward_variable] = other_coef
                obj_coef[rxn.reverse_variable] = other_coef
        self.model.objective = expression_objective
        expression_objective.set_linear_coefficients(obj_coef)
        self.parameters["gfobj"] = self.model.objective