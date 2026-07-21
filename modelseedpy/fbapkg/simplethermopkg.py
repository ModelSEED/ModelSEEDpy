# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from optlang.symbolics import Zero
import re

# Base class for FBA packages
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "simple thermo",
            {"potential": "metabolite", "dgp":"reaction"},
            {"thermo": "reaction", "dgpc": "reaction"},
        )
        self.pkgmgr.addpkgs(["RevBinPkg"])

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {
                "filter": None,
                "max_potential": 10000,
                "exclude_reactions": []
            },
        )
        self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
        for metabolite in self.model.metabolites:
            self.build_potential_variable(metabolite)
        for reaction in self.model.reactions:
            if reaction.id[:3] not in ["EX_", "SK_", "DM_"] and reaction.id not in self.parameters["exclude_reactions"]:
                self.build_dgp_variable(reaction)
                self.build_dgp_constraint(reaction)
                self.build_fluxcoupling_constraint(reaction)

    def build_potential_variable(self, object):
        return self.build_variable(
            "potential",
            -1*self.parameters["max_potential"],
            self.parameters["max_potential"],
            "continuous",
            object,
        )
    
    def build_dgp_variable(self, object):
        return self.build_variable(
            "dgp",
            -1*self.parameters["max_potential"],
            self.parameters["max_potential"],
            "continuous",
            object,
        )

    def build_dgp_constraint(self, object):
        # Gibbs: dg = Sum(st(i,j)*p(j))
        # 0 <= (-1) dgp(i) + Sum(st(i,j)*p(j)) <= 0
        coef = {}
        for metabolite in object.metabolites:
            coef[self.variables["potential"][metabolite.id]] = object.metabolites[
                metabolite
            ]
        coef[self.variables["dgp"][object.id]] = -1

        return self.build_constraint(
            "dgpc", 0, 0, coef, object
        )

    def build_fluxcoupling_constraint(self, object):
        # Gibbs: dg = Sum(st(i,j)*p(j)) 
        #  0 <= dgp(i) + max_energy_magnitude*revbin(i) <= max_energy_magnitude
        coef = {}
        coef[self.variables["dgp"][object.id]] = 1
        coef[self.pkgmgr.getpkg("RevBinPkg").variables["revbin"][object.id]] = self.parameters["max_potential"]+1
        # build the constraint
        return self.build_constraint(
            "thermo", 0, (self.parameters["max_potential"]+1), coef, object
        )
