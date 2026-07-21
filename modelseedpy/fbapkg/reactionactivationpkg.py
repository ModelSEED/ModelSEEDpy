# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging

logger = logging.getLogger(__name__)
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

# Base class for FBA packages
class ReactionActivationPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "ReactionActivation",
            {"fra": "reaction", "rra": "reaction"},
            {
                "fra": "reaction",
                "rra": "reaction"
            }
        )

    def build_package(self, rxn_filter=None,max_value=0.001):
        self.pkgmgr.getpkg("RevBinPkg").build_package(filter=rxn_filter)
        for rxn in self.model.reactions:
            # Checking that reaction passes input filter if one is provided
            if rxn_filter == None:
                self.build_variable(rxn,max_value)
                self.build_constraint(rxn)
            elif rxn.id in rxn_filter:
                self.build_variable(rxn,max_value)
                self.build_constraint(rxn)

    def build_variable(self, cobra_obj,max_value):
        variable = BaseFBAPkg.build_variable(self, "fra", 0,max_value, "continuous", cobra_obj)
        variable = BaseFBAPkg.build_variable(self, "rra", 0,max_value, "continuous", cobra_obj)
        return variable

    def build_constraint(self, cobra_obj):
        constraint = None
        if cobra_obj.id not in self.constraints["fra"]:
            constraint = BaseFBAPkg.build_constraint(
                self,
                "fra",
                None,
                0,
                {
                    self.variables["fra"][cobra_obj.id]: 1,
                    cobra_obj.forward_variable: -1,
                },
                cobra_obj,
            )
        if cobra_obj.id not in self.constraints["rra"]:
            constraint = BaseFBAPkg.build_constraint(
                self,
                "rra",
                None,
                0,
                {
                    self.variables["rra"][cobra_obj.id]: 1,
                    cobra_obj.reverse_variable: -1
                },
                cobra_obj,
            )
        return constraint