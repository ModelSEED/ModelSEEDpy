# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
import json as _json
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.exceptions import FeasibilityError


# Base class for FBA packages
class TotalFluxPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self, model, "totalflux", {"tf": "reaction"}, {"tf": "reaction"}
        )

    def build_package(self, reaction_filter=None, upper_bound=100):
        for reaction in self.model.reactions:
            # Checking that variable has not yet been created
            if reaction.id not in self.variables["tf"]:
                # Checking that reaction passes input filter if one is provided
                if reaction_filter == None or reaction.id in reaction_filter:
                    self.variables["tf"][reaction.id] = self.model.problem.Variable(
                        reaction.id + "_tf", lb=0, ub=upper_bound
                    )
                    self.model.add_cons_vars(self.variables["tf"][reaction.id])
                    self.constraints["tf"][reaction.id] = self.model.problem.Constraint(
                        reaction.forward_variable + reaction.reverse_variable - self.variables["tf"][reaction.id],
                        lb=0,
                        ub=0,
                        name=reaction.id + "_tf",
                    )
                    self.model.add_cons_vars(self.constraints["tf"][reaction.id])
