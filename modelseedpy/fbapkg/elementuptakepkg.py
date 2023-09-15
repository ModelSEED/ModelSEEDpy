# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class ElementUptakePkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "element uptake",
            {"elements": "string"},
            {"elements": "string"},
        )

    def build_package(
        self, element_limits, exception_compounds=[], exception_reactions=[]
    ):
        # Converting exception compounds list into exception reaction list
        self.parameters = {
            "element_limits": element_limits,
            "exception_compounds": exception_compounds,
            "exception_reactions": exception_reactions,
        }
        exchange_hash = self.modelutl.exchange_hash()
        for met in exception_compounds:
            if met in exchange_hash:
                exception_reactions.append(exchange_hash[met])
        # Now building or rebuilding constraints
        for element in element_limits:
            if element not in self.variables["elements"]:
                self.build_variable(element, element_limits[element])
        for element in element_limits:
            # This call will first remove existing constraints then build the new constraint
            self.build_constraint(element, exception_reactions)

    def build_variable(self, element, limit):
        return BaseFBAPkg.build_variable(
            self, "elements", 0, limit, "continuous", element
        )

    def build_constraint(self, element, exception_reactions):
        coef = {self.variables["elements"][element]: -1}
        rxnlist = self.modelutl.exchange_list()
        for reaction in rxnlist:
            if reaction not in exception_reactions:
                total = 0
                for metabolite in reaction.metabolites:
                    elements = metabolite.elements
                    if element in elements:
                        total += elements[element] * reaction.metabolites[metabolite]
                if total < 0:
                    coef[reaction.reverse_variable] = -1 * total
                elif total > 0:
                    coef[reaction.forward_variable] = total
        return BaseFBAPkg.build_constraint(self, "elements", 0, 0, coef, element)
