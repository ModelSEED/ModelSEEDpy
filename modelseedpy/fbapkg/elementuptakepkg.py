# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class ElementUptakePkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"element uptake",{"elements":"string"},{"elements":"string"})
        
    def build_package(self,element_limits):
        for element, limit in element_limits.items():
            if element not in self.variables["elements"]:
                self.build_variable(element, limit)
                self.build_constraint(element)
                   
    def build_variable(self,element,limit):
        return BaseFBAPkg.build_variable(self,"elements",0,limit,"continuous",element)
    
    def build_constraint(self,element):
        coef = {self.variables["elements"][element] : -1}
        for reaction in self.model.reactions:
            if reaction.id[0:3] == "EX_":
                total = 0
                for metabolite in reaction.metabolites:
                    elements = metabolite.elements
                    if element in elements:
                        total += elements[element]*reaction.metabolites[metabolite]
                if total < 0:
                    coef[reaction.reverse_variable] = -1*total
                elif total > 0:
                    coef[reaction.forward_variable] = total
        return BaseFBAPkg.build_constraint(self,"elements",0,0,coef,element)   