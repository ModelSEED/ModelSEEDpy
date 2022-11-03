# -*- coding: utf-8 -*-

from __future__ import absolute_import
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
import logging

class ElementUptakePkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"element uptake",{"elements":"string"},{"elements":"string"})
        
    def build_package(self,element_limits):
        for element, limit in element_limits.items():
            if element not in self.variables["elements"]:
                self.build_variable("elements",0,limit,"continuous",element)
                self.build_constraint(element)
                   
    def build_constraint(self,element):
        coef = {self.variables["elements"][element]: -1}
        for exRXN in self.modelutl.exchange_list():
            totalNumAtoms = sum([met.elements[element]*exRXN.metabolites[met]
                                 for met in exRXN.metabolites if element in met.elements])
            if totalNumAtoms < 0:
                coef[exRXN.reverse_variable] = abs(totalNumAtoms)
            elif totalNumAtoms > 0:
                coef[exRXN.forward_variable] = totalNumAtoms
        return BaseFBAPkg.build_constraint(self,"elements",0,0,coef,element)
