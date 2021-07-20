# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"simple thermo",{"potential":"metabolite"},{"thermo":"reaction"})
        self.pkgmgr.addpkgs(["RevBinPkg"])

    def build_package(self,parameters):
        self.validate_parameters(parameters,[],{
            "filter":None,
            "min_potential":0,
            "max_potential":1000,
            "add_dgbin_variables": False,
        })
        self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
        for metabolite in self.model.metabolites:
            self.build_variable(metabolite)
        for reaction in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if self.parameters["filter"] == None or reaction.id in self.parameters["filter"]:
                self.build_constraint(reaction)
    
    def build_variable(self,object):
        return BaseFBAPkg.build_variable(self,"potential",self.parameters["min_potential"],self.parameters["max_potential"],"continuous",object)

    def build_constraint(self,object):#Gibbs: dg = Sum(st(i,j)*p(j))
        #0 <= 1000*revbin(i) + 1000*dgbin,f - 1000*dgbin,r + Sum(st(i,j)*p(j)) <= 1000
        coef = {self.pkgmgr.getpkg("RevBinPkg").variables["revbin"][object.id] : 1000}
        for metabolite in object.metabolites:
            coef[self.variables["potential"][metabolite.id]] = object.metabolites[metabolite]
        return BaseFBAPkg.build_constraint(self,"thermo",0,1000,coef,object)