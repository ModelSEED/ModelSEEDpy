# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

import re

#Base class for FBA packages
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"simple thermo",{"potential":"metabolite"},{"thermo":"reaction"})
        self.pkgmgr.addpkgs(["RevBinPkg"])

    def build_package(self,parameters):
        reaction_filter = []
        for reaction in self.model.reactions:
            if re.search('EX', reaction.id):
                reaction_filter.append(reaction.id)
        
        self.validate_parameters(parameters,[],{
            "filter":reaction_filter,
            "min_potential":0,
            "max_potential":1000,
            "dgbin": False,
            'reduced_constraints': False
        })
        self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
        for metabolite in self.model.metabolites:
            self.build_variable(metabolite)
        for reaction in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if self.parameters["filter"] == None or reaction.id not in self.parameters["filter"]:                    
                self.build_constraint(reaction)
                
    def build_variable(self,object):
        return BaseFBAPkg.build_variable(self,"potential",self.parameters["min_potential"],self.parameters["max_potential"],"continuous",object)

    def build_constraint(self,object):
        # Gibbs: dg = Sum(st(i,j)*p(j))
        # 0 <= 1000*revbin(i) - 1000*dgbinR + 1000*dgbinF + Sum(st(i,j)*p(j)) <= 1000
        # -1000 <= -1000*revbin + Sum(st(i,j)*p(j)) <= 0
                              
        for metabolite in object.metabolites:
            coef = {self.variables["potential"][metabolite.id]:object.metabolites[metabolite]}
            
        if not self.parameters['reduced_constraints']:
            coef[self.pkgmgr.getpkg("RevBinPkg").variables["revbin"][object.id]] = 1000
            if self.parameters['dgbin']:
                # build the dgbin variables
                BaseFBAPkg.build_variable(self,"dgbinF",0,1,"binary",object)
                BaseFBAPkg.build_variable(self,"dgbinR",0,1,"binary",object)

                # define the dgbin coefficients
                coef[self.variables['dgbinF'][object.id]] = 1000
                coef[self.variables['dgbinR'][object.id]] = -1000            
            
            # build the constraint
            built_constraint = BaseFBAPkg.build_constraint(self,"thermo",0,1000,coef,object)
            
        else:
            built_constraint = None
                
        return built_constraint