# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from optlang.symbolics import Zero

import re

#Base class for FBA packages
class SimpleThermoPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"simple thermo",{"potential":"metabolite", 'dgbinF': 'reaction', 'dgbinR':'reaction'},{"thermo":"reaction"})
        self.pkgmgr.addpkgs(["RevBinPkg"])
        
        
    def build_package(self,parameters):
        self.validate_parameters(parameters,[],{
            "filter":None,
            "min_potential":0,
            "max_potential":1000,
            "dgbin": False,
            'reduced_constraints': False
        })
        self.pkgmgr.getpkg("RevBinPkg").build_package(self.parameters["filter"])
        for metabolite in self.model.metabolites:
            self.build_variable(metabolite)
        for reaction in self.model.reactions:
            if re.search('^EX_', reaction.id) == None and re.search('^SK', reaction.id) == None and re.search('^DM_', reaction.id) == None:                       
                # determine the reaction flux range from the Delta_rG values 
                objective_coefficient = {}
                for metabolite in reaction.metabolites:
                    objective_coefficient[self.variables['potential'][metabolite.id]] = reaction.metabolites[metabolite]

                # define the maximum progression
                self.model.objective = self.model.problem.Objective(Zero, direction='max')
                self.model.objective.set_linear_coefficients(objective_coefficient)
                solution = self.model.optimize()
                max_value = solution.objective_value

                # define the minimum progression
                self.model.objective = self.model.problem.Objective(Zero,direction='min')
                self.model.objective.set_linear_coefficients(objective_coefficient)
                solution = self.model.optimize()
                min_value = solution.objective_value

                # determine the maximum flux
                reaction_flux_range = [min_value, max_value]
                max_flux = max(abs(flux) for flux in reaction_flux_range)

                # build constraints for the filtered reactions
                if self.parameters["filter"] == None or reaction.id in self.parameters["filter"]:                    
                    self.build_constraint(reaction, max_flux)
                   
                
    def build_variable(self,object):
        return BaseFBAPkg.build_variable(self,"potential",self.parameters["min_potential"],self.parameters["max_potential"],"continuous",object)

    
    def build_constraint(self,object, max_flux):
        # Gibbs: dg = Sum(st(i,j)*p(j))
        # 0 <= max_flux*revbin(i) - max_flux*dgbinR + max_flux*dgbinF + Sum(st(i,j)*p(j)) <= max_flux

        coef = {}
        for metabolite in object.metabolites:
            coef[self.variables["potential"][metabolite.id]] = object.metabolites[metabolite]
            
        if not self.parameters['reduced_constraints']:
            coef[self.pkgmgr.getpkg("RevBinPkg").variables["revbin"][object.id]] = max_flux
            if self.parameters['dgbin']:
                # build the dgbin variables
                BaseFBAPkg.build_variable(self,"dgbinF",0,1,"binary",object)
                BaseFBAPkg.build_variable(self,"dgbinR",0,1,"binary",object)

                # define the dgbin coefficients
                coef[self.variables['dgbinF'][object.id]] = max_flux
                coef[self.variables['dgbinR'][object.id]] = -max_flux      
            
            # build the constraint
            built_constraint = BaseFBAPkg.build_constraint(self,"thermo",0,max_flux,coef,object)
            
        else:
            built_constraint = None
                
        return built_constraint