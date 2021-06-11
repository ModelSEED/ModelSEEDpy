# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.fbapkg.revbinpkg import RevBinPkg

#Base class for FBA packages
class FullThermoPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"full thermo",{"logconc":"metabolite"},{"potentialc":"metabolite"})
        self.childpkgs["simple thermo"] = SimpleThermoPkg(model)

    def build_package(self,parameters):
        self.validate_parameters(parameters,[],{
            "default_max_conc":20,#mM
            "default_min_conc":0.001,#mM
            "custom_concentration_constraints":{"cpd000027_e0":["10","10"]},
            "compartment_potential":{"c0":,}
        })
        #Extend the custom concentration for H+, CO2, O2
        self.childpkgs["simple thermo"].build_package(filter)
        for metabolite in self.model.metabolites:
            #Build concentration variable
            self.build_variable(metabolite)
            #Build the potential constraint
            self.build_constraint(reaction)
    
    def build_variable(self,object):
        lb = self.parameters["default_max_conc"]
        if object.id in self.parameters["custom_concentration_constraints"]:
            lb = self.parameters["custom_concentration_constraints"][object.id]
        ub = ?
        return BaseFBAPkg.build_variable(self,"logconc",lb,ub,"continuous",object)
    
    def build_constraint(self,object):
        #potential(i) (KJ/mol) = deltaG(i) (KJ/mol) + R * T(K) * lnconc(i) + charge(i) * compartment_potentil
        constant = ?
        coef = {?}
        return BaseFBAPkg.build_constraint(self,"thermo",constant,constant,coef,object)