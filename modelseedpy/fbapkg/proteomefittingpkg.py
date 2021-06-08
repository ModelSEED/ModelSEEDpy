# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.fbapkg.fluxfittingpkg import FluxFittingPkg
from modelseedpy.fbapkg.revbinpkg import RevBinPkg
from modelseedpy.fbapkg.totalfluxpkg import TotalFluxPkg

#Base class for FBA packages
class ProteomeFittingPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"proteome fitting",{"kapp":"reaction","kvfit":"reaction","kfit":"reaction"},{"vkapp":"reaction","kfitc":"reaction"})
        self.childpkgs["flux fitting"] = FluxFittingPkg(model)
        
    def build_package(self,parameters):
        self.validate_parameters(parameters,["rxn_proteome"],{
            "flux_values":{},
            "kcat_values":{},
            "prot_coef" : 0.1,
            "totalflux" : 1,
            "kcat_coef" : 0.333,
            "obj_kfit":1,
            "obj_kvfit":1,
            "obj_vfit":1,
            "set_objective":1
        })
        objvars = []
        #Adding flux fitting variables and constraints
        self.childpkgs["flux fitting"].build_package({
            "flux_values":self.parameters["flux_values"],
            "totalflux":self.parameters["totalflux"],
            "set_objective":0
        })
        for rxnid in self.childpkgs["flux fitting"].variables["vfit"]:
            objvars.append(objcoef["vfit"] * self.childpkgs["flux fitting"].variables["vfit"][rxnid] ** 2)
        #Adding proteome fitting variables and constraints
        for rxnid in rxn_proteome:
            if rxnid in self.model.reactions:
                rxnobj = self.model.reactions.get_by_id(rxnid)
                self.build_variable(rxnobj,"kapp")
                var = self.build_variable(rxnobj,"kvfit")
                objvars.append(objcoef["kvfit"] * var ** 2)
                const = self.build_constraint(rxnobj,"vkapp") 
        #Adding kcat fitting variables and constraints
        for rxnid in kcat_values:
            if rxnid in self.model.reactions:
                rxnobj = self.model.reactions.get_by_id(rxnid)
                var = self.build_variable(rxnobj,"kfit")
                const = self.build_constraint(rxnobj,"kfitc")
                objvars.append(objcoef["kfit"] * var ** 2)
        #Creating objective function
        if set_objective == 1:
            self.model.objective = self.model.problem.Objective(add(objvars), direction="min", sloppy=True)
                 
    def build_variable(self,object,type):
        if type == "kvfit":
            return BaseFBAPkg.build_variable(self,type,-1000,1000,"continuous",object)
        elif type == "kfit":
            return BaseFBAPkg.build_variable(self,type,-1000000,1000000,"continuous",object)
    
    def build_constraint(self,object,type):
        if type == "vkapp" and object.id in self.parameters["rxn_proteome"]:
            #kvfit(i) = kapp(i)*ProtCoef*Prot(i) - v(i)
            prot = self.parameters["rxn_proteome"][object.id]*self.parameters["prot_coef"]
            coef = {self.variables["kvfit"][object.id]:1,self.variables["kapp"][object.id]:-1*prot}
            if self.parameters["totalflux"] == 1:
                coef[self.childpkgs["flux fitting"].childpkgs["totalflux"].variables["tf"][object.id]] = 1
            else:
                coef[object.forward_variable] = 1
                coef[object.reverse_variable] = 1
            return BaseFBAPkg.build_constraint(self,type,0,0,coef,object)
        elif type == "kfitc" and object.id in self.parameters["kcat_values"]:
            #kfit(i) = kapp(i) - kcoef*kcat(i)
            rhs = -1*self.parameters["kcat_values"][object.id]*self.parameters["kcat_coef"]
            coef = {self.variables["kvfit"][object.id]:1,self.variables["kapp"][object.id]:-1}
            return BaseFBAPkg.build_constraint(self,type,rhs,rhs,coef,object)
