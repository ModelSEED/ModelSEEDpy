# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg import BaseFBAPkg

#Base class for FBA packages
class MetaboFBAPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"metabo fba",{"met":"metabolite","pk":"string"},{"metc":"metabolite","pkc":"string"})
        self.pkgmgr.addpkgs(["SimpleThermoPkg"])

    def build_package(self,parameters):
        self.validate_parameters(parameters,["peaks"],{
            "set_objective":1,
        })
        self.pkgmgr.getpkg("SimpleThermoPkg").build_package()
        peak_hash = {}
        for peak_data in peaks:
            peak_hash[peak_data["id"]] = peak_data
            self.find_metabolites_matching_peak(peak_data)
            self.build_variable(peak_data,"pk")
            for met in peak_data["metabolites"]:
                self.build_variable(met,"met")
                self.build_constraint(met,"metc")
            self.build_constraint(peak_data,"pkc")
        if parameters["set_objective"] == 1:
            metabolite_objective = self.model.problem.Objective(
                Zero,
                direction="max")
            obj_coef = dict()
            for peak_id in self.variables["pk"]:
                if "wieght" in peak_hash[peak_id]:
                    obj_coef[self.variables["pk"][peak_id]] = peak_hash[peak_id]["wieght"]
                else:
                    obj_coef[self.variables["pk"][peak_id]] = 1
            self.model.objective = metabolite_objective
            metabolite_objective.set_linear_coefficients(obj_coef)
    
    def build_variable(self,object,type):
        if type == "met":
            return BaseFBAPkg.build_variable(self,type,0,1,"continuous",object)
        elif type == "pk":
            return BaseFBAPkg.build_variable(self,type,0,1,"continuous",object["id"])
    
    def build_constraint(self,object,type):
        #TODO: need to determine coefficients
        coef = {self.variables["met"][object.id]:1}
        if type == "metc":
            return BaseFBAPkg.build_constraint(self,"metc",0,0,coef,object)
        elif type == "pkc":
            return BaseFBAPkg.build_constraint(self,"pkc",0,0,coef,object["id"])
        
    def find_metabolites_matching_peak(self,data):
        #TODO: need to write this function
        pass
      