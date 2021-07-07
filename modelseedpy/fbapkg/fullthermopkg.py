# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from scipy.constants import physical_constants, calorie, kilo, R
from numpy import log as ln 

from modelseedpy.fbapkg.simplethermopkg import SimpleThermoPkg
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

#Base class for FBA packages
class FullThermoPkg(BaseFBAPkg):
    @staticmethod
    def default_concentration():
        return {
            "cpd00067_c0":[0.0000001,0.0000001],
            "cpd00067_e0":[3.16228E-07,3.16228E-07]
        }

    @staticmethod
    def default_deltaG_error():
        return {
            "cpd00067":0
        }
    
    @staticmethod   
    def default_compartment_potential():
        return {
            "e0":0,
            "c0":-127#mV value for E. coli at delta pH of 0.5
        }
    
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"full thermo",{"logconc":"metabolite","dgerr":"metabolite"},{"potc":"metabolite"})
        self.childpkgs["simple thermo"] = SimpleThermoPkg(model)

    def build_package(self,parameters):
        self.validate_parameters(parameters,["modelseed_path"],{
            "default_max_conc":0.02,#M
            "default_min_conc":0.000001,#M
            "default_max_error":100,#kcal
            "custom_concentrations":{},
            "custom_deltaG_error":{},
            "compartment_potential":{},
            "temperature":298,
            "filter":None
        })
        print("Faraday:"+str(physical_constants['Faraday constant'][0]))
        print("kilo:"+str(kilo))
        print("R:"+str(R))
        print("calorie:"+str(calorie))
        self.parameters["modelseed_api"] = FBAHelper.get_modelseed_db_api(self.parameters["modelseed_path"])
        self.childpkgs["simple thermo"].build_package({
            "filter":self.parameters["filter"],
            "min_potential":-100000,
            "max_potential":100000
        })
        self.parameters["combined_custom_concentrations"] = FullThermoPkg.default_concentration()
        for cpd in self.parameters["custom_concentrations"]:
            self.parameters["combined_custom_concentrations"][cpd] = self.parameters["custom_concentrations"][cpd]
        self.parameters["combined_custom_deltaG_error"] = FullThermoPkg.default_deltaG_error()
        for cpd in self.parameters["custom_deltaG_error"]:
            self.parameters["combined_custom_deltaG_error"][cpd] = self.parameters["custom_deltaG_error"][cpd]
        self.parameters["combined_custom_comp_pot"] = FullThermoPkg.default_compartment_potential()
        for cmp in self.parameters["compartment_potential"]:
            self.parameters["combined_custom_comp_pot"][cmp] = self.parameters["compartment_potential"][cmp]
        for metabolite in self.model.metabolites:
            #Build concentration variable
            self.build_variable(metabolite,"logconc")
            #Build error variable
            self.build_variable(metabolite,"dgerr")
            #Build the potential constraint
            self.build_constraint(metabolite)
    
    def build_variable(self,object,type):
        if type == "logconc":
            lb = ln(self.parameters["default_min_conc"])#Need to log this
            ub = ln(self.parameters["default_max_conc"])#Need to log this
            if object.id in self.parameters["combined_custom_concentrations"]:
                lb = ln(self.parameters["combined_custom_concentrations"][object.id][0])
                ub = ln(self.parameters["combined_custom_concentrations"][object.id][1])
            return BaseFBAPkg.build_variable(self,"logconc",lb,ub,"continuous",object)
        elif type == "dgerr":
            ub = self.parameters["default_max_error"]
            if object.id in self.parameters["combined_custom_deltaG_error"]:
                ub = self.parameters["combined_custom_deltaG_error"][object.id]
            return BaseFBAPkg.build_variable(self,"dgerr",-1*ub,ub,"continuous",object)
    
    def build_constraint(self,object):
        #potential(i) (KJ/mol) = deltaG(i) (KJ/mol) + R * T(K) * lnconc(i) + charge(i) * compartment_potential
        if object.id not in self.childpkgs["simple thermo"].variables["potential"]:
            return None
        msid = FBAHelper.modelseed_id_from_cobra_metabolite(object)
        if msid == None:
            print(object.id+" has no modelseed ID!")
            return None
        mscpd = self.parameters["modelseed_api"].get_seed_compound(msid)
        if mscpd is None:
            print(object.id+" has modelseed ID "+msid+" but cannot be found in ModelSEED DB!")
            return None
        if mscpd.deltag == 10000000:
            print(object.id+" has modelseed ID "+msid+" but does not have a valid deltaG!")
            return None
        Faraday = physical_constants['Faraday constant'][0]#C/mol
        compartment_potential = 0
        if object.compartment in self.parameters["combined_custom_comp_pot"]:
            compartment_potential = self.parameters["combined_custom_comp_pot"][object.compartment]
        constant = mscpd.deltag/calorie + object.charge*Faraday*compartment_potential*kilo
        coef = {
            self.childpkgs["simple thermo"].variables["potential"][object.id]:1,
            self.variables["logconc"][object.id]:-1*R*kilo*self.parameters["temperature"],
            self.variables["dgerr"][object.id]:-1
        }
        return BaseFBAPkg.build_constraint(self,"potc",constant,constant,coef,object)