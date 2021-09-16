# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from scipy.constants import physical_constants, calorie, kilo, R
from numpy import log as ln
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

#Base class for FBA packages
class FullThermoPkg(BaseFBAPkg):
    @staticmethod
    def default_concentration():
        return {
            "cpd00067_c0":[0.0000001,0.0000001],#M H+ - equivalent to pHint = 7
            "cpd00007_c0":[1E-07,8.2E-06],#M O2 instracellular
            "cpd00011_c0":[1E-08,0.0014],#M CO2 instracellular
            "cpd00067_e0":[3.16228E-07,3.16228E-07],#M H+ - equivalent to pHext = 6.5
            "cpd00009_e0":[0.056,0.056],#Extracellular phosphate - overridden by media when media concentration is given
            "cpd00048_e0":[0.0030,0.0030],#Extracellular sulfate - overridden by media when media concentration is given
            "cpd00013_e0":[0.019,0.019],#Extracellular ammonia - overridden by media when media concentration is given
            "cpd00971_e0":[0.16,0.16],#Extracellular sodium - overridden by media when media concentration is given
            "cpd00205_e0":[0.022,0.022],#Extracellular potassium - overridden by media when media concentration is given
            "cpd10515_e0":[0.062,0.062],#Extracellular Fe2+ - overridden by media when media concentration is given
            "cpd00011_e0":[0.00010,0.00010],#Extracellular CO2 - overridden by media when media concentration is given
            "cpd00007_e0":[8.2E-06,8.2E-06],#Extracellular O2 - overridden by media when media concentration is given
            "cpd00027_e0":[0.020,0.020]#Extracellular glucose - overridden by media when media concentration is given
        }

    @staticmethod
    def default_deltaG_error():
        return {
            "cpd00067":0#KJ/mol - H delta G is based on pH and so has no error
        }
    
    @staticmethod   
    def default_compartment_potential():
        return {
            "e0":0,#Extracellular MUST be the zero reference for compartment electrochemical potential (so community modeling works)
            "c0":-160#mV = 0.33 (pHint - pHext) - 143.33 where pHint = 7 and pHext = 6.5
        }
    
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"full thermo",{"logconc":"metabolite","dgerr":"metabolite"},{"potc":"metabolite"})
        self.pkgmgr.addpkgs(["SimpleThermoPkg"])

    def build_package(self,parameters, verbose = True):
        self.validate_parameters(parameters,["modelseed_path"],{
            "default_max_conc":0.02,#M
            "default_min_conc":0.000001,#M
            "default_max_error":5,#KJ/mol
            "custom_concentrations":{},
            "custom_deltaG_error":{},
            "compartment_potential":{},
            "temperature":298,#K
            "filter":None,
            "infeasible_model": False,
            'dgbin':False
        })
        self.parameters["modelseed_api"] = FBAHelper.get_modelseed_db_api(self.parameters["modelseed_path"])
        simple_thermo_parameters = {
                "filter":self.parameters["filter"],
                "min_potential":-100000,#KJ/mol
                "max_potential":100000,#KJ/mol
                'dgbin':self.parameters['dgbin']
            }
        if self.parameters['infeasible_model']:
            simple_thermo_parameters['dgbin'] = True
        self.pkgmgr.getpkg("SimpleThermoPkg").build_package(simple_thermo_parameters)
            
        self.parameters["combined_custom_concentrations"] = FullThermoPkg.default_concentration()
        for cpd in self.parameters["custom_concentrations"]:
            self.parameters["combined_custom_concentrations"][cpd] = self.parameters["custom_concentrations"][cpd]
        self.parameters["combined_custom_deltaG_error"] = FullThermoPkg.default_deltaG_error()
        for cpd in self.parameters["custom_deltaG_error"]:
            self.parameters["combined_custom_deltaG_error"][cpd] = self.parameters["custom_deltaG_error"][cpd]
        self.parameters["combined_custom_comp_pot"] = FullThermoPkg.default_compartment_potential()
        for cmp in self.parameters["compartment_potential"]:
            self.parameters["combined_custom_comp_pot"][cmp] = self.parameters["compartment_potential"][cmp]
        msid_hash = {}
        for metabolite in self.model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid != None:
                if msid not in msid_hash:
                    msid_hash[msid] = {}
                msid_hash[msid][metabolite.id] = metabolite
            #Build concentration variable
            self.build_variable(metabolite,"logconc")
            #Build error variable
            self.build_variable(metabolite,"dgerr")
            #Build the potential constraint
            self.build_constraint(metabolite, verbose)

    def build_variable(self,object,type):
        msid = FBAHelper.modelseed_id_from_cobra_metabolite(object)
        if type == "logconc" and msid != "cpd00001":#Do not make a concentration variable for water
            lb = ln(self.parameters["default_min_conc"])
            ub = ln(self.parameters["default_max_conc"])
            if object.id in self.parameters["combined_custom_concentrations"]:
                lb = ln(self.parameters["combined_custom_concentrations"][object.id][0])
                ub = ln(self.parameters["combined_custom_concentrations"][object.id][1])
            return BaseFBAPkg.build_variable(self,"logconc",lb,ub,"continuous",object)
        elif type == "dgerr":
            ub = self.parameters["default_max_error"]
            if object.id in self.parameters["combined_custom_deltaG_error"]:
                ub = self.parameters["combined_custom_deltaG_error"][object.id]
            return BaseFBAPkg.build_variable(self,"dgerr",-1*ub,ub,"continuous",object)
    
    def build_constraint(self,object, verbose):
        #potential(i) (KJ/mol) = deltaG(i) (KJ/mol) + R * T(K) * lnconc(i) + charge(i) * compartment_potential
        if object.id not in self.pkgmgr.getpkg("SimpleThermoPkg").variables["potential"]:
            return None
        msid = FBAHelper.modelseed_id_from_cobra_metabolite(object)
        if msid == None:
            if verbose:
                print(object.id+" has no modelseed ID!")
            return None
        mscpd = self.parameters["modelseed_api"].get_seed_compound(msid)
        if mscpd is None:
            if verbose:
                print(object.id+" has modelseed ID "+msid+" but cannot be found in ModelSEED DB!")
            return None
        if mscpd.deltag == 10000000:
            if verbose:
                print(object.id+" has modelseed ID "+msid+" but does not have a valid deltaG!")
            return None
        Faraday = physical_constants['Faraday constant'][0]#C/mol
        compartment_potential = 0
        if object.compartment in self.parameters["combined_custom_comp_pot"]:
            compartment_potential = self.parameters["combined_custom_comp_pot"][object.compartment]
        constant = mscpd.deltag/calorie + object.charge*Faraday*compartment_potential/kilo/kilo
        coef = {
            self.pkgmgr.getpkg("SimpleThermoPkg").variables["potential"][object.id]:1,
            self.variables["dgerr"][object.id]:-1
        }
        if msid != "cpd00001":#Water concentration should not contribute to potential
            coef[self.variables["logconc"][object.id]] = -1*R/kilo*self.parameters["temperature"]
        return BaseFBAPkg.build_constraint(self,"potc",constant,constant,coef,object)