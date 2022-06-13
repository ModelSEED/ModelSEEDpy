# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
logger = logging.getLogger("modelseedpy")

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

#Base class for FBA packages
class DrainFluxPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"drain flux",{},{"drain"})
        self.update_drain_fluxes()
        
    def build_package(self,parameters):
        self.validate_parameters(parameters,[],{
            "add_all_intracellular_drains":False,
            "default_uptake":0,
            "default_excretion":100,
            "drain_compounds":{},
            "set_minimal_drain_objective":False,
            "update_drain_fluxes":False
        })
        if self.parameters["update_drain_fluxes"]:
            self.update_drain_fluxes()
        if self.parameters["add_all_intracellular_drains"]:
            for cpd in self.model.metabolites:
                self.add_drain_reaction(cpd,self.parameters["default_uptake"],self.parameters["default_excretion"])
        else:
            for cpd in self.parameters["drain_compounds"]:
                if cpd in self.model.metabolites:
                    cpdobj = self.model.metabolites.get_by_id(cpd)
                    self.add_drain_reaction(cpdobj,self.parameters["drain_compounds"][cpd]["uptake"],self.parameters["drain_compounds"][cpd]["excretion"])

    def add_drain_reaction(self,cpd,uptake,excretion):
        namespace = "DN"
        if cpd.id.split("_")[-1][0:1] == "e" or cpd.compartment == "e":
            namespace = "EX"
        if cpd.id not in self.new_reactions["Existing_"+namespace] and cpd.id not in self.new_reactions["New_"+namespace]:
            rxn = FBAHelper.add_drain_from_metabolite_id(self.model, cpd.id,uptake,excretion)
            if rxn.id in self.model.reactions:
                rxn = self.model.reactions.get_by_id(rxn.id)
            self.new_reactions["New_"+namespace][cpd.id] = rxn    
    
    def update_drain_fluxes(self):
        previous_reactions = self.new_reactions
        self.new_reactions = {"Existing_EX":{},"Existing_DN":{},"New_EX":{},"New_DN":{}}
        for rxn in self.reactions:
            if FBAHelper.is_ex(rxn) and len(rxn.metabolites) == 1:
                cpd = rxn.metabolites.keys()[0]
                if cpd.id.split("_")[-1][0:1] == "e" or cpd.compartment == "e" or rxn.id[0:3] == "EX_":
                    if "Existing_EX" not in previous_reactions or cpd.id in previous_reactions["Existing_EX"]:
                        self.new_reactions["Existing_EX"][cpd.id] = rxn
                    else:
                        self.new_reactions["New_EX"][cpd.id] = rxn
                else:
                    if "Existing_DN" not in previous_reactions or cpd.id in previous_reactions["Existing_DN"]:
                        self.new_reactions["Existing_DN"][cpd.id] = rxn
                    else:
                        self.new_reactions["New_DN"][cpd.id] = rxn
        logger.info("Updated drain fluxes - Exist EX:"+
                    ";".join(self.new_reactions["Existing_EX"].keys())+"|New EX:"+
                    ";".join(self.new_reactions["New_EX"].keys())+"|Existing DN:"+
                    ";".join(self.new_reactions["Existing_DN"].keys())+"|New DN:"+
                    ";".join(self.new_reactions["New_DN"].keys()))
                
    