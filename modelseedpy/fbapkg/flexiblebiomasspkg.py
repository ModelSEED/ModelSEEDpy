# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from cobra import Model, Reaction, Metabolite
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

classes = {
    "rna":{"cpd00052":-1,"cpd00038":-1,"cpd00002":-1,"cpd00062":-1},
    "dna":{"cpd00115":-1,"cpd00356":-1,"cpd00241":-1,"cpd00357":-1},
    "protein":{"cpd00132":-1,"cpd00023":-1,"cpd00053":-1,"cpd00054":-1,"cpd00033":-1,"cpd00039":-1,"cpd00119":-1,"cpd00051":-1,"cpd00041":-1,"cpd00107":-1,"cpd00129":-1,"cpd00322":-1,"cpd00069":-1,"cpd00065":-1,"cpd00060":-1,"cpd00084":-1,"cpd00035":-1,"cpd00161":-1,"cpd00156":-1,"cpd00066":-1},
    "energy":{"cpd00008":1}
}

#Base class for FBA packages
class FlexibleBiomassPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"flexible biomass",{},{"flxbio":"reaction","fflxcpd":"metabolite","rflxcpd":"metabolite","fflxcls":"reaction","rflxcls":"reaction"})
        
    def build_package(self,parameters):
        self.validate_parameters(parameters,["bio_rxn_id"],{
            "flex_coefficient":0.75,
            "use_rna_class":[-0.75,0.75],
            "use_dna_class":[-0.75,0.75],
            "use_protein_class":[-0.75,0.75],
            "use_energy_class":[-0.1,0.1],
        })
        if self.parameters["bio_rxn_id"] not in self.model.reactions:
            raise ValueError(self.parameters["bio_rxn_id"]+" not found in model!")
        self.parameters["bio_rxn"] = self.model.reactions.get_by_id(self.parameters["bio_rxn_id"])
        newrxns = []
        class_coef = {"rna":{},"dna":{},"protein":{},"energy":{}}
        refcpd = {"cpd00001":None,"cpd00009":None,"cpd00012":None,"cpd00067":None,"cpd00002":None}
        for metabolite in self.model.metabolites:
            for msid in refcpd:
                if FBAHelper.modelseed_id_from_cobra_metabolite(metabolite) == msid:
                    refcpd[msid] = metabolite
        for metabolite in self.parameters["bio_rxn"].metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid != "cpd11416":
                met_class = "none"
                if msid != None:
                    for curr_class in classes:
                        if msid in classes[curr_class]:
                            met_class = curr_class
                            class_coef[curr_class][msid] = metabolite
                if (met_class == "none" or self.class_complete(class_coef,met_class) == 0 or self.parameters["use_"+met_class+"_class"] == None) and msid not in refcpd:
                    drain_reaction = FBAHelper.add_drain_from_metabolite_id(self.model,metabolite.id,1000,1000,"FLEX_")
                    if drain_reaction.id not in self.new_reactions:
                        self.new_reactions[drain_reaction.id] = drain_reaction
                        self.model.add_reactions([drain_reaction])
                    self.build_constraint(metabolite,"flxcpd")
        for met_class in class_coef:
            add = 0
            total_coef = 0
            object_stoichiometry = {}
            for msid in class_coef[met_class]:
                if met_class == "rna" and msid == "cpd00002" and "cpd00008" in class_coef["energy"]:
                    object_stoichiometry[class_coef[met_class][msid]] = self.parameters["bio_rxn"].metabolites[class_coef[met_class][msid]] + self.parameters["bio_rxn"].metabolites[class_coef["energy"]["cpd00008"]]        
                else:
                    object_stoichiometry[class_coef[met_class][msid]] = self.parameters["bio_rxn"].metabolites[class_coef[met_class][msid]]
                total_coef += abs(object_stoichiometry[class_coef[met_class][msid]])
            if (met_class == "rna" or met_class == "dna") and refcpd["cpd00012"] != None and refcpd["cpd00001"] != None:
                add = 1
                object_stoichiometry[refcpd["cpd00012"]] = total_coef
                object_stoichiometry[refcpd["cpd00001"]] = total_coef
            if met_class == "protein" and refcpd["cpd00001"] != None:
                add = 1
                object_stoichiometry[refcpd["cpd00001"]] = total_coef
            if met_class == "energy" and refcpd["cpd00001"] != None and refcpd["cpd00002"] != None and refcpd["cpd00067"] != None and refcpd["cpd00009"] != None:
                add = 1
                object_stoichiometry[refcpd["cpd00001"]] = -1*total_coef
                object_stoichiometry[refcpd["cpd00002"]] = -1*total_coef
                object_stoichiometry[refcpd["cpd00009"]] = total_coef
                object_stoichiometry[refcpd["cpd00067"]] = total_coef
            if add == 1:
                if met_class+"_flex" not in self.new_reactions:
                    self.new_reactions[met_class+"_flex"] = Reaction(id=met_class+"_flex", 
                                              name= met_class+"_flex", 
                                              lower_bound=-1000, 
                                              upper_bound=1000)
                    self.new_reactions[met_class+"_flex"].add_metabolites(object_stoichiometry)
                    self.new_reactions[met_class+"_flex"].annotation["sbo"] = 'SBO:0000627'
                    self.model.add_reactions([self.new_reactions[met_class+"_flex"]])
                self.build_constraint(self.new_reactions[met_class+"_flex"],"flxcls")
        self.build_constraint(self.parameters["bio_rxn"],"flxbio")
                   
    def build_variable(self,object,type):
        pass
        
    def build_constraint(self,object,type):
        element_mass = FBAHelper.elemental_mass()
        if type == "flxbio":
            #Sum(MW*(vdrn,for-vdrn,ref)) + Sum(massdiff*(vrxn,for-vrxn,ref)) = 0
            coef = {}
            for metabolite in self.parameters["bio_rxn"].metabolites:
                if "FLEX_"+metabolite.id in self.model.reactions:
                    mw = FBAHelper.metabolite_mw(metabolite)
                    sign = -1
                    if self.parameters["bio_rxn"].metabolites[metabolite] > 0:
                        sign = 1
                    coef[self.model.reactions.get_by_id("FLEX_"+metabolite.id).forward_variable] = sign*mw
                    coef[self.model.reactions.get_by_id("FLEX_"+metabolite.id).reverse_variable] = -1*sign*mw
            for met_class in classes:
                if met_class+"_flex" in self.model.reactions:
                    massdiff = 0
                    rxn = self.model.reactions.get_by_id(met_class+"_flex")
                    for metabolite in rxn.metabolites:
                        mw = FBAHelper.metabolite_mw(metabolite)
                        massdiff += rxn.metabolites[metabolite]*mw
                    if abs(massdiff) > 0.00001:
                        coef[rxn.forward_variable] = massdiff
                        coef[rxn.reverse_variable] = -1*massdiff
            return BaseFBAPkg.build_constraint(self,type,0,0,coef,object)
        elif type == "flxcpd":
            #0.75 * abs(bio_coef) * vbio - vdrn,for >= 0
            #0.75 * abs(bio_coef) * vbio - vdrn,rev >= 0
            coef = self.parameters["flex_coefficient"]*abs(self.parameters["bio_rxn"].metabolites[object])
            if coef > 0.75:
                coef = 0.75
            BaseFBAPkg.build_constraint(self,"f"+type,0,None,{
                self.parameters["bio_rxn"].forward_variable:coef,
                self.model.reactions.get_by_id("FLEX_"+object.id).forward_variable:-1
            },object)
            return BaseFBAPkg.build_constraint(self,"r"+type,0,None,{
                self.parameters["bio_rxn"].forward_variable:coef,
                self.model.reactions.get_by_id("FLEX_"+object.id).reverse_variable:-1
            },object)
        elif type == "flxcls" and object.id[0:-5] != None:
            #0.75 * vbio - vrxn,for >= 0
            #0.75 * vbio - vrxn,rev >= 0
            #First deal with the situation where the flux is locked into a particular value relative to biomass
            const = None
            if self.parameters["use_"+object.id[0:-5]+"_class"][0] == self.parameters["use_"+object.id[0:-5]+"_class"][1]:
                #If the value is positive, lock in the forward variable and set the reverse to zero
                if self.parameters["use_"+object.id[0:-5]+"_class"][0] > 0:
                    const = BaseFBAPkg.build_constraint(self,"f"+type,0,0,{
                        self.parameters["bio_rxn"].forward_variable:self.parameters["use_"+object.id[0:-5]+"_class"][1],
                        object.forward_variable:-1
                    },object)
                    object.lower_bound = 0
                #If the value is negative, lock in the reverse variable and set the forward to zero
                elif self.parameters["use_"+object.id[0:-5]+"_class"][0] < 0:
                    const = BaseFBAPkg.build_constraint(self,"r"+type,0,0,{
                        self.parameters["bio_rxn"].forward_variable:-1*self.parameters["use_"+object.id[0:-5]+"_class"][0],
                        object.reverse_variable:-1
                    },object)
                    object.upper_bound = 0            
                #If the value is zero, lock both variables to zero
                if self.parameters["use_"+object.id[0:-5]+"_class"][0] == 0:
                    object.lower_bound = 0
                    object.upper_bound = 0
            elif self.parameters["use_"+object.id[0:-5]+"_class"][1] >= 0:
                if self.parameters["use_"+object.id[0:-5]+"_class"][0] >= 0:
                    const = BaseFBAPkg.build_constraint(self,"f"+type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:self.parameters["use_"+object.id[0:-5]+"_class"][1],
                        object.forward_variable:-1
                    },object)
                    BaseFBAPkg.build_constraint(self,"r"+type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:-1*self.parameters["use_"+object.id[0:-5]+"_class"][0],
                        object.forward_variable:1
                    },object)
                    object.lower_bound = 0
                else:
                    const = BaseFBAPkg.build_constraint(self,"f"+type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:self.parameters["use_"+object.id[0:-5]+"_class"][1],
                        object.forward_variable:-1
                    },object)
                    BaseFBAPkg.build_constraint(self,"r"+type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:-1*self.parameters["use_"+object.id[0:-5]+"_class"][0],
                        object.reverse_variable:-1
                    },object)
            else:
                const = BaseFBAPkg.build_constraint(self,"f"+type,0,None,{
                    self.parameters["bio_rxn"].forward_variable:self.parameters["use_"+object.id[0:-5]+"_class"][1],
                    object.reverse_variable:1
                },object)
                BaseFBAPkg.build_constraint(self,"r"+type,0,None,{
                    self.parameters["bio_rxn"].forward_variable:-1*self.parameters["use_"+object.id[0:-5]+"_class"][0],
                    object.reverse_variable:-1
                },object)
                object.upper_bound = 0
            return const
            
    def class_complete(self,class_coef,met_class):
        for msid in classes[met_class]:
            if msid not in class_coef[met_class]:
                return 0
        return 1
