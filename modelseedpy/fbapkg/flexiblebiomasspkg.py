# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! Neither import is ever used
from cobra import Model, Reaction, Metabolite  # !!! Model and Metabolite are never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

classes = {
    "rna":{"cpd00052":-1,"cpd00038":-1,"cpd00002":-1,"cpd00062":-1},
    "dna":{"cpd00115":-1,"cpd00356":-1,"cpd00241":-1,"cpd00357":-1},
    "protein":{"cpd00132":-1, "cpd00023":-1, "cpd00053":-1, "cpd00054":-1, "cpd00033":-1, "cpd00039":-1, "cpd00119":-1, "cpd00051":-1, "cpd00041":-1, "cpd00107":-1, "cpd00129":-1, "cpd00322":-1, "cpd00069":-1, "cpd00065":-1, "cpd00060":-1, "cpd00084":-1, "cpd00035":-1, "cpd00161":-1, "cpd00156":-1, "cpd00066":-1},
    "energy":{"cpd00008":1}
}

#Base class for FBA packages
class FlexibleBiomassPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"flexible biomass",{},{
            "flxbio":"reaction", "fflxcpd":"metabolite", "rflxcpd":"metabolite", "fflxcls":"reaction", "rflxcls":"reaction"})
        
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
        newrxns = []  # !!! newrxns is never used
        class_coef = {"rna":{},"dna":{},"protein":{},"energy":{}}
        refcpd = {"cpd00001":None,"cpd00009":None,"cpd00012":None,"cpd00067":None,"cpd00002":None}
        for met in self.model.metabolites:
            for msid in refcpd:
                if FBAHelper.modelseed_id_from_cobra_metabolite(met) == msid:
                    refcpd[msid] = met
        for metabolite in self.parameters["bio_rxn"].metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid != "cpd11416":
                met_class = "none"
                if msid != None:
                    for curr_class, contents in classes.items():
                        if msid in contents:
                            met_class = dict(curr_class)
                            contents[msid] = metabolite
                if any([met_class == "none", self.class_complete(class_coef,met_class), not self.parameters["use_"+met_class+"_class"]]) and msid not in refcpd:
                    drain_reaction = FBAHelper.add_drain_from_metabolite_id(self.model, metabolite.id, 1000, 1000, "FLEX_")
                    if drain_reaction.id not in self.new_reactions:
                        self.new_reactions[drain_reaction.id] = drain_reaction
                        self.model.add_reactions([drain_reaction])
                    self.build_constraint(metabolite,"flxcpd")
        for met_class, content in class_coef.items():
            add = False
            total_coef = 0
            object_stoichiometry = {}
            for msid in content:
                if all([met_class == "rna", msid == "cpd00002", "cpd00008" in class_coef["energy"]]):
                    object_stoichiometry[content[msid]] = self.parameters["bio_rxn"].metabolites[content[msid]] + self.parameters["bio_rxn"].metabolites[class_coef["energy"]["cpd00008"]]        
                else:
                    object_stoichiometry[content[msid]] = self.parameters["bio_rxn"].metabolites[content[msid]]
                total_coef += abs(object_stoichiometry[content[msid]])
            if any([met_class == "rna", met_class == "dna"]) and refcpd["cpd00012"] and refcpd["cpd00001"]:
                add = True
                object_stoichiometry[refcpd["cpd00012"]] = total_coef
                object_stoichiometry[refcpd["cpd00001"]] = total_coef
            if met_class == "protein" and refcpd["cpd00001"]:
                add = True
                object_stoichiometry[refcpd["cpd00001"]] = total_coef
            if met_class == "energy" and all([refcpd["cpd00001"], refcpd["cpd00002"], refcpd["cpd00067"], refcpd["cpd00009"]]):
                add = True
                object_stoichiometry[refcpd["cpd00001"]] = -1*total_coef
                object_stoichiometry[refcpd["cpd00002"]] = -1*total_coef
                object_stoichiometry[refcpd["cpd00009"]] = total_coef
                object_stoichiometry[refcpd["cpd00067"]] = total_coef
            if add:
                if met_class+"_flex" not in self.new_reactions:
                    self.new_reactions[met_class+"_flex"] = Reaction(id=met_class+"_flex", name= met_class+"_flex", 
                                                                     lower_bound=-1000, upper_bound=1000)
                    self.new_reactions[met_class+"_flex"].add_metabolites(object_stoichiometry)
                    self.new_reactions[met_class+"_flex"].annotation["sbo"] = 'SBO:0000627'
                    self.model.add_reactions([self.new_reactions[met_class+"_flex"]])
                self.build_constraint(self.new_reactions[met_class+"_flex"],"flxcls")
        self.build_constraint(self.parameters["bio_rxn"],"flxbio")
        
    def build_variable(self,object,type):  # !!! can the function be removed?
        pass
                   
    def build_constraint(self,cobra_obj,obj_type):
        element_mass = FBAHelper.elemental_mass()  # !!! element_mass is never used
        if obj_type == "flxbio":
            #Sum(MW*(vdrn,for-vdrn,ref)) + Sum(massdiff*(vrxn,for-vrxn,ref)) = 0
            coef = {}
            for metabolite in self.parameters["bio_rxn"].metabolites:
                if "FLEX_"+metabolite.id in self.model.reactions:
                    mw = FBAHelper.metabolite_mw(metabolite)
                    sign = -1
                    if self.parameters["bio_rxn"].metabolites[metabolite] > 0:
                        sign = 1
                    coef[self.model.reactions.get_by_id("FLEX_"+metabolite.id).forward_variable] = sign*mw
                    coef[self.model.reactions.get_by_id("FLEX_"+metabolite.id).reverse_variable] = -sign*mw
            for met_class in classes:
                if met_class+"_flex" in self.model.reactions:
                    massdiff = 0
                    rxn = self.model.reactions.get_by_id(met_class+"_flex")
                    for met in rxn.metabolites:
                        mw = FBAHelper.metabolite_mw(met)
                        massdiff += rxn.metabolites[met]*mw
                    if abs(massdiff) > 0.00001:
                        coef[rxn.forward_variable] = massdiff
                        coef[rxn.reverse_variable] = -massdiff
            return BaseFBAPkg.build_constraint(self,obj_type,0,0,coef,cobra_obj)
        elif obj_type == "flxcpd":
            #0.75 * abs(bio_coef) * vbio - vdrn,for >= 0
            #0.75 * abs(bio_coef) * vbio - vdrn,rev >= 0
            coef = self.parameters["flex_coefficient"]*abs(self.parameters["bio_rxn"].metabolites[cobra_obj])
            if coef > 0.75:
                coef = 0.75
            BaseFBAPkg.build_constraint(self,"f"+obj_type,0,None,{
                self.parameters["bio_rxn"].forward_variable:coef,
                self.model.reactions.get_by_id("FLEX_"+cobra_obj.id).forward_variable:-1
            },cobra_obj)
            return BaseFBAPkg.build_constraint(self,"r"+obj_type,0,None,{
                self.parameters["bio_rxn"].forward_variable:coef,
                self.model.reactions.get_by_id("FLEX_"+cobra_obj.id).reverse_variable:-1
            },cobra_obj)
        elif obj_type == "flxcls" and cobra_obj.id[0:-5] != None:
            #0.75 * vbio - vrxn,for >= 0
            #0.75 * vbio - vrxn,rev >= 0
            #First deal with the situation where the flux is locked into a particular value relative to biomass
            const = None
            first_entry = self.parameters["use_"+cobra_obj.id[0:-5]+"_class"][0]
            second_entry = self.parameters["use_"+cobra_obj.id[0:-5]+"_class"][1]
            if first_entry == second_entry:
                #If the value is positive, lock in the forward variable and set the reverse to zero
                if first_entry > 0:
                    const = BaseFBAPkg.build_constraint(self,"f"+obj_type,0,0,{
                        self.parameters["bio_rxn"].forward_variable:second_entry,
                        cobra_obj.forward_variable:-1
                    },cobra_obj)
                    cobra_obj.lower_bound = 0
                #If the value is negative, lock in the reverse variable and set the forward to zero
                elif first_entry < 0:
                    const = BaseFBAPkg.build_constraint(self,"r"+obj_type,0,0,{
                        self.parameters["bio_rxn"].forward_variable:-first_entry,
                        cobra_obj.reverse_variable:-1
                    },cobra_obj)
                    cobra_obj.upper_bound = 0            
                #If the value is zero, lock both variables to zero
                if first_entry == 0:
                    cobra_obj.lower_bound = 0
                    cobra_obj.upper_bound = 0
            elif second_entry >= 0:
                if self.parameters["use_"+cobra_obj.id[0:-5]+"_class"][0] >= 0:
                    const = BaseFBAPkg.build_constraint(self,"f"+obj_type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:second_entry,
                        cobra_obj.forward_variable:-1
                    },cobra_obj)
                    BaseFBAPkg.build_constraint(self,"r"+obj_type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:-first_entry,
                        cobra_obj.forward_variable:1
                    },cobra_obj)
                    cobra_obj.lower_bound = 0
                else:
                    const = BaseFBAPkg.build_constraint(self,"f"+obj_type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:second_entry,
                        cobra_obj.forward_variable:-1
                    },cobra_obj)
                    BaseFBAPkg.build_constraint(self,"r"+obj_type,0,None,{
                        self.parameters["bio_rxn"].forward_variable:-first_entry,
                        cobra_obj.reverse_variable:-1
                    },cobra_obj)
            else:
                const = BaseFBAPkg.build_constraint(self,"f"+obj_type,0,None,{
                    self.parameters["bio_rxn"].forward_variable:second_entry,
                    cobra_obj.reverse_variable:1
                },cobra_obj)
                BaseFBAPkg.build_constraint(self,"r"+obj_type,0,None,{
                    self.parameters["bio_rxn"].forward_variable:-first_entry,
                    cobra_obj.reverse_variable:-1
                },cobra_obj)
                cobra_obj.upper_bound = 0
            return const
            
    def class_complete(self,class_coef,met_class):
        for msid in classes[met_class]:
            if msid not in class_coef[met_class]:
                return True
        return False
