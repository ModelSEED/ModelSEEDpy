# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

#Base class for FBA packages
class CommKineticPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"community kinetics",{},{"compkin":"string"})

    def build_package(self,kinetic_coef,abundances = None,predict_abundance = False):
        self.validate_parameters({},[],{
            "kinetic_coef":kinetic_coef,
            "abundances":abundances,
            "predict_abundance" : predict_abundance,
            "biohash" : {}
        })
        for reaction in self.model.reactions:
            if re.search('^bio\d+$', reaction.id) != None:
                for metabolite in reaction.metabolites:
                    msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
                    if FBAHelper.modelseed_id_from_cobra_metabolite(metabolite) == "cpd11416" and re.search('[a-z](\d+)', metabolite.compartment) != None:
                        m = re.search('[a-z](\d+)', metabolite.compartment)
                        index = m[1]
                        biorxnid = reaction.id
                        if index not in self.parameters["biohash"]:
                            self.parameters["biohash"][index] = []
                        if index == "0" or reaction.id != "bio1":
                            self.parameters["biohash"][index].append(reaction.id)
        if abundances != None:
            print("Warning:Community biomass reaction is being altered with input abundances")
            if "0" in self.parameters["biohash"]:
                primarybiomass = self.model.reactions.get_by_id(self.parameters["biohash"]["0"][0])
                totalabundance = 0
                for species in abundances:
                    totalabundance += abundances[species]
                for species in abundances:
                    abundances[species] = abundances[species]/totalabundance
                for metabolite in primarybiomass.metabolites:
                    if metabolite.id[0:8] == "cpd11416":
                        index = metabolite.id[10:]
                        if index != "0":
                            if index not in abundances:
                                if int(index) in abundances:
                                    primarybiomass.add_metabolites({metabolite:abundances[index]},combine=False)
                                else:
                                    primarybiomass.metabolites[metabolite] = 0
                            else:
                                primarybiomass.add_metabolites({metabolite:abundances[index]},combine=False)
                self.model.solver.update()
        for index in self.parameters["biohash"]:
            if index != "0" and index in self.parameters["biohash"]:
                self.build_constraint(index)
        
    def build_constraint(self,index):
        coef = dict()
        for bio in self.parameters["biohash"][index]:
            biorxn = self.model.reactions.get_by_id(bio)
            coef[biorxn.forward_variable] = -1*self.parameters["kinetic_coef"]
        for reaction in self.model.reactions:
            comp = reaction.id.split("_").pop()
            if comp[1:] == index:
                coef[reaction.forward_variable] = 1
                coef[reaction.reverse_variable] = 1
        return BaseFBAPkg.build_constraint(self,"compkin",None,0,coef,"Species"+index)