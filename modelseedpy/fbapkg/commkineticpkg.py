# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
import re
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class CommKineticPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"community kinetics",{},{"compkin":"string"})

    def build_package(self,kinetic_coef,abundances = None,predict_abundance = 0):
        biohash = dict()
        for reaction in self.model.reactions:
            if reaction.id != "bio1" and re.search('^bio\d+$', reaction.id) != None:
                for metabolite in reaction.metabolites:
                    if metabolite.id[0:8] == "cpd11416":
                        biorxnid = reaction.id
                        index = metabolite.id[10:]
                        if index not in biohash:
                            biohash[index] = []
                        biohash[index].append(reaction.id)
        if abundances != None:
            print("Warning:Community biomass reaction is being altered with input abundances")
            if "0" in biohash:
                primarybiomass = self.model.reactions.get_by_id(biohash["0"][0])
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
                                    primarybiomass.metabolites[metabolite] = abundances[int(index)]
                                else:
                                    primarybiomass.metabolites[metabolite] = 0
                            else:
                                primarybiomass.metabolites[metabolite] = abundances[index]
                self.model.solver.update()    
        for index in biohash:
            if index != "0":
                self.build_constraint(index,biohash,kinetic_coef)
        
    def build_constraint(self,index,biohash,kinetic_coef):
        indexlen = len(index)
        coef = dict()
        for bio in biohash[index]:
            biorxn = self.model.reactions.get_by_id(bio)
            coef[biorxn.forward_variable] = -1*kinetic_coef
        for reaction in self.model.reactions:
            comp = reaction.id.split("_").pop()
            if comp[1:] == index:
                coef[reaction.forward_variable] = 1
                coef[reaction.reverse_variable] = 1
        return BaseFBAPkg.build_constraint(self,"compkin",None,0,coef,"Species"+index)