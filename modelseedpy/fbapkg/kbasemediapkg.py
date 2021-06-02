# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class KBaseMediaPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"kbase media",{},{})

    def build_package(self,media,default_uptake = None,default_excretion = None):
        if default_uptake == None:
            default_uptake = 0
        if default_excretion == None:
            default_excretion = 100
        if media == None:
            default_uptake = 100
        bound_hash = {}
        if media != None:
            for compound in media.mediacompounds:
                bound_hash[compound.id] = {
                    "lb" : -1 * compound.maxFlux,
                    "ub" : -1 * compound.minFlux,
                }
        for reaction in self.model.reactions:
            if reaction.id[0:3].lower() == "ex_":
                compound = reaction.id[3:]
                if compound[-3:] == "_e0":
                    compound = compound[:-3]
                if compound in bound_hash:
                    reaction.lower_bound = bound_hash[compound]["lb"]
                    reaction.upper_bound = bound_hash[compound]["ub"]
                else:
                    reaction.lower_bound = -1*default_uptake
                    reaction.upper_bound = default_excretion
                reaction.update_variable_bounds()