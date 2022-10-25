# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.core.fbahelper import FBAHelper

#Base class for FBA packages
class CommKineticPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"community kinetics",{},{"commkin":"string"})

    def build_package(self,kinetic_coef,community_model=None):
        self.validate_parameters({},[],{
            "kinetic_coef":kinetic_coef,
            "community":community_model if community_model else MSCommunity(self.model)
        })
        for species in self.parameters["community"].species:
            self.build_constraint(species)

    def build_constraint(self,species):
        bioRXN = species.biomasses[0]
        coef = {bioRXN.forward_variable:-self.parameters["kinetic_coef"]}
        for rxn in self.model.reactions:
            if int(FBAHelper.rxn_compartment(rxn)[1:]) == species.index and rxn != bioRXN:
                coef.update({rxn.forward_variable:1, rxn.reverse_variable:1})
        return BaseFBAPkg.build_constraint(self,"commkin",None,0,coef,"Species"+str(species.index))