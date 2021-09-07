# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper


class KBaseMediaPkg(BaseFBAPkg):
    """
    Base class for FBA packages
    """

    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "kbase media", {}, {})

    def build_package(self, media_or_parameters, default_uptake=None, default_excretion=None):
        if isinstance(media_or_parameters, dict):
            self.validate_parameters(media_or_parameters, [], {
                "default_uptake": 0,
                "default_excretion": 100,
                "media": None
            })
        else:
            self.validate_parameters({}, [], {
                "default_uptake": default_uptake,
                "default_excretion": default_excretion,
                "media": media_or_parameters
            })
            if self.parameters["default_uptake"] is None:
                self.parameters["default_uptake"] = 0
            if self.parameters["default_excretion"] is None:
                self.parameters["default_excretion"] = 100    
        if self.parameters["media"] is None:
            self.parameters["default_uptake"] = 100

        exchange_reactions = {}
        for reaction in self.model.exchanges:
            if reaction.id[:3] == 'EX_':
                print(reaction)
                compound = reaction.id[3:]
                exchange_reactions[compound] = reaction                
                reaction.lower_bound = -1*self.parameters["default_uptake"]
                reaction.upper_bound = self.parameters["default_excretion"]
                reaction.update_variable_bounds()  # FIXME: this seems unnecessary

        if self.parameters["media"]:
            # Searching for media compounds in model
            for compound in self.parameters["media"].mediacompounds:
                mdlcpds = self.find_model_compounds(compound.id)
                for mdlcpd in mdlcpds:
                    if mdlcpd.id in exchange_reactions:
                        print(mdlcpd)
                        exchange_reactions[mdlcpd.id].lower_bound = -1 * compound.maxFlux
                        exchange_reactions[mdlcpd.id].upper_bound = -1 * compound.minFlux
                    if self.pkgmgr != None and "FullThermoPkg" in self.pkgmgr.packages:
                        print(mdlcpd)
                        if mdlcpd.id in self.variables["logconc"] and mdlcpd.compartment == "e0":
                            if compound.concentration != 0.001:
                                self.variables["logconc"][msid_hash[compound.id].id].lb = compound.concentration
                                self.variables["logconc"][msid_hash[compound.id].id].ub = compound.concentration
    
    def find_model_compounds(self, cpd_id):
        if cpd_id in self.model.metabolites:
            return [self.model.metabolites.get_by_id(cpd_id)]
        return [m for m in self.model.metabolites if FBAHelper.modelseed_id_from_cobra_metabolite(m) == cpd_id]
