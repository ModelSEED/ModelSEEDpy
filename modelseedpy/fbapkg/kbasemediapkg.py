# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper   # !!! imported but not used

logger = logging.getLogger(__name__)

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
        if self.parameters["media"] is None and self.parameters["default_uptake"] == 0:
            self.parameters["default_uptake"] = 100
        
        #First initializing all exchanges to default uptake and excretion
        exchange_list = self.modelutl.exchange_list()
        for reaction in exchange_list:
            reaction.lower_bound = -1*self.parameters["default_uptake"]
            reaction.upper_bound = self.parameters["default_excretion"]
        
        #Now constraining exchanges for specific compounds specified in the media
        if self.parameters["media"]:
            exchange_hash = self.modelutl.exchange_hash()
            self.modelutl.build_metabolite_hash()
            for mediacpd in self.parameters["media"].mediacompounds:
                mets = self.modelutl.find_met(mediacpd.id)
                if len(mets) > 0:
                    for met in mets:
                        if met in exchange_hash:
                            exchange_hash[met].lower_bound = -1 * mediacpd.maxFlux
                            exchange_hash[met].upper_bound = -1 * mediacpd.minFlux
                            if self.pkgmgr != None and "FullThermoPkg" in self.pkgmgr.packages:
                                logger.info('FullThermo constrained compound: ', met.id)
                                if met.id in self.variables["logconc"] and met.compartment[0:1] == "e":
                                    if mediacpd.concentration != 0.001:
                                        self.variables["logconc"][met.id].lb = ln(mediacpd.concentration)
                                        self.variables["logconc"][met.id].ub = ln(mediacpd.concentration)
                else:
                    logger.warn('Media compound: ', mediacpd.id,' not found in model.')