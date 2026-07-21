# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging

logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! Zero is never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class FluxFittingPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self, model, "flux fitting", {"vfitf": "reaction","vfitr": "reaction"}, {"vfitc": "reaction"}
        )

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {
                "target_flux": {},
                "totalflux": 0,
                "set_objective": 1,
                "default_rescaling": 0.1,
                #"rescale_vfit_by_flux": True,
                "quadratic_formulation": True,
                "multiply_vfit_by_flux": False,
            },
        )
        if self.parameters["totalflux"] == 1:
            self.pkgmgr.getpkg("RevBinPkg", 1).build_package(
                self.parameters["target_flux"]
            )
            self.pkgmgr.getpkg("TotalFluxPkg", 1).build_package(
                self.parameters["target_flux"]
            )
        objvars = []
        for rxnid in self.parameters["target_flux"]:
            if rxnid in self.model.reactions:
                rxnobj = self.model.reactions.get_by_id(rxnid)
                var = self.build_fluxfit_variable(rxnobj)
                
                self.build_fluxfit_constraint(rxnobj)
        if self.parameters["set_objective"] == 1:
            self.set_fluxfit_objective()

    def set_fluxfit_objective(self):
        objvars = []
        coef = 1
        if self.parameters["quadratic_formulation"]:
            for rxnid in self.variables["vfitf"]:
                if self.parameters["multiply_vfit_by_flux"] and rxnid in self.parameters["target_flux"]:
                    coef = abs(self.parameters["target_flux"][rxnid])
                    if coef == 0:
                        coef = 0.0001
                objvars.append(coef*self.variables["vfitf"][rxnid]**2)
        else:
            for rxnid in self.variables["vfitf"]:
                if self.parameters["multiply_vfit_by_flux"] and rxnid in self.parameters["target_flux"]:
                    coef = abs(self.parameters["target_flux"][rxnid])
                    if coef == 0:
                        coef = 0.0001
                objvars.append(coef*self.variables["vfitf"][rxnid])
            for rxnid in self.variables["vfitr"]:
                if self.parameters["multiply_vfit_by_flux"] and rxnid in self.parameters["target_flux"]:
                    coef = abs(self.parameters["target_flux"][rxnid])
                    if coef == 0:
                        coef = 0.0001
                objvars.append(coef*self.variables["vfitr"][rxnid])
        self.model.objective = self.model.problem.Objective(
            add(objvars), direction="min", sloppy=True
        )

    def build_fluxfit_variable(self, object):
        lower_bound = -100000
        if not self.parameters["quadratic_formulation"]:
            lower_bound = 0
            self.build_variable(
                "vfitr", 0, 100000, "continuous", object
            )
        self.build_variable(
            "vfitf", lower_bound, 100000, "continuous", object
        )

    def build_fluxfit_constraint(self, cobra_obj):
        # vfit(i) = flux(i) - v(i) <OR> vfitf(i)-vfitr(i) = flux(i) + v(i)
        if cobra_obj.id in self.parameters["target_flux"]:
            flux = self.parameters["target_flux"][cobra_obj.id]
            vfitcoef = 1
            # if self.parameters["rescale_vfit_by_flux"] == True:
            #    if flux != None and abs(flux) > 0:
            #        vfitcoef = vfitcoef*flux#Multiply coef by fit flux which rescales by flux
            #    else:
            #        vfitcoef = vfitcoef*self.parameters["default_rescaling"]#Multiply coef by fit flux which rescales by flux
            coef = {self.variables["vfitf"][cobra_obj.id]: vfitcoef}
            if not self.parameters["quadratic_formulation"]:
                coef[self.variables["vfitr"][cobra_obj.id]] = -1*vfitcoef
            if self.parameters["totalflux"] == 0:
                coef[cobra_obj.forward_variable] = -1
                coef[cobra_obj.reverse_variable] = 1
            else:
                coef[
                    self.pkgmgr.getpkg("TotalFluxPkg").variables["tf"][cobra_obj.id]
                ] = 1  # !!! the total flux package does not return anything
                flux = abs(flux)
            return BaseFBAPkg.build_constraint(
                self, "vfitc", flux, flux, coef, cobra_obj
            )
