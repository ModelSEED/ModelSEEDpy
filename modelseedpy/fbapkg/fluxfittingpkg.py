# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

# Base class for FBA packages
class FluxFittingPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self, model, "flux fitting", {"vfit": "reaction"}, {"vfitc": "reaction"}
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
                "rescale_vfit_by_flux": True,
            },
        )
        if self.parameters["totalflux"] == 0:
            self.pkgmgr.getpkg("RevBinPkg", 1).build_package(
                self.parameters["target_flux"]
            )
        else:
            self.pkgmgr.getpkg("TotalFluxPkg", 1).build_package(
                self.parameters["target_flux"]
            )
        objvars = []
        for rxnid in self.parameters["target_flux"]:
            if rxnid in self.model.reactions:
                rxnobj = self.model.reactions.get_by_id(rxnid)
                var = self.build_variable(rxnobj)
                objvars.append(var**2)
                const = self.build_constraint(rxnobj)
        if self.parameters["set_objective"] == 1:
            self.model.objective = self.model.problem.Objective(
                add(objvars), direction="min", sloppy=True
            )

    def build_variable(self, object):
        return BaseFBAPkg.build_variable(
            self, "vfit", -1000, 1000, "continuous", object
        )

    def build_constraint(self, object):
        # vfit(i) = flux(i) - v(i)
        if object.id in self.parameters["target_flux"]:
            flux = self.parameters["target_flux"][object.id]
            vfitcoef = 1
            # if self.parameters["rescale_vfit_by_flux"] == True:
            #    if flux != None and abs(flux) > 0:
            #        vfitcoef = vfitcoef*flux#Multiply coef by fit flux which rescales by flux
            #    else:
            #        vfitcoef = vfitcoef*self.parameters["default_rescaling"]#Multiply coef by fit flux which rescales by flux
            coef = {self.variables["vfit"][object.id]: vfitcoef}
            if self.parameters["totalflux"] == 0:
                coef[object.forward_variable] = 1
                coef[object.reverse_variable] = -1
            else:
                coef[self.pkgmgr.getpkg("TotalFluxPkg").variables["tf"][object.id]] = 1
                flux = abs(flux)
            return BaseFBAPkg.build_constraint(self, "vfitc", flux, flux, coef, object)
