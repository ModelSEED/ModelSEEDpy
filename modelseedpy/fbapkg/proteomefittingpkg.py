# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
logger = logging.getLogger(__name__)

from optlang.symbolics import add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.multiomics.msexpression import GENOME, COLUMN_NORM

# Options for default behavior
LOWEST = 10

# Base class for FBA packages
class ProteomeFittingPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "proteome fitting",
            {"kapp": "reaction", "kvfit": "reaction", "kfit": "reaction"},
            {"vkapp": "reaction", "kfitc": "reaction"},
        )
        self.pkgmgr.addpkgs(["FluxFittingPkg"])

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            ["proteome", "condition"],
            {
                "flux_values": {},
                "kcat_values": {},
                "prot_coef": 0.1,
                "totalflux": 1,
                "kcat_coef": 0.333,
                "obj_kfit": 1,
                "obj_kvfit": 1,
                "obj_vfit": 1,
                "set_objective": 1,
                "rescale_vfit_by_flux": True,
                "default_rescaling": 0.1,
                "default_expression": LOWEST,
            },
        )
        objvars = []
        # Converting genome proteome to reaction proteome if necessary
        if self.parameters["proteome"].type == GENOME:
            self.parameters["proteome"] = self.parameters[
                "proteome"
            ].build_reaction_expression(
                self.model, self.parameters["default_expression"]
            )
        # Checking for condition in proteome and converting condition string to object if necessary
        if isinstance(self.parameters["condition"], str):
            if (
                self.parameters["condition"]
                not in self.parameters["proteome"].conditions
            ):
                logger.critical(
                    "Condition "
                    + self.parameters["condition"]
                    + " not found in proteome!"
                )
            self.parameters["condition"] = self.parameters[
                "proteome"
            ].conditions.get_by_id(self.parameters["condition"])
        # Adding flux fitting variables and constraints
        self.pkgmgr.getpkg("FluxFittingPkg").build_package(
            {
                "target_flux": self.parameters["flux_values"],
                "totalflux": self.parameters["totalflux"],
                "rescale_vfit_by_flux": self.parameters["rescale_vfit_by_flux"],
                "default_rescaling": self.parameters["default_rescaling"],
                "set_objective": 0,
            }
        )
        self.pkgmgr.getpkg("TotalFluxPkg").build_package()
        for rxnid in self.pkgmgr.getpkg("FluxFittingPkg").variables["vfit"]:
            objvars.append(
                self.parameters["obj_vfit"]
                * self.pkgmgr.getpkg("FluxFittingPkg").variables["vfit"][rxnid] ** 2
            )
        # Adding proteome fitting variables and constraints
        for rxnobj in self.model.reactions:
            # Only make constraints and variables if reaction is in the proteome
            if rxnobj.id in self.parameters["proteome"].features:
                self.build_variable(rxnobj, "kapp")
                var = self.build_variable(rxnobj, "kvfit")
                objvars.append(self.parameters["obj_kvfit"] * var**2)
                const = self.build_constraint(rxnobj, "vkapp")
        # Adding kcat fitting variables and constraints
        for rxnid in self.parameters["kcat_values"]:
            for rxnobj in self.model.reactions:
                if rxnid == FBAHelper.modelseed_id_from_cobra_reaction(rxnobj):
                    if rxnobj.id not in self.variables["kapp"]:
                        self.build_variable(rxnobj, "kapp")
                    var = self.build_variable(rxnobj, "kfit")
                    const = self.build_constraint(rxnobj, "kfitc")
                    objvars.append(self.parameters["obj_kfit"] * var**2)
        # Creating objective function
        if self.parameters["set_objective"] == 1:
            self.model.objective = self.model.problem.Objective(
                add(objvars), direction="min", sloppy=True
            )

    def build_variable(self, object, type):
        if type == "kapp":
            return BaseFBAPkg.build_variable(
                self, type, -1000000, 1000000, "continuous", object
            )
        if type == "kvfit":
            return BaseFBAPkg.build_variable(
                self, type, -1000, 1000, "continuous", object
            )
        elif type == "kfit":
            return BaseFBAPkg.build_variable(
                self, type, -1000000, 1000000, "continuous", object
            )

    def build_constraint(self, object, type):
        if type == "vkapp" and object.id in self.parameters["proteome"].features:
            # kvfit(i) = kapp(i)*ProtCoef*Prot(i) - v(i)
            # Pulling expression value for selected condition and reaction
            expval = self.parameters["proteome"].get_value(
                object.id, self.parameters["condition"], COLUMN_NORM
            )
            if expval is None and self.parameters["default_expression"] is not None:
                if self.parameters["default_expression"] == LOWEST:
                    expval = (
                        self.parameters["condition"].lowest
                        / self.parameters["condition"].column_sum
                    )
            if expval is not None:
                prot = expval * self.parameters["prot_coef"]
                if (
                    object.id in self.variables["kvfit"]
                    and object.id in self.variables["kapp"]
                ):
                    coef = {
                        self.variables["kvfit"][object.id]: 1,
                        self.variables["kapp"][object.id]: -1 * prot,
                    }
                    if self.parameters["totalflux"] == 1:
                        coef[
                            self.pkgmgr.getpkg("TotalFluxPkg").variables["tf"][
                                object.id
                            ]
                        ] = 1
                    else:
                        coef[object.forward_variable] = 1
                        coef[object.reverse_variable] = 1
                    return BaseFBAPkg.build_constraint(self, type, 0, 0, coef, object)
        elif type == "kfitc":
            # kfit(i) = kapp(i) - kcoef*kcat(i)
            msid = FBAHelper.modelseed_id_from_cobra_reaction(object)
            rhs = (
                -1 * self.parameters["kcat_values"][msid] * self.parameters["kcat_coef"]
            )
            coef = {
                self.variables["kfit"][object.id]: 1,
                self.variables["kapp"][object.id]: -1,
            }
            return BaseFBAPkg.build_constraint(self, type, rhs, rhs, coef, object)
