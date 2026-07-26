# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from scipy.constants import physical_constants, calorie, kilo, R
from optlang.symbolics import Zero, add
import numpy as np
from numpy import log as ln
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper

# Base class for FBA packages
class FullThermoPkg(BaseFBAPkg):
    @staticmethod
    def default_concentration():
        return {
            "cpd00067_c0": [0.0000001, 0.0000001],  # M H+ - equivalent to pHint = 7
            "cpd00007_c0": [1e-07, 8.2e-06],  # M O2 instracellular
            "cpd00011_c0": [1e-08, 0.0014],  # M CO2 instracellular
            "cpd00067_e0": [
                3.16228e-08,
                3.16228e-07,
            ],  # M H+ - equivalent to pHext = 6.5
            "cpd00009_e0": [
                0.056,
                0.056,
            ],  # Extracellular phosphate - overridden by media when media concentration is given
            "cpd00048_e0": [
                0.0030,
                0.0030,
            ],  # Extracellular sulfate - overridden by media when media concentration is given
            "cpd00013_e0": [
                0.019,
                0.019,
            ],  # Extracellular ammonia - overridden by media when media concentration is given
            "cpd00971_e0": [
                0.16,
                0.16,
            ],  # Extracellular sodium - overridden by media when media concentration is given
            "cpd00205_e0": [
                0.022,
                0.022,
            ],  # Extracellular potassium - overridden by media when media concentration is given
            "cpd10515_e0": [
                0.062,
                0.062,
            ],  # Extracellular Fe2+ - overridden by media when media concentration is given
            "cpd00011_e0": [
                0.00010,
                0.00010,
            ],  # Extracellular CO2 - overridden by media when media concentration is given
            "cpd00007_e0": [
                8.2e-06,
                8.2e-06,
            ],  # Extracellular O2 - overridden by media when media concentration is given
            "cpd00027_e0": [
                0.020,
                0.020,
            ],  # Extracellular glucose - overridden by media when media concentration is given
        }

    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "full thermo",
            {"logconc": "metabolite", "pdgerr": "metabolite","ndgerr": "metabolite", "pconcfit": "metabolite", "nconcfit": "metabolite"},
            {"potc": "metabolite","concfit": "metabolite"},
        )
        self.pkgmgr.addpkgs(["SimpleThermoPkg"])

    def build_package(self, parameters, verbose=True):
        self.validate_parameters(
            parameters,
            [],
            {
                "default_max_conc": 0.02,  # M
                "default_min_conc": 0.000001,  # M
                "default_max_error": 5,  # KJ/mol
                "custom_concentrations": {},
                "custom_deltaG_error": {"cpd00067": 0},
                "compartment_potential": {
                    "e0": 0,    # mV — extracellular MUST be the zero reference for compartment electrochemical potential (so community modeling works)
                    "c0": -160, # mV = 0.33*(pHint - pHext) - 143.33 where pHint=7, pHext=6.5
                },
                "temperature": 298,  # K
                "filter": None,
                "infeasible_model": False,
                "dgbin": False,
                "modelseed_db_path": None,
                "exclude_reactions": [],
                "set_min_error_objective": False,
                "set_concfit_objective": False,
                "concfit_objective_type": "l1",
                "concfit_weights": None,
            },
        )
        simple_thermo_parameters = {
            "filter": self.parameters["filter"],
            "min_potential": -100000,  # KJ/mol
            "max_potential": 100000,  # KJ/mol
            "dgbin": self.parameters["dgbin"],
            "exclude_reactions": self.parameters["exclude_reactions"]
        }
        if self.parameters["infeasible_model"]:
            simple_thermo_parameters["dgbin"] = True
        self.pkgmgr.getpkg("SimpleThermoPkg").build_package(simple_thermo_parameters)

        self.parameters[
            "combined_custom_concentrations"
        ] = FullThermoPkg.default_concentration()
        for cpd in self.parameters["custom_concentrations"]:
            self.parameters["combined_custom_concentrations"][cpd] = self.parameters[
                "custom_concentrations"
            ][cpd]
        msid_hash = {}
        for metabolite in self.model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid != None:
                if msid not in msid_hash:
                    msid_hash[msid] = {}
                msid_hash[msid][metabolite.id] = metabolite
            # Build concentration variable
            self.build_logconc_variable(metabolite)
            # Build error variable
            self.build_dgerr_variables(metabolite)
            # Build the potential constraint
            self.build_gibbs_constraint(metabolite)
        if self.parameters["set_min_error_objective"]:
            self.build_min_error_objective()
        if self.parameters["set_concfit_objective"]:
            if "target_concentrations" not in self.parameters:
                raise ValueError(
                    "build_package: 'target_concentrations' is required in parameters "
                    "when 'set_concfit_objective' is True."
                )
            target_concentrations = self.parameters["target_concentrations"]
            concfit_weights = self.parameters["concfit_weights"]
            objective_type = self.parameters["concfit_objective_type"]
            if objective_type == "l2":
                self.build_concfit_objective_l2(target_concentrations, concfit_weights)
            elif objective_type == "l1":
                self.build_concfit_objective(target_concentrations, concfit_weights)
            else:
                raise ValueError(
                    f"build_package: 'concfit_objective_type' must be 'l1' or 'l2', got {objective_type!r}."
                )

    def build_min_error_objective(self):
        objvars = []
        for metabolite in self.model.metabolites:
            if metabolite.id in self.variables["pdgerr"]:
                objvars.append(1 * self.variables["pdgerr"][metabolite.id])
            if metabolite.id in self.variables["ndgerr"]:
                objvars.append(1 * self.variables["ndgerr"][metabolite.id])
        self.model.objective = self.model.problem.Objective(
            add(objvars), direction="min", sloppy=True
        )

    @staticmethod
    def _validate_concfit_weights(weights):
        """Validate a metabolite-id -> weight dict: every weight must be a
        finite number >= 0. Raises ValueError otherwise. A weight of 0 is
        valid (it drops that metabolite's term from the objective)."""
        for metabolite_id, weight in weights.items():
            if weight is None or not np.isfinite(weight):
                raise ValueError(
                    f"Weight for metabolite '{metabolite_id}' must be a finite number; got {weight!r}."
                )
            if weight < 0:
                raise ValueError(
                    f"Weight for metabolite '{metabolite_id}' must be >= 0; got {weight!r}."
                )

    def _is_qp_capable_solver(self):
        """Detect whether self.model is currently attached to a QP-capable
        optlang backend (gurobi or cplex), by module identity of the solver
        interface. GLPK, SCIPY, and the hybrid interface are not QP-capable
        (see build_concfit_objective_l2 for why this matters)."""
        interface_name = getattr(self.model.solver.interface, "__name__", "")
        return interface_name.startswith("optlang.gurobi_interface") or interface_name.startswith(
            "optlang.cplex_interface"
        )

    def _require_qp_capable_solver(self):
        if not self._is_qp_capable_solver():
            raise NotImplementedError(
                "L2 concentration-fit objective requires a QP-capable solver backend "
                "(gurobi/cplex) - set model.solver before calling build_package"
            )

    def build_concfit_objective(self, target_concentrations, weights=None):
        weights = weights or {}
        self._validate_concfit_weights(weights)
        #Removing all existing concfit constraints
        for metabolite in self.model.metabolites:
            if metabolite.id in self.constraints["concfit"]:
                self.model.remove_cons_vars(self.constraints["concfit"][metabolite.id])
        objvars = []
        for metabolite in self.model.metabolites:
            if metabolite.id in target_concentrations:
                w = weights.get(metabolite.id, 1.0)
                if metabolite.id not in self.variables["pconcfit"]:
                    self.build_variable("pconcfit", 0, 100, "continuous", metabolite)
                if metabolite.id not in self.variables["nconcfit"]:
                    self.build_variable("nconcfit", 0, 100, "continuous", metabolite)
                self.build_constraint(
                    "concfit", target_concentrations[metabolite.id], target_concentrations[metabolite.id], {self.variables["pconcfit"][metabolite.id]:1,self.variables["nconcfit"][metabolite.id]:-1,self.variables["logconc"][metabolite.id]:1}, metabolite
                )
                objvars.append(w * self.variables["pconcfit"][metabolite.id])
                objvars.append(w * self.variables["nconcfit"][metabolite.id])
        self.model.objective = self.model.problem.Objective(
            add(objvars), direction="min", sloppy=True
        )

    def build_concfit_objective_l2(self, target_concentrations, weights=None):
        """Weighted least-squares concentration-fit objective:
        minimize sum_i w_i * (logconc_i - target_i)^2 over the fitted set
        (metabolites present in target_concentrations). Unlike
        build_concfit_objective, this uses no pconcfit/nconcfit slack
        variables — the squared term is expressible directly against the
        logconc variable — so any slacks/constraints left over from a prior
        L1 build for the fitted set are removed first.

        Requires a QP-capable optlang backend (gurobi/cplex); raises
        NotImplementedError early otherwise, since GLPK cannot hold a
        quadratic objective at all.
        """
        weights = weights or {}
        self._require_qp_capable_solver()
        self._validate_concfit_weights(weights)

        # Remove any pre-existing concfit constraints (mirrors build_concfit_objective).
        for metabolite in self.model.metabolites:
            if metabolite.id in self.constraints["concfit"]:
                self.model.remove_cons_vars(self.constraints["concfit"][metabolite.id])
                del self.constraints["concfit"][metabolite.id]

        # Remove any pre-existing pconcfit/nconcfit slack variables for the fitted set —
        # an L2 build after an L1 build must leave no stale slacks behind.
        stale_vars = []
        for metabolite_id in target_concentrations:
            if metabolite_id in self.variables["pconcfit"]:
                stale_vars.append(self.variables["pconcfit"].pop(metabolite_id))
            if metabolite_id in self.variables["nconcfit"]:
                stale_vars.append(self.variables["nconcfit"].pop(metabolite_id))
        if stale_vars:
            self.model.remove_cons_vars(stale_vars)

        objvars = []
        for metabolite in self.model.metabolites:
            if metabolite.id in target_concentrations and metabolite.id in self.variables["logconc"]:
                w = weights.get(metabolite.id, 1.0)
                target = target_concentrations[metabolite.id]
                objvars.append(w * (self.variables["logconc"][metabolite.id] - target) ** 2)
        self.model.objective = self.model.problem.Objective(
            add(objvars), direction="min", sloppy=True
        )

    def build_logconc_variable(self, object):
        msid = self.modelutl.metabolite_msid(object)
        if msid == "cpd00001":
            return None
        lb = ln(self.parameters["default_min_conc"])
        ub = ln(self.parameters["default_max_conc"])
        if object.id in self.parameters["combined_custom_concentrations"]:
            lb = ln(self.parameters["combined_custom_concentrations"][object.id][0])
            ub = ln(self.parameters["combined_custom_concentrations"][object.id][1])
        if msid in self.parameters["combined_custom_concentrations"]:
            lb = ln(self.parameters["combined_custom_concentrations"][msid][0])
            ub = ln(self.parameters["combined_custom_concentrations"][msid][1])
        return self.build_variable(
            "logconc", lb, ub, "continuous", object
        )
    
    def build_dgerr_variables(self, object):
        msid = self.modelutl.metabolite_msid(object)
        ub = self.parameters["default_max_error"]
        if msid in self.parameters["custom_deltaG_error"]:
            ub = self.parameters["custom_deltaG_error"][msid]
        elif object.id in self.parameters["custom_deltaG_error"]:
            ub = self.parameters["custom_deltaG_error"][object.id]
        self.build_variable(
            "pdgerr", 0, ub, "continuous", object
        )
        self.build_variable(
            "ndgerr", 0, ub, "continuous", object
        )
    
    def build_concfit_variables(self, object):
        self.build_variable(
            "pconcfit", 0, 100, "continuous", object
        )
        self.build_variable(
            "nconcfit", 0, 100, "continuous", object
        )

    def build_gibbs_constraint(self, object):
        msid = self.modelutl.metabolite_msid(object)

        # potential(i) (kJ/mol) = deltaG(i) (kJ/mol) + R*T(K)*lnconc(i)
        #                         + charge(i) * Faraday(kC/mol) * compartment_potential(mV) / 1000
        # compartment_potential is in mV throughout (build_package default: e0=0 mV, c0=-160 mV).
        # Dividing by 1000 converts mV -> V so that kC/mol * V = kJ/mol.
        if object.id not in self.pkgmgr.getpkg("SimpleThermoPkg").variables["potential"]:
            return None

        pos_deltagerr_variable = self.variables["pdgerr"][object.id]
        neg_deltagerr_variable = self.variables["ndgerr"][object.id]
        potential_variable = self.pkgmgr.getpkg("SimpleThermoPkg").variables["potential"][object.id]

        #First check if compound has a deltaG value in the model
        deltag = None
        if object.notes and 'deltag' in object.notes:
            deltag = object.notes['deltag']
        elif msid is not None:
            from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
            if self.parameters["modelseed_db_path"] is not None:
                biochem_db =  ModelSEEDBiochem.get(path=self.parameters["modelseed_db_path"])
            else:
                biochem_db =  ModelSEEDBiochem.get()
            mscpd = biochem_db.compounds.get_by_id(msid)
            if mscpd is not None:
                if hasattr(mscpd, 'deltag') and mscpd.deltag is not None:
                    deltag = mscpd.deltag
                elif hasattr(compound, 'delta_g') and compound.delta_g is not None:
                    deltag= mscpd.delta_g
        RValue = 0.008314 # kJ/mol/K
        if deltag is None:
            return None
        if msid == "cpd00067": #Adding pH7 concentration to deltaG so H+ concentration cancels unless it's different from pH7
            deltag += -1 * RValue * self.parameters["temperature"] * np.log(0.0000001)/4.184 #Dividing by 4.184 here because we multiply by this value later
        Faraday = 96.485  # kC/mol
        compartment_potential = 0  # mV
        if object.compartment in self.parameters["compartment_potential"]:
            compartment_potential = self.parameters["compartment_potential"][
                object.compartment
            ]
        # Sanity-check: compartment_potential is expected in mV.
        # Typical biological membrane potentials are well under 300 mV in magnitude.
        # Values above 300 mV are suspicious (possible volts/mV unit slip; e.g. 0.16 V
        # mis-entered as 160 would not trigger this, but -160 mV correctly does not).
        # Values above 1000 mV are egregiously non-physical.
        abs_cp = abs(compartment_potential)
        if abs_cp > 1000:
            logging.warning(
                "compartment_potential for '%s' is %.1f mV, which is egregiously "
                "non-physical (|V| > 1000 mV). Verify units — expected mV.",
                object.compartment, compartment_potential,
            )
        elif abs_cp > 300:
            logging.warning(
                "compartment_potential for '%s' is %.1f mV, which is unusually large "
                "(|V| > 300 mV). Typical biological membrane potentials sit well under "
                "300 mV; verify that mV (not V) was supplied.",
                object.compartment, compartment_potential,
            )
        # Divide compartment_potential by 1000 to convert mV -> V.
        # Faraday (kC/mol) * V = kJ/mol, matching the kJ/mol scale of deltaG.
        constant = deltag * 4.184 + object.charge * Faraday * (compartment_potential / 1000)
        coef = {
            potential_variable: 1,
            pos_deltagerr_variable: -1,
            neg_deltagerr_variable: 1,
        }
        if msid != "cpd00001":  # Water concentration should not contribute to potential
            coef[self.variables["logconc"][object.id]] = -1 * RValue * self.parameters["temperature"]
        return self.build_constraint(
            "potc", constant, constant, coef, object
        )
