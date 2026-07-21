# -*- coding: utf-8 -*-
import logging
import re
import traceback
import cobra
from cobra.flux_analysis import pfba
from cobra.flux_analysis import flux_variability_analysis
from modelseedpy.core.msmodelutl import MSModelUtil

logger = logging.getLogger(__name__)

class MSFBA:
    def __init__(self,model_or_mdlutl,media,objective_reactions={"bio1":1},maximize=True,gene_ko=[],reaction_ko=[],pfba=True,fva=True,clone=True,primary_solution=None,id=None):
        if isinstance(model_or_mdlutl, MSModelUtil):
            model_or_mdlutl = model_or_mdlutl.model
        if clone:
            model_or_mdlutl = cobra.io.json.from_json(cobra.io.json.to_json(model_or_mdlutl))
        self.model = model_or_mdlutl
        self.mdlutl = MSModelUtil.get(model_or_mdlutl)        
        self.media = media
        self.objective_reactions = objective_reactions
        self.maximize = maximize
        self.gene_ko = gene_ko
        self.reaction_ko = reaction_ko
        self.pkgmgr = self.mdlutl.pkgmgr
        self.apply_parameters()
        self.primary_solution = primary_solution
        self.secondary_solutions = None
        self.fva = fva
        self.pfba = pfba
        self.fva_results = None
        self.secondary_fva = None
        if id == None:
            id = self.mdlutl.model.id+".fba"
        self.id = id

    def build_objective(self):
        sense = "max"
        if not self.maximize:
            sense = "min"
        obj = self.model.problem.Objective(0, direction=sense)
        objcoef = {}
        for rxnid in self.objective_reactions:
            if rxnid in self.model.reactions:
                rxn = self.model.reactions.get_by_id(rxnid)
                objcoef[rxn.forward_variable] = self.objective_reactions[rxnid]
                objcoef[rxn.reverse_variable] = -1*self.objective_reactions[rxnid]
            else:
                logger.warning(f"KO reaction {rxnid} not found in model")
        obj.set_linear_coefficients(objcoef)

    def apply_parameters(self):
        self.pkgmgr.getpkg("KBaseMediaPkg").build_package(self.media)
        for gene in self.gene_ko:
            if gene in self.model.genes:
                self.model.genes.get_by_id(gene).knock_out()
            else:
                logger.warning(f"KO gene {gene} not found in model")
        for rxn in self.reaction_ko:
            if rxn in self.model.reactions:
                self.model.reactions.get_by_id(rxn).knock_out()
            else:
                logger.warning(f"KO reaction {rxn} not found in model")
    
    def run(self):
        if self.pfba:
            self.primary_solution = pfba(self.model)
        else:
            self.primary_solution = self.model.optimize()
        if self.fva:
            self.fva_results = flux_variability_analysis(self.model)

    def add_secondary_solution(self,solution,fva=None):
        if self.secondary_solutions == None:
            self.secondary_solutions = []
        self.secondary_solutions.append(solution)
        if fva:
            if self.secondary_fva == None:
                self.secondary_fva = []
            self.secondary_fva.append(fva)

    def get_variable_class(self,variable_min, variable_max):
        variable_class = "Unknown"
        if variable_min is None or variable_max is None:
            return variable_class
        if variable_min == 0 and variable_max == 0:
            variable_class = "Blocked"
        elif variable_min > 0 and variable_max > 0:
            variable_class = "Positive"
        elif variable_min >= 0 and variable_max > 0:
            variable_class = "Positive variable"
        elif variable_min < 0 and variable_max < 0:
            variable_class = "Negative"
        elif variable_min < 0 and variable_max <= 0:
            variable_class = "Negative variable"
        else:
            variable_class = "Variable"
        return variable_class
    
    def generate_kbase_data(self,fbamodel_ref,media_ref):
        output = {
            "FBABiomassVariables": [],
            "FBACompoundBounds": [],
            "FBACompoundVariables": [],
            "FBAConstraints": [],
            "FBADeletionResults": [],
            "FBAMetaboliteProductionResults": [],
            "FBAMinimalMediaResults": [],
            "FBAMinimalReactionsResults": [],
            "FBAPromResults": [],
            "FBAReactionBounds": [],
            "FBAReactionVariables": [],
            "FBATintleResults": [],
            "MFALog": "",
            "PROMKappa": 1,
            "QuantitativeOptimizationSolutions": [],
            "__VERSION__": 1,
            "additionalCpd_refs": [],
            "allReversible": 0,
            "biomassRemovals": {},
            "biomassflux_objterms": {"bio1": 1},
            "calculateReactionKnockoutSensitivity": 0,
            "comboDeletions": 0,
            "compoundflux_objterms": {},
            "decomposeReversibleDrainFlux": 0,
            "decomposeReversibleFlux": 0,
            "defaultMaxDrainFlux": 0,
            "defaultMaxFlux": 1000,
            "defaultMinDrainFlux": -1000,
            "drainfluxUseVariables": 0,
            "fbamodel_ref": fbamodel_ref,
            "findMinimalMedia": 0,
            "fluxMinimization": 1,
            "fluxUseVariables": 0,
            "fva": 0,
            "gapfillingSolutions": [],
            "geneKO_refs": [],
            "id": self.id,
            "inputfiles": {},
            "maximizeActiveReactions": 0,
            "maximizeObjective": 1,
            "media_list_refs": [],
            "media_ref": media_ref,
            "minimizeErrorThermodynamicConstraints": 0,
            "minimize_reaction_costs": {},
            "minimize_reactions": 0,
            "noErrorThermodynamicConstraints": 0,
            "numberOfSolutions": 1,
            "objectiveConstraintFraction": 0.1,
            "objectiveValue": self.primary_solution.objective_value,
            "other_objectives": [],
            "outputfiles": {},
            "parameters": {
                "Auxotrophy metabolite list": "",
                "Beachhead metabolite list": "",
                "minimum_target_flux": "0.01",
                "save phenotype fluxes": "0",
                "suboptimal solutions": "1",
            },
            "quantitativeOptimization": 0,
            "reactionKO_refs": [],
            "reactionflux_objterms": {},
            "simpleThermoConstraints": 0,
            "thermodynamicConstraints": 0,
            "uptakeLimits": {},
        }

        for rxn in self.model.reactions:
            flux = 0
            if rxn.id in self.primary_solution.fluxes:
                flux = self.primary_solution.fluxes[rxn.id]
            min_flux = rxn.lower_bound
            max_flux = rxn.upper_bound
            if self.fva_results and rxn.id in self.fva_results:
                min_flux, max_flux = self.fva_results[rxn.id]
            other_mins= []
            other_maxes = []
            other_fluxes = []
            if self.secondary_solutions:
                for sol in self.secondary_solutions:
                    if rxn.id in sol.fluxes:
                        other_fluxes.append(sol.fluxes[rxn.id])
                    else:
                        other_fluxes.append(0)
            if self.secondary_fva:
                othermin = rxn.lower_bound
                othermax = rxn.upper_bound
                for fva in self.secondary_fva:
                    if rxn.id in fva:
                        othermin, othermax = fva[rxn.id]
                other_mins.append(othermin)
                other_maxes.append(othermax)
            variable_class = self.get_variable_class(min_flux, max_flux)
            variable_data = {
                "class": variable_class,
                "lowerBound": rxn.lower_bound,
                "max": max_flux,
                "min": min_flux,
                "upperBound": rxn.upper_bound,
                "other_max": other_maxes,
                "other_min": other_mins,
                "other_values": other_fluxes,
                "value": flux,
                "variableType": "flux"
            }
            variable_key = "FBAReactionVariables"
            if rxn.id.startswith("EX_"):
                lower = variable_data["lowerBound"]
                variable_data["lowerBound"] = -1 * variable_data["upperBound"]
                variable_data["upperBound"] = -1 * lower
                lower = variable_data["min"]
                variable_data["min"] = -1 * variable_data["max"]
                variable_data["max"] = -1 * lower
                variable_data["value"] = -1 * variable_data["value"]
                variable_data["variableType"] = "drainflux"
                variable_data["modelcompound_ref"] = "~/fbamodel/modelcompounds/id/" + rxn.id[3:]
                variable_key = "FBACompoundVariables"
            elif rxn.id.startswith("bio"):
                variable_data["variableType"] = "biomassflux"
                variable_data["biomass_ref"] = "~/fbamodel/biomasses/id/" + rxn.id
                variable_key = "FBABiomassVariables"
            else:
                variable_data["modelreaction_ref"] = "~/fbamodel/modelreactions/id/" + rxn.id
                variable_data["exp_state"] = "unknown"
                variable_data["biomass_dependencies"] = []
                variable_data["coupled_reactions"] = []
                variable_data["expression"] = 0
                variable_data["scaled_exp"] = 0
            output[variable_key].append(variable_data)    
        return output