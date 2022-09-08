# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper
# from optlang import Objective, Constraint
from optlang.symbolics import Zero
from deepdiff import DeepDiff
import re

logger = logging.getLogger(__name__)

#Base class for FBA packages
class ReactionUsePkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(  # model, name, variables, constraints
            self,model,"reaction use",{"fu":"reaction","ru":"reaction"},
            {"fu":"reaction","ru":"reaction","exclusion":"none","urev":"reaction"}) # minimize ru (influx) for exchanges

    def build_package(self, rxn_filter = None, reversibility = False): # exchange reactions as filter
        for rxn in self.model.reactions:
            #Checking that reaction passes input filter if one is provided
            if rxn_filter == None:
                self.build_variable(rxn,"=")
                self.build_constraint(rxn,reversibility)
            elif rxn.id in rxn_filter:
                self.build_variable(rxn,rxn_filter[rxn.id])
                self.build_constraint(rxn,reversibility)
    
    def build_variable(self,reaction,direction):
        variable = None
        if (direction == ">" or direction == "=") and reaction.upper_bound > 0 and reaction.id not in self.variables["fu"]:
            variable = BaseFBAPkg.build_variable(self,"fu",0,1,"binary",reaction)
        if (direction == "<" or direction == "=") and reaction.lower_bound < 0 and reaction.id not in self.variables["ru"]:
            variable = BaseFBAPkg.build_variable(self,"ru",0,1,"binary",reaction)
        return variable
        
    def build_constraint(self,reaction,reversibility):
        constraint = None
        if reaction.id not in self.constraints["fu"] and reaction.id in self.variables["fu"]:
            constraint = BaseFBAPkg.build_constraint(self,"fu",0,None,{
                self.variables["fu"][reaction.id]:1000, reaction.forward_variable:-1
                },reaction)
        if reaction.id not in self.constraints["ru"] and reaction.id in self.variables["ru"]:
            constraint = BaseFBAPkg.build_constraint(self,"ru",0,None,{
                self.variables["ru"][reaction.id]:1000, reaction.reverse_variable:-1
                },reaction)
        if all([reversibility, reaction.id in self.variables["ru"], reaction.id in self.variables["fu"]]):
            constraint = BaseFBAPkg.build_constraint(self,"urev",None,1,{
                self.variables["ru"][reaction.id]:1, self.variables["fu"][reaction.id]:1
                },reaction)
        return constraint
    
    def build_exclusion_constraint(self, flux_values=None):
        flux_values = flux_values or FBAHelper.compute_flux_values_from_variables(self.model)
        count = len(self.constraints["exclusion"])
        solution_coef = {}
        solution_size = 0 
        for rxnid, flux in flux_values.items():
            if flux > Zero:
                solution_size += 1
                solution_coef[self.variables["fu"][rxnid]] = 1
            elif flux < -1*Zero:
                solution_size += 1
                solution_coef[self.variables["ru"][rxnid]] = 1            
        if len(solution_coef) > 0:
            const_name = "exclusion."+str(count+1)
            self.constraints["exclusion"][const_name] = self.model.problem.Constraint(
                Zero,lb=None,ub=(solution_size-1),name=const_name
            )
            self.model.add_cons_vars(self.constraints["exclusion"][const_name])
            self.model.solver.update()
            self.constraints["exclusion"][const_name].set_linear_coefficients(solution_coef)
            return self.constraints["exclusion"][const_name]
        return None

class MinimalMediaPkg(ReactionUsePkg):
    def __init__(self, model):
        copied_model = model.copy()
        ReactionUsePkg.__init__(self, copied_model)
        self.constraints.update({"knockout": "reaction"})

    def solve_media(self,):
        # define variables and constraints for all 
        for ex_rxn in FBAHelper.exchange_reactions(self.model):
            ReactionUsePkg.build_variable(ex_rxn, '<')
            ReactionUsePkg.build_constraint(ex_rxn, False)
        FBAHelper.add_objective(self.model, sum([var for var in self.variables["ru"].values()]), "min")
            
        # identify the set of minimal media solutions
        solution_dicts = []
        sol = self.model.optimize()
        while sol.status == "optimal":
            sol_dict = FBAHelper.solution_to_variables_dict(sol, self.model)
            solution_dicts.append(sol_dict)
            ## omit the solution from the next search
            BaseFBAPkg.build_constraint(self, "exclusion", len(sol_dict)-1, len(sol_dict)-1,
                                        coef=sol_dict, cobra_obj=f"sol{len(solution_dicts)}")
            sol = self.model.optimize()
        if not solution_dicts:
            logger.error("No simulations were feasible.")
            
        # search the permutation space by omitting previously investigated solution_dicts
        self.interdependencies = {}
        for sol_index, sol_dict in enumerate(solution_dicts):
            sol_ex_cpds = {re.sub("(_[a-z]\d+)", "", rxn.id).replace("EX_", ""): flux 
                           for rxn, flux in sol_dict.items() if "EX_" in rxn.id}
            self.interdependencies[sol_index] = {}
            for ex_cpd in sol_ex_cpds:
                sol_dict_sans_cpd = FBAHelper.solution_to_dict(sol_dict)
                sol_dict_sans_cpd.pop(ex_cpd)
                self.interdependencies[sol_index][ex_cpd] = self._examine_permutations(
                    ex_cpd, sol_dict, sol_index, sol_dict_sans_cpd)
        
    def _examine_permutations(self, ex_cpd, sol_dict, sol_index, sol_dict_sans_cpd):
        interdependencies = {}
        # knockout the selected variables
        coef = {self.variables["ru"][ex_cpd]:0}
        coef.update({self.variables["ru"][ex_cpd2]:1 for ex_cpd2 in sol_dict if ex_cpd != ex_cpd2})
        BaseFBAPkg.build_constraint(self, "knockout", 0.1, None, coef, f"{ex_cpd}-sol{sol_index}")
        new_sol = self.model.optimize()
        
        ## explore permutations after removing the selected variable
        diff = DeepDiff(sol_dict_sans_cpd, FBAHelper.solution_to_dict(new_sol))
        if diff: # either new exchanges or altered fluxes to accommodate the removed compound
            for key, value in diff.items():
                new_mets = {re.search("(?<=[\')(.+)(?=\'])", met):change for met, change in value.items()}
                # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                self.interdependencies[sol_index][ex_cpd].update(new_mets)
            coef = {self.variables["met"][ex_cpd]:0 for cpd in new_mets.keys()}
            coef.update({self.variables["met"][ex_cpd]:1 for ex_cpd in sol_dict if ex_cpd not in new_mets.keys()})
            cpd_name = "_".join(new_mets.keys())
            BaseFBAPkg.build_constraint(self, "met", 0.1, None, coef, f"{cpd_name}-sol{sol_index}")
            new_sol = self.model.optimize()
            if new_sol.status != "optimal":
                return interdependencies
            self._examine_permutations(ex_cpd, new_sol, sol_index, sol_dict_sans_cpd)
        return interdependencies

class MinimalMedia(BaseFBAPkg):
    """A class that determines the minimal media of a model"""  # investigate conversion to a staticmethod or single function for convenience and in-line utility
    def __init__(self, model, min_growth=0.1):
        # define the system
        BaseFBAPkg.__init__(self, model, "Minimal Media", {"met":"metabolite"}, {"met":"string", "obj":"string"})
        
        # define the exchange variables, the minimal objective, and the minimal growth value
        for cpd in FBAHelper.exchange_reactions(self.model):
            BaseFBAPkg.build_variable(self,"met",0,1,"binary",cpd)
        self.model = FBAHelper.add_objective(self.model, sum([var for var in self.variables["met"].values()]), "min")
        BaseFBAPkg.build_constraint(self, "obj", min_growth, None, {
            rxn.forward_variable:1 for rxn in FBAHelper.bio_reactions(self.model)}, "min_growth")
        
        # determine the set of media solutions
        solutions = []
        sol = self.model.optimize()
        while sol.status == "optimal":
            solutions.append(sol)
            sol_dict = FBAHelper.solution_to_variables_dict(sol, model)
            ## omit the solution from the next search
            BaseFBAPkg.build_constraint(self, "met", len(sol_dict)-1, len(sol_dict)-1, 
                                        coef=sol_dict, cobra_obj=f"sol{len(solutions)}")
            sol = self.model.optimize()
        if not solutions:
            logger.error("No simulations were feasible.")
                
        # search the permutation space by omitting previously investigated solutions
        self.interdependencies = {}
        for sol_index, sol in enumerate(solutions):
            self.interdependencies[sol_index] = {}
            for cpd in sol:
                self.interdependencies[sol_index][cpd] = {}
                coef = {self.variables["met"][cpd]:0}
                coef.update({self.variables["met"][cpd2]:1 for cpd2 in sol if cpd != cpd2})
                BaseFBAPkg.build_constraint(self, "met", sol.objective_value, 
                                            sol.objective_value, coef, f"{cpd}-sol{sol_index}")
                new_sol = self.model.optimize()
                diff = DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))
                
                ## track metabolites that fill the void from the removed compound
                while diff:
                    for key, value in diff.items():
                        new_mets = {re.search("(?<=[\')(.+)(?=\'])", met):change for met, change in value.items()}
                        # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                        self.interdependencies[sol_index][cpd].update(new_mets)
                    diff = self.test_compounds(cpd, sol, sol_index, new_mets.keys())
                        
    def test_compounds(self, cpd, sol, sol_index, zero_compounds):
        coef = {self.variables["met"][cpd]:0 for cpd in zero_compounds}
        coef.update({self.variables["met"][cpd]:1 for cpd in sol if cpd not in zero_compounds})
        cpd_name = "_".join(zero_compounds)
        BaseFBAPkg.build_constraint(self, "met", sol.objective_value, 
                                    sol.objective_value, coef, f"{cpd_name}-sol{sol_index}")
        new_sol = self.model.optimize()
        if new_sol.status != "optimal":
            return False
        return DeepDiff(FBAHelper.solution_to_dict(sol), FBAHelper.solution_to_dict(new_sol))