# -*- coding: utf-8 -*-

from __future__ import absolute_import
import logging

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.fbahelper import FBAHelper
from optlang import Objective, Constraint
from optlang.symbolics import Zero
from deepdiff import DeepDiff
import re

#Base class for FBA packages
class ReactionUsePkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(
            self,model,"reaction use",{"fu":"reaction","ru":"reaction"},
            {"fu":"reaction","ru":"reaction","exclusion":"none","urev":"reaction"})

    def build_package(self, rxn_filter = None, reversibility = False):
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

class MinimalMedia(BaseFBAPkg):
    """A class that determines the minimal media of a model while operating like a single function"""
    def __init__(self, model, min_growth=0.1):
        # define the system
        BaseFBAPkg.__init__(self, model, "Minimal Media", {"met":"metabolite"}, {"met":"metabolite", "obj":"reaction"})
        
        # define the variables and objective
        for cpd in FBAHelper.exchange_reactions(self.model):
            BaseFBAPkg.build_variable(self,"met",0,1,"binary",cpd)
        self.model.objective = Objective(sum([var for var in self.variables["met"]]), direction="min")
        BaseFBAPkg.build_constraint(self, "obj", min_growth, None, {
            self.model.reactions.bio1.flux_expression:1}, "min_growth")
        
        # determine the set of media solutions
        solutions = []
        sol = self.model.optimize()
        while sol.status == "optimal":
            solutions.append(sol)
            ## omit the solution from the next search
            obj_val = sol.objective_value
            BaseFBAPkg.build_constraint(self, "met", len(sol)-1, len(sol)-1, coef=FBAHelper.solution_to_dict(sol))
            sol = self.model.optimize()
                
        # search the permutation space by omitting previously investigated solutions
        self.interdependencies = {}
        for sol_index, sol in enumerate(solutions):
            self.interdependencies[sol_index] = {}
            for cpd in sol:
                self.interdependencies[sol_index][cpd] = {}
                coef = {self.variables["met"][cpd]:0}
                coef.update({self.variables["met"][cpd2]:1 for cpd2 in sol if cpd != cpd2})
                obj_val = sol.objective_value
                BaseFBAPkg.build_constraint(self, "met", obj_val, obj_val, coef, f"{cpd}-sol{sol_index}")
                new_sol = FBAHelper.solution_to_dict(self.model.optimize())
                
                ## track metabolites that fill the void from the removed compound
                diff = DeepDiff(FBAHelper.solution_to_dict(sol), new_sol)
                while diff:
                    for key, value in diff.items(): # TODO loop over the replacement compounds
                        new_mets = {re.search("(?<=\[\')(.+)(?=\'\])", met):change for met, change in value.items()}
                        # this dictionary should be parsed into a list of substitute metabolites and a list of functionally coupled reactions
                        self.interdependencies[sol_index][cpd].update(new_mets)
                    diff = self.test_compounds(cpd, sol, sol_index, new_mets.keys())
                        
    def test_compounds(self, cpd, sol, sol_index, zero_compounds):
        coef = {self.variables["met"][cpd]:0 for cpd in zero_compounds}
        coef.update({self.variables["met"][cpd]:1 for cpd in sol if cpd not in zero_compounds})
        obj_val = sol.objective_value
        cpd_name = "_".join(zero_compounds)
        BaseFBAPkg.build_constraint(self, "met", obj_val, obj_val, coef, f"{cpd_name}-sol{sol_index}")
        new_sol = FBAHelper.solution_to_dict(self.model.optimize())
        
        ## track metabolites that fill the void from the removed compound
        diff = DeepDiff(FBAHelper.solution_to_dict(sol), new_sol)
        return diff
                    