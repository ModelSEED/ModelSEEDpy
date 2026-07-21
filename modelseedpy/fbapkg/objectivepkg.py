# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from optlang.symbolics import Zero, add
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.WARNING#INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

class ObjectiveTerm:
    def __init__(self, variable, coefficient,direction=""):
        self.coefficient = coefficient
        self.variable = variable
        self.direction = direction

    @staticmethod
    def from_string(term_string):
        coefficient = 1
        variable = None
        direction = ""
        #Checking for coefficient
        if term_string[0:1] == "(":
            array = term_string.split(")")
            coefficient = float(array[0][1])
            term_string = array[1]
        #Checking for a +/- on term
        if term_string[0:1] == "+" or term_string[0:1] == "-":
            variable = term_string[1:]
            direction = term_string[0:1]
        else:
            variable = term_string
            direction = ""
        return ObjectiveTerm(variable, coefficient, direction)
           
    def to_string(self):
        return "("+str(self.coefficient)+")"+self.direction+self.variable


#Class for defining an objective function in a modelseedpy model.
class ObjectiveData:
    def __init__(self, terms, sense="max"):
        self.sense = sense
        self.terms = terms

    @staticmethod
    def from_string(objective_string):
        sense = "max"
        terms = []
        if objective_string[0:3] == "MAX":
            objective_string = objective_string[4:-1]#Clearing out the directionality MAX{}
        elif objective_string[0:3] == "MIN":
            objective_string = objective_string[4:-1]#Clearing out the directionality MIN{}
            sense = "min"
        term_strings = objective_string.split("|")
        for term_string in term_strings:
            term = ObjectiveTerm.from_string(term_string)
            terms.append(term)
        return ObjectiveData(terms, sense)
        
    def to_string(self):
        objective_string = ""
        if self.sense == "max":
            objective_string += "MAX{"
        else:
            objective_string += "MIN{"
        for term in self.terms:
            objective_string += term.to_string()+"|"
        objective_string = objective_string[:-1] + "}"
        return objective_string
    
    def to_cobrapy_objective(self, model):
        #Creating empty objective
        objective = model.problem.Objective(Zero, direction=self.sense)
        #Parsing the terms
        coefficients = {}
        for term in self.terms:
            if term.variable in model.reactions:
                coef = term.coefficient
                rxnobj = model.reactions.get_by_id(term.variable)
                if term.direction == "+":
                    coefficients[rxnobj.forward_variable] = coef
                elif term.direction == "-":
                    coefficients[rxnobj.reverse_variable] = coef
                else:
                    coefficients[rxnobj.forward_variable] = coef
                    coefficients[rxnobj.reverse_variable] = -1*coef
            else:
                logger.warning("Reaction "+term.variable+" not found in model")
        model.objective = objective
        objective.set_linear_coefficients(coefficients)
        return objective

# Base class for FBA packages
class ObjectivePkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(self, model, "objective builder", {}, {})
        self.original_model_objective = None
        self.objective_name = None
        self.objective_data = None
        self.objective_data_cache = {}

    def build_package(self,objective_or_string,objective_name=None,set_objective=True):
        #Caching the current objective
        self.original_model_objective = self.model.objective
        #check if input is a string or an ObjectiveData object
        if isinstance(objective_or_string, str):
            self.objective_data = ObjectiveData.from_string(objective_or_string)
        elif isinstance(objective_or_string, ObjectiveData):
            self.objective_data = objective_or_string
        else:
            raise TypeError("Input must be a string or an ObjectiveData object")
        #Setting default objective name if not provided
        self.objective_name = objective_name
        if objective_name == None:
            self.objective_name = self.objective_data.to_string()
        #Caching objective with name
        self.objective_data_cache[self.objective_name] = self.objective_data
        #Creating the objective in the model
        if set_objective:
            self.objective_data_cache[self.objective_name].to_cobrapy_objective(self.model)
        return objective_name

    def restore_objective(self,name):
        self.original_model_objective = self.model.objective
        if name in self.objective_data_cache:
            self.model.objective = self.objective_data_cache[name].to_cobrapy_objective(self.model)
        else:
            logger.warning("Objective "+name+" not found in cache")