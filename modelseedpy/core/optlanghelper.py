# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:26:32 2022

@author: Andrew Freiburger
"""
from collections import namedtuple
from optlang import Model
from typing import Iterable
import logging

logger = logging.getLogger(__name__)

Bounds = namedtuple("Bounds", ("lb", "ub"), defaults=(0,1000))
tupVariable = namedtuple("tupVariable", ("name", "bounds", "type"), defaults=("varName", Bounds(), "continuous"))
tupConstraint = namedtuple("tupConstraint", ("name", "bounds", "expr"), defaults=("consName", Bounds(0,0), None))
tupObjective = namedtuple("tupObjective", ("name", "expr", "direction"), defaults=("objectiveName", None, "max"))

def isIterable(term):
    try:
        iter(term)
        if type(term) is not str:
            return True
        return False
    except:
        return False
    
def isnumber(obj):
    try:
        float(obj)
        return True
    except:
        return False
    
def define_term(value):
    if isinstance(value, str):
        return {"type":"Symbol", "name": value}
    elif isinstance(value, (float, int)):
        return {"type":"Number", "value": value}
    print(f"ERROR: The {value} is not known.")
    
def get_expression_template(expr):
    if isinstance(expr, list):
        return {"type": "Add", "args": []}
    return {"type": expr["operation"], "args": []}
    
class OptlangHelper:
    
    @staticmethod
    def add_variables(var_name:str, var_bounds:Iterable, var_type:str="continuous"):
        return {"name": var_name.replace(" ", "_"), "lb": var_bounds[0], "ub": var_bounds[1], "type": var_type}
    
    @staticmethod
    def add_constraint(cons_name:str, cons_bounds:Iterable, cons_expr:dict):
        return {"name": cons_name.replace(" ", "_"),
        "expression": OptlangHelper._define_expression(cons_expr),
         "lb": cons_bounds[0], "ub": cons_bounds[1], "indicator_variable": None, "active_when": 1}
    
    @staticmethod
    def add_objective(obj_name:str, obj_expr:dict, direction:str):
        return {"name": obj_name.replace(" ", "_"),
        "expression": OptlangHelper._define_expression(obj_expr),
        "direction": direction}
    
    @staticmethod
    def define_model(model_name, variables, constraints, objective, optlang=False):
        model = {'name':model_name, 'variables':[], 'constraints':[],
            "objective": OptlangHelper.add_objective(objective[0], objective[1], objective[2])}
        for var in variables:
            if len(var) == 2:
                var.append("continuous")
            model["variables"].append(OptlangHelper.add_variables(var[0], var[1], var[2]))
        for cons in constraints:
            model["constraints"].append(OptlangHelper.add_constraint(cons[0], cons[1], cons[2]))
                                      
        if optlang:
            return Model.from_json(model)
        return model
    
    @staticmethod
    def _define_expression(expr:dict):
        expression = get_expression_template(expr)
        for ele in expr["elements"]:
            if not isnumber(ele) and not isinstance(ele, str):
                arguments = []
                for ele2 in ele["elements"]:
                    if not isnumber(ele2) and not isinstance(ele2, str):
                        print("recursive ele\t\t", type(ele2), ele2)
                        arguments.append(OptlangHelper._define_expression(ele2))
                    else:
                        arguments.append(define_term(ele2))
                expression["args"].append(get_expression_template(ele))
                expression["args"][-1]["args"] = arguments
            else:
                expression["args"].append(define_term(ele))
        return expression      
    
    @staticmethod
    def heuns_dot_product(coefficients, zipped_to_sum):
        expression = {"type": "Add", "args": []}
        for index, (term1, term2) in enumerate(zipped_to_sum):
            expression["args"].append({"type": "Mul", "args": [
                {"type": "Add", "args": [define_term(term1), define_term(term2)]},
                define_term(coefficients[index])
            ]})
        return expression
    
    @staticmethod
    def dot_product(zipped_to_sum):
        expression = {"type": "Add", "args": []}
        for (term1, term2) in zipped_to_sum:
            expression["args"].append({"type": "Mul", "args": [
                {"type": "Add", "args": [define_term(term1), define_term(term2)]}
            ]})
        return expression