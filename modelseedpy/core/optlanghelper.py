# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 10:26:32 2022

@author: Andrew Freiburger
"""
from collections import namedtuple
from optlang import Model
from typing import Iterable, Union
from pprint import pprint
import logging

logger = logging.getLogger(__name__)

Bounds = namedtuple("Bounds", ("lb", "ub"), defaults=(0, 1000))
tupVariable = namedtuple(
    "tupVariable",
    ("name", "bounds", "type"),
    defaults=("varName", Bounds(), "continuous"),
)
tupConstraint = namedtuple(
    "tupConstraint",
    ("name", "bounds", "expr"),
    defaults=("consName", Bounds(0, 0), None),
)
tupObjective = namedtuple(
    "tupObjective",
    ("name", "expr", "direction"),
    defaults=("objectiveName", None, "max"),
)


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
    if isnumber(value):
        return {"type": "Number", "value": value}
    if isinstance(value, str):
        return {"type": "Symbol", "name": value}
    print(f"ERROR: The {value} of type {type(value)} is not known.")


def get_expression_template(expr):
    # print(expr)
    if isinstance(expr, list):
        return {"type": "Add", "args": []}
    return {"type": expr["operation"], "args": []}


class OptlangHelper:

    @staticmethod
    def add_variables(
        var_name: str, var_bounds: (list, tuple), var_type: str = "continuous"
    ):
        return {
            "name": var_name.replace(" ", "_"),
            "lb": var_bounds[0],
            "ub": var_bounds[1],
            "type": var_type,
        }

    @staticmethod
    def add_constraint(cons_name: str, cons_bounds: (list, tuple), cons_expr: dict):
        return {
            "name": cons_name.replace(" ", "_"),
            "expression": OptlangHelper._define_expression(cons_expr),
            "lb": cons_bounds[0],
            "ub": cons_bounds[1],
            "indicator_variable": None,
            "active_when": 1,
        }

    @staticmethod
    def add_objective(obj_name: str, objective_expr: Union[dict, list], direction: str):
        if isinstance(objective_expr, list):
            obj_expr = {
                "type": "Add",
                "args": [
                    OptlangHelper._define_expression(expr) for expr in objective_expr
                ],
            }
        elif isinstance(objective_expr, dict):
            obj_expr = {
                "type": objective_expr["operation"],
                "args": [define_term(term) for term in objective_expr["elements"]],
            }
        return {
            "name": obj_name.replace(" ", "_"),
            "expression": obj_expr,
            "direction": direction,
        }

    @staticmethod
    def define_model(model_name, variables, constraints, objective, optlang=False):
        model = {"name": model_name, "variables": [], "constraints": []}
        # pprint(objective)
        for var in variables:
            if len(var) == 2:
                var.append("continuous")
            model["variables"].append(
                OptlangHelper.add_variables(var[0], var[1], var[2])
            )
        for cons in constraints:
            model["constraints"].append(
                OptlangHelper.add_constraint(cons[0], cons[1], cons[2])
            )
        # if not isinstance(obj, str): # catches a strange error of the objective name as the objective itself
        model["objective"] = OptlangHelper.add_objective(
            objective[0], objective[1], objective[2]
        )
        if optlang:
            return Model.from_json(model)
        return model

    @staticmethod
    def _define_expression(expr: dict):
        expression = get_expression_template(expr)
        level1_coef = 0
        for ele in expr["elements"]:
            if not isnumber(ele) and not isinstance(ele, str):
                # print(expr, ele, end="\r")
                arguments = []
                level2_coef = 0
                for ele2 in ele["elements"]:
                    if not isnumber(ele2) and not isinstance(ele2, str):
                        # print("recursive ele\t\t", type(ele2), ele2)
                        arguments.append(OptlangHelper._define_expression(ele2))
                    elif isinstance(ele2, str):
                        arguments.append(define_term(ele2))
                    else:
                        level2_coef += float(ele2)
                expression["args"].append(get_expression_template(ele))
                if level2_coef != 0:
                    arguments.append(define_term(level2_coef))
                expression["args"][-1]["args"] = arguments
            elif isinstance(ele, str):
                expression["args"].append(define_term(ele))
            else:
                level1_coef += float(ele)
        if level1_coef != 0:
            expression["args"].append(define_term(level1_coef))
        # pprint(expression)
        return expression

    @staticmethod
    def dot_product(zipped_to_sum, heuns_coefs=None):
        # ensure that the lengths are compatible for heun's dot-products
        if heuns_coefs is not None:
            coefs = (
                heuns_coefs
                if isinstance(heuns_coefs, (list, set))
                else heuns_coefs.tolist()
            )
            zipped_length = len(zipped_to_sum)
            coefs_length = len(coefs)
            if zipped_length != coefs_length:
                raise IndexError(
                    f"ERROR: The length of zipped elements {zipped_length}"
                    f" is unequal to that of coefficients {coefs_length}"
                )

        elements = []
        for index, (term1, term2) in enumerate(zipped_to_sum):
            if heuns_coefs is not None:
                elements.extend(
                    [
                        {"operation": "Mul", "elements": [heuns_coefs[index], term1]},
                        {"operation": "Mul", "elements": [heuns_coefs[index], term2]},
                    ]
                )
            else:
                elements.append({"operation": "Mul", "elements": [term1, term2]})
        return elements
