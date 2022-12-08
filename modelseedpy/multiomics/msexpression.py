# -*- coding: utf-8 -*-
import logging

import re
import copy
from cobra.core.dictlist import DictList
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
from ast import And, BitAnd, BitOr, BoolOp, Expression, Name, NodeTransformer, Or
from modelseedpy.core.msgenome import MSGenome, MSFeature

# Types of expression data
GENOME = 10
MODEL = 20

# Types of normalization
COLUMN_NORM = 10

logger = logging.getLogger(__name__)


def compute_gene_score(expr, values, default):
    if isinstance(expr, Expression):
        return compute_gene_score(expr.body, values, default)
    elif isinstance(expr, Name):
        if expr.id in values:
            return values[expr.id]
        else:
            return default
    elif isinstance(expr, BoolOp):
        op = expr.op
        if isinstance(op, Or):
            total = 0
            for subexpr in expr.values:
                total += compute_gene_score(subexpr, values, default)
            return total
        elif isinstance(op, And):
            least = None
            for subexpr in expr.values:
                value = compute_gene_score(subexpr, values, default)
                if least == None or value < least:
                    least = value
            return least
        else:
            raise TypeError("unsupported operation " + op.__class__.__name__)
    elif expr is None:
        return default
    else:
        raise TypeError("unsupported operation  " + repr(expr))


class MSCondition:
    def __init__(self, id):
        self.id = id
        self.column_sum = None
        self.feature_count = None
        self.lowest = None


class MSExpressionFeature:
    def __init__(self, feature, parent):
        self.id = feature.id
        self.feature = feature
        self.values = {}
        self.parent = parent

    def add_value(self, condition, value):
        if condition in self.values:
            condition.feature_count += -1
            condition.column_sum += -1 * value
            logger.warning(
                "Overwriting value "
                + str(self.values[condition])
                + " with "
                + str(value)
                + " in feature "
                + self.feature.id
            )
        if condition.lowest is None or condition.lowest > value:
            condition.lowest = value
        condition.feature_count += 1
        condition.column_sum += value
        self.values[condition] = value

    def get_value(self, condition, normalization=None):
        if isinstance(condition, str):
            if condition not in self.parent.conditions:
                logger.warning(
                    "Condition " + condition + " not found in expression object!"
                )
                return None
            condition = self.parent.conditions.get_by_id(condition)
        if condition not in self.values:
            logger.info(
                "Condition " + condition.id + " has no value in " + self.feature.id
            )
            return None
        if normalization == COLUMN_NORM:
            return self.values[condition] / condition.column_sum
        return self.values[condition]


class MSExpression:
    def __init__(self, type):
        self.type = type
        self.object = None
        self.features = DictList()
        self.conditions = DictList()

    @staticmethod
    def from_gene_feature_file(filename, genome=None, create_missing_features=False):
        expression = MSExpression(GENOME)
        if genome == None:
            expression.object = MSGenome()
            create_missing_features = True
        else:
            expression.object = genome
        data = ""
        with open(filename, "r") as file:
            data = file.read()
        lines = data.split("\n")
        conditions = None
        for line in lines:
            if conditions == None:
                conditions = []
                headers = line.split("\t")
                for i in range(1, len(headers)):
                    if headers[i] not in expression.conditions:
                        conditions.append(MSCondition(headers[i]))
                        expression.conditions.append(conditions[i - 1])
                    else:
                        conditions.append(self.conditions.get_by_id(headers[i]))
                    conditions[i - 1].column_sum = 0
                    conditions[i - 1].feature_count = 0
            else:
                array = line.split("\t")
                protfeature = expression.add_feature(array[0], create_missing_features)
                if protfeature != None:
                    for i in range(1, len(array)):
                        protfeature.add_value(conditions[i - 1], float(array[i]))
        return expression

    def add_feature(self, id, create_gene_if_missing=False):
        if id in self.features:
            return self.features.get_by_id(id)
        feature = None
        if self.type == GENOME:
            if self.object.search_for_gene(id) == None:
                if create_gene_if_missing:
                    self.object.features.append(MSFeature(id, ""))
            feature = self.object.search_for_gene(id)
        else:
            if id in self.object.reactions:
                feature = self.object.reactions.get_by_id(id)
        if feature == None:
            logger.warning(
                "Feature referred by expression " + id + " not found in genome object!"
            )
            return None
        if feature.id in self.features:
            return self.features.get_by_id(feature.id)
        protfeature = MSExpressionFeature(feature, self)
        self.features.append(protfeature)
        return protfeature

    def get_value(self, feature, condition, normalization=None):
        if isinstance(feature, str):
            if feature not in self.features:
                logger.warning(
                    "Feature " + feature + " not found in expression object!"
                )
                return None
            feature = self.features.get_by_id(feature)
        return feature.get_value(condition, normalization)

    def build_reaction_expression(self, model, default):
        if self.type == MODEL:
            logger.critical(
                "Cannot build a reaction expression from a model-based expression object!"
            )
        # Creating the expression and features
        rxnexpression = MSExpression(MODEL)
        rxnexpression.object = model
        for rxn in model.reactions:
            if len(rxn.genes) > 0:
                rxnexpression.add_feature(rxn.id)
        for condition in self.conditions:
            rxnexpression.conditions.append(condition)
        # Pulling the gene values from the current expression
        values = {}
        logger.warning("TESTING!")
        for gene in model.genes:
            feature = self.object.search_for_gene(gene.id)
            if feature == None:
                logger.warning(
                    "Model gene " + gene.id + " not found in genome of expression"
                )
            elif feature.id not in self.features:
                logger.warning(
                    "Model gene " + gene.id + " in genome but not in expression"
                )
            else:
                feature = self.features.get_by_id(feature.id)
                for condition in self.conditions:
                    if condition.id not in values:
                        values[condition.id] = {}
                    if condition in feature.values:
                        values[condition.id][gene.id] = feature.values[condition]
        # Computing the reaction level values
        for condition in rxnexpression.conditions:
            for feature in rxnexpression.features:
                tree = parse_gpr(feature.feature.gene_reaction_rule)[0]
                feature.add_value(
                    condition, compute_gene_score(tree, values[condition.id], default)
                )
        return rxnexpression
