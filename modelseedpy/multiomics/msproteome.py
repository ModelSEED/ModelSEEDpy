import logging

import re
import copy
from cobra.core.dictlist import DictList
from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
from ast import And, BitAnd, BitOr, BoolOp, Expression, Name, NodeTransformer, Or
from modelseedpy.core.msgenome import MSGenome, MSFeature

logger = logging.getLogger(__name__)

def compute_gene_score(expr,values,default):
    if isinstance(expr, Expression):
        return compute_gene_score(expr.body,values,default)
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
                total += compute_gene_score(subexpr,values,default)
            return total
        elif isinstance(op, And):
            least = None
            for subexpr in expr.values:
                value = compute_gene_score(subexpr,values,default)
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
        
class MSProteomeFeature:

    def __init__(self, feature):
        self.id = feature.id
        self.feature = feature
        self.values = {}
        
    def add_value(self,condition,value):
        if condition in self.values:
            logger.warning("Overwriting value "+str(self.values[condition])+" with "+str(value)+" in feature "+self.feature.id)
        self.values[condition] = value

class MSProteome:

    def __init__(self,type):
        self.type = type
        self.object = None
        self.features = DictList()
        self.conditions = DictList()

    @staticmethod
    def from_gene_feature_file(filename, genome = None,create_missing_features = False):
        proteome = MSProteome("genome")
        if genome == None:
            proteome.object = MSGenome()
            create_missing_features = True
        else:
            proteome.object = genome
        data = ""
        with open(filename, 'r') as file:
            data = file.read()
        lines = data.split("\n")
        conditions = None
        for line in lines:
            if conditions == None:
                conditions = []
                headers = line.split("\t")
                for i in range(1,len(headers)):
                    if headers[i] not in proteome.conditions:
                        conditions.append(MSCondition(headers[i]))
                        proteome.conditions.append(conditions[i-1])
                    else:
                        conditions.append(self.conditions.get_by_id(headers[i]))
            else:
                array = line.split("\t")
                protfeature = proteome.add_feature(array[0], create_missing_features)
                if protfeature != None:
                    for i in range(1,len(array)):
                        protfeature.add_value(conditions[i-1], float(array[i]))
        return proteome

    def add_feature(self,id,create_gene_if_missing = False):
        if id in self.features:
            return self.features.get_by_id(id)
        feature = None
        if self.type == "genome": 
            if self.object.search_for_gene(id) == None:
                if create_gene_if_missing:
                    self.object.features.append(MSFeature(id,""))
            feature = self.object.search_for_gene(id)
        else:
            if id in self.object.reactions:
                feature = self.object.reactions.get_by_id(id)
        if feature == None:
            logger.warning("Gene referred by proteome "+id+" not found in genome object!")
            return None
        if feature.id in self.features:
            return self.features.get_by_id(feature.id)
        protfeature = MSProteomeFeature(feature)
        self.features.append(protfeature)
        return protfeature
            
    def get_value(self,feature,condition):
        if feature not in self.features:
            logger.warning("Feature "+feature+" not found in proteome object!")
            return None
        if condition not in self.conditions:
            logger.warning("Condition "+condition+" not found in proteome object!")
            return None
        feature = self.features.get_by_id(feature)
        condition = self.conditions.get_by_id(condition)
        if condition not in feature.values:
            logger.warning("Feature "+feature.id+" has no data for condition "+condition.id)
            return None
        return feature.values[condition]
            
    def build_reaction_proteome(self,model,default):
        if self.type == "model":
            logger.critical("Cannot build a reaction proteome from a model-based proteome object!")
        #Creating the proteome and features
        rxnproteome = MSProteome("model")
        rxnproteome.object = model
        for rxn in model.reactions:
            rxnproteome.add_feature(rxn.id)
        for condition in self.conditions:
            rxnproteome.conditions.append(condition)
        #Pulling the gene values from the current proteome
        values = {}
        logger.warning("TESTING!")
        for gene in model.genes:
            feature = self.object.search_for_gene(gene.id)
            if feature == None:
                logger.warning("Model gene "+gene.id+" not found in genome of proteome")
            elif feature.id not in self.features:
                logger.warning("Model gene "+gene.id+" in genome but not in proteome")
            else:
                feature = self.features.get_by_id(feature.id)
                for condition in self.conditions:
                    if condition.id not in values:
                        values[condition.id] = {}
                    if condition in feature.values:
                        values[condition.id][gene.id] = feature.values[condition]
        #Computing the reaction level values
        for condition in rxnproteome.conditions:    
            for feature in rxnproteome.features:
                tree = parse_gpr(feature.feature.gene_reaction_rule)[0]
                feature.add_value(condition,compute_gene_score(tree,values[condition.id],default))