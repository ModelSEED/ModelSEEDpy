# -*- coding: utf-8 -*-

from __future__ import absolute_import

from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

#Base class for FBA packages
class BilevelPkg(BaseFBAPkg):
    def __init__(self,model):
        BaseFBAPkg.__init__(self,model,"reaction use",{"dualconst":"string","dualub":"string","duallb":"string"},{"dualvar":"string","objective":"string"})
    
    def build_package(self,filter = None,binary_variable_count = 0):
        self.validate_parameters(parameters, [], {
            "binary_variable_count":binary_variable_count
        });
        coefficients = {}
        obj_coef = {}
        obj = self.model.solver.objective
        #Using the JSON calls because get_linear_coefficients is REALLY slow
        #Creating new objective coefficient and bound variables
        bound_variables = {}
        reactions = self.model.reactions
        if binary_variable_count > 0:
            for reaction in self.model.reactions: 
                var = self.build_variable("fluxcomp",reaction,None)
                const = self.build_constraint("fluxcomp",reaction,None,None,None)
        #Retreiving model data with componenent flux variables
        mdldata = self.model.solver.to_json()
        consthash = {}
        for const in mdldata["constraints"]:
            consthash[const["name"]] = const
        constraints = list(self.model.solver.constraints)
        variables = list(self.model.solver.variables)
        objterms = obj.get_linear_coefficients(variables)
        #Adding binary variables and constraints which should not be included in dual formulation
        if binary_variable_count > 0:
            for variable in self.variables["fluxcomp"]: 
                var = self.build_variable("bfluxcomp",variable,None)
                const = self.build_constraint("bfluxcomp",variable,None,None,None)                
        #Now implementing dual variables and constraints
        varhash = {}
        for var in variables:
            varhash[var.name] = var
        for const in constraints:
            var = self.build_variable("dualconst",const,obj_coef)
            if var != None and const.name in consthash and "expression" in consthash[const.name] and "args" in consthash[const.name]["expression"]:
                for item in consthash[const.name]["expression"]["args"]:
                    if "args" in item and len(item["args"]) >= 2 and item["args"][1]["name"] in varhash:
                        if varhash[item["args"][1]["name"]] not in coefficients:
                            coefficients[varhash[item["args"][1]["name"]]] = {}
                        coefficients[varhash[item["args"][1]["name"]]][var] = item["args"][0]["value"]
        for var in variables:
            if var.type == "continuous":
                dvar = self.build_variable("duallb",var,obj_coef)
                if dvar != None:
                    if var not in coefficients:
                        coefficients[var] = {}
                    coefficients[var][dvar] = 1
                dvar = self.build_variable("dualub",var,obj_coef)
                if dvar != None:
                    if var not in coefficients:
                        coefficients[var] = {}
                    coefficients[var][dvar] = 1
                self.build_constraint("dualvar",var,obj,objterms,coefficients)
        self.build_constraint("objective",None,obj,objterms,obj_coef)
        #Add binary variables and constaints which should not be part of the formulation
    
    def build_variable(self,type,object,obj_coef):
        if type == "dualconst":
            lb = -1000000
            ub = 1000000
            coef = 0
            if object.lb == None:
                lb = 0
                coef = object.ub
            if object.ub == None:
                ub = 0
                coef = object.lb
            var = BaseFBAPkg.build_variable(self,type,lb,ub,"continuous",object.name)
            obj_coef[var] = coef
            return var
        if type == "dualub":
            var = BaseFBAPkg.build_variable(self,type,0,1000000,"continuous",object.name)
            obj_coef[var] = object.ub
            return var
        if type == "duallb":
            var = BaseFBAPkg.build_variable(self,type,-1000000,0,"continuous",object.name)
            obj_coef[var] = object.lb
            return var
        if type == "flxcmp" and self.parameters["binary_variable_count"] > 0:
            denominator = 2**self.parameters["binary_variable_count"]-1
            coefs = [{},{}]
            for i in range(0,self.parameters["binary_variable_count"]):
                value = 2**i
                if object.lower_bound < 0:
                    var = BaseFBAPkg.build_variable(self,"rflxcmp"+str(i),0,-1*value*object.lower_bound/denominator,"continuous",object.id)
                    coefs[0][var] = -1
                if object.upper_bound > 0:
                    var = BaseFBAPkg.build_variable(self,"fflxcmp"+str(i),0,value*object.upper_bound/denominator,"continuous",object.id)
                    coefs[1][var] = -1
            if object.lower_bound < 0:
                #flux - flux_comp_0 - flux_comp_n = 0 - restriction of reverse fluxes by component fluxes
                coefs[0][object.reverse_variable] = 1
                BaseFBAPkg.build_constraint(self,"rflxcmpc",0,0,coefs[0],object.id)
            if object.upper_bound > 0:
                #flux - flux_comp_0 - flux_comp_n = 0 - restriction of forward fluxes by component fluxes
                coefs[1][object.forward_variable] = 1
                BaseFBAPkg.build_constraint(self,"fflxcmpc",0,0,coefs[1],object.id)
            return None
        if type == "bflxcmp" and self.parameters["binary_variable_count"] > 0:
            for i in range(0,self.parameters["binary_variable_count"]):
                if object.lower_bound < 0:
                    var = BaseFBAPkg.build_variable(self,"brflxcmp"+str(i),0,1,"binary",object.id)
                    othervar = self.variables["rflxcmp"+str(i)][object.id]
                    BaseFBAPkg.build_constraint(self,"brflxcmpc"+str(i),None,0,{othervar:1,var:-1000},object.id)
                if object.upper_bound > 0:
                    var = BaseFBAPkg.build_variable(self,"bfflxcmp"+str(i),0,1,"binary",object.id)
                    othervar = self.variables["fflxcmp"+str(i)][object.id]
                    BaseFBAPkg.build_constraint(self,"bfflxcmpc"+str(i),None,0,{othervar:1,var:-1000},object.id)
            return None
                    
    def build_constraint(self,type,object,objective,objterms,coefficients):
        if type == "dualvar":
            coef = {}
            lb = 0
            ub = 0
            objsign = 1
            if objective.direction == "min":
                objsign = -1
            if object in objterms:
                lb = objterms[object]
                ub = objterms[object]
            if object in coefficients:
                for var in coefficients[object]:
                    coef[var] = coefficients[object][var]
            if object.lb == 0:
                ub = None
            elif object.ub == 0:
                lb = None
            return BaseFBAPkg.build_constraint(self,type,lb,ub,coef,object.name)
        elif type == "objective":
            coef = {}
            objsign = 1
            if objective.direction == "min":
                objsign = -1
            for var in objterms:
                coef[var] = objsign*objterms[var]
            for dvar in coefficients:
                coef[dvar] = -1*coefficients[dvar]
            return BaseFBAPkg.build_constraint(self,type,0,0,coef,"dualobjconst")