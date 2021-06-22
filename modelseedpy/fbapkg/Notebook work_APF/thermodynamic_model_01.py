# -*- coding: utf-8 -*-

from __future__ import absolute_import

#import logging
import cobra
bigg_model_path = '..\COBRA function scripts\e_coli_core metabolism from BiGG.json'
model = cobra.io.load_json_model(bigg_model_path)



#Adding a few exception classes to handle different types of errors
class FeasibilityError(Exception):
    """Error in FBA formulation
    """
    if Exception == '':
        pass
    
    if Exception == '':
        pass
    

#Base class for FBA packages
class BaseFBAPkg:
    def __init__(self, model, name, object, object_type = 'reaction', variable_types = {'concentration': 'float', 'lnconc': 'float'}, constraint_types = {'concentration': 'float'}, parent=None):
        '''Intantiate the model
        'model' (COBRA object): The COBRApy model object
        'name' (Python object, string): The name of the model
        'object' (COBRA object): The name of a COBRA reaction or metabolite, although, the former is the essential intention of the API  
        'object_type' (Python object, string): A description of the COBRA object, which is used to apply the pertinent code for the passed object
        'variable_types' (Python object, dictionary): The types and variables examples for the model variables
        'constraint_types' (Python object, dictionary): The names and values of the model constraints
        'parent' (Python object, boolean): The categorical description of the model 
        '''
        import re
        
        self.model = model
        self.name = name
        self.childpkgs = dict()
        self.parentpkg = parent
        self.constraints = dict()
        self.variables = dict()
        self.parameters = dict()
        self.variable_types = variable_types
        self.constraint_types = constraint_types
        for kind in variable_types:
            self.variables[kind] = dict()
        for kind in constraint_types:
            self.constraints[kind] = dict()
    

        from scipy.constants import calorie
        import pandas
        import re

        #import modelseedpy      
        #modelseed_path = '..\..\..\Biofilm growth code\GSWL code\ModelSEEDDatabase'
        #modelseed = modelseedpy.biochem.from_local(modelseed_path)
        
        # import the Gibbs data from the TMFA supplementary excel file
        reactions_data = pandas.read_excel('Supplementary file.xls', sheet_name = 'Reactions')
        reactions_data.columns = reactions_data.iloc[0]
        reactions_data.drop(reactions_data.index[0], axis = 0, inplace = True)
        reactions_data.head()
        reactions_dict = {}
        for index, row in reactions_data.iterrows():
            reaction_abbreviation = reactions_data.at[index, 'iJR904 Abbreviation']
            reaction_name = reactions_data.at[index, 'Reaction name']
            reaction_gibbs = reactions_data.at[index, 'Estimated Gibbs free energy change of reaction (kcal/mol)']
            
            reactions_dict[reaction_abbreviation] = {'name': reaction_name, 'gibbs': reaction_gibbs}
        
        self.thermo_reactions = reactions_dict
        
        compounds_data = pandas.read_excel('Supplementary file.xls', sheet_name = 'Compounds')
        compounds_dict = {}
        for index, row in compounds_data.iterrows():
            compound_abbreviation = compounds_data.at[index, 'iJR904 Abbreviation']
            compound_name = compounds_data.at[index, 'Compound Name']
            compound_charge = compounds_data.at[index, 'Charge at pH 7*']
            try:
                compound_gibbs = float(compounds_data.at[index, 'Estimated Gibbs free energy of formation (kcal/mol)']) / calorie
            except:
                compound_gibbs = 0
            
            compounds_dict[compound_abbreviation] = {'name': compound_name, 'gibbs': compound_gibbs, 'charge': compound_charge}
        
        self.thermo_compounds = compounds_dict      
        
        # expand the metabolic thermodynamic data of the reaction in the argument 
        if object_type == 'reaction':
            for metabolite in object.metabolites:
                try:
                    tmfa_name = re.sub('(_\w$)', '', metabolite.id)
                    self.thermo_compounds[tmfa_name]['stoichiometry'] = object.metabolites[metabolite]
                except:
                    print('ERROR: The metabolite {} is unrepresented in the thermodynamic database.'. format(metabolite.id))

        elif object_type == 'metabolite':
            tmfa_name = re.sub('(_\w)', '', metabolite.id)
            self.thermo_compounds[tmfa_name] = object
        else:
            print('ERROR: The object_type is not compatible with this API.')
            
                        
    def validate_parameters(self, params, required, defaults):
        '''Validate the passed parameters with the required parameters and update the instance parameters with the default parameters
        'params' (Python object, dictionary): The parameters and values of the model
        'required' (Python object, list): The required parameters for the model 
        'defaults' (Python object, dictionary): The default parameters and values for the model 
        '''
        missing_parameters = set(required) - set(params)
        missing_string = ', '.join(list(missing_parameters))
        if len(missing_parameters) > 1:
            raise ValueError('The required parameters < {} > are missing!'.format(missing_string))
        elif len(missing_parameters) == 1:
            raise ValueError('The required parameter < {} > is missing!'.format(missing_string))
        
        self.parameters = params
        for key in defaults:
            if key not in params:
                self.parameters[key] = defaults[key]
        
        
    def build_variable(self, kind, lower_bound, upper_bound, vartype, object, object_type = 'reaction'):
        '''Create variables of the specified type in the COBRA model
        'kind' (Python object, string): The variable type within which variables will be created
        'lower_bound' (Python object, float): The lower bound value for the added variable
        'upper_bound' (Python object, float): The upper bound value for the added variable
        'vartype' (Python object, string): The variable type as either 'continuous', 'integer', or 'binary' 
        'object' (COBRA object): The COBRA entity into which a variable will be build
        'object_type' (Python object): The variable type of the COBRA object 
        '''
        # assign a variable name based upon the passed arguments
        self.variable_types[kind] = vartype  
        if object_type == "none":
            count = len(self.variables[type])
            name = str(count + 1)
        elif object_type == "string":
            name = object
        elif object_type in ['reaction', 'metabolite']:
            name = object.id
            
        # add an optlang variable, when the variable is undefined
        #raise ValueError('The object name {} is not recognized by your model'.format(missing_string))
        if name in self.variables[kind]:
            print('ERROR: The {} constraint already exists with a value of {}.'.format(name, self.variables[kind][name]))
            variable_definition = self.variables[kind][name] 
        
        elif name not in self.variables[kind]:
            variable_name = name + "_" + kind
            variable_definition = self.model.problem.Variable(name = variable_name, lb = lower_bound, ub = upper_bound, type = vartype)
            self.model.add_cons_vars(variable_definition)
            self.variables[kind][name] = variable_definition 
            
        return variable_definition
        
        
    def build_constraint(self, constraint_expression, kind, lower_bound, upper_bound, coef, object, object_type = 'reaction'):
        '''Create constraints for the COBRA model
        'kind' (Python object, string): The type of the constraint that will be created 
        'lower_bound' (Python object, float): The lower bound value for the added constraint
        'upper_bound' (Python object, float): The upper bound value for the added constraint
        'coef' (Python object, dictionary): The set of coefficients that define the COBRA model
        'object' (Python object, string): The variable name when the name is defined
        'object_type' (Python object, string): The variable type of the COBRA object 
        '''
        from optlang.symbolics import Zero
        
        # assign a constraint name based upon the passed arguments
        if object_type == "none":
            count = len(self.constraints[type])
            name = str(count + 1)
        elif object_type == "string":
            name = object
        elif object_type in ['reaction', 'metabolite']:
            self.constraint_types[kind] = coef
            name = object.id
                   
        # add an optlang constraint, when the constraint is undefined 
        if name in self.constraints[kind]:
            print('ERROR: The {} constraint has already been added with a value of {}.'.format(name, self.constraints[kind][name]))
            constraint_definition = None
        elif name not in self.constraints[kind]:
            constraint_name = '{}_{}'.format(name, kind)
            constraint_definition = self.model.problem.Constraint(expression = Zero, lb = lower_bound, ub = upper_bound, name = constraint_name)
            self.model.add_cons_vars(constraint_definition)
            self.model.solver.update()
            self.constraints[kind][name] = constraint_definition
            if len(coef) > 0:
                self.constraints[kind][name].set_linear_coefficients(coef)
                
            self.model.solver.update()

        
        return constraint_definition
    
    
    def all_variables(self):
        '''Quantify the variables in the class
        '''
        vars = {}
        for child in self.childpkgs:
            for kind in child.variables:
                vars[kind] = child.variables[kind]

        for kind in self.variables:
            vars[kind] = self.variables[kind]
        
        return vars
    
    
    def all_constraints(self):
        '''Quantify the constraints in the class
        '''
        consts = {}
        for child in self.childpkgs:
            for kind in child.constraints:
                consts[kind] = child.constraints[kind]

        for kind in self.constraints:
            consts[kind] = self.constraints[kind]
            
        return consts
    
    
    def revert_to_original(self, cobra_model_path):
        '''Remove changed components of the COBRA model
        'model_path' (Python object, string): The path string of the COBRA model
        '''               
        global model
        # remove added variables and constants from the model by re-uploading the COBRA model  
        model = cobra.io.load_json_model(cobra_model_path)
        
        return model
    

    def write_lp_file(self, model, export_filename = 'test'):
        '''Export the LP file of the COBRA model
        'model' (COBRA object): The COBRA model that is expanded through this API
        'export_filename' (Python object, string): The string of the lp file that will be exported
        '''
        from datetime import date
        import os
        
        import_iteration = 0
        lp_filename = '{}_{}_{}.lp'.format(date.today(), export_filename, import_iteration)
        while os.path.exists(lp_filename):
            import_iteration += 1
            lp_filename = '{}_{}_{}.lp'.format(date.today(), export_filename, import_iteration)
            
        with open(lp_filename, 'w') as lp_file:
            lp_file.write(str(model.solver))
         
        
# ---------------------------------------------- Revbin -------------------------------------------------

# The base class for FBA packages is inherited
class RevBinPkg(BaseFBAPkg):
    def __init__(self, model, object):
        '''Redefining the inherited __init__ function
        'model' (COBRA object): The COBRApy FBA model
        '''
        BaseFBAPkg.__init__(self, model = model, object = object, object_type = 'reaction', name = "reversible binary", variable_types = {"revbin":"reaction", "forv":"reaction", "revv":"reaction"}, constraint_types = {"revbinF":"reaction", "revbinR":"reaction"})

        # define the variables that are used in the constraints
        BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", object = object)
        
        
    def build_constraint(self, object):
        '''Build constraints through the inherited function and the calculated coefficient fluxes
        'object' (Python object, string): The variable name when the name is defined
        '''
        from optlang.symbolics import Zero
        
        # define the constraints of the system
        coef = {self.variables['revbin'][object.id]:-1000, object.forward_variable:1}
        built_forward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinF", lower_bound = None, upper_bound = 0, coef = coef, object = object)
        
        coef = {self.variables["revbin"][object.id]:1000, object.reverse_variable:1}
        built_backward_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = "revbinR", lower_bound = None, upper_bound = 1000, coef = coef, object = object)
        
        return built_backward_constraint
    
    
    def build_package(self, filter):
        '''Build variables and constraints through the inherited function
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        '''
        added_constraints = []
        for reaction in model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id not in filter:
                BaseFBAPkg.build_variable(self, kind = "revbin", lower_bound = 0, upper_bound = 1, vartype = "binary", object = reaction)
                self.build_constraint(object = reaction)
                added_constraints.append(reaction.id)
                
        if len(added_constraints) > 0:
            constraints = ', '.join(added_constraints)
            print('The constraints for {} were added to the model.'.format(constraints))
        else:
            print('ERROR: No reactions were added to the model.')
        

                
# ------------------------------------------ Simple Thermo package ------------------------------------------                

# The base class for FBA packages is inherited
class SimpleThermoPkg(RevBinPkg):
    def __init__(self, model, object):
        '''Redefining the inherited __init__ function
        'model' (COBRA object): The COBRApy FBA model
        '''        
        from scipy.constants import physical_constants, kilo, R
        from numpy import log as ln 
        import re

        F = physical_constants['Faraday constant'][0]
                
        # execute the parent __init__ and arbitrarily assign potentials to each metabolite 
        BaseFBAPkg.__init__(self, model = model, name = "simple thermo", object = object, object_type = 'reaction', variable_types = {"potential":"metabolite", 'concentration_potential': 'metabolite', "revbin":"reaction"}, constraint_types = {"thermo":"reaction"})
        
        # parameterize the initial chemical concentrations 
        self.validate_parameters(params = self.parameters, required = [], defaults = {
            "default_conc_range": [0.001, 20],  # an arbitrary range
            "custom_concentration_constraints": {"glc-D": [10,10],  #cpd00027,  
                                                'co2': [20, 24], #cpd00011, as bicarbonate, E. B. Brown and Richard L. Clancy, 1964
                                                'h': [0.0053, 0.0053], #cpd00067, Battaglia, Hellegers, & Seeds, 1965
                                                'o2': [0.672, 0.672] #cpd00007, 0.3 mL / dL serum O2 concentration
                                                 }
        })
        
        # calculate the total energy of a reaction based upon the values for each constitutent metabolite in the reaction 
        delta_g = 0
        sum_concentration_potential = 0
        sum_electro_potential = 0
        temperature = 25 # degrees kelvin  
        self.variables['concentration_potential'] = {}
        self.variables['electro_potential'] = {}
        self.parameters['activity_coefficient'] = {}
        self.variables['total_energy'] = {}
        for metabolite in object.metabolites:
            tmfa_name = re.sub('(_.)', '', metabolite.id)
            
            # calculate the concentration range for the chemical specie
            if tmfa_name not in self.parameters['custom_concentration_constraints']:
                concentration_range = self.parameters['default_conc_range']
                self.thermo_compounds[tmfa_name]['concentration'] = sum(concentration_range) / len(concentration_range)
            elif tmfa_name in self.parameters['custom_concentration_constraints']:
                concentration_range = self.parameters['custom_concentration_constraints'][tmfa_name]
                self.thermo_compounds[tmfa_name]['concentration'] = sum(concentration_range) / len(concentration_range)
                
                
            # calculate the electrochemical potential term
            if metabolite.compartment == 'c':
                ph_gradient = 2 # arbitrary value
            elif metabolite.compartment == 'e':
                ph_gradient = 0 # by defintion of a zero ph gradient between the metabolite compartment and the extracellular compartment
                
            psi_electro_potential = 33.33 * ph_gradient - 143.33  # millivolts, equation 9 from the TMFA paper 
            self.variables['electro_potential'][metabolite.id] = psi_electro_potential * F * self.thermo_compounds[tmfa_name]['charge'] * kilo

            
            # calculate the concentration potential term from the average concentration
            if metabolite.id not in self.parameters['activity_coefficient']:
                activity_coefficient = 1
            elif metabolite.id in self.parameters['activity_coefficient']:
                activity_coefficient = self.parameters['activity_coefficient'][metabolite.id]
            
            self.variables['concentration_potential'][metabolite.id] = self.thermo_compounds[tmfa_name]['stoichiometry'] * ln(self.thermo_compounds[tmfa_name]['concentration'] * activity_coefficient)
            
            
            # sum the all energetic descriptions of each metabolite in a reaction
            sum_concentration_potential += self.variables['concentration_potential'][metabolite.id]
            sum_electro_potential += self.variables['electro_potential'][metabolite.id]
            delta_g += self.thermo_compounds[tmfa_name]['gibbs']                
            

        # calculation of the total energetic potential of the reaction based upon the metabolite calculations
        self.variables['total_energy'][object.id] = delta_g + R * temperature * sum_concentration_potential + sum_electro_potential  
        
        # inherit the RevBinPkg instance variables
        self.childpkgs["reverse_binary"] = RevBinPkg(model = model, object = model.reactions.get_by_id('PFK'))
            
            
    def build_constraint(self, object):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'object' (Python object, string): The variable name when the name is defined
        '''
        from optlang.symbolics import Zero
        import re
        
        delta_g = 0 #sum(st(i,j)*p(j))
        
        for metabolite in object.metabolites:
            tmfa_name = re.sub('(_.)', '', metabolite.id)
            delta_g += self.thermo_compounds[tmfa_name]['stoichiometry'] * (self.thermo_compounds[tmfa_name]['gibbs'] + self.variables['concentration_potential'][metabolite.id] + self.variables['electro_potential'][metabolite.id])
            
        if object.reversibility:
            binary = 0
        elif not object.reversibility:
            binary = 1
        else:
            print('ERROR: The reaction object possesses unpredictable data structure.')

        # create thermodynamic constraints and the associated variables
        for metabolite in self.model.metabolites:
            BaseFBAPkg.build_variable(self, kind = "potential", lower_bound = 0, upper_bound = 1000, vartype = "continuous", object = metabolite)
            
        '''different_constraints = set(model.constraints) - set(self.constraints)
        for constraint in different_constraints:
            self.constraints[constraint.name] = constraint.expression'''

        kind = 'thermo'
        constraint_name = '{}_{}'.format(object.id, kind)
        if object.id not in self.constraints[kind]:
            coef = {self.childpkgs["reverse_binary"].variables["revbin"][object.id] : 1000}
            BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = kind, lower_bound = 0, upper_bound = 1000, coef = coef, object = object)
        elif object.id not in self.constraints[kind]:
            print('ERROR: The {} constraint already exists in the model with the a value of {}.'.format(constraint_name, self.constraints[constraint_name]))


        
    def build_package(self, filter = model.reactions):
        '''Build variables and constraints through the inherited function
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        '''
        RevBinPkg.build_package(self, filter = filter)
                
                
# ------------------------------------------ Full Thermo package ------------------------------------------

# The base class for FBA packages is inherited
class FullThermoPkg(SimpleThermoPkg):
    def __init__(self, model, object):
        from numpy import log as ln

        '''Redefining the inherited __init__ function and importing thermodynamic data
        'model' (COBRA object): The COBRApy FBA model
        '''                 
        # execute the base __init__ file
        BaseFBAPkg.__init__(self, model = model, object = object, name = "full thermo", variable_types = {"lnconc":"metabolite", 'thermo': 'reactions', "revbin":"reaction"}, constraint_types = {"potentialc":"metabolite", 'thermo': 'metabolite'})    
        
    '''def build_variable(self, object): #, modelseed = modelseed):   
        #flux variability analysis? 
        
        from numpy import log as ln
        if object.id in self.parameters["custom_concentration_constraints"]:
            lb = ln(self.parameters["custom_concentration_constraints"][object.id][0])
            ub = ln(self.parameters["custom_concentration_constraints"][object.id][1])
        else:
            lb = ln(self.parameters["default_min_conc"])
            ub = ln(self.parameters["default_max_conc"])
            
        return BaseFBAPkg.build_variable(self, type = "lnconc", lower_bound = lb, upper_bound = ub, vartype = "continuous", object = object)'''
    
    
    def build_constraint(self, model, object, kind = 'thermo'): #, modelseed = modelseed):
        '''Build constraints through the inherited function and the calculated variable coeffiients 
        'object' (Python object, string): The variable name when the name is defined
        Notes - Equation 14 in the TMFA paper, with the addition of the (charge * compartment_potential) term?
        '''
        from optlang.symbolics import Zero
        
        # calculate the parameters for the constraint expression calculation               
        constant = 20  # arbitrary value
        
        constraint_names = []
        for constraint in model.constraints:
            constraint_names.append(constraint.name)
        
        constraint_name = '{}_{}'.format(object.id, kind)
        if constraint_name not in constraint_names:
            self.childpkgs["reverse_binary"] = RevBinPkg(model = model, object = object)
            coef = {self.childpkgs["reverse_binary"].variables["revbin"][object.id] : 1000}
            built_constraint = BaseFBAPkg.build_constraint(self, constraint_expression = Zero, kind = kind, lower_bound = constant, upper_bound = constant, coef = coef, object = object)
            
        else:
            print('ERROR: The {} contraint already exists in the model: {}, ub: {}, lb: {}'.format(constraint_name, model.constraints[constraint_name], object.upper_bound, object.lower_bound))
            object.upper_bound = constant
            object.lower_bound = constant
            #model.constraints[constraint_name] = constraint_expression
            print('The {} contraint is updated in the model: {}, ub: {}, lb: {}'.format(constraint_name, model.constraints[constraint_name], object.upper_bound, object.lower_bound))
            built_constraint = None
            
        return built_constraint

        
    def build_package(self):
        '''Create the final model package
        'filter' (Python object, list): The accepted list of reactions that will be built into the model
        Notes - The concentrations are expressed in millimolar
        '''           
        filter = self.constraints
        SimpleThermoPkg.build_package(filter)
        
        # The concentration variable and potential constraint are built
        for metabolite in self.model.metabolites:
            BaseFBAPkg.build_variable(self, kind = "thermo", lower_bound = 0, upper_bound = 1000, vartype = "continuous", object = metabolite)
            
        for reaction in self.model.reactions:
            # Unfiltered reactions are constructed through the aforementioned functions
            if reaction.id not in filter:
                self.build_constraint(model, object = reaction)