import logging

import re
from cobra.core import Gene, Metabolite, Model, Reaction

logger = logging.getLogger(__name__)

#Source:https://stackoverflow.com/questions/16699180/how-to-get-molecular-weight-of-a-compound-in-python/45557858
elementmass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
  'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
  'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
  'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
  'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
  'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
  'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
  'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
  'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
  'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
  'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
  'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
  'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
  'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
  'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
  'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
  'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
  'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,'R':0,
  'ZERO': 0} 

class FBAHelper:
    @staticmethod
    def add_drain_from_metabolite_id(model,cpd_id,uptake,excretion,prefix = 'EX_', prefix_name = 'Exchange for '):
        if cpd_id in model.metabolites:
            id = prefix + cpd_id
            cobra_metabolite = model.metabolites.get_by_id(cpd_id)
            object_stoichiometry = {cobra_metabolite : -1}
            drain_reaction = Reaction(id=id, 
                                      name= prefix_name + cobra_metabolite.name, 
                                      lower_bound=-1*uptake, 
                                      upper_bound=excretion)
            drain_reaction.add_metabolites(object_stoichiometry)
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            return drain_reaction
        return None

    @staticmethod
    def set_reaction_bounds_from_direction(reaction,direction,add=0):
        if direction == "<":
            reaction.lower_bound = -100
            if add == 0:
                reaction.upper_bound = 0
        if direction == ">":
            reaction.upper_bound = 100
            if add == 0:
                reaction.lower_bound = 0
        reaction.update_variable_bounds()

    @staticmethod
    def set_objective_from_target_reaction(model,target_reaction,maximize = 1):
        target_reaction = model.reactions.get_by_id(target_reaction)
        sense = "max"
        if maximize == 0:
            sense = "min"
        target_objective = model.problem.Objective(
            1 * target_reaction.flux_expression,
            direction=sense)
        model.objective = target_objective
        return target_reaction

    @staticmethod
    def compute_flux_values_from_variables(model):
        flux_values = {}
        for rxnobj in model.reactions:
            flux_values[rxnobj.id] = {}
            flux_values[rxnobj.id]["reverse"] = rxnobj.reverse_variable.primal
            flux_values[rxnobj.id]["forward"] = rxnobj.forward_variable.primal
        return flux_values
    
    @staticmethod
    def modelseed_id_from_cobra_metabolite(metabolite):
        if re.search('^(cpd\d+)', metabolite.id) != None:
            m = re.search('^(cpd\d+)', metabolite.id)
            return m[1]
        #TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None
        
    def modelseed_id_from_cobra_reaction(reaction):
        if re.search('^(rxn\d+)', reaction.id) != None:
            m = re.search('^(rxn\d+)', reaction.id)
            return m[1]
        #TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None
    
    @staticmethod    
    def metabolite_mw(metabolite):
        mw = 0
        elements = metabolite.elements
        for element in elements:
            if element not in elementmass:
                print("Missing mass for element "+element+" in compound "+metabolite.id+". Element will be ignored when computing MW")
            else:
                mw += elements[element]*elementmass[element]
        return mw
    
    @staticmethod    
    def elemental_mass():
        #Source:https://stackoverflow.com/questions/16699180/how-to-get-molecular-weight-of-a-compound-in-python/45557858
        return elementmass    