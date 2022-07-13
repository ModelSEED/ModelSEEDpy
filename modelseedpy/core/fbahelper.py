from __future__ import absolute_import

import logging
from chemicals import periodic_table
import re
import time
from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util import solver as sutil
from modelseedpy.biochem import from_local
from scipy.odr.odrpack import Output
#from Carbon.Aliases import false

logger = logging.getLogger(__name__)

elementmass = {}
for element in periodic_table:
    elementmass[element.symbol] = element.MW


class FBAHelper:

    @staticmethod
    def add_autodrain_reactions_to_community_model(model,auto_sink = ["cpd02701", "cpd11416", "cpd15302"]):
        #Adding missing drains in the base model
        drain_reactions = []
        for metabolite in model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(metabolite)
            if msid in auto_sink:
                if msid != "cpd11416" or metabolite.compartment == "c0":
                    met_id = metabolite.id
                    if all(rxn not in model.reactions for rxn in [f"EX_{met_id}", f"DM_{met_id}", f"SK_{met_id}"]):
                        drain_reaction = FBAHelper.add_drain_from_metabolite_id(model,metabolite.id,0,100,"DM_")
                        if drain_reaction != None:
                            logger.info("Adding "+met_id+" DM")
                            drain_reactions.append(drain_reaction)
        model.add_reactions(drain_reactions)
        
    @staticmethod
    def add_drain_from_metabolite_id(model, cpd_id, uptake, excretion, prefix='EX_', prefix_name='Exchange for '):
        """
        :param model:
        :param cpd_id:
        :param uptake:
        :param excretion:
        :param prefix:
        :param prefix_name:
        :return:
        """
        if cpd_id in model.metabolites:
            cobra_metabolite = model.metabolites.get_by_id(cpd_id)
            drain_reaction = Reaction(id=f'{prefix}{cpd_id}',
                                      name=prefix_name + cobra_metabolite.name,
                                      lower_bound=-1*uptake, 
                                      upper_bound=excretion)
            drain_reaction.add_metabolites({cobra_metabolite : -1})
            drain_reaction.annotation["sbo"] = 'SBO:0000627'    
            #model.add_reactions([drain_reaction])
            return drain_reaction
        return None
    
    @staticmethod
    def set_reaction_bounds_from_direction(reaction, direction, add=0):
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
    def set_objective_from_target_reaction(model,target_reaction,minimize = False):
        target_reaction = model.reactions.get_by_id(target_reaction)
        sense = "max"
        if minimize:
            sense = "min"
        model.objective = model.problem.Objective(
            1 * target_reaction.flux_expression,
            direction=sense)
        return target_reaction
    
    @staticmethod
    def modelseed_id_from_cobra_metabolite(metabolite):
        if re.search('^(cpd\d+)', metabolite.id):
            m = re.search('^(cpd\d+)', metabolite.id)
            return m[1]
        #TODO: should check to see if ModelSEED ID is in the annotations for the compound
        else:
            return None
        
    @staticmethod
    def modelseed_id_from_cobra_reaction(reaction):
        if re.search('^(rxn\d+)', reaction.id):
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
        return elementmass
    
    @staticmethod
    def get_modelseed_db_api(modelseed_path):
        return from_local(modelseed_path)
    
    @staticmethod
    def is_ex(reaction):
        # TODO: check for SBO
        if len(reaction.id) > 3 and reaction.id[0:3] in ["EX_", "DM_", "SK_"]:
            return True
        return False

    @staticmethod
    def is_biomass(reaction):
        # TODO: check for SBO
        return reaction.id[0:3] == "bio"
    
    @staticmethod
    def exchange_hash(model):
        exchange_hash = {}
        for reaction in model.reactions:
            if len(reaction.metabolites) == 1:
                for metabolite in reaction.metabolites:
                    (base,comp,index) = FBAHelper.parse_id(metabolite)
                    #exchange_hash[base][comp]

    @staticmethod
    def find_reaction(model,stoichiometry):
        output = FBAHelper.stoichiometry_to_string(stoichiometry)
        atpstring = output[0]
        rxn_hash = FBAHelper.rxn_hash(model)
        if atpstring in rxn_hash:
            return rxn_hash[atpstring]
        return None
    
    @staticmethod
    def msid_hash(model): 
        output = {}
        for cpd in model.metabolites:
            msid = FBAHelper.modelseed_id_from_cobra_metabolite(cpd)
            if msid != None:
                if msid not in output:
                    output[msid] = []
                output[msid].append(cpd)
        return output
    
    @staticmethod
    def rxn_hash(model): 
        output = {}
        for rxn in model.reactions:
            strings = FBAHelper.stoichiometry_to_string(rxn.metabolites)
            output[strings[0]] = [rxn,1]
            output[strings[1]] = [rxn,-1]
        return output
    
    @staticmethod
    def rxn_compartment(reaction): 
        compartments = list(reaction.compartments)
        if len(compartments) == 1:
            return compartments[0]
        cytosol = None
        othercomp = None
        for comp in compartments:
            if comp[0:1] != "e":
                if comp[0:1] == "c":
                    cytosol = comp
                else:
                    othercomp = comp
        if othercomp is not None:
            return othercomp
        return cytosol
    
    @staticmethod
    def stoichiometry_to_string(stoichiometry):
        reactants = []
        products = []
        for met in stoichiometry:
            coef = stoichiometry[met]
            if not isinstance(met, str):
                if FBAHelper.modelseed_id_from_cobra_metabolite(met) == "cpd00067":
                    met = None
                else:
                    met = met.id
            if met != None:
                if coef < 0:
                    reactants.append(met)
                else:
                    products.append(met)
        reactants.sort()
        products.sort()
        return ["+".join(reactants)+"="+"+".join(products),"+".join(products)+"="+"+".join(reactants)]
    
    @staticmethod
    def add_atp_hydrolysis(model,compartment):
        #Searching for ATP hydrolysis compounds
        coefs = {"cpd00002":[-1,compartment],"cpd00001":[-1,compartment],"cpd00008":[1,compartment],"cpd00009":[1,compartment],"cpd00067":[1,compartment]}
        msids = ["cpd00002","cpd00001","cpd00008","cpd00009","cpd00067"]
        stoichiometry = {}
        id_hash = FBAHelper.msid_hash(model)
        for msid in msids:
            if msid not in id_hash:
                logger.warning("Compound "+msid+" not found in model!")
                return None
            else:
                for cpd in id_hash[msid]:
                    if cpd.compartment == coefs[msid][1]:
                        stoichiometry[cpd] = coefs[msid][0]
        output = FBAHelper.find_reaction(model,stoichiometry)
        if output != None and output[1] == ">":
            return {"reaction":output[0],"direction":">","new":False}
        cobra_reaction = Reaction("rxn00062_"+compartment, 
                                  name="ATP hydrolysis", 
                                  lower_bound=0, 
                                  upper_bound=1000)
        cobra_reaction.annotation["sbo"] = "SBO:0000176" #biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = "rxn00062"
        cobra_reaction.add_metabolites(stoichiometry)
        model.add_reactions([cobra_reaction])
        return {"reaction":cobra_reaction,"direction":">","new":True}
        
    @staticmethod
    def parse_id(object):
        if re.search('(.+)_([a-z])(\d+)$', object.id) != None:
            m = re.search('(.+)_([a-z])(\d+)$', object.id)
            return (m[1],m[2],int(m[3]))
        return None
    
    @staticmethod
    def medianame(media):
        if media == None:
            return "Complete"
        return media.id
    
    @staticmethod
    def validate_dictionary(dictionary,required_keys,optional_keys):
        for item in required_keys:
            if item not in dictionary:
                raise ValueError('Required key '+item+' is missing!')
        for key in optional_keys:
            if key not in dictionary:
                dictionary[key] = defaults[key]
        return dictionary
