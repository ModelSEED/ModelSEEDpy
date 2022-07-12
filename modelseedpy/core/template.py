import logging
logger = logging.getLogger(__name__)

import re  # !!! import never used
import copy  # !!! import never used
from cobra.core.dictlist import DictList
from cobra.core.model import Metabolite, Reaction

class Template():
    def __init__(self):
        self.compounds, self.compcompounds, self.reactions = DictList(), DictList(), DictList()
        
    def convert_template_compound(self,cpdid,index):
        comp_compound = self.compcompounds.get_by_id(cpdid)
        base_id = cpdid.split("_")[0]
        base_compound = self.compounds.get_by_id(base_id)
        compartment = comp_compound.templatecompartment_ref.split("/").pop() + str(index)
        new_id = base_compound.id + str(index)
        met = Metabolite(new_id, formula=base_compound.formula, name=base_compound.name, 
            charge=comp_compound.charge, compartment=compartment)
        met.annotation["sbo"] = "SBO:0000247" #simple chemical - Simple, non-repetitive chemical entity.
        met.annotation["seed.compound"] = base_id
        return met
    
    def convert_template_reaction(self,model,rxnid,index,for_gapfilling=True):   
        template_reaction = self.reactions.get_by_id(rxnid)
        base_id = template_reaction["id"].split("_")[0]
        new_id = template_reaction["id"] + str(index)

        lower_bound = template_reaction["maxrevflux"]
        upper_bound = template_reaction["maxforflux"]

        direction = template_reaction["GapfillDirection"]
        if not for_gapfilling:
            direction = template_reaction["direction"]
        if direction == ">":
            lower_bound = 0
        elif direction == "<":
            upper_bound = 0

        cobra_reaction = Reaction(
            new_id, 
            name=template_reaction["name"], 
            lower_bound=lower_bound, 
            upper_bound=upper_bound)

        object_stoichiometry = {}
        for item in template_reaction["templateReactionReagents"]:
            metabolite_id = item["templatecompcompound_ref"].split("/").pop()
            template_compound = self.compcompounds.get_by_id(metabolite_id)
            compartment = template_compound["templatecompartment_ref"].split("/").pop()
            if compartment == "e":
                metabolite_id = metabolite_id + "0"
            else:
                metabolite_id = metabolite_id + str(index)

            metabolite = model.metabolites.get_by_id(metabolite_id)
            object_stoichiometry[metabolite] = item["coefficient"]

        cobra_reaction.add_metabolites(object_stoichiometry)

        cobra_reaction.annotation["sbo"] = "SBO:0000176" #biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = base_id

        return cobra_reaction
