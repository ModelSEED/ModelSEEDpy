import logging

import re
import copy
import cobra
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg import GapfillingPkg, KBaseMediaPkg

#from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)

class MSBuilder:
    @staticmethod
    def gapfill_model(original_mdl,target_reaction,template,media):
        FBAHelper.set_objective_from_target_reaction(original_mdl,target_reaction)
        model = cobra.io.json.from_json(cobra.io.json.to_json(original_mdl))
        gfp = GapfillingPkg(model)
        gfp.build_package({
            "default_gapfill_templates":[template],
            "gapfill_all_indecies_with_default_templates":1,
            "minimum_obj":0.01,
            "set_objective":1
        })
        kmp = KBaseMediaPkg(model)
        kmp.build_package(media)
        sol=model.optimize()
        gfresults = gfp.compute_gapfilled_solution()
        for rxnid in gfresults["reversed"]:
            rxn = original_mdl.reactions.get_by_id(rxnid)
            if gfresults["reversed"][rxnid] == ">":
                rxn.upper_bound = 100
            else:
                rxn.lower_bound = -100
        
        for rxnid in gfresults["new"]:
            rxn = model.reactions.get_by_id(rxnid)
            rxn = rxn.copy()
            original_mdl.add_reactions([rxn])
            if gfresults["new"][rxnid] == ">":
                rxn.upper_bound = 100
                rxn.lower_bound = 0 
            else:
                rxn.upper_bound = 0
                rxn.lower_bound = -100
        return original_mdl