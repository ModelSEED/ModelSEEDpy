import logging
import itertools
import cobra
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.fbapkg.gapfillingpkg import default_blacklist

logger = logging.getLogger(__name__)


#Adding a few exception classes to handle different types of errors
class GapfillingError(Exception):
    """Error in model gapfilling"""
    pass


class MSGapfill:

    def __init__(self, model, default_gapfill_templates=[], default_gapfill_models=[],
                 test_conditions=[], reaction_scores={}, blacklist=[]):
        self.auto_sink = ["cpd02701", "cpd11416", "cpd15302"]
        self.model = model
        self.gfmodel = None
        self.model_penalty = 1
        self.default_gapfill_models = default_gapfill_models
        self.default_gapfill_templates = default_gapfill_templates
        self.gapfill_templates_by_index = {}
        self.gapfill_models_by_index = {}
        self.gapfill_all_indecies_with_default_templates = True
        self.gapfill_all_indecies_with_default_models = True
        self.blacklist = default_blacklist
        for rxnid in blacklist:
            if rxnid not in self.blacklist:
                self.blacklist.append(rxnid)
        self.lp_filename = None
        self.test_condition_iteration_limit = 10
        self.test_conditions = test_conditions
        self.reaction_scores = reaction_scores
        self.solutions = {}
        
    def run_gapfilling(self, media=None, target="bio1", minimum_obj=0.01, binary_check=False):
        self.model.objective = self.model.problem.Objective(
            self.model.reactions.get_by_id(target).flux_expression, 
            direction='max'
        )
        self.gfmodel = cobra.io.json.from_json(cobra.io.json.to_json(self.model))
        pkgmgr = MSPackageManager.get_pkg_mgr(self.gfmodel)
        pkgmgr.getpkg("GapfillingPkg").build_package({
            "auto_sink": self.auto_sink,
            "model_penalty": self.model_penalty,
            "default_gapfill_models": self.default_gapfill_models,
            "default_gapfill_templates": self.default_gapfill_templates,
            "gapfill_templates_by_index": self.gapfill_templates_by_index,
            "gapfill_models_by_index": self.gapfill_models_by_index,
            "gapfill_all_indecies_with_default_templates": self.gapfill_all_indecies_with_default_templates,
            "gapfill_all_indecies_with_default_models": self.gapfill_all_indecies_with_default_models,
            "default_excretion":100,
            "default_uptake":100,
            "minimum_obj": minimum_obj,
            "blacklist": self.blacklist,
            "reaction_scores": self.reaction_scores,
            "set_objective": 1
        })
        pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        if self.lp_filename:
            with open(self.lp_filename, 'w') as out:
                out.write(str(self.gfmodel.solver))
        sol = self.gfmodel.optimize()
        logger.debug('gapfill solution objective value %f (%s) for media %s', sol.objective_value, sol.status, media)

        if sol.status != 'optimal':
            logger.warning("No solution found for %s", media)
            return None

        if media not in self.solutions:
            self.solutions[media] = {}
        self.solutions[media][target] = pkgmgr.getpkg("GapfillingPkg").compute_gapfilled_solution()
        if self.test_conditions:
            self.solutions[media][target] = pkgmgr.getpkg("GapfillingPkg").run_test_conditions(self.test_conditions,self.solutions[media][target],self.test_condition_iteration_limit)
            if self.solutions[media][target] is None:
                logger.warning("no solution could be found that satisfied all specified test conditions in specified iterations!")
                return None
        if binary_check:
            return pkgmgr.getpkg("GapfillingPkg").binary_check_gapfilling_solution()
        return self.solutions[media][target] 
    
    def integrate_gapfill_solution(self, solution):
        for rxn_id in solution["reversed"]:
            rxn = self.model.reactions.get_by_id(rxn_id)
            if solution["reversed"][rxn_id] == ">":
                rxn.upper_bound = 100
            else:
                rxn.lower_bound = -100
        for rxn_id in solution["new"]:
            rxn = self.gfmodel.reactions.get_by_id(rxn_id)
            rxn = rxn.copy()
            self.model.add_reactions([rxn])
            if solution["new"][rxn_id] == ">":
                rxn.upper_bound = 100
                rxn.lower_bound = 0 
            else:
                rxn.upper_bound = 0
                rxn.lower_bound = -100
        return self.model
    
    @staticmethod
    def gapfill(model,media = None,target_reaction = "bio1",default_gapfill_templates = [],default_gapfill_models = [],test_conditions = [],reaction_scores = {},blacklist = []):
        gapfiller = MSGapfill(model,default_gapfill_templates,default_gapfill_models,test_conditions,reaction_scores,blacklist)
        gfresults = gapfiller.run_gapfilling(media,target_reaction)
        return gapfiller.integrate_gapfill_solution(gfresults)
