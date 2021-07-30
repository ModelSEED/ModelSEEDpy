import logging
import itertools
import cobra
from modelseedpy.core.rast_client import RastClient
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msmodel import get_gpr_string, get_reaction_constraints_from_direction
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core import FBAHelper, MSGapfill
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

logger = logging.getLogger(__name__)

class MSATPCorrection:

    def __init__(self,template):
        self.max_gapfilling = None
        self.gapfilling_delta = 0
        self.default_coretemplate = template
        self.msgapfill = MSGapfill()
        self.original_bounds = {}
    
    #This function disables all noncore reactions in the model
    def disable_noncore_reactions(self,model,coretemplate = None):
        if coretemplate == None:
            coretemplate = self.default_coretemplate
        original_bounds = {}
        noncore_reactions = []
        for reaction in model.reactions:
            if reaction.id not in coretemplate.reactions and not FBAHelper.is_ex(reaction):
                original_bounds[reaction.id] = [reaction.lower_bound,reaction.upper_bound]
                if reaction.lower_bound < 0:
                    noncore_reactions.append([reaction,"<"])
                if reaction.upper_bound > 0:
                    noncore_reactions.append([reaction,">"])
                reaction.lower_bound = 0
                reaction.upper_bound = 0
                reaction.update_variable_bounds()
        return (noncore_reactions,original_bounds)
    
    #This function determines how much gapfilling each input test media requires to make ATP
    def evaluate_growth_media(self,model,test_medias = None,atp_objective = "bio2",genome_class = None):
        if test_medias == None:
            test_medias = self.compute_default_medias(genome_class)
        media_gapfill_stats = {}
        with model:
            FBAHelper.set_objective_from_target_reaction(model,atp_objective)
            pkgmgr = MSPackageManager.get_pkg_mgr(model)
            for media in test_medias:
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                solution = model.optimize()
                media_gapfill_stats[media] = None
                if solution.objective_value == 0:
                    media_gapfill_stats[media] = self.msgapfill.run_gapfilling(model,media,atp_objective)
        return media_gapfill_stats
    
    #This function decides which of the test media to use as growth conditions for this model
    def determine_growth_media(self,model,media_gapfill_stats,max_gapfilling = None,gapfilling_delta = 0):
        media_list = []
        best_score = None
        for media in media_gapfill_stats:
            gfscore = 0
            if media_gapfill_stats[media] != None:
                gfscore = len(media_gapfill_stats[media]["new"].keys()) + 0.5*len(media_gapfill_stats[media]["reversed"].keys())
            if best_score == None or gfscore < best_score:
                best_score = gfscore
        if max_gapfilling == None:
            max_gapfilling = best_score
        for media in media_gapfill_stats:
            gfscore = 0
            if media_gapfill_stats[media] != None:
                gfscore = len(media_gapfill_stats[media]["new"].keys()) + 0.5*len(media_gapfill_stats[media]["reversed"].keys())
            if gfscore <= max_gapfilling and gfscore <= (best_score+gapfilling_delta):
                media_list.append(media)
        return media_list
        
    #This function applies the gapfilling to all selected growth media
    def apply_growth_media_gapfilling(self,model,media_list,media_gapfill_stats):
        for media in media_list:
            if media in media_gapfill_stats and media_gapfill_stats[media] != None:
                model = self.msgapfill.integrate_gapfill_solution(model,media_gapfill_stats[media])
        return model
    
    #This function expands the model to genome-scale while preventing ATP overproduction
    def expand_model_to_genome_scale(self,model,media_list,noncore_reactions,atp_objective):
        gapfilling_tests = []
        pkgmgr = MSPackageManager.get_pkg_mgr(model)
        for media in media_list:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(media)    
            solution = model.optimize()
            gapfilling_tests.append({"media":media,"is_max_threshold": True,"threshold":1.2*solution.objective_value,"objective":atp_objective})
        #Extending model with noncore reactions while retaining ATP accuracy
        filtered = FBAHelper.reaction_expansion_test(model,noncore_reactions,gapfilling_tests)
        #Removing filtered reactions
        for item in filtered:
            if item[1] == ">":
                item[0].upper_bound = 0
            else:
                item[0].lower_bound = 0
            reaction.update_variable_bounds()
            if item[0].lower_bound == 0 and item[0].upper_bound == 0:
                model.remove_reactions(item[0])
        return model
    
    #This function restores the bounds on all noncore reactions
    def restore_noncore_reactions(self,model,noncore_reactions,original_bounds):
        #Restoring original reaction bounds
        for reaction in noncore_reactions:
            if reaction.id in original_bounds and reaction in model.reactions:
                reaction.lower_bound = original_bounds[reaction.id][0]
                reaction.upper_bound = original_bounds[reaction.id][1]
                reaction.update_variable_bounds()
        return model
    
    #This function runs the entire ATP method    
    def run_atp_correction(self,model,atp_medias = None,atp_objective = "bio2",genome_class = None):
        if atp_medias == None:
            pass
        coretemplate = self.default_coretemplate
        (noncore_reactions,original_bounds) = self.disable_noncore_reactions(model,coretemplate)
        media_gapfill_stats = self.evaluate_growth_media(model,atp_medias,atp_objective,genome_class)
        media_list = self.determine_growth_media(model,media_gapfill_stats,self.max_gapfilling,self.gapfilling_delta)
        model = self.apply_growth_media_gapfilling(model,media_list,media_gapfill_stats)
        model = self.expand_model_to_genome_scale(model,media_list,noncore_reactions,atp_objective)
        model = self.restore_noncore_reactions(model,noncore_reactions,original_bounds)
        
    @staticmethod
    def atp_correction(model,coretemplate,atp_medias = None,atp_objective = "bio2",genome_class = None):
        return MSATPCorrection(coretemplate).run_atp_correction(model,atp_medias,atp_objective,genome_class)
            
        