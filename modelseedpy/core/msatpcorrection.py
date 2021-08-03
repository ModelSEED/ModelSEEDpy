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

    def __init__(self,model,coretemplate,atp_medias,atp_objective = "bio2",max_gapfilling = None,gapfilling_delta = 0):
        self.atp_medias = atp_medias
        self.atp_objective = atp_objective
        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta
        self.coretemplate = coretemplate
        self.msgapfill = MSGapfill(model)
        self.original_bounds = {}
        self.noncore_reactions = []
        self.model = model
        self.media_gapfill_stats = {}
        self.gapfilling_tests = []
        self.selected_media = []
        self.filtered_noncore = []
        self.lp_filename = None
    
    #This function disables all noncore reactions in the model
    def disable_noncore_reactions(self):
        self.original_bounds = {}
        self.noncore_reactions = []
        for reaction in self.model.reactions:
            if reaction.id[0:-1] not in self.coretemplate.reactions and not FBAHelper.is_ex(reaction) and reaction.id[0:3] != "bio":
                self.original_bounds[reaction.id] = [reaction.lower_bound,reaction.upper_bound]
                if reaction.lower_bound < 0:
                    self.noncore_reactions.append([reaction,"<"])
                if reaction.upper_bound > 0:
                    self.noncore_reactions.append([reaction,">"])
                reaction.lower_bound = 0
                reaction.upper_bound = 0
                reaction.update_variable_bounds()
            elif reaction.id[0:-1] in self.coretemplate.reactions:
                self.original_bounds[reaction.id] = [reaction.lower_bound,reaction.upper_bound]
                if reaction.lower_bound < 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).lower_bound >= 0:
                    self.noncore_reactions.append([reaction,"<"])
                    reaction.lower_bound = 0
                    reaction.update_variable_bounds()
                if reaction.upper_bound > 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).upper_bound <= 0:
                    self.noncore_reactions.append([reaction,">"])
                    reaction.upper_bound = 0
                    reaction.update_variable_bounds()
    
    #This function determines how much gapfilling each input test media requires to make ATP
    def evaluate_growth_media(self):
        self.media_gapfill_stats = {}
        self.msgapfill.default_gapfill_templates = [self.coretemplate]
        if self.lp_filename != None:
            self.msgapfill.lp_filename = self.lp_filename
        with self.model:
            self.model.objective = self.atp_objective
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            for media in self.atp_medias:
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                solution = self.model.optimize()
                self.media_gapfill_stats[media] = None
                if solution.objective_value == 0:
                    self.media_gapfill_stats[media] = self.msgapfill.run_gapfilling(media,self.atp_objective)
                    print(media.id+":"+str(len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())))
    
    #This function decides which of the test media to use as growth conditions for this model
    def determine_growth_media(self):
        self.selected_media = []
        best_score = None
        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media] != None:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            if best_score == None or gfscore < best_score:
                best_score = gfscore
        if self.max_gapfilling == None:
            self.max_gapfilling = best_score
        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media] != None:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            if gfscore <= self.max_gapfilling and gfscore <= (best_score+self.gapfilling_delta):
                print("Keeping:"+media.id)
                self.selected_media.append(media)
        
    #This function applies the gapfilling to all selected growth media
    def apply_growth_media_gapfilling(self):
        for media in self.selected_media:
            if media in self.media_gapfill_stats and self.media_gapfill_stats[media] != None:
                print("Merging gapfill for media "+media.id)
                self.model = self.msgapfill.integrate_gapfill_solution(self.model,self.media_gapfill_stats[media])
    
    #This function expands the model to genome-scale while preventing ATP overproduction
    def expand_model_to_genome_scale(self):
        self.gapfilling_tests = []
        self.filtered_noncore = []
        self.model.objective = self.atp_objective
        pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
        for media in self.selected_media:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(media)    
            solution = self.model.optimize()
            print("Objective in "+media.id+":"+str(solution.objective_value))
            self.gapfilling_tests.append({"media":media,"is_max_threshold": True,"threshold":1.2*solution.objective_value,"objective":self.atp_objective})
        #Extending model with noncore reactions while retaining ATP accuracy
        self.filtered_noncore = FBAHelper.reaction_expansion_test(self.model,self.noncore_reactions,self.gapfilling_tests)
        #Removing filtered reactions
        for item in self.filtered_noncore:
            if item[1] == ">":
                item[0].upper_bound = 0
            else:
                item[0].lower_bound = 0
            reaction.update_variable_bounds()
            if item[0].lower_bound == 0 and item[0].upper_bound == 0:
                self.model.remove_reactions(item[0])
    
    #This function restores the bounds on all noncore reactions
    def restore_noncore_reactions(self):
        #Restoring original reaction bounds
        for item in self.noncore_reactions:
            reaction = item[0]
            if reaction.id in self.original_bounds and reaction in self.model.reactions:
                reaction.lower_bound = self.original_bounds[reaction.id][0]
                reaction.upper_bound = self.original_bounds[reaction.id][1]
                reaction.update_variable_bounds()
    
    #This function runs the entire ATP method    
    def run_atp_correction(self):
        self.disable_noncore_reactions()
        self.evaluate_growth_media()
        self.determine_growth_media()
        self.apply_growth_media_gapfilling()
        self.expand_model_to_genome_scale()
        self.restore_noncore_reactions()
        
    @staticmethod
    def atp_correction(model,coretemplate,atp_medias = None,atp_objective = "bio2",max_gapfilling = None,gapfilling_delta = 0):
        msatpobj = MSATPCorrection(model,coretemplate,atp_medias,atp_objective,max_gapfilling,gapfilling_delta)
        msatpobj.run_atp_correction()
        return msatpobj
    
            
        