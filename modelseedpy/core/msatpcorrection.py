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

    def __init__(self, model, core_template, atp_medias,compartment="c0", max_gapfilling=None, gapfilling_delta=0):
        """

        :param model:
        :param core_template:
        :param atp_medias:
        :param atp_objective:
        :param max_gapfilling:
        :param gapfilling_delta:
        """
        self.model = model
        self.compartment = compartment
        output = FBAHelper.add_atp_hydrolysis(self.model,compartment)
        self.atp_hydrolysis = output["reaction"]
        self.atp_medias = atp_medias
        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta
        self.coretemplate = core_template
        self.msgapfill = MSGapfill(model,default_gapfill_templates=core_template)
        self.original_bounds = {}
        self.noncore_reactions = []
        
        self.media_gapfill_stats = {}
        self.gapfilling_tests = []
        self.selected_media = []
        self.filtered_noncore = []
        self.lp_filename = None

    def disable_noncore_reactions(self):
        """
        Disables all noncore reactions in the model
        :return:
        """
        self.original_bounds = {}
        self.noncore_reactions = []
        for reaction in self.model.reactions:
            if FBAHelper.is_ex(reaction):
                continue
            if FBAHelper.is_biomass(reaction):
                continue
            msid = FBAHelper.modelseed_id_from_cobra_reaction(reaction)
            msid += "_"+self.compartment[0:1]
            if msid in self.coretemplate.reactions and FBAHelper.rxn_compartment(reaction) == self.compartment:
                self.original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
                if reaction.lower_bound < 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).lower_bound >= 0:
                    self.noncore_reactions.append([reaction, "<"])
                    reaction.lower_bound = 0
                if reaction.upper_bound > 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).upper_bound <= 0:
                    self.noncore_reactions.append([reaction, ">"])
                    reaction.upper_bound = 0
            else:
                self.original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
                if reaction.lower_bound < 0:
                    self.noncore_reactions.append([reaction, "<"])
                if reaction.upper_bound > 0:
                    self.noncore_reactions.append([reaction, ">"])
                reaction.lower_bound = 0
                reaction.upper_bound = 0

    def evaluate_growth_media(self):
        """
        Determines how much gap filling each input test media requires to make ATP

        :return:
        """
        self.media_gapfill_stats = {}
        self.msgapfill.default_gapfill_templates = [self.coretemplate]
        if self.lp_filename:
            self.msgapfill.lp_filename = self.lp_filename
        with self.model:
            self.model.objective = self.model.problem.Objective(Zero,direction="max")
            self.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            for media in self.atp_medias:
                logger.debug('evaluate media %s', media)
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                solution = self.model.optimize()
                logger.debug('evaluate media %s - %f (%s)', media, solution.objective_value, solution.status)
                self.media_gapfill_stats[media] = None
                if solution.objective_value == 0 or solution.status != 'optimal':
                    self.media_gapfill_stats[media] = self.msgapfill.run_gapfilling(media, self.atp_objective)
                    #IF gapfilling fails - need to activate and penalize the noncore and try again

    def determine_growth_media(self):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        self.selected_media = []
        best_score = None
        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media]:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            if best_score is None or gfscore < best_score:
                best_score = gfscore
        if self.max_gapfilling is None:
            self.max_gapfilling = best_score
        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media]:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            if gfscore <= self.max_gapfilling and gfscore <= (best_score+self.gapfilling_delta):
                self.selected_media.append(media)

    def determine_growth_media2(self, max_gapfilling=None):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        def scoring_function(media):
            return len(self.media_gapfill_stats[media]["new"].keys()) + 0.5 * \
                   len(self.media_gapfill_stats[media]["reversed"].keys())

        self.selected_media = []
        media_scores = dict(
            (media, scoring_function(media)) for media in self.media_gapfill_stats if self.media_gapfill_stats[media])
        best_score = min(media_scores.values())
        if max_gapfilling is None:
            max_gapfilling = best_score
        for media in media_scores:
            score = media_scores[media]
            print(score, best_score, max_gapfilling)
            if score <= max_gapfilling and score <= (best_score + self.gapfilling_delta):
                self.selected_media.append(media)

    def apply_growth_media_gapfilling(self):
        """
        Applies the gapfilling to all selected growth media
        :return:
        """
        for media in self.selected_media:
            if media in self.media_gapfill_stats and self.media_gapfill_stats[media]:
                self.model = self.msgapfill.integrate_gapfill_solution(self.model, self.media_gapfill_stats[media])

    def expand_model_to_genome_scale(self):
        """
        Expands the model to genome-scale while preventing ATP overproduction
        :return:
        """
        self.gapfilling_tests = []
        self.filtered_noncore = []
        self.model.objective = self.atp_objective
        pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
        for media in self.selected_media:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            solution = self.model.optimize()
            self.gapfilling_tests.append({"media":media,"is_max_threshold": True,"threshold":1.2*solution.objective_value,"objective":self.atp_objective})
        # Extending model with noncore reactions while retaining ATP accuracy
        self.filtered_noncore = FBAHelper.reaction_expansion_test(self.model,self.noncore_reactions,self.gapfilling_tests)
        # Removing filtered reactions
        for item in self.filtered_noncore:
            if item[1] == ">":
                item[0].upper_bound = 0
            else:
                item[0].lower_bound = 0
            # reaction.update_variable_bounds()
            if item[0].lower_bound == 0 and item[0].upper_bound == 0:
                self.model.remove_reactions(item[0])

    def restore_noncore_reactions(self):
        """
        Restores the bounds on all noncore reactions
        :return:
        """
        # Restoring original reaction bounds
        for item in self.noncore_reactions:
            reaction = item[0]
            if reaction.id in self.original_bounds and reaction in self.model.reactions:
                reaction.lower_bound = self.original_bounds[reaction.id][0]
                reaction.upper_bound = self.original_bounds[reaction.id][1]

    def run_atp_correction(self):
        """
        Runs the entire ATP method
        :return:
        """
        #Ensure all specified media work
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
