import logging
logger = logging.getLogger(__name__)
import itertools  # !!! import never used
import cobra  # !!! import never used
from optlang.symbolics import Zero, add  # !!! neither Zero nor Add are ever used
from modelseedpy.core.rast_client import RastClient  # !!! import never used
from modelseedpy.core.msgenome import normalize_role  # !!! import never used
from modelseedpy.core.msmodel import get_gpr_string, get_reaction_constraints_from_direction  # !!! neither import is ever used
from cobra.core import Gene, Metabolite, Model, Reaction  # !!! none of these imports used
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper

class MSATPCorrection:
    def __init__(self, model, core_template, atp_medias, compartment="c0",
                 max_gapfilling=None, gapfilling_delta=0, atp_hydrolysis_id=None):
        """
        :param model:
        :param core_template:
        :param atp_medias:
        :param atp_objective:
        :param max_gapfilling:
        :param gapfilling_delta:
        :param atp_hydrolysis_id: ATP Hydrolysis reaction ID, if None it will perform a SEED reaction search
        """
        self.model = model
        self.compartment = compartment
        if atp_hydrolysis_id and atp_hydrolysis_id in self.model.reactions:
            self.atp_hydrolysis = self.model.reactions.get_by_id(atp_hydrolysis_id)
        else:
            output = FBAHelper.add_atp_hydrolysis(self.model, compartment)
            self.atp_hydrolysis = output["reaction"]
        self.atp_medias = atp_medias
        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta
        self.coretemplate = core_template
        self.msgapfill = MSGapfill(model, default_gapfill_templates=core_template)
        
        self.original_bounds, self.media_gapfill_stats = {}, {}
        self.noncore_reactions, self.other_compartments, self.gapfilling_tests = [], [], []
        self.selected_media, self.filtered_noncore = [], []
        self.lp_filename = None

    def disable_noncore_reactions(self):
        """
        Disables all noncore reactions in the model
        :return:
        """
        #Must restore reactions before disabling to ensure bounds are not overwritten
        if len(self.noncore_reactions) > 0:
            self.restor_noncore_reactions(noncore = True,othercompartment = True)
        #Now clearing the existing noncore datastructures
        self.original_bounds = {}
        self.noncore_reactions, self.other_compartments = [], []
        #Iterating through reactions and disabling
        for reaction in self.model.reactions:
            if not any([reaction.id == self.atp_hydrolysis.id, FBAHelper.is_ex(reaction), FBAHelper.is_biomass(reaction)]):
                msid = FBAHelper.modelseed_id_from_cobra_reaction(reaction)
                if msid:
                    msid += "_"+self.compartment[0:1]
                if msid in self.coretemplate.reactions:
                    self.original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
                    logger.debug(reaction.id+" core")
                    if reaction.lower_bound < 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).lower_bound >= 0:
                        logger.debug(reaction.id+" core but reversible")
                        self.noncore_reactions.append([reaction, "<"])
                        reaction.lower_bound = 0
                    if reaction.upper_bound > 0 and self.coretemplate.reactions.get_by_id(reaction.id[0:-1]).upper_bound <= 0:
                        logger.debug(reaction.id+" core but reversible")
                        self.noncore_reactions.append([reaction, ">"])
                        reaction.upper_bound = 0
                else:
                    logger.debug(reaction.id+" noncore")
                    self.original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)
                    if reaction.lower_bound < 0:
                        if FBAHelper.rxn_compartment(reaction) != self.compartment:
                            self.other_compartments.append([reaction, "<"])
                        else:
                            self.noncore_reactions.append([reaction, "<"])
                    if reaction.upper_bound > 0:
                        if FBAHelper.rxn_compartment(reaction) != self.compartment:
                            self.other_compartments.append([reaction, ">"])
                        else:
                            self.noncore_reactions.append([reaction, ">"])
                    reaction.lower_bound = reaction.upper_bound = 0

    def evaluate_growth_media(self):
        """
        Determines how much gap filling each input test media requires to make ATP
        
        :return:
        """
        self.disable_noncore_reactions()
        self.media_gapfill_stats = {}
        self.msgapfill.default_gapfill_templates = [self.coretemplate]
        if self.lp_filename:
            self.msgapfill.lp_filename = self.lp_filename
        output = {}
        with self.model:
            self.model.objective = self.atp_hydrolysis.id
            #self.model.objective = self.model.problem.Objective(Zero,direction="max")
            #self.atp_hydrolysis.update_variable_bounds()
            logger.debug(f"ATP bounds: {self.atp_hydrolysis.lower_bound} : {self.atp_hydrolysis.upper_bound}")
            #self.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            for media in self.atp_medias:
                logger.debug('evaluate media %s', media)
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                solution = self.model.optimize()
                logger.debug('evaluate media %s - %f (%s)', media.id, solution.objective_value, solution.status)
                self.media_gapfill_stats[media] = None
                output[media.id] = solution.objective_value
                if solution.objective_value == 0 or solution.status != 'optimal':
                    self.media_gapfill_stats[media] = self.msgapfill.run_gapfilling(media, self.atp_hydrolysis.id)
                    #IF gapfilling fails - need to activate and penalize the noncore and try again
                elif solution.objective_value > 0 or solution.status == 'optimal':
                    self.media_gapfill_stats[media] = {'reversed': {}, 'new': {}}
        return output

    def determine_growth_media(self):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        from math import inf 
        
        self.selected_media = []
        best_score = inf
        for media in self.media_gapfill_stats:
            if self.media_gapfill_stats[media]:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            best_score = min(best_score, gfscore)
        if self.max_gapfilling is None:
            self.max_gapfilling = best_score
        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media]:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())
            if gfscore <= self.max_gapfilling and gfscore <= (best_score+self.gapfilling_delta):
                self.selected_media.append(media)

    def determine_growth_media2(self, max_gapfilling=None):  #!!! unused function
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        def scoring_function(media):
            return len(self.media_gapfill_stats[media]["new"].keys()) + 0.5 * \
                    len(self.media_gapfill_stats[media]["reversed"].keys())

        self.selected_media = []
        media_scores = {media:scoring_function(media) for media in self.media_gapfill_stats if self.media_gapfill_stats[media]}
        best_score = min(media_scores.values())
        if max_gapfilling is None:
            max_gapfilling = best_score
        for media, score in media_scores.items():
            logger.debug(score, best_score, max_gapfilling)
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
        self.gapfilling_tests, self.filtered_noncore = [], []
        self.model.objective = self.atp_hydrolysis.id
        pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
        for media in self.selected_media:
            pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            solution = self.model.optimize()
            logger.debug(media.name+" = "+str(solution.objective_value))
            self.gapfilling_tests.append({"media":media,"is_max_threshold": True,"threshold":1.2*solution.objective_value,"objective":self.atp_hydrolysis.id})
        #Must restore noncore reactions and NOT other compartment reactions before running this function - it is not detrimental to run this twice
        self.restore_noncore_reactions(noncore = True,othercompartment = False)
        # Extending model with noncore reactions while retaining ATP accuracy
        self.filtered_noncore = FBAHelper.reaction_expansion_test(self.model,self.noncore_reactions,self.gapfilling_tests,pkgmgr)
        # Removing filtered reactions
        for item in self.filtered_noncore:
            print("Removing "+item[0].id+" "+item[1])
            if item[1] == ">":
                item[0].upper_bound = 0
            else:
                item[0].lower_bound = 0
            # reaction.update_variable_bounds()
            if item[0].lower_bound == 0 and item[0].upper_bound == 0:
                self.model.remove_reactions([item[0]])
        #Restoring other compartment reactions but not the core because this would undo reaction filtering
        self.restore_noncore_reactions(noncore = False,othercompartment = True)

    def restore_noncore_reactions(self,noncore = True,othercompartment = True):
        """
        Restores the bounds on all noncore reactions
        :return:
        """
        # Restoring original reaction bounds
        if noncore:
            for item in self.noncore_reactions:
                reaction = item[0]
                if reaction.id in self.original_bounds and reaction in self.model.reactions:
                    reaction.lower_bound = self.original_bounds[reaction.id][0]
                    reaction.upper_bound = self.original_bounds[reaction.id][1]
        if othercompartment:
            for item in self.other_compartments:
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
        self.evaluate_growth_media()
        self.determine_growth_media()
        self.apply_growth_media_gapfilling()
        self.expand_model_to_genome_scale()

    @staticmethod
    def atp_correction(model, coretemplate, atp_medias=None, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
        msatpobj = MSATPCorrection(model, coretemplate, atp_medias, atp_objective, max_gapfilling, gapfilling_delta)
        msatpobj.run_atp_correction()
        return msatpobj

