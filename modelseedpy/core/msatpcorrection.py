import logging
import itertools  # !!! import never used
import cobra  # !!! import never used
import json
import time
from optlang.symbolics import Zero, add  # !!! neither Zero nor Add are ever used
from modelseedpy.core.rast_client import RastClient  # !!! import never used
from modelseedpy.core.msgenome import normalize_role  # !!! import never used
from modelseedpy.core.msmodel import get_gpr_string, get_reaction_constraints_from_direction  # !!! neither import is ever used
from cobra.core import Gene, Metabolite, Model, Reaction  # !!! none of these imports used
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core import FBAHelper, MSGapfill, MSMedia
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

logger = logging.getLogger(__name__)

class MSATPCorrection:

    DEBUG = False

    def __init__(self, model, core_template, atp_medias: list, compartment="c0",
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
        if isinstance(model, MSModelUtil):
            self.model = model.model
            self.modelutl = model
        else:
            self.model = model
            self.modelutl = MSModelUtil(model)
        self.compartment = compartment
        if atp_hydrolysis_id and atp_hydrolysis_id in self.model.reactions:
            self.atp_hydrolysis = self.model.reactions.get_by_id(atp_hydrolysis_id)
        else:
            output = self.modelutl.add_atp_hydrolysis(compartment)
            self.atp_hydrolysis = output["reaction"]
        self.atp_medias = []
        for media in atp_medias:
            if isinstance(media, MSMedia):
                self.atp_medias.append([media, 0.01])
            else:
                self.atp_medias.append(media)        
        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta
        self.coretemplate = core_template
        self.msgapfill = MSGapfill(self.modelutl, default_gapfill_templates=core_template)
        self.original_bounds, self.media_gapfill_stats = {}, {}
        self.noncore_reactions, self.other_compartments = [], []
        self.selected_media, self.filtered_noncore = [], []
        self.lp_filename = None
        self.multiplier = 1.2

    @staticmethod
    def find_reaction_in_template(model_reaction, template, compartment):
        template_reaction = None  # we save lookup result here
        if model_reaction.id in template.reactions:
            template_reaction = template.reactions.get_by_id(model_reaction.id)
        else:
            msid = FBAHelper.modelseed_id_from_cobra_reaction(model_reaction)
            if msid is not None:
                msid += "_" + compartment
            if msid in template.reactions:
                template_reaction = template.reactions.get_by_id(model_reaction.id[0:-1])
            else:
                # will leave this here for now
                def split_id_from_index(s):
                    """
                    Extracts the last digits of a string example: rxn12345, returns rxn 12345

                    @param s: any string
                    @return: string split into head (remaining) + tail (digits)
                    """
                    str_pos = len(s) - 1
                    while str_pos >= 0:
                        if not s[str_pos].isdigit():
                            break
                        str_pos -= 1

                    return s[:str_pos + 1], s[str_pos + 1:]

                rxn_id, index = split_id_from_index(model_reaction.id)
                if rxn_id in template.reactions:
                    template_reaction = template.reactions.get_by_id(rxn_id)

        return template_reaction

    def disable_noncore_reactions(self):
        """
        Disables all non core reactions in the model
        :return:
        """
        # Must restore reactions before disabling to ensure bounds are not overwritten
        if len(self.noncore_reactions) > 0:
            self.restore_noncore_reactions(noncore=True, othercompartment=True)
        # Now clearing the existing noncore data structures
        self.original_bounds = {}
        self.noncore_reactions, self.other_compartments = [], []
        #Iterating through reactions and disabling
        self.noncore_reactions, self.other_compartments = [], []
        # Iterating through reactions and disabling
        for reaction in self.model.reactions:
            if any([reaction.id == self.atp_hydrolysis.id, FBAHelper.is_ex(reaction), FBAHelper.is_biomass(reaction)]):
                continue
            self.original_bounds[reaction.id] = (reaction.lower_bound, reaction.upper_bound)

            # check if reaction is in core template
            template_reaction = self.find_reaction_in_template(reaction, self.coretemplate, self.compartment[0:1])

            # update bounds to reaction
            if template_reaction is not None:
                logger.debug(f"{reaction.id} core")
                if reaction.lower_bound < 0 and template_reaction.lower_bound >= 0:
                    logger.debug(reaction.id + " core but reversible")
                    self.noncore_reactions.append([reaction, "<"])
                    reaction.lower_bound = 0
                if reaction.upper_bound > 0 and template_reaction.upper_bound <= 0:
                    logger.debug(reaction.id + " core but reversible")
                    self.noncore_reactions.append([reaction, ">"])
                    reaction.upper_bound = 0
            else:
                logger.debug(f"{reaction.id} non core")
                if FBAHelper.rxn_compartment(reaction) != self.compartment:
                    if reaction.lower_bound < 0:
                        self.other_compartments.append([reaction, "<"])
                    if reaction.upper_bound > 0:
                        self.other_compartments.append([reaction, ">"])
                else:
                    if reaction.lower_bound < 0:
                        self.noncore_reactions.append([reaction, "<"])
                    if reaction.upper_bound > 0:
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

            logger.debug(f'ATP bounds: ({self.atp_hydrolysis.lower_bound}, {self.atp_hydrolysis.upper_bound})')
            #self.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            for media, minimum_obj in self.atp_medias:
                logger.debug('evaluate media %s', media)
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                logger.debug('model.medium %s', self.model.medium)
                solution = self.model.optimize()
                logger.debug('evaluate media %s - %f (%s)', media.id, solution.objective_value, solution.status)
                self.media_gapfill_stats[media] = None
                output[media.id] = solution.objective_value
                if solution.objective_value < minimum_obj or solution.status != 'optimal':
                    self.media_gapfill_stats[media] = self.msgapfill.run_gapfilling(media,
                                                                                    self.atp_hydrolysis.id,
                                                                                    minimum_obj)
                    #IF gapfilling fails - need to activate and penalize the noncore and try again
                elif solution.objective_value >= minimum_obj:
                    self.media_gapfill_stats[media] = {'reversed': {}, 'new': {}}
                logger.debug('gapfilling stats: %s', json.dumps(self.media_gapfill_stats[media], indent=2))

        if MSATPCorrection.DEBUG:
            with open('debug.json', 'w') as outfile:
                json.dump(self.media_gapfill_stats[media], outfile)

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

        logger.debug(f'max_gapfilling: {self.max_gapfilling}, best_score: {best_score}')

        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media]:
                gfscore = len(self.media_gapfill_stats[media]["new"].keys()) + 0.5*len(self.media_gapfill_stats[media]["reversed"].keys())

            logger.debug(f'media gapfilling score: {media.id}: {gfscore}')
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
                self.model = self.msgapfill.integrate_gapfill_solution(self.media_gapfill_stats[media])

    def expand_model_to_genome_scale(self):
        """Restores noncore reactions to model while filtering out reactions that break ATP
        Parameters
        ----------
        Returns
        -------
        Raises
        ------
        """
        self.filtered_noncore = []
        tests = self.build_tests()
        # Must restore non core reactions and NOT other compartment reactions before running this function
        # it is not detrimental to run this twice
        self.restore_noncore_reactions(noncore=True, othercompartment=False)
        # Extending model with non core reactions while retaining ATP accuracy
        self.filtered_noncore = self.modelutl.reaction_expansion_test(self.noncore_reactions,tests)        
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
        # Restoring other compartment reactions but not the core because this would undo reaction filtering
        self.restore_noncore_reactions(noncore=False, othercompartment=True)

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

    def build_tests(self, multiplier=None):
        """Build tests based on ATP media evaluations
        
        Parameters
        ----------
        multiplier : float
            Override for multiplier to ATP threshold dictating when tests will pass and fail
                    
        Returns
        -------
        list<{"media":obj media,"is_max_threshold":bool,"threshold":float,"objective":string}>
            List of test specifications
            
        Raises
        ------
        """
        if multiplier is None:
            multiplier = self.multiplier
        tests = []
        self.model.objective = self.atp_hydrolysis.id
        for media in self.selected_media:
            self.modelutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            obj_value = self.model.slim_optimize()
            logger.debug(f'{media.name} = {obj_value}')
            tests.append({
                "media": media,
                "is_max_threshold": True,
                "threshold": multiplier*obj_value,
                "objective": self.atp_hydrolysis.id
            })
        return tests

    def run_atp_correction(self):
        """
        Runs the entire ATP method
        :return:
        """
        # Ensure all specified media work
        self.evaluate_growth_media()
        self.determine_growth_media()
        self.apply_growth_media_gapfilling()
        self.evaluate_growth_media()
        self.expand_model_to_genome_scale()
        return self.build_tests()
    
    @staticmethod
    def atp_correction(model, core_template, atp_medias=None, atp_objective="bio2", max_gapfilling=None, gapfilling_delta=0):
          atp_correction = MSATPCorrection(model, core_template, atp_medias, max_gapfilling=max_gapfilling, 
                                           gapfilling_delta=gapfilling_delta, atp_hydrolysis_id=atp_objective)
          return atp_correction.run_atp_correction()
