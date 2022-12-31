# -*- coding: utf-8 -*-
import logging
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.fbahelper import FBAHelper
import itertools
import cobra
import json
import time
import pandas as pd
from os.path import abspath as _abspath
from os.path import dirname as _dirname
from optlang.symbolics import Zero, add
from modelseedpy.core.rast_client import RastClient
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msmodel import (
    get_gpr_string,
    get_reaction_constraints_from_direction,
)
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.msmedia import MSMedia
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.helpers import get_template

logger = logging.getLogger(__name__)

_path = _dirname(_abspath(__file__))

min_gap = {
    "Glc/O2": 5,
    "Etho/O2": 0.01,
    "Ac/O2": 1,
    "Pyr/O2": 3,
    "Glyc/O2": 2,
    "Fum/O2": 3,
    "Succ/O2": 2,
    "Akg/O2": 2,
    "LLac/O2": 2,
    "Dlac/O2": 2,
    "For/O2": 2,
    "For/NO3": 1.5,
    "Pyr/NO": 2.5,
    "Pyr/NO2": 2.5,
    "Pyr/NO3": 2.5,
    "Pyr/SO4": 2.5,
}


class MSATPCorrection:

    DEBUG = False

    def __init__(
        self,
        model_or_mdlutl,
        core_template=None,
        atp_medias=[],
        compartment="c0",
        max_gapfilling=10,
        gapfilling_delta=0,
        atp_hydrolysis_id=None,
        load_default_medias=True,
        forced_media=[],
        default_media_path=None,
    ):
        """
        :param model:
        :param core_template:
        :param atp_medias: list<MSMedia> : list of additional medias to test
        :param load_default_medias: Bool : load default media set
        :param forced_media: list<string> : name of medias in which ATP production should be forced
        :param compartment: string : ID of compartment to test ATP in
        :param max_gapfilling: string : maximum gapfilling allowed in accepted media
        :param gapfilling_delta: string : difference between lowest gapfilling and current gapfilling where media will be accepted
        :param atp_hydrolysis_id: string : ATP Hydrolysis reaction ID, if None it will perform a SEED reaction search
        """
        # Discerning input is model or mdlutl and setting internal links
        if isinstance(model_or_mdlutl, MSModelUtil):
            self.model = model_or_mdlutl.model
            self.modelutl = model_or_mdlutl
        else:
            self.model = model_or_mdlutl
            self.modelutl = MSModelUtil.get(model_or_mdlutl)
        # Setting atpcorrection attribute in model utl so link is bidirectional
        self.modelutl.atputl = self

        if default_media_path:
            self.default_media_path = default_media_path
        else:
            self.default_media_path = _path + "/../data/atp_medias.tsv"

        self.compartment = compartment

        if atp_hydrolysis_id and atp_hydrolysis_id in self.model.reactions:
            self.atp_hydrolysis = self.model.reactions.get_by_id(atp_hydrolysis_id)
        else:
            output = self.modelutl.add_atp_hydrolysis(compartment)
            self.atp_hydrolysis = output["reaction"]

        self.atp_medias = []
        if load_default_medias:
            self.load_default_medias()
        for media in atp_medias:
            if isinstance(media, list):
                self.atp_medias.append(media)
            else:
                self.atp_medias.append([media, 0.01])

        self.forced_media = []
        for media_id in forced_media:
            for media in self.atp_medias:
                if media.id == media_id:
                    self.forced_media.append(media)
                    break

        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta

        if not core_template:
            self.load_default_template()
        else:
            self.coretemplate = core_template

        self.msgapfill = MSGapfill(
            self.modelutl, default_gapfill_templates=core_template
        )
        # These should stay as None until atp correction is actually run
        self.cumulative_core_gapfilling = None
        self.selected_media = None
        self.original_bounds = {}
        self.noncore_reactions = []
        self.other_compartments = []
        self.media_gapfill_stats = {}
        self.filtered_noncore = []
        self.lp_filename = None
        self.multiplier = 1.2

    def load_default_template(self):
        self.coretemplate = MSTemplateBuilder.from_dict(
            get_template("template_core"), None
        ).build()

    def load_default_medias(self):
        filename = self.default_media_path
        medias = pd.read_csv(filename, sep="\t", index_col=0).to_dict()
        for media_id in medias:
            media_d = {}
            for exchange, v in medias[media_id].items():
                if v > 0:
                    k = exchange.split("_")[1]
                    media_d[k] = v
            media_d["cpd00001"] = 1000
            media_d["cpd00067"] = 1000
            media = MSMedia.from_dict(media_d)
            media.id = media_id
            media.name = media_id
            min_obj = 0.01
            if media_id in min_gap:
                min_obj = min_gap[media_id]
            self.atp_medias.append([media, min_obj])

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
                template_reaction = template.reactions.get_by_id(
                    model_reaction.id[0:-1]
                )
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

                    return s[: str_pos + 1], s[str_pos + 1 :]

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
        # Iterating through reactions and disabling
        for reaction in self.model.reactions:
            if reaction.id == self.atp_hydrolysis.id:
                continue
            if FBAHelper.is_ex(reaction):
                continue
            if FBAHelper.is_biomass(reaction):
                continue

            self.original_bounds[reaction.id] = (
                reaction.lower_bound,
                reaction.upper_bound,
            )

            # check if reaction is in core template
            template_reaction = self.find_reaction_in_template(
                reaction, self.coretemplate, self.compartment[0:1]
            )

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
                reaction.lower_bound = 0
                reaction.upper_bound = 0

    def evaluate_growth_media(self):
        """Determines how much gap filling each input test media requires to make ATP"""
        self.disable_noncore_reactions()
        self.media_gapfill_stats = {}
        self.msgapfill.default_gapfill_templates = [self.coretemplate]
        if self.lp_filename:
            self.msgapfill.lp_filename = self.lp_filename
        results = {}
        with self.model:
            self.model.objective = self.atp_hydrolysis.id
            # self.model.objective = self.model.problem.Objective(Zero,direction="max")

            logger.debug(
                f"ATP bounds: ({self.atp_hydrolysis.lower_bound}, {self.atp_hydrolysis.upper_bound})"
            )
            # self.model.objective.set_linear_coefficients({self.atp_hydrolysis.forward_variable:1})
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            for media, minimum_obj in self.atp_medias:
                logger.debug("evaluate media %s", media)
                pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                logger.debug("model.medium %s", self.model.medium)
                solution = self.model.optimize()
                logger.debug(
                    "evaluate media %s - %f (%s)",
                    media.id,
                    solution.objective_value,
                    solution.status,
                )
                self.media_gapfill_stats[media] = None
                output[media.id] = solution.objective_value

                if (
                    solution.objective_value < minimum_obj
                    or solution.status != "optimal"
                ):
                    self.media_gapfill_stats[media] = self.msgapfill.run_gapfilling(
                        media, self.atp_hydrolysis.id, minimum_obj
                    )
                    # IF gapfilling fails - need to activate and penalize the noncore and try again
                elif solution.objective_value >= minimum_obj:
                    self.media_gapfill_stats[media] = {"reversed": {}, "new": {}}
                logger.debug(
                    "gapfilling stats: %s",
                    json.dumps(self.media_gapfill_stats[media], indent=2, default=vars),
                )

        if MSATPCorrection.DEBUG:
            with open("debug.json", "w") as outfile:
                json.dump(self.media_gapfill_stats[media], outfile)

        return output

    def determine_growth_media(self, max_gapfilling=None):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        self.selected_media = []
        best_score = inf
        for media in self.media_gapfill_stats:
            if self.media_gapfill_stats[media]:
                gfscore = len(
                    self.media_gapfill_stats[media]["new"].keys()
                ) + 0.5 * len(self.media_gapfill_stats[media]["reversed"].keys())
            best_score = min(best_score, gfscore)
        if self.max_gapfilling is None:
            self.max_gapfilling = best_score

        logger.debug(f"max_gapfilling: {self.max_gapfilling}, best_score: {best_score}")

        for media in self.media_gapfill_stats:
            gfscore = 0
            if self.media_gapfill_stats[media]:
                gfscore = len(
                    self.media_gapfill_stats[media]["new"].keys()
                ) + 0.5 * len(self.media_gapfill_stats[media]["reversed"].keys())

            logger.debug(f"media gapfilling score: {media.id}: {gfscore}")
            if gfscore <= self.max_gapfilling and gfscore <= (
                best_score + self.gapfilling_delta
            ):
                self.selected_media.append(media)

    def determine_growth_media2(self, max_gapfilling=None):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """

        def scoring_function(media):
            return len(self.media_gapfill_stats[media]["new"].keys()) + 0.5 * len(
                self.media_gapfill_stats[media]["reversed"].keys()
            )

        if not max_gapfilling:
            max_gapfilling = self.max_gapfilling
        self.selected_media = []
        media_scores = dict(
            (media, scoring_function(media))
            for media in self.media_gapfill_stats
            if self.media_gapfill_stats[media]
        )
        best_score = min(media_scores.values())
        if max_gapfilling is None or max_gapfilling > (
            best_score + self.gapfilling_delta
        ):
            max_gapfilling = best_score + self.gapfilling_delta
        for media in media_scores:
            score = media_scores[media]
            logger.debug(score, best_score, max_gapfilling)
            if score <= max_gapfilling:
                self.selected_media.append(media)

    def apply_growth_media_gapfilling(self):
        """
        Applies the gapfilling to all selected growth media
        :return:
        """
        self.cumulative_core_gapfilling = (
            []
        )  # TODO: In case someone runs ATP correction twice with different parameters, before resetting this, maybe check if any of these reactions are already in the model and remove them so we're starting fresh???
        for media in self.selected_media:
            if (
                media in self.media_gapfill_stats
                and self.media_gapfill_stats[media]
                and MSGapfill.gapfill_count(self.media_gapfill_stats[media]) > 0
            ):
                self.msgapfill.integrate_gapfill_solution(
                    self.media_gapfill_stats[media], self.cumulative_core_gapfilling
                )

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
        self.filtered_noncore = self.modelutl.reaction_expansion_test(
            self.noncore_reactions, tests
        )
        # Removing filtered reactions
        for item in self.filtered_noncore:
            print("Removing " + item[0].id + " " + item[1])
            if item[1] == ">":
                item[0].upper_bound = 0
            else:
                item[0].lower_bound = 0
            # reaction.update_variable_bounds()
            if item[0].lower_bound == 0 and item[0].upper_bound == 0:
                self.model.remove_reactions([item[0]])
        # Restoring other compartment reactions but not the core because this would undo reaction filtering
        self.restore_noncore_reactions(noncore=False, othercompartment=True)

    def restore_noncore_reactions(self, noncore=True, othercompartment=True):
        """
        Restores the bounds on all noncore reactions
        :return:
        """
        # Restoring original reaction bounds
        if noncore:
            for item in self.noncore_reactions:
                reaction = item[0]
                if (
                    reaction.id in self.original_bounds
                    and reaction in self.model.reactions
                ):
                    reaction.lower_bound = self.original_bounds[reaction.id][0]
                    reaction.upper_bound = self.original_bounds[reaction.id][1]
        if othercompartment:
            for item in self.other_compartments:
                reaction = item[0]
                if (
                    reaction.id in self.original_bounds
                    and reaction in self.model.reactions
                ):
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
            logger.debug(f"{media.name} = {obj_value}")
            tests.append(
                {
                    "media": media,
                    "is_max_threshold": True,
                    "threshold": multiplier * obj_value,
                    "objective": self.atp_hydrolysis.id,
                }
            )
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
    def atp_correction(
        model,
        core_template,
        atp_medias=None,
        atp_objective="bio2",
        max_gapfilling=None,
        gapfilling_delta=0,
    ):
        atp_correction = MSATPCorrection(
            model,
            core_template,
            atp_medias,
            atp_hydrolysis_id=atp_objective,
            max_gapfilling=max_gapfilling,
            gapfilling_delta=gapfilling_delta,
        )
        return atp_correction.run_atp_correction()
