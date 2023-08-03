# -*- coding: utf-8 -*-
import logging
import cobra
import copy
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
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core import FBAHelper, MSGapfill, MSMedia
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.helpers import get_template

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
)  # When debugging - set this to INFO then change needed messages below from DEBUG to INFO

_path = _dirname(_abspath(__file__))

min_gap = {
    "Glc.O2": 5,
    "Etho.O2": 0.01,
    "Ac.O2": 1,
    "Pyr.O2": 3,
    "Glyc.O2": 2,
    "Fum.O2": 3,
    "Succ.O2": 2,
    "Akg.O2": 2,
    "LLac.O2": 2,
    "Dlac.O2": 2,
    "For.O2": 1.875,
    "For.NO3": 1.5,
    "Pyr.NO": 2.5,
    "Pyr.NO2": 2.5,
    "Pyr.NO3": 2.5,
    "Pyr.SO4": 2.5,
}

default_threshold_multipiers = {
    "Glc": 2,
    "default":1.2,
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

        self.media_hash = {}
        self.atp_medias = []
        if load_default_medias:
            self.load_default_medias()
        for media in atp_medias:
            if isinstance(media, list):
                self.atp_medias.append(media)
            else:
                self.atp_medias.append([media, 0.01])
            self.media_hash[media.id] = media
        if "empty" not in self.media_hash:
            media = MSMedia.from_dict({})
            media.id = "empty"
            media.name = "empty"
            self.media_hash[media.id] = media
        
        self.forced_media = []
        for media_id in forced_media:
            for item in self.atp_medias:
                if item[0].id == media_id:
                    self.forced_media.append(item[0])
                    break

        self.max_gapfilling = max_gapfilling
        self.gapfilling_delta = gapfilling_delta

        if not core_template:
            self.load_default_template()
        else:
            self.coretemplate = core_template

        self.msgapfill = MSGapfill(
            self.modelutl,
            default_gapfill_templates=[core_template],
            default_target=self.atp_hydrolysis.id,
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
        self.noncore_reactions = []
        self.other_compartments = []
        # Iterating through reactions and disabling
        for reaction in self.model.reactions:
            gfrxn = self.msgapfill.gfmodel.reactions.get_by_id(reaction.id)
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
                    gfrxn.lower_bound = 0
                if reaction.upper_bound > 0 and template_reaction.upper_bound <= 0:
                    logger.debug(reaction.id + " core but reversible")
                    self.noncore_reactions.append([reaction, ">"])
                    reaction.upper_bound = 0
                    gfrxn.upper_bound = 0
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
                gfrxn.lower_bound = 0
                gfrxn.upper_bound = 0

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
            pkgmgr = MSPackageManager.get_pkg_mgr(self.model)
            # First prescreening model for ATP production without gapfilling
            media_list = []
            min_objectives = {}
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
                    media_list.append(media)
                    min_objectives[media] = minimum_obj
                elif solution.objective_value >= minimum_obj:
                    self.media_gapfill_stats[media] = {"reversed": {}, "new": {}}

            # Now running gapfilling on all conditions where initially there was no growth
            all_solutions = self.msgapfill.run_multi_gapfill(
                media_list,
                self.atp_hydrolysis.id,
                min_objectives,
                check_for_growth=False,
            )

            # Adding the new solutions to the media gapfill stats
            for media in all_solutions:
                self.media_gapfill_stats[media] = all_solutions[media]

        if MSATPCorrection.DEBUG:
            export_data = {}
            for media in self.media_gapfill_stats:
                export_data[media.id] = self.media_gapfill_stats[media]
            with open("debug.json", "w") as outfile:
                json.dump(export_data, outfile)

        return output

    def determine_growth_media(self, max_gapfilling=None):
        """
        Decides which of the test media to use as growth conditions for this model
        :return:
        """
        atp_att = {"tests": {}, "selected_media": {}, "core_atp_gapfilling": {}}
        self.selected_media = []
        best_score = None
        for media in self.media_gapfill_stats:
            atp_att["core_atp_gapfilling"][media.id] = {
                "score": 0,
                "new": {},
                "reversed": {},
            }
            if self.media_gapfill_stats[media]:
                atp_att["core_atp_gapfilling"][media.id]["score"] = len(
                    self.media_gapfill_stats[media]["new"].keys()
                ) + 0.5 * len(self.media_gapfill_stats[media]["reversed"].keys())
                atp_att["core_atp_gapfilling"][media.id][
                    "new"
                ] = self.media_gapfill_stats[media]["new"]
                atp_att["core_atp_gapfilling"][media.id][
                    "reversed"
                ] = self.media_gapfill_stats[media]["reversed"]
            else:
                atp_att["core_atp_gapfilling"][media.id] = {
                    "score": 1000,
                    "failed": True,
                }
            if (
                best_score is None
                or atp_att["core_atp_gapfilling"][media.id]["score"] < best_score
            ):
                best_score = atp_att["core_atp_gapfilling"][media.id]["score"]

        if self.max_gapfilling is None:
            self.max_gapfilling = best_score

        logger.info(f"max_gapfilling: {self.max_gapfilling}, best_score: {best_score}")

        for media in self.media_gapfill_stats:
            if atp_att["core_atp_gapfilling"][media.id][
                "score"
            ] <= self.max_gapfilling and atp_att["core_atp_gapfilling"][media.id][
                "score"
            ] <= (
                best_score + self.gapfilling_delta
            ):
                self.selected_media.append(media)
                atp_att["selected_media"][media.id] = 0

        self.modelutl.save_attributes(atp_att, "ATP_analysis")

    def apply_growth_media_gapfilling(self):
        """
        Applies the gapfilling to all selected growth media
        :return:
        """
        self.cumulative_core_gapfilling = []
        #  TODO: In case someone runs ATP correction twice with different parameters,
        #   before resetting this, maybe check if any of these reactions are already in
        #   the model and remove them so we're starting fresh???
        for media in self.selected_media:
            stats = self.media_gapfill_stats.get(media, None)
            if (
                stats is not None
                and MSGapfill.gapfill_count(self.media_gapfill_stats[media]) > 0
            ):
                self.msgapfill.integrate_gapfill_solution(
                    stats,
                    self.cumulative_core_gapfilling,
                    link_gaps_to_objective=False
                )
                #Adding reactions to gapfilling sensitivity structure so we can track all gapfilled reactions
                gf_sensitivity = self.modelutl.get_attributes("gf_sensitivity", {})
                if media.id not in gf_sensitivity:
                    gf_sensitivity[media.id] = {}
                if self.atp_hydrolysis.id not in gf_sensitivity[media.id]:
                    gf_sensitivity[media.id][self.atp_hydrolysis.id] = {}
                gf_sensitivity[media.id][self.atp_hydrolysis.id]["success"] = {}
                for item in stats["new"]:
                    gf_sensitivity[media.id][self.atp_hydrolysis.id]["success"][item] = {
                        stats["new"][item] : []
                    }
                for item in stats["reversed"]:
                    gf_sensitivity[media.id][self.atp_hydrolysis.id]["success"][item] = {
                        stats["reversed"][item] : []
                    }
                self.modelutl.save_attributes(gf_sensitivity, "gf_sensitivity")  
        self.modelutl.save_attributes(len(self.cumulative_core_gapfilling), "total_core_gapfilling")

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
            self.noncore_reactions, tests,atp_expansion=True
        )
        # Removing filtered reactions
        for item in self.filtered_noncore:
            logger.info("Removing " + item[0].id + " " + item[1])
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

    def build_tests(self,multiplier_hash_override={}):
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
        #Applying threshold multiplier
        for key in default_threshold_multipiers:
            if key not in multiplier_hash_override:
                multiplier_hash_override[key] = default_threshold_multipiers[key]
        #Initialzing atp test attributes
        atp_att = self.modelutl.get_attributes(
            "ATP_analysis",
            {"tests": {}, "selected_media": {}, "core_atp_gapfilling": {}},
        )
        #Initializing tests and adding empty media every time
        tests = []
        if "empty" in self.media_hash:
            tests.append(
                {
                    "media": self.media_hash["empty"],
                    "is_max_threshold": True,
                    "threshold": 0.00001,
                    "objective": self.atp_hydrolysis.id,
                }
            )
            atp_att["tests"]["empty"] = {
                "threshold": 0.00001,
                "objective": self.atp_hydrolysis.id,
            }
        #Setting objective to ATP hydrolysis
        self.model.objective = self.atp_hydrolysis.id
        for media in self.selected_media:
            #Setting multiplier for test threshold
            multiplier = multiplier_hash_override["default"]
            if media.id in multiplier_hash_override:
                 multiplier = multiplier_hash_override[media.id]
            #Constraining model exchanges for media
            self.modelutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
            #Computing core ATP production
            obj_value = self.model.slim_optimize()
            logger.debug(f"{media.name} = {obj_value};{multiplier}")
            threshold = multiplier * obj_value
            if threshold == 0:
                threshold += 0.00001
            tests.append(
                {
                    "media": media,
                    "is_max_threshold": True,
                    "threshold": threshold,
                    "objective": self.atp_hydrolysis.id,
                }
            )
            atp_att["selected_media"][media.id] = obj_value
            atp_att["tests"][media.id] = {
                "threshold": multiplier * obj_value,
                "objective": self.atp_hydrolysis.id,
            }
        #Saving test attributes to the model
        self.modelutl.save_attributes(atp_att, "ATP_analysis")
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
        # self.evaluate_growth_media()
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
