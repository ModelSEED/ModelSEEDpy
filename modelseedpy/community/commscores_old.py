from modelseedpy.core.exceptions import ObjectiveError, ParameterError
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.community.mscompatibility import MSCompatibility
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msgapfill import MSGapfill
from itertools import combinations, permutations, chain
from optlang import Variable, Constraint, Objective
from numpy import array, unique, ndarray, where, sort, array_split, nan
from collections import Counter
from deepdiff import DeepDiff  # (old, new)
from typing import Iterable, Union
from pprint import pprint
from numpy.random import shuffle
from multiprocess import current_process
from math import inf
import sigfig

# from icecream import ic
import re

# from math import prod

# silence deprecation warnings from DeepDiff parsing the syntrophy
import warnings

warnings.simplefilter("ignore", category=DeprecationWarning)

rm_comp = FBAHelper.remove_compartment


def _compatibilize(member_models: Iterable, printing=False):
    # return member_models
    models = MSCompatibility.standardize(
        member_models, conflicts_file_name="exchanges_conflicts.json", printing=printing
    )
    if not isinstance(member_models, (set, list, tuple)):
        return models[0]
    return models


def _load_models(
    member_models: Iterable, com_model=None, compatibilize=True, printing=False
):
    # ic(member_models, com_model, compatibilize)
    if not com_model and member_models:
        model = build_from_species_models(member_models, name="SMETANA_pair")
        return member_models, model  # (model, names=names, abundances=abundances)
    # models = PARSING_FUNCTION(community_model)  # TODO the individual models of a community model can be parsed
    if compatibilize:
        return (
            _compatibilize(member_models, printing),
            _compatibilize([com_model], printing)[0],
        )
    return member_models, com_model


def _get_media(
    media=None,
    com_model=None,
    model_s_=None,
    min_growth=None,
    environment=None,
    interacting=True,
    printing=False,
    minimization_method="minFlux",
    skip_bad_media=False,
):
    # ic(media, com_model, model_s_)
    if com_model is None and model_s_ is None:
        raise TypeError("< com_model > or < model_s_ > must be parameterized.")
    if media is not None:
        if model_s_ is not None and not isinstance(model_s_, (list, set, tuple)):
            return media["members"][model_s_.id]["media"]
        elif com_model is not None:
            return media["community_media"]
        return media
    # model_s_ is either a singular model or a list of models
    if com_model is not None:
        try:
            com_media, media_sol = MSMinimalMedia.determine_min_media(
                com_model,
                minimization_method,
                min_growth,
                None,
                interacting,
                5,
                printing,
            )
        except Exception as e:
            if skip_bad_media:
                com_media, media_sol = None, None
            else:
                print(e)
    if model_s_ is not None:
        if not isinstance(model_s_, (list, set, tuple, ndarray)):
            try:
                return MSMinimalMedia.determine_min_media(
                    model_s_,
                    minimization_method,
                    min_growth,
                    environment,
                    interacting,
                    printing,
                )
            except Exception as e:
                if not skip_bad_media:
                    print(e)
                return None
        members_media = {}
        for model in model_s_:
            try:
                members_media[model.id] = {
                    "media": MSMinimalMedia.determine_min_media(
                        model,
                        minimization_method,
                        min_growth,
                        environment,
                        interacting,
                        printing,
                    )[0]
                }
                continue
            except Exception as e:
                if skip_bad_media:
                    continue
                else:
                    print(e)
        # print(members_media)
        if com_model is None:
            return members_media
    else:
        return com_media, media_sol
    return {"community_media": com_media, "members": members_media}


def _sigfig_check(value, sigfigs, default):
    if str(value) in ["inf", "nan"]:
        value = ""
    if FBAHelper.isnumber(value):
        return sigfig.round(value, sigfigs)
    else:
        return default


def nanFilter(value, string=True):
    if isinstance(value, str) or value is None:
        if string:
            return value
        else:
            return nan
    if any([value < 0, value > 1e5]):
        return "" if string else nan
    return value


class CommScores:
    def __init__(
        self,
        member_models,
        min_growth=0.1,
        n_solutions=100,
        environment=None,
        abstol=1e-3,
        media_dict=None,
        printing=True,
        raw_content=False,
        antismash_json_path: str = None,
        antismash_zip_path: str = None,
        minimal_media_method="minFlux",
    ):
        self.min_growth = min_growth
        self.abstol = abstol
        self.n_solutions = n_solutions
        self.printing = printing
        self.raw_content = raw_content
        self.antismash_json_path = antismash_json_path
        self.antismash_zip_path = antismash_zip_path

        # process the models
        self.models = _compatibilize(member_models)
        self.community = MSModelUtil(build_from_species_models(self.models))
        ## define the environment
        if environment:
            if hasattr(environment, "get_media_constraints"):
                ### standardize modelseed media into COBRApy media
                environment = {
                    "EX_" + exID: -bound[0]
                    for exID, bound in environment.get_media_constraints().items()
                }
            self.community.add_medium(environment)
        self.environment = environment
        ## test growth
        for model in self.models:
            if model.slim_optimize() == 0:
                raise ObjectiveError(
                    f"The model {model.id} possesses an objective value of 0 in complete media, "
                    "which is incompatible with minimal media computations and hence SMETANA."
                )
        if self.community.model.slim_optimize() == 0:
            raise ObjectiveError(
                f"The community model {self.community.model.id} possesses an objective "
                "value of 0 in complete media, which is incompatible with minimal "
                "media computations and hence SMETANA."
            )
        ## determine the minimal media for each model, including the community
        self.media = (
            media_dict
            if media_dict
            else MSMinimalMedia.comm_media_est(
                member_models,
                self.community.model,
                minimal_media_method,
                min_growth,
                self.environment,
                True,
                n_solutions,
                printing,
            )
        )

    def all_scores(
        self,
        mp_score=True,
        kbase_obj=None,
        cobrakbase_path: str = None,
        kbase_token_path: str = None,
        annotated_genomes: dict = None,
    ):
        mro = self.mro_score()
        mip = self.mip_score(interacting_media=self.media)
        mp = None if not mp_score else self.mp_score()
        mu = None  # self.mu_score()
        sc = None  # self.sc_score()
        smetana = None  # self.smetana_score()
        gyd = self.gyd_score()
        fs = (
            self.fs_score()
            if any(
                [
                    kbase_obj is not None,
                    annotated_genomes != [],
                    cobrakbase_path is not None and kbase_token_path is not None,
                ]
            )
            else None
        )
        return {
            "mro": mro,
            "mip": mip,
            "mp": mp,
            "mu": mu,
            "sc": sc,
            "smetana": smetana,
            "gyd": gyd,
            "fs": fs,
        }

    def mro_score(self):
        self.mro_val = CommScores.mro(
            self.models,
            self.media["members"],
            self.min_growth,
            self.media,
            self.raw_content,
            self.environment,
            self.printing,
            True,
        )
        if not self.printing:
            return self.mro_val
        if self.raw_content:
            for pair, (interaction, media) in self.mro_val.items():
                newcomer, established = pair.split("---")
                print(
                    f"\n(MRO) The {newcomer} media {media} possesses {interaction} shared "
                    f"requirements with the {established} established member."
                )
                return self.mro_val
        for pair, mro in self.mro_val.items():
            newcomer, established = pair.split("---")
            print(
                f"\nThe {newcomer} on {established} MRO score: {mro[0]} ({mro[0]*100:.2f}%). "
                f"This is the percent of nutritional requirements in {newcomer} "
                f"that overlap with {established} ({mro[1]}/{mro[2]})."
            )
        return self.mro_val

    def mip_score(
        self, interacting_media: dict = None, noninteracting_media: dict = None
    ):
        interacting_media = interacting_media or self.media or None
        diff, self.mip_val = CommScores.mip(
            self.models,
            self.community.model,
            self.min_growth,
            interacting_media,
            noninteracting_media,
            self.environment,
            self.printing,
            True,
        )
        if not self.printing:
            return self.mip_val
        print(
            f"\nMIP score: {self.mip_val}\t\t\t{self.mip_val} required compound(s) can be sourced via syntrophy:"
        )
        if self.raw_content:
            pprint(diff)
        return self.mip_val

    def gyd_score(self, coculture_growth=False):
        self.gyd_val = CommScores.gyd(
            self.models, environment=self.environment, coculture_growth=coculture_growth
        )
        if not self.printing:
            return self.gyd
        growth_type = "monocultural" if not coculture_growth else "cocultural"
        for pair, score in self.gyd_val.items():
            print(
                f"\nGYD score: The {growth_type} growth difference between the {pair} member models"
                f" is {score} times greater than the growth of the slower member."
            )
        return self.gyd

    def fs_score(
        self,
        kbase_obj=None,
        cobrakbase_path: str = None,
        kbase_token_path: str = None,
        annotated_genomes: dict = None,
    ):
        self.fs_val = CommScores.fs(
            self.models, kbase_obj, cobrakbase_path, kbase_token_path, annotated_genomes
        )
        if not self.printing:
            return self.fs
        for pair, score in self.fs_val.items():
            print(
                f"\nFS Score: The similarity of RAST functional SSO ontology "
                f"terms between the {pair} members is {score}."
            )
        return self.fs

    def mp_score(self):
        print("executing MP")
        self.mp_val = CommScores.mp(
            self.models,
            self.environment,
            self.community.model,
            None,
            self.abstol,
            self.printing,
        )
        if not self.printing:
            return self.mp_val
        if self.raw_content:
            print(
                "\n(MP) The possible contributions of each member in the member media include:\n"
            )
            pprint(self.mp_val)
        else:
            print(
                "\nMP score:\t\t\tEach member can possibly contribute the following to the community:\n"
            )
            for member, contributions in self.mp_val.items():
                print(member, "\t", len(contributions))
        return self.mp_val

    def mu_score(self):
        member_excreta = self.mp_score() if not hasattr(self, "mp_val") else self.mp_val
        self.mu_val = CommScores.mu(
            self.models,
            self.environment,
            member_excreta,
            self.n_solutions,
            self.abstol,
            True,
            self.printing,
        )
        if not self.printing:
            return self.mu_val
        print(
            "\nMU score:\t\t\tThe fraction of solutions in which each member is the "
            "syntrophic receiver that contain a respective metabolite:\n"
        )
        pprint(self.mu_val)
        return self.mu_val

    def sc_score(self):
        self.sc_val = CommScores.sc(
            self.models,
            self.community.model,
            self.min_growth,
            self.n_solutions,
            self.abstol,
            True,
            self.printing,
        )
        if not self.printing:
            return self.sc_val
        print(
            "\nSC score:\t\t\tThe fraction of community members who syntrophically contribute to each species:\n"
        )
        pprint(self.sc_val)
        return self.sc_val

    def smetana_score(self):
        if not hasattr(self, "sc_val"):
            self.sc_val = self.sc_score()
        sc_coupling = all(array(list(self.sc.values())) is not None)
        if not hasattr(self, "mu_val"):
            self.mu_val = self.mu_score()
        if not hasattr(self, "mp_val"):
            self.mp_val = self.mp_score()

        self.smetana = CommScores.smetana(
            self.models,
            self.community.model,
            self.min_growth,
            self.n_solutions,
            self.abstol,
            (self.sc_val, self.mu_val, self.mp_val),
            True,
            sc_coupling,
            self.printing,
        )
        if self.printing:
            print("\nsmetana score:\n")
            pprint(self.smetana)
        return self.smetana

    def antiSMASH_scores(self, antismash_json_path=None):
        self.antismash = CommScores.antiSMASH(
            antismash_json_path or self.antismash_json_path
        )
        if not self.printing:
            return self.antismash
        if self.raw_content:
            print(
                "\n(antismash) The biosynthetic_areas, BGCs, protein_annotations, clusterBlast, and "
                "num_clusterBlast from the provided antiSMASH results:\n"
            )
            print(
                "The 'areas' that antiSMASH determines produce biosynthetic products:"
            )
            pprint(self.antismash[0])
            print("The set of biosynthetic gene clusters:")
            pprint(self.antismash[1])
            print("The set of clusterblast protein annotations:")
            pprint(self.antismash[2])
            print("Resistance information from clusterblast")
            pprint(self.antismash[3])
            print("The number of proteins associated with resistance")
            pprint(self.antismash[4])
            return self.antismash
        print("\nantiSMASH scores:\n")
        print(
            "The community exhibited:"
            f"- {len(self.antismash[0])}'areas' that antiSMASH determines produce biosynthetic products."
            f"- {len(self.antismash[1])} biosynthetic gene clusters."
            f"- {len(self.antismash[2])} clusterblast protein annotations."
            f"- {len(self.antismash[3])} parcels of resistance information from clusterblast."
            f"- {self.antismash[4]} proteins associated with resistance."
        )
        return list(map(len, self.antismash[:4])) + [self.antismash[4]]

    ###### STATIC METHODS OF THE SMETANA SCORES, WHICH ARE APPLIED IN THE ABOVE CLASS OBJECT ######

    @staticmethod
    def _check_model(model_util, media, model_str, skip_bad_media):
        default_media = model_util.model.medium
        if media is not None:
            model_util.add_medium(media)
        obj_val = model_util.model.slim_optimize()
        if obj_val == 0 or not FBAHelper.isnumber(obj_val):
            print(
                f"The {model_str} model input does not yield an operational model, and will therefore be gapfilled."
            )
            # if not skip_bad_media:  return MSGapfill.gapfill(model_util.model, media)
        model_util.add_medium(default_media)
        return model_util.model

    @staticmethod
    def _load(model, kbase_obj):
        model_str = model
        if len(model) == 2:
            model = kbase_obj.get_from_ws(*model)
        else:
            model = kbase_obj.get_from_ws(model)
        return model, model_str

    @staticmethod
    def _determine_growths(modelUtils):
        return [util.model.slim_optimize() for util in modelUtils]

    @staticmethod
    def calculate_scores(
        pairs,
        models_media=None,
        environments=None,
        annotated_genomes=True,
        lazy_load=False,
        kbase_obj=None,
        cip_score=True,
        costless=True,
        skip_bad_media=False,
        anme_comm=False,
        print_progress=False,
    ):
        from pandas import Series

        if isinstance(pairs, list):
            (
                pairs,
                models_media,
                environments,
                annotated_genomes,
                lazy_load,
                kbase_obj,
            ) = pairs
        series, mets = [], []
        if not isinstance(environments, (list, tuple)):
            environments = [environments]
        if isinstance(environments, (list, tuple)) and hasattr(environments[0], "name"):
            environments = {
                m.name: FBAHelper.convert_kbase_media(m, 1000) for m in environments
            }
        elif not isinstance(environments, dict):
            environments = {f"media{i}": m for i, m in enumerate(environments)}
        pid = current_process().name
        model_utils = {}
        count = 0
        for model1, models in pairs.items():
            if model1.id == "":
                model1.id = "model1"
            if lazy_load:
                model1, model1_str = CommScores._load(model1, kbase_obj)
            else:
                model1_str = model1.id
            if model1.id not in models_media:
                models_media[model1.id] = {
                    "media": _get_media(model_s_=model1, skip_bad_media=skip_bad_media)
                }
                if models_media[model1.id] is None:
                    continue
            if model1.id not in model_utils:
                model_utils[model1.id] = MSModelUtil(model1)
            # print(pid, model1)
            for model_index, model2 in enumerate(models):
                if model2.id == "":
                    model2.id = "model2"
                if lazy_load:
                    model2, model2_str = CommScores._load(model2, kbase_obj)
                else:
                    model2_str = model2.id
                if model2.id not in models_media:
                    models_media[model2.id] = {
                        "media": _get_media(
                            model_s_=model2, skip_bad_media=skip_bad_media
                        )
                    }
                    if models_media[model2.id] is None:
                        continue
                if model2.id not in model_utils:
                    model_utils[model2.id] = MSModelUtil(model2)
                grouping = [model1, model2]
                grouping_utils = [model_utils[model1.id], model_utils[model2.id]]
                modelIDs = [model.id for model in grouping]
                comm_model = build_from_species_models(grouping)
                community = MSCommunity(comm_model, ids=modelIDs)
                comm_sol = comm_model.optimize()
                print(f"{pid}~~{count}\t{modelIDs}")
                for environName, environ in environments.items():
                    if print_progress:
                        print(f"\tEnvironment\t{environName}", end="\t")
                    if not anme_comm:
                        model1 = CommScores._check_model(
                            model_utils[model1.id], environ, model1_str, skip_bad_media
                        )
                        model2 = CommScores._check_model(
                            model_utils[model2.id], environ, model2_str, skip_bad_media
                        )
                    # initiate the KBase output
                    report_dic = {
                        f"model{i+1}": modelID for i, modelID in enumerate(modelIDs)
                    }
                    g1, g2, comm = CommScores._determine_growths(
                        [model_utils[model1.id], model_utils[model2.id], community.util]
                    )
                    g1, g2, comm = (
                        _sigfig_check(g1, 5, ""),
                        _sigfig_check(g2, 5, ""),
                        _sigfig_check(comm, 5, ""),
                    )
                    report_dic.update(
                        {
                            "media": environName,
                            "model1 growth": g1,
                            "model2 growth": g2,
                            "community growth": comm,
                        }
                    )
                    coculture_growths = {
                        mem.id: comm_sol.fluxes[mem.primary_biomass.id]
                        for mem in community.members
                    }
                    report_dic.update(
                        {
                            f"coculture growth model{modelIDs.index(memID)}": growth
                            for memID, growth in coculture_growths.items()
                        }
                    )
                    # define the MRO content
                    mro_values = CommScores.mro(
                        grouping, models_media, raw_content=True, environment=environ
                    )
                    report_dic.update(
                        {
                            f"MRO_model{modelIDs.index(models_string.split('--')[0])+1}": f"{100*len(intersection)/len(memMedia):.3f}% ({len(intersection)}/{len(memMedia)})"
                            for models_string, (
                                intersection,
                                memMedia,
                            ) in mro_values.items()
                        }
                    )
                    mets.append({"MRO metabolites": list(mro_values.values())[0][0]})
                    if print_progress:
                        print("MRO done", end="\t")
                    # define the CIP content
                    if cip_score:
                        cip_values = CommScores.cip(
                            modelutils=[model_utils[mem.id] for mem in grouping]
                        )
                        report_dic.update({"CIP": cip_values[1]})
                        mets[-1].update({"CIP metabolites": list(cip_values[0])})
                        if print_progress:
                            print("CIP done", end="\t")
                    # define the MIP content
                    mip_values = CommScores.mip(
                        grouping,
                        comm_model,
                        0.1,
                        None,
                        None,
                        environ,
                        print_progress,
                        True,
                        costless,
                        costless,
                        skip_bad_media,
                    )
                    # print(mip_values)
                    if mip_values is not None:
                        report_dic.update(
                            {
                                f"MIP_model{modelIDs.index(models_name)+1}": str(
                                    len(received)
                                )
                                for models_name, received in mip_values[0].items()
                            }
                        )
                        mets[-1].update(
                            {
                                "MIP model1 metabolites": list(mip_values[0].values())[
                                    0
                                ],
                                "MIP model2 metabolites": list(mip_values[0].values())[
                                    1
                                ],
                            }
                        )
                        if costless:
                            for models_name, received in mip_values[1].items():
                                report_dic[
                                    f"MIP_model{modelIDs.index(models_name)+1} (costless)"
                                ] = (
                                    report_dic[
                                        f"MIP_model{modelIDs.index(models_name)+1}"
                                    ]
                                    + f" ({len(received)})"
                                )
                                del report_dic[
                                    f"MIP_model{modelIDs.index(models_name)+1}"
                                ]
                            if print_progress:
                                print("costless_MIP  done", end="\t")
                    else:
                        report_dic.update(
                            {f"MIP_model1 (costless)": "", f"MIP_model2 (costless)": ""}
                        )
                        mets[-1].update(
                            {
                                "MIP model1 metabolites": [None],
                                "MIP model2 metabolites": [None],
                            }
                        )
                    if print_progress:
                        print("MIP done", end="\t")
                    # define the BSS content
                    bss_values = CommScores.bss(
                        grouping,
                        grouping_utils,
                        environments,
                        models_media,
                        skip_bad_media,
                    )
                    report_dic.update(
                        {
                            f"BSS_model{modelIDs.index(name.split(' supporting ')[0])+1}": f"{_sigfig_check(100*val, 5, '')}%"
                            for name, (mets, val) in bss_values.items()
                        }
                    )
                    mets[-1].update(
                        {
                            "BSS model1 metabolites": [
                                met_set for met_set, val in bss_values.values()
                            ][0],
                            "BSS model2 metabolites": [
                                met_set for met_set, val in bss_values.values()
                            ][1],
                        }
                    )
                    # mets[-1].update({"bss_mets": list(bss_values[0].values())})
                    if print_progress:
                        print("BSS done", end="\t")
                    # define the PC content
                    pc_values = CommScores.pc(
                        grouping,
                        grouping_utils,
                        comm_model,
                        None,
                        comm_sol,
                        environ,
                        True,
                        community,
                    )
                    report_dic.update(
                        {
                            "PC_comm": _sigfig_check(pc_values[0], 5, ""),
                            "PC_model1": _sigfig_check(
                                list(pc_values[1].values())[0], 5, ""
                            ),
                            "PC_model2": _sigfig_check(
                                list(pc_values[1].values())[1], 5, ""
                            ),
                            "BIT": pc_values[3],
                        }
                    )
                    if print_progress:
                        print("PC  done\tBIT done", end="\t")
                    # print([mem.slim_optimize() for mem in grouping])
                    # define the GYD content
                    gyd1, gyd2, g1, g2 = list(
                        CommScores.gyd(
                            grouping,
                            grouping_utils,
                            environ,
                            False,
                            community,
                            anme_comm,
                        ).values()
                    )[0]
                    report_dic.update(
                        {
                            "GYD1": _sigfig_check(gyd1, 5, ""),
                            "GYD2": _sigfig_check(gyd2, 5, ""),
                        }
                    )
                    if print_progress:
                        print("GYD done\t\t", end="\t" if annotated_genomes else "\n")
                    # define the FS content
                    if kbase_obj is not None and annotated_genomes and not anme_comm:
                        fs_values = list(
                            CommScores.fs(
                                grouping, kbase_obj, annotated_genomes=annotated_genomes
                            ).values()
                        )[0]
                        print(
                            len(fs_values[0]) if fs_values[0] is not None else "NaN",
                            fs_values[1],
                        )
                        report_dic.update({"FS": sigfig.round(fs_values[1], 5)})
                        if fs_values is not None:
                            mets[-1].update({"FS features": fs_values[0]})
                        if print_progress:
                            print("FS done\t\t")
                    # return a pandas Series, which can be easily aggregated with other results into a DataFrame
                    series.append(Series(report_dic))
                count += 1
        return series, mets

    @staticmethod
    def html_report(
        df, mets, export_html_path="commscores_report.html", msdb_path=None
    ):
        from modelseedpy.core.report import commscores_report

        return commscores_report(df, mets, export_html_path, msdb_path)

    @staticmethod
    def report_generation(
        all_models: iter = None,  # a list of distinct lists is provided for specifying exclusive groups
        pairs: dict = None,
        mem_media: dict = None,
        pair_limit: int = None,
        exclude_pairs: list = None,
        kbase_obj=None,
        annotated_genomes: dict = True,  # True triggers internal acquisition of the genomes, where None skips
        see_media=True,
        environments: iter = None,  # a collection of environment dicts or KBase media objects
        pool_size: int = None,
        cip_score=True,
        costless=True,
        skip_bad_media=False,
        anme_comm=False,
        print_progress=False,
    ):
        from pandas import concat

        if pairs:
            model_pairs = unique(
                [
                    {model1, model2}
                    for model1, models in pairs.items()
                    for model2 in models
                ]
            )
        elif all_models is not None:
            if not isinstance(all_models[0], list):
                all_models = list(set(all_models))
                model_pairs = array(list(combinations(all_models, 2)))
            else:
                model_pairs = []
                for models1, models2 in combinations(all_models, 2):
                    models1 = set(models1)
                    models2 = set(models2)
                    if len(models1) > len(models2):
                        larger_list = models1
                        smaller_list = models2
                    else:
                        larger_list = models2
                        smaller_list = models1
                    model_pairs.append(
                        [
                            list(zip(combin, smaller_list))
                            for combin in permutations(larger_list, len(smaller_list))
                        ]
                    )
                # flatten the assembled pairs and filter duplicates
                model_pairs = array(
                    [
                        x
                        for x in set(
                            tuple(x)
                            for x in [
                                i
                                for y in list(chain.from_iterable(model_pairs))
                                for i in y
                            ]
                        )
                    ]
                )
                all_models = list(chain.from_iterable(all_models))
            if pair_limit is not None:
                shuffle(model_pairs)
                new_pairs = []
                for index, pair in enumerate(model_pairs):
                    if set(pair) not in exclude_pairs and index < pair_limit:
                        new_pairs.append(pair)
                    elif index >= pair_limit:
                        break
                model_pairs = array(new_pairs)
            if isinstance(model_pairs[0], str):
                model_pairs = unique(sort(model_pairs, axis=1))
            pairs = {
                first: model_pairs[where(model_pairs[:, 0] == first)][:, 1]
                for first in model_pairs[:, 0]
            }
        else:
            raise ValueError(
                "Either < all_models > or < pairs > must be defined to simulate interactions."
            )
        if not all_models:
            all_models = list(
                chain(*[list(values) for values in pairs.values()])
            ) + list(pairs.keys())
        lazy_load = len(model_pairs) > 10000  # all_models[0], (list,set,tuple))
        if lazy_load and not kbase_obj:
            ValueError(
                "The < kbase_obj > argument must be provided to lazy load models."
            )
        new_models = []
        for index, model in enumerate(all_models):
            if model.id == "":
                model.id = f"model_index{index}"
            new_models.append(model)
        all_models = new_models[:]
        if not mem_media:
            models_media = _get_media(
                model_s_=all_models, skip_bad_media=skip_bad_media
            )
        else:
            models_media = mem_media.copy()
            missing_models = set()
            missing_modelID = []
            for model in all_models:
                if model is not None and model.id not in models_media:
                    missing_models.add(model)
                    missing_modelID.append(
                        model if not hasattr(model, "id") else model.id
                    )
            if missing_models != set():
                print(
                    f"Media of the {missing_modelID} models are not defined, and will be calculated separately."
                )
                models_media.update(
                    _get_media(model_s_=missing_models), skip_bad_media=skip_bad_media
                )
        if see_media:
            print(f"The minimal media of all members:\n{models_media}")
        print(f"\nExamining the {len(list(model_pairs))} model pairs")
        if pool_size is not None:
            from datetime import datetime
            from multiprocess import Pool

            print(
                f"Loading {int(pool_size)} workers and computing the scores",
                datetime.now(),
            )
            pool = Pool(
                int(pool_size)
            )  # .map(calculate_scores, [{k: v} for k,v in pairs.items()])
            args = [
                [
                    dict([pair]),
                    models_media,
                    environments,
                    annotated_genomes,
                    lazy_load,
                    kbase_obj,
                ]
                for pair in list(pairs.items())
            ]
            output = pool.map(CommScores.calculate_scores, args)
            series = chain.from_iterable([ele[0] for ele in output])
            mets = chain.from_iterable([ele[1] for ele in output])
        else:
            series, mets = CommScores.calculate_scores(
                pairs,
                models_media,
                environments,
                annotated_genomes,
                lazy_load,
                kbase_obj,
                cip_score,
                costless,
                skip_bad_media,
                anme_comm,
                print_progress,
            )
        return concat(series, axis=1).T, mets

    @staticmethod
    def mro(
        member_models: Iterable = None,
        mem_media: dict = None,
        min_growth=0.1,
        media_dict=None,
        raw_content=False,
        environment=None,
        skip_bad_media=False,
        printing=False,
        compatibilized=False,
    ):
        """Determine the overlap of nutritional requirements (minimal media) between member organisms."""
        # determine the member minimal media if they are not parameterized
        if not mem_media:
            if not member_models:
                raise ParameterError(
                    "The either member_models or minimal_media parameter must be defined."
                )
            member_models = (
                member_models
                if compatibilized
                else _compatibilize(member_models, printing)
            )
            mem_media = _get_media(
                media_dict,
                None,
                member_models,
                min_growth,
                environment,
                printing=printing,
                skip_bad_media=skip_bad_media,
            )
            if "community_media" in mem_media:
                mem_media = mem_media["members"]
        # MROs = array(list(map(len, pairs.values()))) / array(list(map(len, mem_media.values())))
        mro_values = {}
        for model1, model2 in combinations(member_models, 2):
            intersection = set(mem_media[model1.id]["media"].keys()) & set(
                mem_media[model2.id]["media"].keys()
            )
            inter = [ex.replace("EX_", "").replace("_e0", "") for ex in intersection]
            m1_media = mem_media[model1.id]["media"]
            m2_media = mem_media[model2.id]["media"]
            if raw_content:
                mro_values.update(
                    {
                        f"{model1.id}---{model2.id})": (inter, m1_media),
                        f"{model2.id}---{model1.id})": (inter, m2_media),
                    }
                )
            else:
                mro_values.update(
                    {
                        f"{model1.id}---{model2.id})": 100
                        * (len(inter) / len(m1_media), len(inter), len(m1_media)),
                        f"{model2.id}---{model1.id})": 100
                        * (len(inter) / len(m2_media), len(inter), len(m2_media)),
                        "mets": inter,
                    }
                )
        return mro_values
        # return mean(list(map(len, pairs.values()))) / mean(list(map(len, mem_media.values())))

    @staticmethod
    def mip(
        member_models: Iterable,
        com_model=None,
        min_growth=0.1,
        interacting_media_dict=None,
        noninteracting_media_dict=None,
        environment=None,
        printing=False,
        compatibilized=False,
        costless=False,
        multi_output=False,
        skip_bad_media=False,
    ):
        """Determine the quantity of nutrients that can be potentially sourced through syntrophy"""
        member_models, community = _load_models(
            member_models, com_model, not compatibilized, printing=printing
        )
        # determine the interacting and non-interacting media for the specified community  .util.model
        noninteracting_medium, noninteracting_sol = _get_media(
            noninteracting_media_dict,
            community,
            None,
            min_growth,
            environment,
            False,
            skip_bad_media=skip_bad_media,
        )
        if noninteracting_medium is None:
            return None
        if "community_media" in noninteracting_medium:
            noninteracting_medium = noninteracting_medium["community_media"]
        interacting_medium, interacting_sol = _get_media(
            interacting_media_dict,
            community,
            None,
            min_growth,
            environment,
            True,
            skip_bad_media=skip_bad_media,
        )
        if interacting_medium is None:
            return None
        if "community_media" in interacting_medium:
            interacting_medium = interacting_medium["community_media"]
        interact_diff = DeepDiff(noninteracting_medium, interacting_medium)
        if "dictionary_item_removed" not in interact_diff:
            return None
        cross_fed_exIDs = [
            re.sub("(root\['|'\])", "", x)
            for x in interact_diff["dictionary_item_removed"]
        ]
        # Determine each direction of the MIP score interactions
        comm_util = MSModelUtil(community)
        cross_fed_metIDs = [
            ex.replace("EX_", "").replace("_e0", "") for ex in cross_fed_exIDs
        ]
        cross_fed_copy = cross_fed_metIDs[:]
        directionalMIP = {mem.id: [] for mem in member_models}
        for rxn in comm_util.transport_list():
            # print(rxn.reaction, "\t", [met.id for met in rxn.metabolites if "_e0" in met.id])
            metIDs = list(
                set([met.id.split("_")[0] for met in rxn.reactants]).intersection(
                    set([met.id.split("_")[0] for met in rxn.products])
                )
            )
            if len(metIDs) == 1:
                metID = metIDs[0]
            else:
                if "cpd00067" in metIDs:
                    metIDs.remove("cpd00067")
                metID = metIDs[0]
            if metID not in cross_fed_metIDs:
                continue
            rxn_index = FBAHelper.compartment_index(rxn.id.split("_")[-1])
            if rxn_index == 0:
                continue
            mets = [met for met in rxn.metabolites if met.id == f"{metID}_c{rxn_index}"]
            if mets == []:
                print(f"The {metID}_c{rxn_index} is missing in {rxn.reaction}.")
                continue
            rxn_model = member_models[rxn_index - 1]
            # comm_trans[metID] = comm_trans.get(f"{metID}_c{rxn_index}", {})
            if (
                rxn.metabolites[mets[0]] > 0
                and interacting_sol.fluxes[rxn.id] > 0
                or rxn.metabolites[mets[0]] < 0
                and interacting_sol.fluxes[rxn.id] < 0
            ):  # donor
                directionalMIP[rxn_model.id].append(metID)
                if metID in cross_fed_copy:
                    cross_fed_copy.remove(metID)
                    continue
            # if printing:  print(f"{mets[0]} in {rxn.id} ({rxn.reaction}) is not assigned a receiving member.")
        if cross_fed_copy != [] and printing:
            print(f"Missing directions for the {cross_fed_copy} cross-fed metabolites")
        outputs = [directionalMIP]
        # TODO categorize all of the cross-fed substrates to examine potential associations of specific compounds
        if costless:
            costless_mets, numExs = CommScores.cip(member_models=member_models)
            # print(list(directionalMIP.values()), costless_mets)
            costlessDirectionalMIP = {
                member_name: set(receive_mets).intersection(costless_mets)
                for member_name, receive_mets in directionalMIP.items()
            }
            if not multi_output:
                return costlessDirectionalMIP
            outputs.append(costlessDirectionalMIP)
        return outputs

    @staticmethod
    def cip(modelutils=None, member_models=None):  # costless interaction potential
        if not modelutils:
            modelutils = {MSModelUtil(model) for model in member_models}
        costless_mets = set(
            chain.from_iterable(
                [modelutil.costless_excreta() for modelutil in modelutils]
            )
        )
        return costless_mets, len(costless_mets)

    @staticmethod
    def contributions(org_possible_contributions, scores, model_util, abstol):
        # identify and log excreta from the solution
        model_util.add_objective(
            sum(ex_rxn.flux_expression for ex_rxn in org_possible_contributions)
        )
        sol = model_util.model.optimize()
        if sol.status != "optimal":
            # exit the while loop by returning the original possible_contributions,
            ## hence DeepDiff == {} and the while loop terminates
            return scores, org_possible_contributions
        # identify and log excreta from the solution
        possible_contributions = org_possible_contributions[:]
        for ex in org_possible_contributions:
            if ex.id in sol.fluxes.keys() and sol.fluxes[ex.id] >= abstol:
                possible_contributions.remove(ex)
                scores[model_util.model.id].update([met.id for met in ex.metabolites])
        return scores, possible_contributions

    @staticmethod
    def mp(
        member_models: Iterable,
        environment,
        com_model=None,
        minimal_media=None,
        abstol=1e-3,
        printing=False,
    ):
        """Discover the metabolites that each species can contribute to a community"""
        community = (
            _compatibilize(com_model)
            if com_model
            else build_from_species_models(member_models, standardize=True)
        )
        community.medium = minimal_media or MSMinimalMedia.minimize_flux(community)
        scores = {}
        for (
            org_model
        ) in (
            member_models
        ):  # TODO support parsing the individual members through the MSCommunity object
            model_util = MSModelUtil(org_model)
            model_util.compatibilize(printing=printing)
            if environment:
                model_util.add_medium(environment)
            scores[model_util.model.id] = set()
            # determines possible member contributions in the community environment, where the excretion of media compounds is irrelevant
            org_possible_contr = [
                ex_rxn
                for ex_rxn in model_util.exchange_list()
                if (ex_rxn.id not in community.medium and ex_rxn.upper_bound > 0)
            ]
            # ic(org_possible_contributions, len(model_util.exchange_list()), len(community.medium))
            scores, possible_contr = CommScores.contributions(
                org_possible_contr, scores, model_util, abstol
            )
            while DeepDiff(org_possible_contr, possible_contr):
                print("remaining possible_contributions", len(possible_contr), end="\r")
                ## optimize the sum of the remaining exchanges that have not surpassed the abstol
                org_possible_contr = possible_contr[:]
                scores, possible_contr = CommScores.contributions(
                    org_possible_contr, scores, model_util, abstol
                )

            ## individually checks the remaining possible contributions
            for ex_rxn in possible_contr:
                model_util.model.objective = Objective(ex_rxn.flux_expression)
                sol = model_util.model.optimize()
                if sol.status == "optimal" or sol.objective_value > abstol:
                    for met in ex_rxn.metabolites:
                        if met.id in scores[model_util.model.id]:
                            scores[model_util.model.id].remove(met.id)
                            print("removing", met.id)
        return scores

    @staticmethod
    def mu(
        member_models: Iterable,
        environment=None,
        member_excreta=None,
        n_solutions=100,
        abstol=1e-3,
        compatibilized=False,
        printing=True,
    ):
        """the fractional frequency of each received metabolite amongst all possible alternative syntrophic solutions"""
        # member_solutions = member_solutions if member_solutions else {model.id: model.optimize() for model in member_models}
        scores = {}
        member_models = (
            member_models if compatibilized else _compatibilize(member_models, printing)
        )
        if member_excreta:
            missing_members = [
                model for model in member_models if model.id not in member_excreta
            ]
            if missing_members:
                print(
                    f"The {','.join(missing_members)} members are missing from the defined "
                    f"excreta list and will therefore be determined through an additional MP simulation."
                )
                member_excreta.update(CommScores.mp(missing_members, environment))
        else:
            member_excreta = CommScores.mp(
                member_models, environment, None, abstol, printing
            )
        for org_model in member_models:
            other_excreta = set(
                chain.from_iterable(
                    [
                        excreta
                        for model, excreta in member_excreta.items()
                        if model != org_model.id
                    ]
                )
            )
            print(f"\n{org_model.id}\tOther Excreta", other_excreta)
            model_util = MSModelUtil(org_model, True)
            if environment:
                model_util.add_medium(environment)
            ex_rxns = {
                ex_rxn: list(ex_rxn.metabolites)[0]
                for ex_rxn in model_util.exchange_list()
            }
            print(f"\n{org_model.id}\tExtracellular reactions", ex_rxns)
            variables = {
                ex_rxn.id: Variable(
                    "___".join([model_util.model.id, ex_rxn.id]),
                    lb=0,
                    ub=1,
                    type="binary",
                )
                for ex_rxn in ex_rxns
            }
            model_util.add_cons_vars(list(variables.values()))
            media, solutions = [], []
            sol = model_util.model.optimize()
            while sol.status == "optimal" and len(solutions) < n_solutions:
                solutions.append(sol)
                medium = set(
                    [
                        ex
                        for ex in ex_rxns
                        if sol.fluxes[ex.id] < -abstol and ex in other_excreta
                    ]
                )
                model_util.create_constraint(
                    Constraint(
                        sum([variables[ex.id] for ex in medium]),
                        ub=len(medium) - 1,
                        name=f"iteration_{len(solutions)}",
                    )
                )
                media.append(medium)
                sol = model_util.model.optimize()
            counter = Counter(chain(*media))
            scores[model_util.model.id] = {
                met.id: counter[ex] / len(media)
                for ex, met in ex_rxns.items()
                if counter[ex] > 0
            }
        return scores

    @staticmethod
    def sc(
        member_models: Iterable = None,
        com_model=None,
        min_growth=0.1,
        n_solutions=100,
        abstol=1e-6,
        compatibilized=True,
        printing=False,
    ):
        """Calculate the frequency of interspecies dependency in a community"""
        member_models, community = _load_models(
            member_models, com_model, not compatibilized, printing=printing
        )
        for rxn in com_model.reactions:
            rxn.lower_bound = 0 if "bio" in rxn.id else rxn.lower_bound

        # c_{rxn.id}_lb: rxn < 1000*y_{species_id}
        # c_{rxn.id}_ub: rxn > -1000*y_{species_id}
        variables = {}
        constraints = []
        # TODO this can be converted to an MSCommunity object by looping through each index
        # leverage CommKinetics
        for org_model in member_models:
            model_util = MSModelUtil(org_model, True)
            variables[model_util.model.id] = Variable(
                name=f"y_{model_util.model.id}", lb=0, ub=1, type="binary"
            )
            model_util.add_cons_vars([variables[model_util.model.id]])
            for rxn in model_util.model.reactions:
                if "bio" not in rxn.id:
                    # print(rxn.flux_expression)
                    lb = Constraint(
                        rxn.flux_expression + 1000 * variables[model_util.model.id],
                        name="_".join(["c", model_util.model.id, rxn.id, "lb"]),
                        lb=0,
                    )
                    ub = Constraint(
                        rxn.flux_expression - 1000 * variables[model_util.model.id],
                        name="_".join(["c", model_util.model.id, rxn.id, "ub"]),
                        ub=0,
                    )
                    constraints.extend([lb, ub])

        # calculate the SCS
        scores = {}
        for model in member_models:
            com_model_util = MSModelUtil(com_model)
            com_model_util.add_cons_vars(constraints, sloppy=True)
            # model growth is guaranteed while minimizing the growing members of the community
            ## SMETANA_Biomass: {biomass_reactions} > {min_growth}
            com_model_util.create_constraint(
                Constraint(
                    sum(
                        rxn.flux_expression
                        for rxn in model.reactions
                        if "bio" in rxn.id
                    ),
                    name="SMETANA_Biomass",
                    lb=min_growth,
                )
            )  # sloppy = True)
            other_members = [other for other in member_models if other.id != model.id]
            com_model_util.add_objective(
                sum([variables[other.id] for other in other_members]), "min"
            )
            previous_constraints, donors_list = [], []
            for i in range(n_solutions):
                sol = com_model.optimize()  # FIXME The solution is not optimal
                if sol.status != "optimal":
                    scores[model.id] = None
                    break
                donors = [
                    o
                    for o in other_members
                    if com_model.solver.primal_values[f"y_{o.id}"] > abstol
                ]
                donors_list.append(donors)
                previous_con = f"iteration_{i}"
                previous_constraints.append(previous_con)
                com_model_util.add_cons_vars(
                    [
                        Constraint(
                            sum(variables[o.id] for o in donors),
                            name=previous_con,
                            ub=len(previous_constraints) - 1,
                        )
                    ],
                    sloppy=True,
                )
            if i != 0:
                donors_counter = Counter(chain(*donors_list))
                scores[model.id] = {
                    o.id: donors_counter[o] / len(donors_list) for o in other_members
                }
        return scores

    @staticmethod
    def gyd(
        member_models: Iterable = None,
        model_utils: Iterable = None,
        environment=None,
        coculture_growth=False,
        community=None,
        anme_comm=False,
    ):
        gyds = {}
        for combination in combinations(model_utils or member_models, 2):
            if model_utils is None:
                model1_util = MSModelUtil(combination[0], True)
                model2_util = MSModelUtil(combination[1], True)
                print(
                    f"{model1_util.model.id} ++ {model2_util.model.id}",
                    model1_util.model.slim_optimize(),
                    model2_util.model.slim_optimize(),
                )
                if environment and not anme_comm:
                    model1_util.add_medium(environment)
                    model2_util.add_medium(environment)
            else:
                model1_util = combination[0]
                model2_util = combination[1]
            if not coculture_growth:
                G_m1, G_m2 = CommScores._determine_growths([model1_util, model2_util])
                G_m1, G_m2 = G_m1 if FBAHelper.isnumber(str(G_m1)) else 0, (
                    G_m2 if FBAHelper.isnumber(str(G_m2)) else 0
                )
            else:
                community = community or MSCommunity(
                    member_models=[model1_util.model, model2_util.model],
                    ids=[mem.id for mem in member_models],
                )
                community.run_fba()
                member_growths = community.parse_member_growths()
                G_m1, G_m2 = (
                    member_growths[model1_util.model.id],
                    member_growths[model2_util.model.id],
                )
            if G_m2 <= 0 or G_m1 <= 0:
                gyds[f"{model1_util.model.id} ++ {model2_util.model.id}"] = (
                    "",
                    "",
                    G_m1,
                    G_m2,
                )
                continue
            gyds[f"{model1_util.model.id} ++ {model2_util.model.id}"] = (
                abs(G_m1 - G_m2) / G_m1,
                abs(G_m2 - G_m1) / G_m2,
                G_m1,
                G_m2,
            )
        return gyds

    @staticmethod
    def pc(
        member_models=None,
        modelutils=None,
        com_model=None,
        isolate_growths=None,
        comm_sol=None,
        environment=None,
        comm_effects=True,
        community=None,
        interaction_threshold=0.1,
        compatibilized=False,
    ):
        assert member_models or modelutils or community, (
            "Members must be defined through either < member_models >"
            "or < modelutils > or < community >."
        )
        member_models = (
            member_models or [mem.model for mem in modelutils] or community.members
        )
        if com_model is None:
            member_models, com_model = _load_models(
                member_models, None, not compatibilized, printing=False
            )
        community = community or MSCommunity(com_model, member_models)
        if comm_sol is None:
            community.util.add_medium(environment)
            comm_sol = community.util.model.optimize()
        model_utils = modelutils or [MSModelUtil(mem, True) for mem in member_models]
        modelutils = []
        for mem in model_utils:
            mem.add_medium(environment)
            modelutils.append(mem)
        if isolate_growths is None:
            isolate_growths = {mem.id: mem.model.slim_optimize() for mem in modelutils}
        pc_score = comm_sol.objective_value / sum(list(isolate_growths.values()))
        if not comm_effects:
            return pc_score

        comm_member_growths = {
            mem.id: comm_sol.fluxes[mem.primary_biomass.id] for mem in community.members
        }
        comm_growth_effect = {
            memID: nanFilter(comm_environ / isolate_growths[memID])
            for memID, comm_environ in comm_member_growths.items()
        }
        growth_diffs = array(
            [nanFilter(x, False) for x in list(comm_growth_effect.values())]
        )
        th_pos, th_neg = 1 + interaction_threshold, 1 - interaction_threshold
        if all(growth_diffs > th_pos):
            bit = "mutualism"
        elif all(growth_diffs < th_neg):
            bit = "competitive"
        elif ((th_pos > growth_diffs) & (growth_diffs > th_neg)).all():
            bit = "neutral"
        elif all(growth_diffs > th_neg) and any(growth_diffs > th_pos):
            bit = "commensalism"
        elif all(growth_diffs < th_pos) and any(growth_diffs < th_neg):
            bit = "amensalism"
        elif any(growth_diffs > th_pos) and any(growth_diffs < th_neg):
            bit = "parasitism"
        else:
            print(
                f"The relative growths {comm_growth_effect} from {comm_member_growths} coculture and"
                f" {isolate_growths} monoculture are not captured."
            )
            bit = ""
        return (pc_score, comm_growth_effect, comm_member_growths, bit)

    @staticmethod
    def bss(
        member_models: Iterable = None,
        model_utils: Iterable = None,
        environments=None,
        minMedia=None,
        skip_bad_media=False,
    ):
        def compute_score(minMedia, environment=None, index=0):
            minMedia = minMedia or _get_media(
                model_s_=[modelUtil.model for modelUtil in model_utils],
                environment=environment,
                skip_bad_media=skip_bad_media,
            )
            model1_media = set(
                [
                    re.sub(r"(\_\w\d+$)", "", rxnID.replace("EX_", ""))
                    for rxnID in minMedia[model1_util.id]["media"].keys()
                ]
            )
            model2_media = set(
                [
                    re.sub(r"(\_\w\d+$)", "", rxnID.replace("EX_", ""))
                    for rxnID in minMedia[model2_util.id]["media"].keys()
                ]
            )
            model1_internal = {
                rm_comp(met.id)
                for rxn in model1_util.internal_list()
                for met in rxn.products
            }
            model2_internal = {
                rm_comp(met.id)
                for rxn in model2_util.internal_list()
                for met in rxn.products
            }
            bss_scores[
                f"{model1_util.id} supporting {model2_util.id} in media{index}"
            ] = (
                model1_internal,
                len(model2_media.intersection(model1_internal)) / len(model2_media),
            )
            bss_scores[
                f"{model2_util.id} supporting {model1_util.id} in media{index}"
            ] = (
                model2_internal,
                len(model1_media.intersection(model2_internal)) / len(model1_media),
            )

        bss_scores = {}
        for combination in combinations(model_utils or member_models, 2):
            if model_utils is None:
                model1_util = MSModelUtil(combination[0], True)
                model2_util = MSModelUtil(combination[1], True)
                model_utils = [model1_util, model2_util]
            else:
                model1_util = combination[0]
                model2_util = combination[1]
            if environments:
                for index, environment in enumerate(environments):
                    compute_score(minMedia, environment, index)
            else:
                compute_score(minMedia)
        return bss_scores

    @staticmethod
    def mqs():
        pass

    @staticmethod
    def _calculate_jaccard_score(set1, set2):
        if set1 == set2:
            print(f"The sets are identical, with a length of {len(set1)}.")
        if len(set1.union(set2)) == 0:
            return (None, None)
        return (
            set1.intersection(set2),
            len(set1.intersection(set2)) / len(set1.union(set2)),
        )

    @staticmethod
    def get_all_genomes_from_ws(
        ws_id,
        kbase_object=None,
        cobrakbase_repo_path: str = None,
        kbase_token_path: str = None,
    ):
        def get_genome(genome_name):
            return kbase_object.ws_client.get_objects2(
                {"objects": [{"ref": f"{ws_id}/{genome_name}"}]}
            )["data"][0]["data"]

        # load the kbase client instance
        if not kbase_object:
            import os

            os.environ["HOME"] = cobrakbase_repo_path
            import cobrakbase

            with open(kbase_token_path) as token_file:
                kbase_object = cobrakbase.KBaseAPI(token_file.readline())

        # calculate the complementarity
        genome_list = kbase_object.ws_client.list_objects(
            {
                "ids": [ws_id],
                "type": "KBaseGenomes.Genome",
                "minObjectID": 0,
                "maxObjectID": 10000,
            }
        )
        genome_names = [g[1] for g in genome_list if g[1].endswith("RAST")]
        return {
            genome_name: set(
                [
                    sso
                    for j in get_genome(genome_name)["cdss"]
                    for sso in j["ontology_terms"]["SSO"].keys()
                ]
            )
            for genome_name in genome_names
        }

    @staticmethod
    def fs(
        models: Iterable = None,
        kbase_object=None,
        cobrakbase_repo_path: str = None,
        kbase_token_path: str = None,
        annotated_genomes: dict = None,
        printing=False,
    ):
        if not isinstance(annotated_genomes, dict):
            if not kbase_object:
                import os

                os.environ["HOME"] = cobrakbase_repo_path
                import cobrakbase

                with open(kbase_token_path) as token_file:
                    kbase_object = cobrakbase.KBaseAPI(token_file.readline())
            annotated_genomes = {
                model.id: kbase_object.get_from_ws(model.genome_ref)
                for model in models
                if hasattr(model, "genome_ref")
            }
        elif isinstance(annotated_genomes, list):
            annotated_genomes = dict(
                zip([model.id for model in models], annotated_genomes)
            )
        elif models is not None:
            annotated_genomes = {
                k: v
                for k, v in annotated_genomes.items()
                if k in [model.id for model in models]
            }
        genome_combinations = list(combinations(annotated_genomes.keys(), 2))
        if printing:
            print(
                f"The Functionality Score (FS) will be calculated for {len(genome_combinations)} pairs."
            )
        if not isinstance(list(annotated_genomes.values())[0], dict):
            genome1_set, genome2_set = set(), set()
            distances = {}
            for genome1, genome2 in genome_combinations:
                for j in annotated_genomes[genome1].features:
                    for key, val in j.ontology_terms.items():
                        if key == "SSO":
                            genome1_set.update(val)
                for j in annotated_genomes[genome2].features:
                    for key, val in j.ontology_terms.items():
                        if key == "SSO":
                            genome2_set.update(val)
                distances[f"{genome1} ++ {genome2}"] = (
                    CommScores._calculate_jaccard_score(genome1_set, genome2_set)
                )
        else:
            distances = {
                f"{genome1} ++ {genome2}": CommScores._calculate_jaccard_score(
                    set(
                        list(content["SSO"].keys())[0]
                        for dic in annotated_genomes[genome1]["cdss"]
                        for x, content in dic.items()
                        if x == "ontology_terms" and len(content["SSO"].keys()) > 0
                    ),
                    set(
                        list(content["SSO"].keys())[0]
                        for dic in annotated_genomes[genome2]["cdss"]
                        for x, content in dic.items()
                        if x == "ontology_terms" and len(content["SSO"].keys()) > 0
                    ),
                )
                for genome1, genome2 in combinations(annotated_genomes.keys(), 2)
            }
        return distances

    @staticmethod
    def smetana(
        member_models: Iterable,
        environment,
        com_model=None,
        min_growth=0.1,
        n_solutions=100,
        abstol=1e-6,
        prior_values=None,
        compatibilized=False,
        sc_coupling=False,
        printing=False,
    ):
        """Quantifies the extent of syntrophy as the sum of all exchanges in a given nutritional environment"""
        member_models, community = _load_models(
            member_models, com_model, compatibilized == False, printing=printing
        )
        sc = None
        if not prior_values:
            mp = CommScores.mp(member_models, environment, com_model, abstol)
            mu = CommScores.mu(
                member_models, environment, mp, n_solutions, abstol, compatibilized
            )
            if sc_coupling:
                sc = CommScores.sc(
                    member_models,
                    com_model,
                    min_growth,
                    n_solutions,
                    abstol,
                    compatibilized,
                )
        elif len(prior_values) == 3:
            sc, mu, mp = prior_values
        else:
            mu, mp = prior_values

        smetana_scores = {}
        for pairs in combinations(member_models, 2):
            for model1, model2 in permutations(pairs):
                if model1.id not in smetana_scores:
                    smetana_scores[model1.id] = {}
                if not any([not mu[model1.id], not mp[model1.id]]):
                    sc_score = 1 if not sc_coupling else sc[model1.id][model2.id]
                    models_mets = list(model1.metabolites) + list(model2.metabolites)
                    unique_mets = set([met.id for met in models_mets])
                    smetana_scores[model1.id][model2.id] = 0
                    for met in models_mets:
                        if met.id in unique_mets:
                            mp_score = 0 if met.id not in mp[model1.id] else 1
                            smetana_scores[model1.id][model2.id] += (
                                mu[model1.id].get(met.id, 0) * sc_score * mp_score
                            )
        return smetana_scores

    @staticmethod
    def antiSMASH(json_path=None, zip_path=None):
        # TODO Scores 2, 4, and 5 are being explored for relevance to community formation and reveal specific member interactions/targets
        # load the antiSMASH report from either the JSON or the raw ZIP, or both
        from os import mkdir, listdir, path
        from zipfile import ZipFile
        from json import load

        if json_path:
            cwd_files = listdir()
            if json_path not in cwd_files and zip_path:
                with ZipFile(zip_path, "r") as zip_file:
                    zip_file.extract(json_path)
            with open(json_path, "r") as json_file:
                data = load(json_file)
        elif zip_path:
            mkdir("extracted_antiSMASH")
            with ZipFile(zip_path, "r") as zip_file:
                zip_file.extractall("extracted_antiSMASH")
            json_files = [
                x for x in listdir("extracted_antiSMASH") if x.endswith("json")
            ]
            if len(json_files) > 1:
                print(
                    f"The antiSMASH report describes {len(json_files)} JSON files, the first of which is selected "
                    f"{json_files[0]} for analysis, otherwise explicitly identify the desired JSON file in the json_path parameter."
                )
            with open(
                path.join("extracted_antiSMASH", json_files[0]), "r"
            ) as json_file:
                data = load(json_file)
        else:
            raise ParameterError(
                "Either the json_path or zip_path from the antiSMASH analysis must be provided,"
                " for these scores to be determined."
            )
        # Parse data and scores from the antiSMASH report
        biosynthetic_areas = data["records"][0]["areas"]
        BGCs = set(
            array(
                [
                    data["records"][0]["areas"][i]["products"]
                    for i in range(biosynthetic_areas)
                ]
            ).flatten()
        )
        len_proteins = len(
            data["records"][0]["modules"]["antismash.modules.clusterblast"][
                "knowncluster"
            ]["proteins"]
        )
        protein_annotations = [
            data["records"][0]["modules"]["antismash.modules.clusterblast"][
                "knowncluster"
            ]["proteins"][i]["annotations"]
            for i in range(len_proteins)
        ]
        clusterBlast = [s for s in protein_annotations if "resistance" in s]
        num_clusterBlast = sum(
            [item.count("resistance") for item in protein_annotations]
        )

        return (
            biosynthetic_areas,
            BGCs,
            protein_annotations,
            clusterBlast,
            num_clusterBlast,
        )
