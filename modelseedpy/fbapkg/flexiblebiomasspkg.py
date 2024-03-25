# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging

logger = logging.getLogger(__name__)
from optlang.symbolics import Zero, add  # !!! Neither import is ever used
from cobra import Model, Reaction, Metabolite  # !!! Model and Metabolite are never used
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper

classes = {
    "rna": {"cpd00052": -1, "cpd00038": -1, "cpd00002": -1, "cpd00062": -1},
    "dna": {"cpd00115": -1, "cpd00356": -1, "cpd00241": -1, "cpd00357": -1},
    "protein": {
        "cpd00132": -1,
        "cpd00023": -1,
        "cpd00053": -1,
        "cpd00054": -1,
        "cpd00033": -1,
        "cpd00039": -1,
        "cpd00119": -1,
        "cpd00051": -1,
        "cpd00041": -1,
        "cpd00107": -1,
        "cpd00129": -1,
        "cpd00322": -1,
        "cpd00069": -1,
        "cpd00065": -1,
        "cpd00060": -1,
        "cpd00084": -1,
        "cpd00035": -1,
        "cpd00161": -1,
        "cpd00156": -1,
        "cpd00066": -1,
    },
    "energy": {"cpd00008": 1},
}

# Base class for FBA packages
class FlexibleBiomassPkg(BaseFBAPkg):
    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "flexible biomass",
            {},
            {
                "flxbio": "reaction",
                "fflxcpd": "metabolite",
                "rflxcpd": "metabolite",
                "fflxcls": "reaction",
                "rflxcls": "reaction",
            },
        )

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            ["bio_rxn_id"],
            {
                "flex_coefficient": [-0.75, 0.75],
                "use_rna_class": [-0.75, 0.75],
                "use_dna_class": [-0.75, 0.75],
                "use_protein_class": [-0.75, 0.75],
                "use_energy_class": [-0.1, 0.1],
                "add_total_biomass_constraint": True,
            },
        )
        if self.parameters["bio_rxn_id"] not in self.model.reactions:
            raise ValueError(self.parameters["bio_rxn_id"] + " not found in model!")
        self.parameters["bio_rxn"] = self.model.reactions.get_by_id(
            self.parameters["bio_rxn_id"]
        )
        newrxns = []
        class_coef = {"rna": {}, "dna": {}, "protein": {}, "energy": {}}
        refcpd = {
            "cpd00001": None,
            "cpd00009": None,
            "cpd00012": None,
            "cpd00067": None,
            "cpd00002": None,
        }
        # Finding all reference compounds in the model
        msid_hash = self.modelutl.msid_hash()
        for msid in refcpd:
            if msid in msid_hash:
                refcpd[msid] = msid_hash[msid][0]
        met_class = {}
        # Determining class for each metabolite in biomass reaction
        for metabolite in self.parameters["bio_rxn"].metabolites:
            met_class[metabolite] = None
            msid = MSModelUtil.metabolite_msid(metabolite)
            if (
                msid != "cpd11416"
                and msid != "cpd11463"
                and msid != "cpd11462"
                and msid != "cpd11461"
                and msid != None
            ):
                if msid in refcpd:
                    met_class[metabolite] = "refcpd"
                else:
                    for curr_class in classes:
                        if (
                            self.parameters["use_" + curr_class + "_class"]
                            and msid in classes[curr_class]
                        ):
                            met_class[metabolite] = curr_class
                            class_coef[curr_class][msid] = metabolite
        # Eliminating any classes that are incomplete
        for curr_class in classes:
            for msid in classes[curr_class]:
                if msid not in class_coef[curr_class]:
                    self.parameters["use_" + curr_class + "_class"] = None
                    break
        # Creating FLEX reactions and constraints for unclassified compounds
        flexcpds = {}
        for metabolite in self.parameters["bio_rxn"].metabolites:
            if not met_class[metabolite]:
                flexcpds[metabolite] = self.parameters["bio_rxn"].metabolites[
                    metabolite
                ]
            elif (
                met_class[metabolite] != "refcpd"
                and not self.parameters["use_" + met_class[metabolite] + "_class"]
            ):
                flexcpds[metabolite] = self.parameters["bio_rxn"].metabolites[
                    metabolite
                ]
        self.modelutl.add_exchanges_for_metabolites(
            flexcpds,
            uptake=1000,
            excretion=1000,
            prefix="FLEX_" + self.parameters["bio_rxn"].id + "_",
            prefix_name="Biomass flex for ",
        )
        for metabolite in flexcpds:
            self.build_constraint(metabolite, "flxcpd")
        # Creating metabolite class constraints
        for met_class in classes:
            if self.parameters["use_" + met_class + "_class"]:
                add = 0
                total_coef = 0
                object_stoichiometry = {}
                for msid in class_coef[met_class]:
                    if (
                        met_class == "rna"
                        and msid == "cpd00002"
                        and "cpd00008" in class_coef["energy"]
                    ):
                        object_stoichiometry[class_coef[met_class][msid]] = (
                            self.parameters["bio_rxn"].metabolites[
                                class_coef[met_class][msid]
                            ]
                            + self.parameters["bio_rxn"].metabolites[
                                class_coef["energy"]["cpd00008"]
                            ]
                        )
                    else:
                        object_stoichiometry[
                            class_coef[met_class][msid]
                        ] = self.parameters["bio_rxn"].metabolites[
                            class_coef[met_class][msid]
                        ]
                    total_coef += abs(object_stoichiometry[class_coef[met_class][msid]])
                if (
                    (met_class == "rna" or met_class == "dna")
                    and refcpd["cpd00012"] != None
                    and refcpd["cpd00001"] != None
                ):
                    add = 1
                    object_stoichiometry[refcpd["cpd00012"]] = total_coef
                    object_stoichiometry[refcpd["cpd00001"]] = total_coef
                if met_class == "protein" and refcpd["cpd00001"] != None:
                    add = 1
                    object_stoichiometry[refcpd["cpd00001"]] = total_coef
                if (
                    met_class == "energy"
                    and refcpd["cpd00001"] != None
                    and refcpd["cpd00002"] != None
                    and refcpd["cpd00067"] != None
                    and refcpd["cpd00009"] != None
                ):
                    add = 1
                    object_stoichiometry[refcpd["cpd00001"]] = -1 * total_coef
                    object_stoichiometry[refcpd["cpd00002"]] = -1 * total_coef
                    object_stoichiometry[refcpd["cpd00009"]] = total_coef
                    object_stoichiometry[refcpd["cpd00067"]] = total_coef
                if add == 1:
                    if met_class + "_flex" not in self.new_reactions:
                        self.new_reactions[met_class + "_flex"] = Reaction(
                            id=met_class + "_flex",
                            name=met_class + "_flex",
                            lower_bound=-1000,
                            upper_bound=1000,
                        )
                        self.new_reactions[met_class + "_flex"].add_metabolites(
                            object_stoichiometry
                        )
                        self.new_reactions[met_class + "_flex"].annotation[
                            "sbo"
                        ] = "SBO:0000627"
                        self.model.add_reactions(
                            [self.new_reactions[met_class + "_flex"]]
                        )
                    self.build_constraint(
                        self.new_reactions[met_class + "_flex"], "flxcls"
                    )
        if parameters["add_total_biomass_constraint"]:
            self.build_constraint(self.parameters["bio_rxn"], "flxbio")

    def build_variable(self, object, type):  # !!! can the function be removed?
        pass

    def build_constraint(self, cobra_obj, obj_type):
        if obj_type == "flxbio":
            # Sum(MW*(vdrn,for-vdrn,ref)) + Sum(massdiff*(vrxn,for-vrxn,ref)) = 0
            coef = {}
            for metabolite in self.parameters["bio_rxn"].metabolites:
                if (
                    "FLEX_" + self.parameters["bio_rxn"].id + "_" + metabolite.id
                    in self.model.reactions
                ):
                    mw = FBAHelper.metabolite_mw(metabolite)
                    sign = -1
                    if self.parameters["bio_rxn"].metabolites[metabolite] > 0:
                        sign = 1
                    coef[
                        self.model.reactions.get_by_id(
                            "FLEX_"
                            + self.parameters["bio_rxn"].id
                            + "_"
                            + metabolite.id
                        ).forward_variable
                    ] = (sign * mw)
                    coef[
                        self.model.reactions.get_by_id(
                            "FLEX_"
                            + self.parameters["bio_rxn"].id
                            + "_"
                            + metabolite.id
                        ).reverse_variable
                    ] = (-1 * sign * mw)
            for met_class in classes:
                if met_class + "_flex" in self.model.reactions:
                    massdiff = 0
                    rxn = self.model.reactions.get_by_id(met_class + "_flex")
                    for met in rxn.metabolites:
                        mw = FBAHelper.metabolite_mw(met)
                        massdiff += rxn.metabolites[met] * mw
                    if abs(massdiff) > 0.00001:
                        coef[rxn.forward_variable] = massdiff
                        coef[rxn.reverse_variable] = -massdiff
            return BaseFBAPkg.build_constraint(self, obj_type, 0, 0, coef, cobra_obj)
        elif obj_type == "flxcpd" or obj_type == "flxcls":
            first_entry = None
            second_entry = None
            product = False
            biovar = self.parameters["bio_rxn"].forward_variable
            object = None
            const = None
            if obj_type == "flxcpd":
                # 0.75 * abs(bio_coef) * vbio - vdrn,for >= 0
                # 0.75 * abs(bio_coef) * vbio - vdrn,rev >= 0
                first_entry = self.parameters["flex_coefficient"][0] * abs(
                    self.parameters["bio_rxn"].metabolites[cobra_obj]
                )
                second_entry = self.parameters["flex_coefficient"][1] * abs(
                    self.parameters["bio_rxn"].metabolites[cobra_obj]
                )
                if self.parameters["bio_rxn"].metabolites[cobra_obj] > 0:
                    product = True
                object = self.model.reactions.get_by_id(
                    "FLEX_" + self.parameters["bio_rxn"].id + "_" + cobra_obj.id
                )
            elif (
                cobra_obj.id[0:-5] == None
                or not self.parameters["use_" + cobra_obj.id[0:-5] + "_class"]
            ):
                return None
            else:
                # 0.75 * vbio - vrxn,for >= 0
                # 0.75 * vbio - vrxn,rev >= 0
                first_entry = self.parameters["use_" + cobra_obj.id[0:-5] + "_class"][0]
                second_entry = self.parameters["use_" + cobra_obj.id[0:-5] + "_class"][
                    1
                ]
                object = cobra_obj
            if first_entry == second_entry:
                # If the value is positive, lock in the forward variable and set the reverse to zero
                if first_entry > 0:
                    if product:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            0,
                            {biovar: second_entry, object.forward_variable: -1},
                            cobra_obj,
                        )
                        object.lower_bound = 0
                    else:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            0,
                            {biovar: second_entry, object.reverse_variable: -1},
                            cobra_obj,
                        )
                        object.upper_bound = 0
                # If the value is negative, lock in the reverse variable and set the forward to zero
                elif first_entry < 0:
                    if product:
                        const = self.build_constraint(
                            "r" + obj_type,
                            0,
                            0,
                            {biovar: -first_entry, object.reverse_variable: -1},
                            cobra_obj,
                        )
                        object.upper_bound = 0
                    else:
                        const = self.build_constraint(
                            "r" + obj_type,
                            0,
                            0,
                            {biovar: -first_entry, object.forward_variable: -1},
                            cobra_obj,
                        )
                        object.lower_bound = 0
                # If the value is zero, lock both variables to zero
                if first_entry == 0:
                    object.lower_bound = 0
                    object.upper_bound = 0
            elif second_entry >= 0:
                if first_entry >= 0:
                    if product:
                        const = BaseFBAPkg.build_constraint(
                            self,
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.forward_variable: -1},
                            cobra_obj,
                        )
                        object.lower_bound = 0
                        if first_entry > 0:
                            BaseFBAPkg.build_constraint(
                                self,
                                "r" + obj_type,
                                0,
                                None,
                                {biovar: -first_entry, object.forward_variable: 1},
                                cobra_obj,
                            )
                    else:
                        const = BaseFBAPkg.build_constraint(
                            self,
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.reverse_variable: -1},
                            cobra_obj,
                        )
                        object.upper_bound = 0
                        if first_entry > 0:
                            BaseFBAPkg.build_constraint(
                                self,
                                "r" + obj_type,
                                0,
                                None,
                                {biovar: -first_entry, object.reverse_variable: 1},
                                cobra_obj,
                            )
                else:
                    if product:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.forward_variable: -1},
                            cobra_obj,
                        )
                        self.build_constraint(
                            "r" + obj_type,
                            0,
                            None,
                            {biovar: -first_entry, object.reverse_variable: -1},
                            cobra_obj,
                        )
                    else:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.reverse_variable: -1},
                            cobra_obj,
                        )
                        self.build_constraint(
                            "r" + obj_type,
                            0,
                            None,
                            {biovar: -first_entry, object.forward_variable: -1},
                            cobra_obj,
                        )
            else:
                if second_entry < 0:
                    if product:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.reverse_variable: 1},
                            cobra_obj,
                        )
                    else:
                        const = self.build_constraint(
                            "f" + obj_type,
                            0,
                            None,
                            {biovar: second_entry, object.forward_variable: 1},
                            cobra_obj,
                        )
                if product:
                    self.build_constraint(
                        "r" + obj_type,
                        0,
                        None,
                        {biovar: -first_entry, object.reverse_variable: -1},
                        cobra_obj,
                    )
                    object.lower_bound = 0
                else:
                    self.build_constraint(
                        "r" + obj_type,
                        0,
                        None,
                        {biovar: -first_entry, object.forward_variable: -1},
                        cobra_obj,
                    )
                    object.upper_bound = 0
            return const
