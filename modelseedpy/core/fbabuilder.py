import logging

import re
import copy
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from cobrakbase.core.kbaseobject import AttrDict
from cobrakbase.annotation_ontology_api.annotation_ontology_apiServiceClient import (
    annotation_ontology_api,
)
import modelseedpy.core.fbahelper

logger = logging.getLogger(__name__)


def build_cpd_id(str):
    if str.startswith("M_"):
        str = str[2:]
    elif str.startswith("M-"):
        str = str[2:]
    str_fix = str
    if "-" in str_fix:
        str_fix = str_fix.replace("-", "__DASH__")
    if not str == str_fix:
        logger.debug("[Species] rename: [%s] -> [%s]", str, str_fix)
    return str


def build_rxn_id(str):
    if str.startswith("R_"):
        str = str[2:]
    elif str.startswith("R-"):
        str = str[2:]
    str_fix = str
    if "-" in str_fix:
        str_fix = str_fix.replace("-", "__DASH__")
    if not str == str_fix:
        logger.debug("[Reaction] rename: [%s] -> [%s]", str, str_fix)
    return str_fix


# Adding a few exception classes to handle different types of errors
class ObjectError(Exception):
    """Error in the construction of a base KBase object"""

    pass


class FeasibilityError(Exception):
    """Error in FBA formulation"""

    pass


# New class to store functions to building and tracking new constraints and variables related to our own custom FBA formulations
class KBaseFBAUtilities:
    def __init__(
        self,
        cobramodel,
        fbamodel,
        kbapi,
        media=None,
        default_uptake=100,
        default_excretion=100,
        blacklist=[],
        auto_sink=["cpd02701_c", "cpd11416_c0", "cpd15302_c"],
    ):
        self.cobramodel = cobramodel
        self.SBO_ANNOTATION = "sbo"
        self.metabolites_remap = {}
        self.solution_exclusion_constraints = []
        self.kbapi = kbapi
        self.potential_variables = dict()
        self.reversibility_binary = dict()
        self.reversibility_binary_constraints = dict()
        self.binary_flux_variables = dict()
        self.total_flux_variables = dict()
        self.total_flux_constraints = dict()
        self.binary_flux_constraints = dict()
        self.simple_thermo_constraints = dict()
        self.metabolomics_peak_variables = dict()
        self.metabolomics_peak_constraints = dict()
        self.compound_flux_variables = dict()
        self.compound_flux_constraints = dict()
        self.metabolomics_constraints = dict()
        self.media = None
        self.default_uptake = default_uptake
        self.default_excretion = default_excretion
        self.apply_media_to_model(media, self.default_uptake, self.default_excretion)
        self.blacklist = [
            "rxn12985",
            "rxn00238",
            "rxn07058",
            "rxn05305",
            "rxn00154",
            "rxn09037",
            "rxn10643",
            "rxn11317",
            "rxn05254",
            "rxn05257",
            "rxn05258",
            "rxn05259",
            "rxn05264",
            "rxn05268",
            "rxn05269",
            "rxn05270",
            "rxn05271",
            "rxn05272",
            "rxn05273",
            "rxn05274",
            "rxn05275",
            "rxn05276",
            "rxn05277",
            "rxn05278",
            "rxn05279",
            "rxn05280",
            "rxn05281",
            "rxn05282",
            "rxn05283",
            "rxn05284",
            "rxn05285",
            "rxn05286",
            "rxn05963",
            "rxn05964",
            "rxn05971",
            "rxn05989",
            "rxn05990",
            "rxn06041",
            "rxn06042",
            "rxn06043",
            "rxn06044",
            "rxn06045",
            "rxn06046",
            "rxn06079",
            "rxn06080",
            "rxn06081",
            "rxn06086",
            "rxn06087",
            "rxn06088",
            "rxn06089",
            "rxn06090",
            "rxn06091",
            "rxn06092",
            "rxn06138",
            "rxn06139",
            "rxn06140",
            "rxn06141",
            "rxn06145",
            "rxn06217",
            "rxn06218",
            "rxn06219",
            "rxn06220",
            "rxn06221",
            "rxn06222",
            "rxn06223",
            "rxn06235",
            "rxn06362",
            "rxn06368",
            "rxn06378",
            "rxn06474",
            "rxn06475",
            "rxn06502",
            "rxn06562",
            "rxn06569",
            "rxn06604",
            "rxn06702",
            "rxn06706",
            "rxn06715",
            "rxn06803",
            "rxn06811",
            "rxn06812",
            "rxn06850",
            "rxn06901",
            "rxn06971",
            "rxn06999",
            "rxn07123",
            "rxn07172",
            "rxn07254",
            "rxn07255",
            "rxn07269",
            "rxn07451",
            "rxn09037",
            "rxn10018",
            "rxn10077",
            "rxn10096",
            "rxn10097",
            "rxn10098",
            "rxn10099",
            "rxn10101",
            "rxn10102",
            "rxn10103",
            "rxn10104",
            "rxn10105",
            "rxn10106",
            "rxn10107",
            "rxn10109",
            "rxn10111",
            "rxn10403",
            "rxn10410",
            "rxn10416",
            "rxn11313",
            "rxn11316",
            "rxn11318",
            "rxn11353",
            "rxn05224",
            "rxn05795",
            "rxn05796",
            "rxn05797",
            "rxn05798",
            "rxn05799",
            "rxn05801",
            "rxn05802",
            "rxn05803",
            "rxn05804",
            "rxn05805",
            "rxn05806",
            "rxn05808",
            "rxn05812",
            "rxn05815",
            "rxn05832",
            "rxn05836",
            "rxn05851",
            "rxn05857",
            "rxn05869",
            "rxn05870",
            "rxn05884",
            "rxn05888",
            "rxn05896",
            "rxn05898",
            "rxn05900",
            "rxn05903",
            "rxn05904",
            "rxn05905",
            "rxn05911",
            "rxn05921",
            "rxn05925",
            "rxn05936",
            "rxn05947",
            "rxn05956",
            "rxn05959",
            "rxn05960",
            "rxn05980",
            "rxn05991",
            "rxn05992",
            "rxn05999",
            "rxn06001",
            "rxn06014",
            "rxn06017",
            "rxn06021",
            "rxn06026",
            "rxn06027",
            "rxn06034",
            "rxn06048",
            "rxn06052",
            "rxn06053",
            "rxn06054",
            "rxn06057",
            "rxn06059",
            "rxn06061",
            "rxn06102",
            "rxn06103",
            "rxn06127",
            "rxn06128",
            "rxn06129",
            "rxn06130",
            "rxn06131",
            "rxn06132",
            "rxn06137",
            "rxn06146",
            "rxn06161",
            "rxn06167",
            "rxn06172",
            "rxn06174",
            "rxn06175",
            "rxn06187",
            "rxn06189",
            "rxn06203",
            "rxn06204",
            "rxn06246",
            "rxn06261",
            "rxn06265",
            "rxn06266",
            "rxn06286",
            "rxn06291",
            "rxn06294",
            "rxn06310",
            "rxn06320",
            "rxn06327",
            "rxn06334",
            "rxn06337",
            "rxn06339",
            "rxn06342",
            "rxn06343",
            "rxn06350",
            "rxn06352",
            "rxn06358",
            "rxn06361",
            "rxn06369",
            "rxn06380",
            "rxn06395",
            "rxn06415",
            "rxn06419",
            "rxn06420",
            "rxn06421",
            "rxn06423",
            "rxn06450",
            "rxn06457",
            "rxn06463",
            "rxn06464",
            "rxn06466",
            "rxn06471",
            "rxn06482",
            "rxn06483",
            "rxn06486",
            "rxn06492",
            "rxn06497",
            "rxn06498",
            "rxn06501",
            "rxn06505",
            "rxn06506",
            "rxn06521",
            "rxn06534",
            "rxn06580",
            "rxn06585",
            "rxn06593",
            "rxn06609",
            "rxn06613",
            "rxn06654",
            "rxn06667",
            "rxn06676",
            "rxn06693",
            "rxn06730",
            "rxn06746",
            "rxn06762",
            "rxn06779",
            "rxn06790",
            "rxn06791",
            "rxn06792",
            "rxn06793",
            "rxn06794",
            "rxn06795",
            "rxn06796",
            "rxn06797",
            "rxn06821",
            "rxn06826",
            "rxn06827",
            "rxn06829",
            "rxn06839",
            "rxn06841",
            "rxn06842",
            "rxn06851",
            "rxn06866",
            "rxn06867",
            "rxn06873",
            "rxn06885",
            "rxn06891",
            "rxn06892",
            "rxn06896",
            "rxn06938",
            "rxn06939",
            "rxn06944",
            "rxn06951",
            "rxn06952",
            "rxn06955",
            "rxn06957",
            "rxn06960",
            "rxn06964",
            "rxn06965",
            "rxn07086",
            "rxn07097",
            "rxn07103",
            "rxn07104",
            "rxn07105",
            "rxn07106",
            "rxn07107",
            "rxn07109",
            "rxn07119",
            "rxn07179",
            "rxn07186",
            "rxn07187",
            "rxn07188",
            "rxn07195",
            "rxn07196",
            "rxn07197",
            "rxn07198",
            "rxn07201",
            "rxn07205",
            "rxn07206",
            "rxn07210",
            "rxn07244",
            "rxn07245",
            "rxn07253",
            "rxn07275",
            "rxn07299",
            "rxn07302",
            "rxn07651",
            "rxn07723",
            "rxn07736",
            "rxn07878",
            "rxn11417",
            "rxn11582",
            "rxn11593",
            "rxn11597",
            "rxn11615",
            "rxn11617",
            "rxn11619",
            "rxn11620",
            "rxn11624",
            "rxn11626",
            "rxn11638",
            "rxn11648",
            "rxn11651",
            "rxn11665",
            "rxn11666",
            "rxn11667",
            "rxn11698",
            "rxn11983",
            "rxn11986",
            "rxn11994",
            "rxn12006",
            "rxn12007",
            "rxn12014",
            "rxn12017",
            "rxn12022",
            "rxn12160",
            "rxn12161",
            "rxn01267",
            "rxn05294",
            "rxn04656",
        ]
        for item in blacklist:
            if item not in self.blacklist:
                self.blacklist.append(item)
        self.auto_sink = []
        full_id = re.compile("\d+$")
        for id in auto_sink:
            if full_id.search(id):
                self.auto_sink.append(id)
            else:
                for i in range(0, 100):
                    newid = id + str(i)
                    self.auto_sink.append(newid)

        self.auto_exchange = "e0"
        self.sink_compounds = set()
        self.demand_compounds = set()
        self.exchange_compounds = set()
        self.COBRA_0_BOUND = 0
        self.COBRA_DEFAULT_LB = -1000
        self.COBRA_DEFAULT_UB = 1000

    def media_const_hash(self):
        bound_hash = dict()
        if not self.media == None:
            for compound in self.media.mediacompounds:
                bound_hash[compound.id] = {
                    "lb": -1 * compound.maxFlux,
                    "ub": -1 * compound.minFlux,
                }
        return bound_hash

    def apply_media_to_model(
        self, media=None, default_uptake=None, default_excretion=None
    ):
        self.media = media
        if default_uptake == None:
            default_uptake = self.default_uptake
        if default_excretion == None:
            default_excretion = self.default_excretion

        bound_hash = self.media_const_hash()
        for reaction in self.cobramodel.reactions:
            if reaction.id[0:3].lower() == "ex_":
                compound = reaction.id[3:]
                if compound[-3:] == "_e0":
                    compound = compound[:-3]
                if compound in bound_hash:
                    reaction.lower_bound = bound_hash[compound]["lb"]
                    reaction.upper_bound = bound_hash[compound]["ub"]
                else:
                    reaction.lower_bound = -1 * default_uptake
                    reaction.upper_bound = default_excretion
                reaction.update_variable_bounds()

    def add_total_flux_constraints(self, reaction_filter=None):
        for reaction in self.cobramodel.reactions:
            if reaction_filter == None or reaction.id in reaction_filter:
                self.total_flux_variables[
                    reaction.id
                ] = self.cobramodel.problem.Variable(
                    reaction.id + "_tot", lb=0, ub=self.COBRA_DEFAULT_UB
                )
                self.cobramodel.add_cons_vars(self.total_flux_variables[reaction.id])
                self.total_flux_constraints[
                    reaction.id
                ] = self.cobramodel.problem.Constraint(
                    reaction.forward_variable
                    + reaction.reverse_variable
                    - self.total_flux_variables[reaction.id],
                    lb=0,
                    ub=0,
                    name=reaction.id + "_tot",
                )
                self.cobramodel.add_cons_vars(self.total_flux_constraints[reaction.id])

    def add_reversibility_binary_constraints(self, reaction_filter=None):
        # Adding thermodynamic constraints
        for reaction in self.cobramodel.reactions:
            if reaction.id not in self.reversibility_binary and (
                reaction_filter == None or reaction.id in reaction_filter
            ):
                self.reversibility_binary[
                    reaction.id
                ] = self.cobramodel.problem.Variable(
                    reaction.id + "_rb", lb=0, ub=1, type="binary"
                )
                self.cobramodel.add_cons_vars(self.reversibility_binary[reaction.id])
                self.reversibility_binary_constraints[reaction.id] = dict()
                self.reversibility_binary_constraints[reaction.id][
                    "ff"
                ] = self.cobramodel.problem.Constraint(
                    1000 * self.reversibility_binary[reaction.id]
                    - reaction.forward_variable,
                    lb=0,
                    ub=None,
                    name=reaction.id + "_FB",
                )
                self.cobramodel.add_cons_vars(
                    self.reversibility_binary_constraints[reaction.id]["ff"]
                )
                self.reversibility_binary_constraints[reaction.id][
                    "rf"
                ] = self.cobramodel.problem.Constraint(
                    -1000 * self.reversibility_binary[reaction.id]
                    - reaction.reverse_variable,
                    lb=-1000,
                    ub=None,
                    name=reaction.id + "_RB",
                )
                self.cobramodel.add_cons_vars(
                    self.reversibility_binary_constraints[reaction.id]["rf"]
                )

    def set_objective_from_target_reaction(self, target_reaction, maximize=1):
        target_reaction = self.cobramodel.reactions.get_by_id(target_reaction)
        sense = "max"
        if maximize == 0:
            sense = "min"
        target_objective = self.cobramodel.problem.Objective(
            1 * target_reaction.flux_expression, direction=sense
        )
        self.cobramodel.objective = target_objective
        return target_reaction

    def add_simple_thermo_constraints(self):
        # Creating potential variables for all compounds
        for metabolite in self.cobramodel.metabolites:
            if metabolite.id not in self.potential_variables:
                self.potential_variables[
                    metabolite.id
                ] = self.cobramodel.problem.Variable(
                    metabolite.id + "_u", lb=0, ub=1000
                )
                self.cobramodel.add_cons_vars(self.potential_variables[metabolite.id])
        # Adding thermodynamic constraints
        for reaction in self.cobramodel.reactions:
            if (
                reaction.id not in self.simple_thermo_constraints
                and reaction.id[0:3].lower() != "ex_"
                and reaction.id[0:3].lower() != "dm_"
            ):
                if reaction.id not in self.reversibility_binary:
                    self.reversibility_binary[
                        reaction.id
                    ] = self.cobramodel.problem.Variable(
                        reaction.id + "_rb", lb=0, ub=1, type="binary"
                    )
                    self.cobramodel.add_cons_vars(
                        self.reversibility_binary[reaction.id]
                    )
                    self.reversibility_binary_constraints[reaction.id] = dict()
                    self.reversibility_binary_constraints[reaction.id][
                        "ff"
                    ] = self.cobramodel.problem.Constraint(
                        1000 * self.reversibility_binary[reaction.id]
                        - reaction.forward_variable,
                        lb=0,
                        ub=None,
                        name=reaction.id + "_FB",
                    )
                    self.cobramodel.add_cons_vars(
                        self.reversibility_binary_constraints[reaction.id]["ff"]
                    )
                    self.reversibility_binary_constraints[reaction.id][
                        "rf"
                    ] = self.cobramodel.problem.Constraint(
                        -1000 * self.reversibility_binary[reaction.id]
                        - reaction.reverse_variable,
                        lb=-1000,
                        ub=None,
                        name=reaction.id + "_RB",
                    )
                    self.cobramodel.add_cons_vars(
                        self.reversibility_binary_constraints[reaction.id]["rf"]
                    )
                self.simple_thermo_constraints[
                    reaction.id
                ] = self.cobramodel.problem.Constraint(
                    Zero, lb=0, ub=1000, name=reaction.id + "_therm"
                )
                self.cobramodel.add_cons_vars(
                    self.simple_thermo_constraints[reaction.id]
                )
                self.cobramodel.solver.update()
                const_coef = {self.reversibility_binary[reaction.id]: 1000}
                for metabolite in reaction.metabolites:
                    const_coef[
                        self.potential_variables[metabolite.id]
                    ] = reaction.metabolites[metabolite]
                self.simple_thermo_constraints[reaction.id].set_linear_coefficients(
                    const_coef
                )
        # Updating solver one final time
        self.cobramodel.solver.update()

    def add_intracellular_metabolomics_constraints(
        self, peakstring, relevant_peaks=None
    ):
        drain_fluxes = list()
        peak_array = peakstring.split(";")
        compound_reactions = dict()
        reaction_hash = dict()
        for reaction in self.cobramodel.reactions:
            reaction_hash[reaction.id] = 1
            for compound in reaction.metabolites:
                if compound.id not in compound_reactions:
                    compound_reactions[compound.id] = dict()
                compound_reactions[compound.id][reaction.id] = reaction.metabolites[
                    compound
                ]
        compartment_tag = re.compile("_[a-z]\d+$")
        for peak in peak_array:
            sub_array = peak.split(":")
            if len(sub_array) > 2:
                peakid = sub_array[0]
                if relevant_peaks == None or peakid in relevant_peaks:
                    coef = sub_array[1]
                    peak_coef = dict()
                    pfound = 0
                    for i in range(2, len(sub_array)):
                        compound_list = []
                        compound = sub_array[i]
                        if compartment_tag.search(compound):
                            compound_list = [compound]
                        else:
                            for i in range(0, 1000):
                                compound_list.append(compound + "_c" + str(i))
                        for compound in compound_list:
                            if compound in compound_reactions:
                                cfound = 0
                                compound_coef = dict()
                                for reaction in compound_reactions[compound]:
                                    if (
                                        reaction[0:3].lower() != "ex_"
                                        and reaction[0:3].lower() != "dm_"
                                    ):
                                        cfound = 1
                                        rxnobj = self.cobramodel.reactions.get_by_id(
                                            reaction
                                        )
                                        compound_coef[rxnobj.forward_variable] = 1000
                                        compound_coef[rxnobj.reverse_variable] = 1000
                                if cfound == 1:
                                    if compound not in self.compound_flux_variables:
                                        self.compound_flux_variables[
                                            compound
                                        ] = self.cobramodel.problem.Variable(
                                            compound + "_f", lb=0, ub=1
                                        )
                                        self.cobramodel.add_cons_vars(
                                            self.compound_flux_variables[compound]
                                        )
                                        self.compound_flux_constraints[
                                            compound
                                        ] = self.cobramodel.problem.Constraint(
                                            Zero, lb=0, ub=None, name=compound + "_flux"
                                        )
                                        self.cobramodel.add_cons_vars(
                                            self.compound_flux_constraints[compound]
                                        )
                                    compound_coef[
                                        self.compound_flux_variables[compound]
                                    ] = -1
                                    self.cobramodel.solver.update()
                                    self.compound_flux_constraints[
                                        compound
                                    ].set_linear_coefficients(compound_coef)
                                    peak_coef[
                                        self.compound_flux_variables[compound]
                                    ] = 1
                                    pfound = 1
                                    drain_reaction = (
                                        self.helper.add_drain_from_metabolite_id(
                                            self.cobramodel, compound
                                        )
                                    )
                                    if (
                                        drain_reaction.id
                                        not in self.cobramodel.reactions
                                    ):
                                        self.cobramodel.add_reactions([drain_reaction])
                    if pfound == 1:
                        if peakid not in self.metabolomics_peak_variables:
                            self.metabolomics_peak_variables[
                                peakid
                            ] = self.cobramodel.problem.Variable(peakid, lb=0, ub=1)
                            self.cobramodel.add_cons_vars(
                                self.metabolomics_peak_variables[peakid]
                            )
                            self.metabolomics_peak_constraints[
                                peakid
                            ] = self.cobramodel.problem.Constraint(
                                Zero, lb=0, ub=None, name=peakid
                            )
                            self.cobramodel.add_cons_vars(
                                self.metabolomics_peak_constraints[peakid]
                            )
                        peak_coef[self.metabolomics_peak_variables[peakid]] = -1
                        self.cobramodel.solver.update()
                        self.metabolomics_peak_constraints[
                            peakid
                        ].set_linear_coefficients(peak_coef)

        return drain_fluxes

    def convert_template_compound(self, template_compound, index, template):
        base_id = template_compound.id.split("_")[0]
        base_compound = template.compounds.get_by_id(base_id)
        new_id = template_compound.id
        new_id += str(index)
        compartment = template_compound.templatecompartment_ref.split("/").pop()
        compartment += str(index)

        met = Metabolite(
            new_id,
            formula=base_compound.formula,
            name=base_compound.name,
            charge=template_compound.charge,
            compartment=compartment,
        )

        met.annotation[
            "sbo"
        ] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        met.annotation["seed.compound"] = base_id
        return met

    def convert_template_reaction(
        self, template_reaction, index, template, for_gapfilling=1
    ):
        array = template_reaction.id.split("_")
        base_id = array[0]
        new_id = template_reaction.id
        new_id += str(index)

        lower_bound = template_reaction.maxrevflux
        upper_bound = template_reaction.maxforflux

        direction = template_reaction.GapfillDirection
        if for_gapfilling == 0:
            direction = template_reaction.direction

        if direction == ">":
            lower_bound = 0
        elif direction == "<":
            upper_bound = 0

        cobra_reaction = Reaction(
            new_id,
            name=template_reaction.name,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )

        object_stoichiometry = {}
        for item in template_reaction.templateReactionReagents:
            metabolite_id = item["templatecompcompound_ref"].split("/").pop()
            template_compound = template.compcompounds.get_by_id(metabolite_id)
            compartment = template_compound.templatecompartment_ref.split("/").pop()
            if compartment == "e":
                metabolite_id = metabolite_id + "0"
            else:
                metabolite_id = metabolite_id + str(index)

            metabolite = self.cobramodel.metabolites.get_by_id(metabolite_id)
            object_stoichiometry[metabolite] = item["coefficient"]

        cobra_reaction.add_metabolites(object_stoichiometry)

        cobra_reaction.annotation["sbo"] = "SBO:0000176"  # biochemical reaction
        cobra_reaction.annotation["seed.reaction"] = base_id

        return cobra_reaction

    def build_model_extended_for_gapfilling(
        self,
        extend_with_template=1,
        source_models=[],
        input_templates=[],
        model_penalty=1,
        reaction_scores={},
    ):
        model_id = self.fbamodel["id"] + ".gf"

        # Determine all indecies that should be gapfilled
        indexlist = [0] * 1000
        compounds = self.fbamodel["modelcompounds"]
        for compound in compounds:
            compartment = compound["modelcompartment_ref"].split("/").pop()
            basecomp = compartment[0:1]
            if not basecomp == "e":
                index = compartment[1:]
                index = int(index)
                indexlist[index] += 1

        # Iterating over all indecies with more than 10 intracellular compounds:
        gapfilling_penalties = dict()
        for i in range(0, 1000):
            if indexlist[i] > 10:
                if extend_with_template == 1:
                    new_penalties = self.temp_extend_model_index_for_gapfilling(
                        i, input_templates
                    )
                    gapfilling_penalties.update(new_penalties)
                if i < len(source_models) and source_models[i] != None:
                    new_penalties = self.mdl_extend_model_index_for_gapfilling(
                        i, source_models[i], model_penalty
                    )
                    gapfilling_penalties.update(new_penalties)
        # Rescaling penalties by reaction scores and saving genes
        for reaction in gapfilling_penalties:
            array = reaction.split("_")
            rxnid = array[0]
            if rxnid in reaction_scores:
                highest_score = 0
                for gene in reaction_scores[rxnid]:
                    if highest_score < reaction_scores[rxnid][gene]:
                        highest_score = reaction_scores[rxnid][gene]
                factor = 1 - 0.9 * highest_score
                if "reverse" in gapfilling_penalties[reaction]:
                    penalties[reaction.id]["reverse"] = (
                        factor * penalties[reaction.id]["reverse"]
                    )
                if "forward" in gapfilling_penalties[reaction]:
                    penalties[reaction.id]["forward"] = (
                        factor * penalties[reaction.id]["forward"]
                    )
        self.cobramodel.solver.update()
        return gapfilling_penalties

    # Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def mdl_extend_model_index_for_gapfilling(self, index, source_model, model_penalty):
        new_metabolites = {}
        new_reactions = {}
        new_exchange = []
        new_demand = []
        new_penalties = dict()
        local_remap = {}

        comp = re.compile("(.*_*)(.)\d+$")
        for modelcompound in source_model.metabolites:
            cobra_metabolite = self.convert_modelcompound(modelcompound)
            original_id = cobra_metabolite.id
            groups = comp.match(cobra_metabolite.compartment)
            if groups[2] == "e":
                cobra_metabolite.compartment = groups[1] + groups[2] + "0"
                groups = comp.match(cobra_metabolite.id)
                cobra_metabolite.id = groups[1] + groups[2] + "0"
            else:
                cobra_metabolite.compartment = groups[1] + groups[2] + str(index)
                groups = comp.match(cobra_metabolite.id)
                cobra_metabolite.id = groups[1] + groups[2] + str(index)
            if (
                cobra_metabolite.id not in self.cobramodel.metabolites
                and cobra_metabolite.id not in new_metabolites
            ):
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                if original_id in self.auto_sink:
                    self.demand_compounds.add(cobra_metabolite.id)
                    new_demand.append(cobra_metabolite)
                if cobra_metabolite.compartment == self.auto_exchange:
                    self.exchange_compounds.add(cobra_metabolite.id)
                    new_exchange.append(cobra_metabolite)
            if cobra_metabolite.id in self.cobramodel.metabolites:
                cobra_metabolite = self.cobramodel.metabolites.get_by_id(
                    cobra_metabolite.id
                )
            else:  # Just in case the same compound is added twice - we want to switch the metabolite to the first new version
                cobra_metabolite = new_metabolites[cobra_metabolite.id]
            local_remap[original_id] = cobra_metabolite
        # Adding all metabolites to model prior to adding reactions
        self.cobramodel.add_metabolites(new_metabolites.values())

        for modelreaction in source_model.reactions:
            if modelreaction.id.split("_")[0] in self.blacklist:
                next
            # cobra_reaction = self.convert_modelreaction(modelreaction)
            cobra_reaction = modelreaction.copy()
            groups = comp.match(cobra_reaction.id)
            cobra_reaction.id = groups[1] + groups[2] + str(index)
            new_penalties[cobra_reaction.id] = dict()
            # Updating metabolites in reaction to new model
            metabolites = cobra_reaction.metabolites
            new_stoichiometry = {}
            for metabolite in metabolites:
                # Adding new coefficient:
                new_stoichiometry[local_remap[metabolite.id]] = metabolites[metabolite]
                # Zeroing out current coefficients
                if local_remap[metabolite.id] != metabolite:
                    new_stoichiometry[metabolite] = 0
            cobra_reaction.add_metabolites(new_stoichiometry, combine=False)
            if (
                cobra_reaction.id not in self.cobramodel.reactions
                and cobra_reaction.id not in new_reactions
            ):
                new_reactions[cobra_reaction.id] = cobra_reaction
                new_penalties[cobra_reaction.id]["added"] = 1
                if cobra_reaction.lower_bound < 0:
                    new_penalties[cobra_reaction.id]["reverse"] = model_penalty
                if cobra_reaction.upper_bound > 0:
                    new_penalties[cobra_reaction.id]["forward"] = model_penalty
            elif (
                cobra_reaction.lower_bound < 0
                and self.cobramodel.reactions.get_by_id(cobra_reaction.id).lower_bound
                == 0
            ):
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).lower_bound = cobra_reaction.lower_bound
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).update_variable_bounds()
                new_penalties[cobra_reaction.id]["reverse"] = model_penalty
                new_penalties[cobra_reaction.id]["reversed"] = 1
            elif (
                cobra_reaction.upper_bound > 0
                and self.cobramodel.reactions.get_by_id(cobra_reaction.id).upper_bound
                == 0
            ):
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).upper_bound = cobra_reaction.upper_bound
                self.cobramodel.reactions.get_by_id(
                    cobra_reaction.id
                ).update_variable_bounds()
                new_penalties[cobra_reaction.id]["forward"] = model_penalty
                new_penalties[cobra_reaction.id]["reversed"] = 1

        # Only run this on new exchanges so we don't readd for all exchanges
        for cpd in new_exchange:
            drain_reaction = self.helper.add_drain_from_metabolite_id(cpd.id)
            if (
                drain_reaction.id not in self.cobramodel.reactions
                and drain_reaction.id not in new_reactions
            ):
                new_reactions[drain_reaction.id] = drain_reaction

        # Only run this on new demands so we don't readd for all exchanges
        for cpd_id in new_demand:
            drain_reaction = self.helper.add_drain_from_metabolite_id(
                cpd_id,
                lower_bound=self.COBRA_0_BOUND,
                upper_bound=self.COBRA_DEFAULT_UB,
                prefix="DM_",
                prefix_name="Demand for ",
                sbo="SBO:0000627",
            )
            if (
                drain_reaction.id not in self.cobramodel.reactions
                and drain_reaction.id not in new_reactions
            ):
                new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.cobramodel.add_reactions(new_reactions.values())
        return new_penalties

    # Possible new function to add to the KBaseFBAModelToCobraBuilder to extend a model with a template for gapfilling for a specific index
    def temp_extend_model_index_for_gapfilling(self, index, input_templates=[]):
        new_metabolites = {}
        new_reactions = {}
        new_exchange = []
        new_demand = []
        new_penalties = dict()
        template = None
        if index < len(input_templates):
            template = input_templates[index]
        elif index in self.fbamodel["template_refs"]:
            template = self.kbapi.get_from_ws(self.fbamodel["template_refs"][index])
        else:
            template = self.kbapi.get_from_ws(self.fbamodel["template_ref"])

        if template.info.type != "KBaseFBA.NewModelTemplate":
            raise ObjectError(
                template.info.type + " loaded when KBaseFBA.NewModelTemplate expected"
            )

        for template_compound in template.compcompounds:
            tempindex = index
            compartment = template_compound.templatecompartment_ref.split("/").pop()
            if compartment == "e":
                tempindex = 0

            cobra_metabolite = self.convert_template_compound(
                template_compound, tempindex, template
            )
            if (
                cobra_metabolite.id not in self.cobramodel.metabolites
                and cobra_metabolite.id not in new_metabolites
            ):
                new_metabolites[cobra_metabolite.id] = cobra_metabolite
                self.cobramodel.add_metabolites([cobra_metabolite])
                if cobra_metabolite.id in self.auto_sink:
                    self.demand_compounds.add(cobra_metabolite.id)
                    new_demand.append(cobra_metabolite.id)
                if cobra_metabolite.compartment == self.auto_exchange:
                    new_exchange.append(cobra_metabolite.id)
                    self.exchange_compounds.add(cobra_metabolite.id)
        # Adding all metabolites to model prior to adding reactions
        self.cobramodel.add_metabolites(new_metabolites.values())

        for template_reaction in template.reactions:
            if template_reaction.id.split("_")[0] in self.blacklist:
                continue
            cobra_reaction = self.convert_template_reaction(
                template_reaction, index, template, 1
            )
            new_penalties[cobra_reaction.id] = dict()
            if (
                cobra_reaction.id not in self.cobramodel.reactions
                and cobra_reaction.id not in new_reactions
            ):
                # Adding any template reactions missing from the present model
                new_reactions[cobra_reaction.id] = cobra_reaction
                if cobra_reaction.lower_bound < 0:
                    new_penalties[cobra_reaction.id]["reverse"] = (
                        template_reaction.base_cost + template_reaction.reverse_penalty
                    )
                if cobra_reaction.upper_bound > 0:
                    new_penalties[cobra_reaction.id]["forward"] = (
                        template_reaction.base_cost + template_reaction.forward_penalty
                    )
                new_penalties[cobra_reaction.id]["added"] = 1
            elif template_reaction.GapfillDirection == "=":
                # Adjusting directionality as needed for existing reactions
                new_penalties[cobra_reaction.id]["reversed"] = 1
                if (
                    self.cobramodel.reactions.get_by_id(cobra_reaction.id).lower_bound
                    == 0
                ):
                    self.cobramodel.reactions.get_by_id(
                        cobra_reaction.id
                    ).lower_bound = template_reaction.maxrevflux
                    self.cobramodel.reactions.get_by_id(
                        cobra_reaction.id
                    ).update_variable_bounds()
                    new_penalties[cobra_reaction.id]["reverse"] = (
                        template_reaction.base_cost + template_reaction.reverse_penalty
                    )
                if (
                    self.cobramodel.reactions.get_by_id(cobra_reaction.id).upper_bound
                    == 0
                ):
                    self.cobramodel.reactions.get_by_id(
                        cobra_reaction.id
                    ).upper_bound = template_reaction.maxforflux
                    self.cobramodel.reactions.get_by_id(
                        cobra_reaction.id
                    ).update_variable_bounds()
                    new_penalties[cobra_reaction.id]["forward"] = (
                        template_reaction.base_cost + template_reaction.forward_penalty
                    )

        # Only run this on new exchanges so we don't readd for all exchanges
        for cpd_id in new_exchange:
            drain_reaction = self.helper.add_drain_from_metabolite_id(cpd_id)
            if drain_reaction != None and drain_reaction.id not in new_reactions:
                new_reactions[drain_reaction.id] = drain_reaction

        # Only run this on new demands so we don't readd for all exchanges
        for cpd_id in new_demand:
            drain_reaction = self.helper.add_drain_from_metabolite_id(
                cpd_id, self.COBRA_0_BOUND, self.COBRA_DEFAULT_UB, "DM_", "Demand for "
            )
            if drain_reaction != None and drain_reaction.id not in new_reactions:
                new_reactions[drain_reaction.id] = drain_reaction

        # Adding all new reactions to the model at once (much faster than one at a time)
        self.cobramodel.add_reactions(new_reactions.values())
        return new_penalties

    def convert_modelreaction(self, reaction, bigg=False):
        mr_id = reaction.id
        name = reaction.name
        annotation = reaction.annotation
        lower_bound, upper_bound = reaction.get_reaction_constraints()

        id = build_rxn_id(mr_id)
        if bigg and "bigg.reaction" in annotation:
            id = annotation["bigg.reaction"]

        gpr = reaction.get_gpr()

        cobra_reaction = Reaction(
            id, name=name, lower_bound=lower_bound, upper_bound=upper_bound
        )
        cobra_reaction.annotation[
            self.SBO_ANNOTATION
        ] = "SBO:0000176"  # biochemical reaction
        cobra_reaction.annotation.update(annotation)

        if id.startswith("rxn"):
            cobra_reaction.annotation["seed.reaction"] = id.split("_")[0]

        cobra_reaction.add_metabolites(
            self.convert_modelreaction_stoichiometry(reaction)
        )

        cobra_reaction.gene_reaction_rule = reaction.gene_reaction_rule

        for genes in gpr:
            for gene in genes:
                if not gene in self.genes:
                    self.genes[gene] = gene

        return cobra_reaction

    def convert_modelcompound(self, metabolite, bigg=False):
        formula = metabolite.formula
        name = metabolite.name
        charge = metabolite.charge
        mc_id = metabolite.id
        compartment = metabolite.compartment
        annotation = metabolite.annotation

        id = build_cpd_id(mc_id)

        if bigg and "bigg.metabolite" in annotation:
            id = annotation["bigg.metabolite"] + "_" + compartment
            # print(id)

        met = Metabolite(
            id, formula=formula, name=name, charge=charge, compartment=compartment
        )

        met.annotation[
            self.SBO_ANNOTATION
        ] = "SBO:0000247"  # simple chemical - Simple, non-repetitive chemical entity.
        if id.startswith("cpd"):
            met.annotation["seed.compound"] = id.split("_")[0]
        met.annotation.update(annotation)
        return met

    def convert_modelreaction_stoichiometry(self, reaction):
        object_stoichiometry = {}
        s = reaction.stoichiometry
        for metabolite_id in s:
            if metabolite_id in self.metabolites_remap:
                object_stoichiometry[
                    self.cobramodel.metabolites.get_by_id(
                        self.metabolites_remap[metabolite_id]
                    )
                ] = s[metabolite_id]
        return object_stoichiometry

    def create_binary_variables(self, rxnobj, forward=1, reverse=1):
        if rxnobj.id not in self.binary_flux_variables:
            self.binary_flux_variables[rxnobj.id] = dict()
            self.binary_flux_constraints[rxnobj.id] = dict()
        if (
            forward == 1
            and rxnobj.upper_bound > 0
            and "forward" not in self.binary_flux_variables[rxnobj.id]
        ):
            self.binary_flux_variables[rxnobj.id][
                "forward"
            ] = self.cobramodel.problem.Variable(
                rxnobj.id + "_fb", lb=0, ub=1, type="binary"
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_variables[rxnobj.id]["forward"]
            )
            self.binary_flux_constraints[rxnobj.id][
                "forward"
            ] = self.cobramodel.problem.Constraint(
                1000 * self.binary_flux_variables[rxnobj.id]["forward"]
                - rxnobj.forward_variable,
                lb=0,
                ub=None,
                name=rxnobj.id + "_fb",
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_constraints[rxnobj.id]["forward"]
            )
        if (
            reverse == 1
            and rxnobj.lower_bound < 0
            and "reverse" not in self.binary_flux_variables[rxnobj.id]
        ):
            self.binary_flux_variables[rxnobj.id][
                "reverse"
            ] = self.cobramodel.problem.Variable(
                rxnobj.id + "_bb", lb=0, ub=1, type="binary"
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_variables[rxnobj.id]["reverse"]
            )
            self.binary_flux_constraints[rxnobj.id][
                "reverse"
            ] = self.cobramodel.problem.Constraint(
                1000 * self.binary_flux_variables[rxnobj.id]["reverse"]
                - rxnobj.forward_variable,
                lb=0,
                ub=None,
                name=rxnobj.id + "_bb",
            )
            self.cobramodel.add_cons_vars(
                self.binary_flux_constraints[rxnobj.id]["reverse"]
            )

    def binary_check_gapfilling_solution(
        self, gapfilling_penalties, add_solution_exclusion_constraint
    ):
        objcoef = {}
        flux_values = self.compute_flux_values_from_variables()
        for rxnobj in self.cobramodel.reactions:
            if rxnobj.id in gapfilling_penalties:
                if (
                    "reverse" in gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["reverse"] > Zero
                ):
                    self.create_binary_variables(rxnobj, 0, 1)
                    objcoef[self.binary_flux_variables[rxnobj.id]["reverse"]] = 1
                if (
                    "forward" in gapfilling_penalties[rxnobj.id]
                    and flux_values[rxnobj.id]["forward"] > Zero
                ):
                    self.create_binary_variables(rxnobj, 1, 0)
                    objcoef[self.binary_flux_variables[rxnobj.id]["forward"]] = 1
        with self.cobramodel:
            # Setting all gapfilled reactions not in the solution to zero
            min_reaction_objective = self.cobramodel.problem.Objective(
                Zero, direction="min"
            )
            for rxnobj in self.cobramodel.reactions:
                if rxnobj.id in gapfilling_penalties:
                    if (
                        "reverse" in gapfilling_penalties[rxnobj.id]
                        and flux_values[rxnobj.id]["reverse"] <= Zero
                    ):
                        rxnobj.lower_bound = 0
                    if (
                        "forward" in gapfilling_penalties[rxnobj.id]
                        and flux_values[rxnobj.id]["forward"] <= Zero
                    ):
                        rxnobj.upper_bound = 0
                    rxnobj.update_variable_bounds()
            # Setting the objective to be minimization of sum of binary variables
            self.cobramodel.objective = min_reaction_objective
            min_reaction_objective.set_linear_coefficients(objcoef)
            with open("GapfillBinary.lp", "w") as out:
                out.write(str(self.cobramodel.solver))
            self.cobramodel.optimize()
            flux_values = self.compute_flux_values_from_variables()
        if add_solution_exclusion_constraint == 1:
            self.add_binary_solution_exclusion_constraint(flux_values)
        return flux_values

    # Adds a constraint that eliminates a gapfilled solution from feasibility so a new solution can be obtained
    def add_binary_solution_exclusion_constraint(self, flux_values):
        count = len(self.solution_exclusion_constraints)
        solution_coef = {}
        solution_size = 0
        for reaction in self.binary_flux_variables:
            for direction in self.binary_flux_variables[reaction]:
                if flux_values[reaction][direction] > Zero:
                    solution_size += 1
                    solution_coef[self.binary_flux_variables[reaction][direction]] = 1
        if len(solution_coef) > 0:
            new_exclusion_constraint = self.cobramodel.problem.Constraint(
                Zero,
                lb=None,
                ub=(solution_size - 1),
                name="exclusion." + str(count + 1),
            )
            self.cobramodel.add_cons_vars(new_exclusion_constraint)
            self.cobramodel.solver.update()
            new_exclusion_constraint.set_linear_coefficients(solution_coef)
            self.solution_exclusion_constraints.append(new_exclusion_constraint)
            return new_exclusion_constraint
        return None

    # Takes gapfilled penalties and creates and objective function minimizing gapfilled reactions
    def create_minimal_reaction_objective(self, penalty_hash, default_penalty=0):
        reaction_objective = self.cobramodel.problem.Objective(Zero, direction="min")
        obj_coef = dict()
        for reaction in self.cobramodel.reactions:
            if reaction.id in penalty_hash:
                # Minimizing gapfilled reactions
                if "reverse" in penalty_hash[reaction.id]:
                    obj_coef[reaction.reverse_variable] = abs(
                        penalty_hash[reaction.id]["reverse"]
                    )
                elif default_penalty != 0:
                    obj_coef[reaction.reverse_variable] = default_penalty
                if "forward" in penalty_hash[reaction.id]:
                    obj_coef[reaction.forward_variable] = abs(
                        penalty_hash[reaction.id]["forward"]
                    )
                elif default_penalty != 0:
                    obj_coef[reaction.forward_variable] = default_penalty
            else:
                obj_coef[reaction.forward_variable] = default_penalty
                obj_coef[reaction.reverse_variable] = default_penalty

        self.cobramodel.objective = reaction_objective
        reaction_objective.set_linear_coefficients(obj_coef)

    # Required this function to add gapfilled compounds to a KBase model for saving gapfilled model
    def convert_cobra_compound_to_kbcompound(self, cpd, kbmodel, add_to_model=1):
        refid = "cpd00000"
        if re.search("cpd\d+_[a-z]+", cpd.id):
            refid = cpd.id
            refid = re.sub("_[a-z]\d+$", "", refid)
        cpd_data = {
            "aliases": [],
            "charge": cpd.charge,
            "compound_ref": "~/template/compounds/id/" + refid,
            "dblinks": {},
            "formula": cpd.formula,
            "id": cpd.id,
            "inchikey": "ALYNCZNDIQEVRV-UHFFFAOYSA-M",
            "modelcompartment_ref": "~/modelcompartments/id/" + cpd.id.split("_").pop(),
            "name": cpd.name(),
            "numerical_attributes": {},
            "string_attributes": {},
        }
        cpd_data = AttrDict(cpd_data)
        if add_to_model == 1:
            kbmodel.modelcompounds.append(cpd_data)
        return cpd_data

    # Required this function to add gapfilled reactions to a KBase model for saving gapfilled model
    def convert_cobra_reaction_to_kbreaction(
        self, rxn, kbmodel, direction="=", add_to_model=1
    ):
        rxnref = "~/template/reactions/id/rxn00000_c"
        if re.search("rxn\d+_[a-z]+", rxn.id):
            rxnref = "~/template/reactions/id/" + rxn.id
            rxnref = re.sub("\d+$", "", rxnref)
        rxn_data = {
            "id": rxn.id,
            "aliases": [],
            "dblinks": {},
            "direction": direction,
            "edits": {},
            "gapfill_data": {},
            "maxforflux": 1000000,
            "maxrevflux": 1000000,
            "modelReactionProteins": [],
            "modelReactionReagents": [],
            "modelcompartment_ref": "~/modelcompartments/id/" + rxn.id.split("_").pop(),
            "name": rxn.name,
            "numerical_attributes": {},
            "probability": 0,
            "protons": 0,
            "reaction_ref": rxnref,
            "string_attributes": {},
        }
        rxn_data = AttrDict(rxn_data)
        for cpd in rxn.metabolites:
            if cpd.id not in kbmodel.modelcompounds:
                convert_cobra_compound_to_kbcompound(cpd, kbmodel, 1)
            rxn_data.modelReactionReagents.append(
                {
                    "coefficient": rxn.metabolites[cpd],
                    "modelcompound_ref": "~/modelcompounds/id/" + cpd.id,
                }
            )
        if add_to_model == 1:
            kbmodel.modelreactions.append(rxn_data)
        return rxn_data

    def convert_objective_to_constraint(self, lower_bound, upper_bound):
        old_obj_variable = self.cobramodel.problem.Variable(
            name="old_objective_variable", lb=lower_bound, ub=upper_bound
        )
        old_obj_constraint = self.cobramodel.problem.Constraint(
            self.cobramodel.solver.objective.expression - old_obj_variable,
            lb=0,
            ub=0,
            name="old_objective_constraint",
        )
        self.cobramodel.add_cons_vars([old_obj_variable, old_obj_constraint])

    def compute_flux_values_from_variables(self):
        flux_values = {}
        for rxnobj in self.cobramodel.reactions:
            flux_values[rxnobj.id] = {}
            flux_values[rxnobj.id]["reverse"] = rxnobj.reverse_variable.primal
            flux_values[rxnobj.id]["forward"] = rxnobj.forward_variable.primal
        return flux_values

    def compute_gapfilled_solution(self, penalties, flux_values=None):
        if flux_values == None:
            flux_values = self.compute_flux_values_from_variables()
        output = {"reversed": {}, "new": {}}
        for reaction in self.cobramodel.reactions:
            if reaction.id in penalties:
                if (
                    flux_values[reaction.id]["forward"] > Zero
                    and "forward" in penalties[reaction.id]
                ):
                    if "added" in penalties[reaction.id]:
                        output["new"][reaction.id] = ">"
                    else:
                        output["reversed"][reaction.id] = ">"
                elif (
                    flux_values[reaction.id]["reverse"] > Zero
                    and "reverse" in penalties[reaction.id]
                ):
                    if "added" in penalties[reaction.id]:
                        output["new"][reaction.id] = "<"
                    else:
                        output["reversed"][reaction.id] = "<"
        return output

    def add_gapfilling_solution_to_kbase_model(self, newmodel, penalties, media_ref):
        gfid = None
        if gfid == None:
            largest_index = 0
            for gapfilling in newmodel.gapfillings:
                current_index = gapfilling.id.split(".").pop()
                if largest_index == 0 or largest_index < current_index:
                    largest_index = current_index
            gfid = "gf." + str(largest_index + 1)
        newmodel.gapfillings.append(
            {
                "gapfill_id": newmodel.id + "." + gfid,
                "id": gfid,
                "integrated": 1,
                "integrated_solution": "0",
                "media_ref": media_ref,
            }
        )
        for reaction in self.cobramodel.reactions:
            if reaction.id in penalties:
                if (
                    reaction.forward_variable.primal > Zero
                    and "forward" in penalties[reaction.id]
                ):
                    if reaction.id not in newmodel.modelreactions:
                        self.convert_cobra_reaction_to_kbreaction(
                            reaction, newmodel, ">", 1
                        )
                    gfrxn = newmodel.modelreactions.get_by_id(reaction.id)
                    gfrxn.gapfill_data[gfid] = dict()
                    gfrxn.gapfill_data[gfid]["0"] = [">", 1, []]
                elif (
                    reaction.forward_variable.primal > Zero
                    and "reverse" in penalties[reaction.id]
                ):
                    if reaction.id not in newmodel.modelreactions:
                        self.convert_cobra_reaction_to_kbreaction(
                            reaction, newmodel, "<", 1
                        )
                    gfrxn = newmodel.modelreactions.get_by_id(reaction.id)
                    gfrxn.gapfill_data[gfid] = dict()
                    gfrxn.gapfill_data[gfid]["0"] = ["<", 1, []]

    def compute_reaction_scores(self, weigh_all_events_equally=1, weights=None):
        reaction_genes = {}
        if "genome_ref" in self.fbamodel:
            anno_api = annotation_ontology_api()
            events = anno_api.get_annotation_ontology_events(
                {
                    "input_ref": self.fbamodel["genome_ref"],
                }
            )
            for event in events:
                for gene in event["ontology_terms"]:
                    if "modelseed_ids" in event["ontology_terms"][gene]:
                        for rxn in event["ontology_terms"][gene]["modelseed_ids"]:
                            newrxn = re.sub("^MSRXN:", "", rxn)
                            if newrxn not in reaction_genes:
                                reaction_genes[newrxn] = {}
                            if gene not in reaction_genes[newrxn]:
                                reaction_genes[newrxn][gene] = 0
                            if weigh_all_events_equally == 1 or weights == None:
                                reaction_genes[newrxn][gene] += 1
                            elif event["description"] in weights:
                                reaction_genes[newrxn][gene] += weights[
                                    event["description"]
                                ]
                            elif event["event_id"] in weights:
                                reaction_genes[newrxn][gene] += weights[
                                    event["event_id"]
                                ]
                            elif event["id"] in weights:
                                reaction_genes[newrxn][gene] += weights[event["id"]]
        return reaction_genes

    def replicate_model(self, count):
        newmodel = Model(self.cobramodel.id + "_rep" + str(count))
        utilities = KBaseFBAUtilities(
            newmodel,
            newmodel,
            self.kbapi,
            self.media,
            default_uptake=self.default_uptake,
            default_excretion=self.default_excretion,
            blacklist=self.blacklist,
        )
        metabolites = []
        reactions = []
        metabolite_hash = {}
        for i in range(0, count):
            for metabolite in self.cobramodel.metabolites:
                metabolite = metabolite.copy()
                metabolite.id = metabolite.id + "__" + str(i)
                metabolite_hash[metabolite.id] = metabolite
                metabolites.append(metabolite)
            for reaction in self.cobramodel.reactions:
                reaction = reaction.copy()
                reaction.id = reaction.id + "__" + str(i)
                input_metabolites = {}
                for metabolite in reaction.metabolites:
                    newid = metabolite.id + "__" + str(i)
                    input_metabolites[metabolite_hash[newid]] = reaction.metabolites[
                        metabolite
                    ]
                reaction.add_metabolites(input_metabolites, combine=False)
                reactions.append(reaction)
        newmodel.add_metabolites(metabolites)
        newmodel.add_reactions(reactions)
        return utilities
