# -*- coding: utf-8 -*-
import logging
import itertools
from enum import Enum
import cobra
from modelseedpy.core.exceptions import ModelSEEDError
from modelseedpy.core.rast_client import RastClient
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.mstemplate import TemplateReactionType
from modelseedpy.core.msmodel import (
    get_gpr_string,
    get_reaction_constraints_from_direction,
)
from cobra.core import Gene, Metabolite, Model, Reaction, Group
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager
from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
from modelseedpy.biochem.modelseed_to_cobra import modelseed_to_cobra_reaction

SBO_ANNOTATION = "sbo"

DEFAULT_SINKS = {
    "cpd02701_c": 1000,  # S-Adenosyl-4-methylthio-2-oxobutanoate
    "cpd11416_c": 1000,  # Biomass
    "cpd15302_c": 1000,  # glycogen(n-1)
    "cpd03091_c": 1000,  # 5'-Deoxyadenosine
    "cpd01042_c": 1000,  # p-Cresol
}

logger = logging.getLogger(__name__)

### temp stuff ###
core_biomass = {
    "cpd00032_c": -1.7867,
    "cpd00005_c": -1.8225,
    "cpd00169_c": -1.496,
    "cpd11416_c": 1,
    "cpd00003_c": -3.547,
    "cpd00008_c": 41.257,
    "cpd00024_c": -1.0789,
    "cpd00009_c": 41.257,
    "cpd00102_c": -0.129,
    "cpd00101_c": -0.8977,
    "cpd00236_c": -0.8977,
    "cpd00002_c": -41.257,
    "cpd00022_c": -3.7478,
    "cpd00020_c": -2.8328,
    "cpd00006_c": 1.8225,
    "cpd00001_c": -41.257,
    "cpd00072_c": -0.0709,
    "cpd00010_c": 3.7478,
    "cpd00004_c": 3.547,
    "cpd00061_c": -0.5191,
    "cpd00067_c": 46.6265,
    "cpd00079_c": -0.205,
}

core_atp2 = {
    "cpd00067_c": 46.6265,
    "cpd00002_c": -41.257,
    "cpd00008_c": 41.257,
    "cpd00001_c": -41.257,
    "cpd00009_c": 41.257,
}

core_atp = {
    "cpd00067_c": 1,
    "cpd00002_c": -1,
    "cpd00008_c": 1,
    "cpd00001_c": -1,
    "cpd00009_c": 1,
}

gramneg = {
    "cpd11463_c0": 1,
    "cpd00008_c0": 40,
    "cpd00067_c0": 40,
    "cpd00009_c0": 40,
    "cpd00001_c0": -40,
    "cpd00002_c0": -40,
    "cpd00166_c0": -0.00280615915959131,
    "cpd00087_c0": -0.00280615915959131,
    "cpd15560_c0": -0.00280615915959131,
    "cpd00028_c0": -0.00280615915959131,
    "cpd10515_c0": -0.00280615915959131,
    "cpd15665_c0": -0.0250105977108944,
    "cpd12370_c0": 0.00280615915959131,
    "cpd15500_c0": -0.00280615915959131,
    "cpd00220_c0": -0.00280615915959131,
    "cpd00003_c0": -0.00280615915959131,
    "cpd00557_c0": -0.00280615915959131,
    "cpd00050_c0": -0.00280615915959131,
    "cpd00264_c0": -0.00280615915959131,
    "cpd00010_c0": -0.00280615915959131,
    "cpd15533_c0": -0.0311453449430676,
    "cpd15540_c0": -0.0311453449430676,
    "cpd00048_c0": -0.00280615915959131,
    "cpd00030_c0": -0.00280615915959131,
    "cpd00034_c0": -0.00280615915959131,
    "cpd00201_c0": -0.00280615915959131,
    "cpd00016_c0": -0.00280615915959131,
    "cpd00104_c0": -0.00280615915959131,
    "cpd11493_c0": -0.00280615915959131,
    "cpd00017_c0": -0.00280615915959131,
    "cpd00042_c0": -0.00280615915959131,
    "cpd00149_c0": -0.00280615915959131,
    "cpd00058_c0": -0.00280615915959131,
    "cpd03736_c0": -0.0250105977108944,
    "cpd15352_c0": -0.00280615915959131,
    "cpd00006_c0": -0.00280615915959131,
    "cpd00345_c0": -0.00280615915959131,
    "cpd00063_c0": -0.00280615915959131,
    "cpd00118_c0": -0.00280615915959131,
    "cpd15793_c0": -0.0311453449430676,
    "cpd01997_c0": 0.00280615915959131,
    "cpd00037_c0": -0.00280615915959131,
    "cpd00099_c0": -0.00280615915959131,
    "cpd02229_c0": -0.0250105977108944,
    "cpd15666_c0": 0.0250105977108944,
    "cpd10516_c0": -0.00280615915959131,
    "cpd00056_c0": -0.00280615915959131,
    "cpd03422_c0": 0.00280615915959131,
    "cpd00015_c0": -0.00280615915959131,
    "cpd00205_c0": -0.00280615915959131,
    "cpd00254_c0": -0.00280615915959131,
    "cpd11463_c0": -0.5,
    "cpd11461_c0": -0.1,
    "cpd11462_c0": -0.2,
}

grampos = {
    "cpd11416_c0": 1,
    "cpd00001_c0": -40,
    "cpd00009_c0": 40,
    "cpd00008_c0": 40,
    "cpd00002_c0": -40,
    "cpd00067_c0": 40,
    "cpd00017_c0": -0.00719527989638797,
    "cpd03422_c0": 0.00719527989638797,
    "cpd00557_c0": -0.00719527989638797,
    "cpd00166_c0": -0.00719527989638797,
    "cpd00030_c0": -0.00719527989638797,
    "cpd00087_c0": -0.00719527989638797,
    "cpd00015_c0": -0.00719527989638797,
    "cpd00038_c0": -0.0424026391792249,
    "cpd15667_c0": -0.00309563020839783,
    "cpd15540_c0": -0.0251172136637642,
    "cpd15748_c0": -0.00309563020839783,
    "cpd00048_c0": -0.00719527989638797,
    "cpd12370_c0": 0.00719527989638797,
    "cpd00052_c0": -0.0261242266150642,
    "cpd15757_c0": -0.00309563020839783,
    "cpd15533_c0": -0.0251172136637642,
    "cpd00006_c0": -0.00719527989638797,
    "cpd10515_c0": -0.00719527989638797,
    "cpd00104_c0": -0.00719527989638797,
    "cpd00028_c0": -0.00719527989638797,
    "cpd00118_c0": -0.00719527989638797,
    "cpd00037_c0": -0.00719527989638797,
    "cpd15793_c0": -0.0251172136637642,
    "cpd00010_c0": -0.00719527989638797,
    "cpd11493_c0": -0.00719527989638797,
    "cpd00264_c0": -0.00719527989638797,
    "cpd15766_c0": -0.00309563020839783,
    "cpd00056_c0": -0.00719527989638797,
    "cpd01997_c0": 0.00719527989638797,
    "cpd15668_c0": -0.00309563020839783,
    "cpd00254_c0": -0.00719527989638797,
    "cpd02229_c0": -0.00309563020839783,
    "cpd00003_c0": -0.00719527989638797,
    "cpd15665_c0": -0.00309563020839783,
    "cpd11459_c0": -0.00309563020839783,
    "cpd15666_c0": 0.0123825208335913,
    "cpd00050_c0": -0.00719527989638797,
    "cpd00063_c0": -0.00719527989638797,
    "cpd00205_c0": -0.00719527989638797,
    "cpd00042_c0": -0.00719527989638797,
    "cpd00034_c0": -0.00719527989638797,
    "cpd15500_c0": -0.00719527989638797,
    "cpd15775_c0": -0.00309563020839783,
    "cpd00099_c0": -0.00719527989638797,
    "cpd00220_c0": -0.00719527989638797,
    "cpd00062_c0": -0.0282246669459237,
    "cpd00345_c0": -0.00719527989638797,
    "cpd10516_c0": -0.00719527989638797,
    "cpd00201_c0": -0.00719527989638797,
    "cpd15669_c0": -0.00309563020839783,
    "cpd15560_c0": -0.00719527989638797,
    "cpd00149_c0": -0.00719527989638797,
    "cpd00058_c0": -0.00719527989638797,
    "cpd00016_c0": -0.00719527989638797,
    "cpd15352_c0": -0.00719527989638797,
    "cpd11463_c0": -0.5,
    "cpd11461_c0": -0.1,
    "cpd11462_c0": -0.2,
}


class MSGenomeClass(Enum):
    P = "Gram Positive"
    N = "Gram Negative"
    C = "Cyano"
    A = "Archaea"


def build_biomass(rxn_id, cobra_model, template, biomass_compounds, index="0"):
    bio_rxn = Reaction(rxn_id, "biomass", "", 0, 1000)
    metabolites = {}
    for cpd_id in biomass_compounds:
        if cpd_id in cobra_model.metabolites:
            cpd = cobra_model.metabolites.get_by_id(cpd_id)
            metabolites[cpd] = biomass_compounds[cpd_id]
        else:
            cpd = template.compcompounds.get_by_id(cpd_id[:-1])
            compartment = f"{cpd.compartment}{index}"
            name = f"{cpd.name}_{compartment}"
            cpd = Metabolite(cpd_id, cpd.formula, name, cpd.charge, compartment)
            metabolites[cpd] = biomass_compounds[cpd_id]
    bio_rxn.add_metabolites(metabolites)
    bio_rxn.annotation[SBO_ANNOTATION] = "SBO:0000629"
    return bio_rxn


def _aaaa(genome, ontology_term):
    search_name_to_genes = {}
    search_name_to_orginal = {}
    for f in genome.features:
        if ontology_term in f.ontology_terms:
            functions = f.ontology_terms[ontology_term]
            for function in functions:
                f_norm = normalize_role(function)
                if f_norm not in search_name_to_genes:
                    search_name_to_genes[f_norm] = set()
                    search_name_to_orginal[f_norm] = set()
                search_name_to_orginal[f_norm].add(function)
                search_name_to_genes[f_norm].add(f.id)
    return search_name_to_genes, search_name_to_orginal


def aux_template(template):
    rxn_roles = {}
    roles = dict(map(lambda x: (x["id"], x), template.roles))
    for r in template.reactions:
        rxn_roles[r.id] = set()
        complex_roles = r.get_complex_roles()
        if len(complex_roles) > 0:
            for cpx_id in complex_roles:
                for role_id in complex_roles[cpx_id]:
                    rxn_roles[r.id].add(normalize_role(roles[role_id]["name"]))
                    # print(role_id, normalize_role(roles[role_id]['name']))
    return rxn_roles


def build_gpr2(cpx_sets):
    list_of_ors = []
    for cpx in cpx_sets:
        list_of_ands = []
        for role_id in cpx_sets[cpx]:
            gene_ors = cpx_sets[cpx][role_id]
            if len(gene_ors) > 1:
                list_of_ands.append("(" + " or ".join(gene_ors) + ")")
            else:
                list_of_ands.append(list(gene_ors)[0])
        list_of_ors.append("(" + " and ".join(list_of_ands) + ")")
    if len(list_of_ors) > 1:
        return " or ".join(list_of_ors)
    return list_of_ors[0]


def build_gpr(cpx_gene_role):
    """
    example input:
     {'sdh': [{'b0721': 'sdhC', 'b0722': 'sdhD', 'b0723': 'sdhA', 'b0724': 'sdhB'}]}

     (b0721 and b0722 and b0724 and b0723)

     {'cpx1': [{'g1': 'role1', 'g3': 'role2'}, {'g2': 'role1', 'g3': 'role2'}]}

     (g1 and g3) or (g2 and g3)

    :param cpx_gene_role:
    :return:
    """
    # save cpx_id and role_id in reaction annotation for KBase object
    gpr_or_ll = []
    for cpx_id in cpx_gene_role:
        # print(cpx_id)
        for complex_set in cpx_gene_role[cpx_id]:
            # print(complex_set)
            gpr_or_ll.append("({})".format(" and ".join(set(complex_set))))

    return " or ".join(gpr_or_ll)


class MSBuilder:
    def __init__(
        self, genome, template=None, name=None, ontology_term="RAST", index="0"
    ):
        """

        @param genome: MSGenome
        @param template: MSTemplate
        @param name:
        @param ontology_term:
        """
        if index is None or type(index) != str:
            raise TypeError("index must be str")
        if ontology_term is None or type(ontology_term) != str:
            raise TypeError("ontology_term must be str")
        self.name = name
        self.genome = genome
        self.template = template
        self.genome_class = None
        self.search_name_to_genes, self.search_name_to_original = _aaaa(
            genome, ontology_term
        )
        self.template_species_to_model_species = None
        self.reaction_to_complex_sets = None
        self.compartments = None
        self.base_model = None
        self.compartments_index = None  # TODO: implement custom index by compartment
        self.index = index

    def build_drains(self):
        if self.template_species_to_model_species is None:
            logger.warning("cannot build model drains without generating model species")
            return None
        if self.template.drains:
            sinks = self.build_sinks()
            demands = self.build_demands()
            return sinks + demands
        else:
            # template without drain specification we build only default sinks
            return self.build_sinks()

    def build_sinks(self):
        if self.template_species_to_model_species is None:
            logger.warning("cannot build model sinks without generating model species")
            return None
        if self.template.drains:
            sinks = {
                x.id: t[1]
                for x, t in self.template.drains.items()
                if t[1] > 0 and x.id in self.template_species_to_model_species
            }
            return [self.build_sink_reaction(x, v) for x, v in sinks.items()]
        else:
            # template without drain specification we build only default sinks
            in_model = {
                k: v
                for k, v in DEFAULT_SINKS.items()
                if k in self.template_species_to_model_species
            }
            return [self.build_sink_reaction(x, v) for x, v in in_model.items()]

    def build_demands(self):
        if self.template_species_to_model_species is None:
            logger.warning("cannot build model sinks without generating model species")
            return None
        if self.template.drains:
            demands = {
                x.id: t[0]
                for x, t in self.template.drains.items()
                if t[0] < 0 and x.id in self.template_species_to_model_species
            }
            return [self.build_demand_reaction(x, v) for x, v in demands.items()]
        else:
            return []

    def build_drain_reaction(
        self,
        template_cpd_id,
        prefix="EX_",
        name_prefix="Exchange for ",
        subsystem="exchanges",
        lower_bound=0,
        upper_bound=1000,
        sbo_term="SBO:0000627",
    ):
        """
        SK_ for sink (SBO_0000632) DM_ for demand (SBO_0000628) EX_ for exchange (SBO_0000627)
        @param template_cpd_id:
        @param prefix:
        @param name_prefix:
        @param subsystem:
        @param lower_bound:
        @param upper_bound:
        @param sbo_term:
        @return:
        """

        if self.template_species_to_model_species is None:
            logger.warning("cannot build model drains without generating model species")
            return None
        else:
            m = self.template_species_to_model_species[template_cpd_id]
            drain = Reaction(
                f"{prefix}{m.id}",
                f"{name_prefix}{m.name}",
                subsystem,
                lower_bound,
                upper_bound,
            )
            drain.add_metabolites({m: -1})
            drain.annotation[SBO_ANNOTATION] = sbo_term
            return drain

    def build_sink_reaction(self, template_cpd_id, upper_bound):
        if upper_bound <= 0:
            raise ModelSEEDError("Sink reactions must have upper bound > 0")
        return self.build_drain_reaction(
            template_cpd_id,
            "SK_",
            "Sink for ",
            "exchanges",
            0,
            upper_bound,
            "SBO:0000632",
        )

    def build_demand_reaction(self, template_cpd_id, lower_bound):
        if lower_bound >= 0:
            raise ModelSEEDError("Demand reactions must have lower bound < 0")
        return self.build_drain_reaction(
            template_cpd_id,
            "DM_",
            "Demand for ",
            "exchanges",
            lower_bound,
            0,
            "SBO:0000628",
        )

    def _get_template_reaction_complexes(self, template_reaction):
        """

        :param template_reaction:
        :return:
        """
        template_reaction_complexes = {}
        for cpx in template_reaction.get_complexes():
            template_reaction_complexes[cpx.id] = {}
            for role, (triggering, optional) in cpx.roles.items():
                sn = normalize_role(role.name)
                template_reaction_complexes[cpx.id][role.id] = [
                    sn,
                    triggering,
                    optional,
                    set()
                    if sn not in self.search_name_to_genes
                    else set(self.search_name_to_genes[sn]),
                ]
        return template_reaction_complexes

    @staticmethod
    def _build_reaction_complex_gpr_sets2(
        match_complex, allow_incomplete_complexes=True
    ):
        complexes = {}
        for cpx_id in match_complex:
            complete = True
            roles = set()
            role_genes = {}
            for role_id in match_complex[cpx_id]:
                t = match_complex[cpx_id][role_id]
                complete &= len(t[3]) > 0 or not t[1] or t[2]
                if len(t[3]) > 0:
                    roles.add(role_id)
                    role_genes[role_id] = t[3]
            # print(cpx_id, complete, roles)
            if len(roles) > 0 and (allow_incomplete_complexes or complete):
                complexes[cpx_id] = {}
                for role_id in role_genes:
                    complexes[cpx_id][role_id] = role_genes[role_id]
                    # print(role_id, role_genes[role_id])
        return complexes

    @staticmethod
    def _build_reaction_complex_gpr_sets(
        match_complex, allow_incomplete_complexes=True
    ):
        complexes = {}
        for cpx_id in match_complex:
            complete = True
            roles = set()
            role_genes = {}
            for role_id in match_complex[cpx_id]:
                t = match_complex[cpx_id][role_id]
                complete &= (
                    len(t[3]) > 0 or not t[1] or t[2]
                )  # true if has genes or is not triggering or is optional
                # print(t[3])
                if len(t[3]) > 0:
                    roles.add(role_id)
                    role_genes[role_id] = t[3]
                # print(t)
            # it is never complete if has no genes, only needed if assuming a complex can have all
            # roles be either non triggering or optional
            logger.debug("[%s] maps to %s and complete: %s", cpx_id, roles, complete)
            if len(roles) > 0 and (allow_incomplete_complexes or complete):
                logger.debug("[%s] role_genes: %s", cpx_id, role_genes)
                ll = []
                for role_id in role_genes:
                    role_gene_set = []
                    for gene_id in role_genes[role_id]:
                        role_gene_set.append({gene_id: role_id})
                    ll.append(role_gene_set)
                logger.debug("[%s] complex lists: %s", cpx_id, ll)
                complexes[cpx_id] = list(
                    map(
                        lambda x: dict(map(lambda o: list(o.items())[0], x)),
                        itertools.product(*ll),
                    )
                )
        return complexes

    def get_gpr_from_template_reaction(
        self, template_reaction, allow_incomplete_complexes=True
    ):
        template_reaction_complexes = self._get_template_reaction_complexes(
            template_reaction
        )
        if len(template_reaction_complexes) == 0:
            return None

        # self.map_gene(template_reaction_complexes)
        # print(template_reaction_complexes)
        gpr_set = self._build_reaction_complex_gpr_sets2(
            template_reaction_complexes, allow_incomplete_complexes
        )
        return gpr_set

    @staticmethod
    def add_exchanges_to_model(model, extra_cell="e0"):
        """
        Build exchange reactions for the "extra_cell" compartment
        :param model: Cobra Model
        :param extra_cell: compartment representing extracellular
        :return:
        """
        reactions_exchanges = []
        for m in model.metabolites:
            if m.compartment == extra_cell:
                rxn_exchange_id = "EX_" + m.id
                if rxn_exchange_id not in model.reactions:
                    rxn_exchange = Reaction(
                        rxn_exchange_id,
                        "Exchange for " + m.name,
                        "exchanges",
                        -1000,
                        1000,
                    )
                    rxn_exchange.add_metabolites({m: -1})
                    rxn_exchange.annotation[SBO_ANNOTATION] = "SBO:0000627"
                    reactions_exchanges.append(rxn_exchange)
        model.add_reactions(reactions_exchanges)

        return reactions_exchanges

    @staticmethod
    def get_or_create_metabolite(
        model, template, cpd_base_id, compartment="c", index=0
    ):
        if isinstance(index, int):
            index = str(index)
        full_id = cpd_base_id + "_" + compartment + index
        if full_id not in model.metabolites:
            pass
        return model.metabolites.get_by_id(full_id)

    def build_static_biomasses(self, model, template):
        res = []
        if template.name.startswith("CoreModel"):
            res.append(self.build_biomass("bio1", model, template, core_biomass))
            res.append(self.build_biomass("bio2", model, template, core_atp))
        if template.name.startswith("GramNeg"):
            res.append(self.build_biomass("bio1", model, template, gramneg))
        if template.name.startswith("GramPos"):
            res.append(self.build_biomass("bio1", model, template, grampos))
        return res

    def auto_select_template(self):
        """

        :return: genome class
        """
        from modelseedpy.helpers import get_template, get_classifier
        from modelseedpy.core.mstemplate import MSTemplateBuilder

        genome_classifier = get_classifier("knn_ACNP_RAST_filter_01_17_2023")
        self.genome_class = genome_classifier.classify(self.genome)

        # TODO: update with enum MSGenomeClass
        template_genome_scale_map = {
            "A": "template_gram_neg",
            "C": "template_gram_neg",
            "N": "template_gram_neg",
            "P": "template_gram_pos",
        }
        template_core_map = {
            "A": "template_core",
            "C": "template_core",
            "N": "template_core",
            "P": "template_core",
        }

        if (
            self.genome_class in template_genome_scale_map
            and self.genome_class in template_core_map
        ):
            self.template = MSTemplateBuilder.from_dict(
                get_template(template_genome_scale_map[self.genome_class])
            ).build()
        elif self.template is None:
            raise Exception(f"unable to select template for {self.genome_class}")

        return self.genome_class

    def generate_reaction_complex_sets(self, allow_incomplete_complexes=True):
        self.reaction_to_complex_sets = {}
        for template_reaction in self.template.reactions:
            gpr_set = self.get_gpr_from_template_reaction(
                template_reaction, allow_incomplete_complexes
            )
            if gpr_set:
                self.reaction_to_complex_sets[template_reaction.id] = gpr_set
                logger.debug("[%s] gpr set: %s", template_reaction.id, gpr_set)

        return self.reaction_to_complex_sets

    """
    def _build_reaction(self, reaction_id, gpr_set, template, index="0", sbo=None):
        template_reaction = template.reactions.get_by_id(reaction_id)

        reaction_compartment = template_reaction.compartment
        metabolites = {}

        for cpd, value in template_reaction.metabolites.items():
            compartment = f"{cpd.compartment}{index}"
            name = f"{cpd.name}_{compartment}"
            cpd = Metabolite(
                cpd.id + str(index), cpd.formula, name, cpd.charge, compartment
            )
            metabolites[cpd] = value

        reaction = Reaction(
            "{}{}".format(template_reaction.id, index),
            "{}_{}{}".format(template_reaction.name, reaction_compartment, index),
            "",
            template_reaction.lower_bound,
            template_reaction.upper_bound,
        )

        gpr_str = build_gpr2(gpr_set) if gpr_set else ""
        reaction.add_metabolites(metabolites)
        if gpr_str and len(gpr_str) > 0:
            reaction.gene_reaction_rule = gpr_str  # get_gpr_string(gpr_ll)

        reaction.annotation["seed.reaction"] = template_reaction.reference_id
        if sbo:
            reaction.annotation[SBO_ANNOTATION] = sbo
        return reaction
    """

    def build_complex_groups(self, complex_sets):
        """
        Builds complex Group from complex sets computed from template and genome
        Example: {'cpx00700': {'ftr01608': {'b3177'}}, 'cpx01370': {'ftr01607': {'b0142'}}}
        @param complex_sets:
        @return:
        """
        group_complexes = {}
        for complex_set in complex_sets:
            for complex_id in complex_set:
                if (
                    complex_id not in group_complexes
                    and complex_id in self.template.complexes
                ):
                    cpx = self.template.complexes.get_by_id(complex_id)
                    g = Group(complex_id)
                    g.notes["complex_source"] = cpx.source
                    for role, (t, o) in cpx.roles.items():
                        if role.id in complex_set[complex_id]:
                            g.notes[f"complex_subunit_note_{role.id}"] = role.name
                            g.notes[f"complex_subunit_optional_{role.id}"] = (
                                1 if o else 0
                            )
                            g.notes[f"complex_subunit_triggering_{role.id}"] = (
                                1 if t else 0
                            )
                            g.notes[f"complex_subunit_features_{role.id}"] = ";".join(
                                sorted(list(complex_set[complex_id][role.id]))
                            )
                    group_complexes[g.id] = g

        return group_complexes

    def build_metabolic_reactions(self):
        if self.base_model is None:
            raise ModelSEEDError(
                "unable to generate metabolic reactions without base model"
            )
        if self.reaction_to_complex_sets is None:
            raise ModelSEEDError(
                "unable to generate metabolic reactions without generate complex sets"
            )

        if self.template_species_to_model_species is None:
            self.template_species_to_model_species = {}
        if self.compartments is None:
            self.compartments = {}

        reactions = []
        for rxn_id, complex_set in self.reaction_to_complex_sets.items():
            template_reaction = self.template.reactions.get_by_id(rxn_id)
            for m in template_reaction.metabolites:
                if m.compartment not in self.compartments:
                    self.compartments[
                        m.compartment
                    ] = self.template.compartments.get_by_id(m.compartment)
                if m.id not in self.template_species_to_model_species:
                    model_metabolite = m.to_metabolite(self.index)
                    self.template_species_to_model_species[m.id] = model_metabolite
                    self.base_model.add_metabolites([model_metabolite])
            reaction = template_reaction.to_reaction(self.base_model, self.index)
            gpr_str = build_gpr2(complex_set) if complex_set else ""
            if gpr_str and len(gpr_str) > 0:
                reaction.gene_reaction_rule = gpr_str
            reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"
            reaction.notes["modelseed_complex"] = ";".join(sorted(list(complex_set)))
            reactions.append(reaction)

        return reactions

    def build_from_annotaton_ontology(
        self,
        model_or_id,
        anno_ont,
        index="0",
        allow_all_non_grp_reactions=False,
        annotate_with_rast=False,
        biomass_classic=False,
        biomass_gc=0.5,
        add_non_template_reactions=True,
        prioritized_event_list=None,
        ontologies=None,
        merge_all=True,
        convert_to_sso=True,
    ):
        # Build base model without annotation
        self.search_name_to_orginal = {}
        self.search_name_to_genes = {}
        gene_term_hash = anno_ont.get_gene_term_hash(
            prioritized_event_list, ontologies, merge_all, convert_to_sso
        )
        residual_reaction_gene_hash = {}
        for gene in gene_term_hash:
            for term in gene_term_hash[gene]:
                if term.ontology.id == "SSO":
                    name = anno_ont.get_term_name(term)
                    f_norm = normalize_role(name)
                    if f_norm not in self.search_name_to_genes:
                        self.search_name_to_genes[f_norm] = set()
                        self.search_name_to_orginal[f_norm] = set()
                    self.search_name_to_orginal[f_norm].add(name)
                    self.search_name_to_genes[f_norm].add(gene.id)
                else:
                    for rxn_id in term.msrxns:
                        if rxn_id not in residual_reaction_gene_hash:
                            residual_reaction_gene_hash[rxn_id] = {}
                        if gene not in residual_reaction_gene_hash[rxn_id]:
                            residual_reaction_gene_hash[rxn_id][gene] = []
                        residual_reaction_gene_hash[rxn_id][gene] = gene_term_hash[
                            gene
                        ][term]

        model_or_id = self.build(
            model_or_id,
            index,
            allow_all_non_grp_reactions,
            annotate_with_rast,
            biomass_classic,
            biomass_gc,
        )
        for rxn in model_or_id.reactions:
            probability = None
            for gene in rxn.genes():
                annoont_gene = anno_ont.get_feature(gene.id)
                if annoont_gene and annoont_gene in gene_term_hash:
                    for term in gene_term_hash[annoont_gene]:
                        if rxn.id[0:-3] in term.msrxns:
                            for item in gene_term_hash[gene][term]:
                                if "probability" in item.scores:
                                    if (
                                        not probability
                                        or item.scores["probability"] > probability
                                    ):
                                        probability = item.scores["probability"]
            if hasattr(rxn, "probability"):
                rxn.probability = probability

        reactions = []
        modelseeddb = ModelSEEDBiochem.get()
        for rxn_id in residual_reaction_gene_hash:
            if rxn_id + "_c0" not in model_or_id.reactions:
                reaction = None
                template_reaction = None
                if rxn_id + "_c" in self.template.reactions:
                    template_reaction = self.template.reactions.get_by_id(rxn_id + "_c")
                elif rxn_id in modelseeddb.reactions:
                    msrxn = modelseeddb.reactions.get_by_id(rxn_id)
                    template_reaction = msrxn.to_template_reaction({0: "c", 1: "e"})
                if template_reaction:
                    for m in template_reaction.metabolites:
                        if m.compartment not in self.compartments:
                            self.compartments[
                                m.compartment
                            ] = self.template.compartments.get_by_id(m.compartment)
                        if m.id not in self.template_species_to_model_species:
                            model_metabolite = m.to_metabolite(self.index)
                            self.template_species_to_model_species[
                                m.id
                            ] = model_metabolite
                            self.base_model.add_metabolites([model_metabolite])
                    reaction = template_reaction.to_reaction(
                        self.base_model, self.index
                    )
                    gpr = ""
                    probability = None
                    for gene in residual_reaction_gene_hash[rxn_id]:
                        for item in residual_reaction_gene_hash[rxn_id][gene]:
                            if "probability" in item["scores"]:
                                if (
                                    not probability
                                    or item["scores"]["probability"] > probability
                                ):
                                    probability = item["scores"]["probability"]
                        if len(gpr) > 0:
                            gpr += " or "
                        gpr += gene.id
                    if hasattr(rxn, "probability"):
                        reaction.probability = probability
                    reaction.gene_reaction_rule = gpr
                    reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"
                    reactions.append(reaction)
                if not reaction:
                    print("Reaction ", rxn_id, " not found in template or database!")

        model_or_id.add_reactions(reactions)
        return model_or_id

    def build_non_metabolite_reactions(
        self, cobra_model, allow_all_non_grp_reactions=False
    ):
        if self.base_model is None:
            raise ModelSEEDError(
                "unable to generate metabolic reactions without base model"
            )
        if self.reaction_to_complex_sets is None:
            raise ModelSEEDError(
                "unable to generate metabolic reactions without generate complex sets"
            )

        if self.template_species_to_model_species is None:
            self.template_species_to_model_species = {}
        if self.compartments is None:
            self.compartments = {}

        reactions = []
        for template_reaction in self.template.reactions:
            rxn_type = template_reaction.type
            if (
                rxn_type == "universal"
                or rxn_type == "spontaneous"
                or rxn_type == TemplateReactionType.UNIVERSAL
                or rxn_type == TemplateReactionType.SPONTANEOUS
            ):
                reaction_metabolite_ids = {m.id for m in template_reaction.metabolites}
                if (
                    len(
                        set(self.template_species_to_model_species)
                        & reaction_metabolite_ids
                    )
                    > 0
                    or allow_all_non_grp_reactions
                ):
                    for m in template_reaction.metabolites:
                        if m.compartment not in self.compartments:
                            self.compartments[
                                m.compartment
                            ] = self.template.compartments.get_by_id(m.compartment)
                        if m.id not in self.template_species_to_model_species:
                            model_metabolite = m.to_metabolite(self.index)
                            self.template_species_to_model_species[
                                m.id
                            ] = model_metabolite
                            self.base_model.add_metabolites([model_metabolite])

                    reaction = template_reaction.to_reaction(
                        self.base_model, self.index
                    )
                    reaction.annotation[SBO_ANNOTATION] = "SBO:0000672"
                    # if template_reaction.type == "spontaneous":
                    #    reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"

                    if reaction.id not in cobra_model.reactions:
                        reactions.append(reaction)

        return reactions

    def build_biomass(self, rxn_id, cobra_model, template, biomass_compounds):
        bio_rxn = Reaction(rxn_id, "biomass", "", 0, 1000)
        metabolites = {}
        for template_cpd_id in biomass_compounds:
            if template_cpd_id in self.template_species_to_model_species:
                model_species_id = self.template_species_to_model_species[
                    template_cpd_id
                ].id
                cpd = cobra_model.metabolites.get_by_id(model_species_id)
                metabolites[cpd] = biomass_compounds[template_cpd_id]
            else:
                template_cpd = template.compcompounds.get_by_id(template_cpd_id[:-1])
                m = template_cpd.to_metabolite(self.index)
                metabolites[m] = biomass_compounds[template_cpd_id]
                self.template_species_to_model_species[template_cpd_id] = m
                cobra_model.add_metabolites([m])
        bio_rxn.add_metabolites(metabolites)
        bio_rxn.annotation[SBO_ANNOTATION] = "SBO:0000629"
        return bio_rxn

    def build(
        self,
        model_or_id,
        index="0",
        allow_all_non_grp_reactions=False,
        annotate_with_rast=True,
        biomass_classic=False,
        biomass_gc=0.5,
        add_reaction_from_rast_annotation=True,
    ):
        """

        @param model_or_id: a string ID to build from cobra.core.Model otherwise a type of cobra.core.Model
        as Base Model
        @param index:
        @param allow_all_non_grp_reactions:
        @param annotate_with_rast:
        @param biomass_classic:
        @param biomass_gc:
        @return:
        """
        self.index = index

        if annotate_with_rast:
            rast = RastClient()
            res = rast.annotate_genome(self.genome)
            self.search_name_to_genes, self.search_name_to_original = _aaaa(
                self.genome, "RAST"
            )

        # rxn_roles = aux_template(self.template)  # needs to be fixed to actually reflect template GPR rules
        if self.template is None:
            self.auto_select_template()

        cobra_model = model_or_id
        if type(model_or_id) == str:
            from cobra.core import Model

            cobra_model = Model(model_or_id)

        self.base_model = cobra_model

        self.generate_reaction_complex_sets()
        complex_groups = self.build_complex_groups(
            self.reaction_to_complex_sets.values()
        )

        if add_reaction_from_rast_annotation:
            metabolic_reactions = self.build_metabolic_reactions()
            cobra_model.add_reactions(metabolic_reactions)

        non_metabolic_reactions = self.build_non_metabolite_reactions(
            cobra_model, allow_all_non_grp_reactions
        )
        cobra_model.add_reactions(non_metabolic_reactions)
        cobra_model.add_groups(list(complex_groups.values()))
        self.add_exchanges_to_model(cobra_model)

        biomass_reactions = []
        for rxn_biomass in self.template.biomasses:
            reaction = rxn_biomass.build_biomass(
                cobra_model, index, biomass_classic, biomass_gc
            )
            for m in reaction.metabolites:
                if "modelseed_template_id" in m.notes:
                    self.template_species_to_model_species[
                        m.notes["modelseed_template_id"]
                    ] = m
            biomass_reactions.append(reaction)

        if len(biomass_reactions) > 0:
            for rxn in biomass_reactions:
                if rxn.id not in cobra_model.reactions:
                    cobra_model.add_reactions([rxn])
            cobra_model.objective = biomass_reactions[0].id

        """
        if (
            self.template.name.startswith("CoreModel")
            or self.template.name.startswith("GramNeg")
            or self.template.name.startswith("GramPos")
        ):
            gc = 0.5
            if hasattr(self.genome,"info"):
                gc = float(self.genome.info.metadata["GC content"])
                print("Genome custom GC:",gc)
            for bio in self.template.biomasses:
                bio.build_biomass(cobra_model, index, classic=False, GC=gc,add_to_model=True)
            cobra_model.objective = "bio1"
        """

        reactions_sinks = self.build_drains()
        cobra_model.add_reactions(reactions_sinks)

        compartment_data = {}
        for cmp_id, data in self.compartments.items():
            cmp_index_id = f"{cmp_id}{self.index}"
            compartment_data[cmp_index_id] = data.name
            kbase_compartment_data_key = f"kbase_compartment_data_{cmp_index_id}"
            kbase_compartment_data = {
                "pH": data.ph,
                "potential": 0,
                "compartmentIndex": self.index,
            }
            cobra_model.notes[kbase_compartment_data_key] = kbase_compartment_data

        cobra_model.compartments = compartment_data

        """
        for cpd_id in ["cpd02701_c0", "cpd11416_c0", "cpd15302_c0"]:
           if cpd_id in cobra_model.metabolites:
                m = cobra_model.metabolites.get_by_id(cpd_id)
                rxn_exchange = Reaction(
                    "SK_" + m.id, "Sink for " + m.name, "exchanges", 0, 1000
                )
                rxn_exchange.add_metabolites({m: -1})
                rxn_exchange.annotation[SBO_ANNOTATION] = "SBO:0000627"
                reactions_sinks.append(rxn_exchange)
        """

        return cobra_model

    @staticmethod
    def build_full_template_model(template, model_id=None, index="0"):
        """

        :param template:
        :param model_id: ID for the model otherwise template.id
        :param index: index for the metabolites
        :return:
        """
        model = MSModel(model_id if model_id else template.id, template=template)
        all_reactions = []
        for rxn in template.reactions:
            reaction = rxn.to_reaction(model, index)
            reaction.annotation["seed.reaction"] = rxn.id
            all_reactions.append(reaction)
        model.add_reactions(all_reactions)
        MSBuilder.add_exchanges_to_model(model)

        if template.name.startswith("CoreModel"):
            bio_rxn1 = build_biomass("bio1", model, template, core_biomass, index)
            bio_rxn2 = build_biomass("bio2", model, template, core_atp, index)
            model.add_reactions([bio_rxn1, bio_rxn2])
            model.objective = "bio1"
        else:
            for bio in template.biomasses:
                bio.build_biomass(
                    model, index, classic=False, GC=0.5, add_to_model=True
                )
            model.objective = "bio1"

        reactions_sinks = []
        for cpd_id in ["cpd02701_c0", "cpd11416_c0", "cpd15302_c0", "cpd03091_c0"]:
            if cpd_id in model.metabolites:
                m = model.metabolites.get_by_id(cpd_id)
                rxn_exchange = Reaction(
                    "SK_" + m.id, "Sink for " + m.name, "exchanges", 0, 1000
                )
                rxn_exchange.add_metabolites({m: -1})
                rxn_exchange.annotation[SBO_ANNOTATION] = "SBO:0000627"
                reactions_sinks.append(rxn_exchange)
        model.add_reactions(reactions_sinks)
        return model

    @staticmethod
    def build_metabolic_model(
        model_id,
        genome,
        gapfill_media=None,
        template=None,
        index="0",
        allow_all_non_grp_reactions=False,
        annotate_with_rast=True,
        gapfill_model=True,
        classic_biomass=False,
    ):
        builder = MSBuilder(genome, template)
        model = builder.build(
            model_id,
            index,
            allow_all_non_grp_reactions,
            annotate_with_rast,
            classic_biomass,
        )
        # Gapfilling model
        if gapfill_model:
            model = MSBuilder.gapfill_model(
                model, "bio1", builder.template, gapfill_media
            )
        return model

    @staticmethod
    def gapfill_model(original_mdl, target_reaction, template, media):
        FBAHelper.set_objective_from_target_reaction(original_mdl, target_reaction)
        model = cobra.io.json.from_json(cobra.io.json.to_json(original_mdl))
        pkgmgr = MSPackageManager.get_pkg_mgr(model)
        pkgmgr.getpkg("GapfillingPkg").build_package(
            {
                "default_gapfill_templates": [template],
                "gapfill_all_indecies_with_default_templates": 1,
                "minimum_obj": 0.01,
                "set_objective": 1,
            }
        )
        pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
        # with open('Gapfilling.lp', 'w') as out:
        #    out.write(str(model.solver))
        sol = model.optimize()
        gfresults = pkgmgr.getpkg("GapfillingPkg").compute_gapfilled_solution()
        for rxnid in gfresults["reversed"]:
            rxn = original_mdl.reactions.get_by_id(rxnid)
            if gfresults["reversed"][rxnid] == ">":
                rxn.upper_bound = 100
            else:
                rxn.lower_bound = -100
        for rxnid in gfresults["new"]:
            rxn = model.reactions.get_by_id(rxnid)
            rxn = rxn.copy()
            original_mdl.add_reactions([rxn])
            if gfresults["new"][rxnid] == ">":
                rxn.upper_bound = 100
                rxn.lower_bound = 0
            else:
                rxn.upper_bound = 0
                rxn.lower_bound = -100
        return original_mdl


def build_metabolic_model(
    genome,
    media=None,
    atp_test_medias=None,
    core_template=None,
    genome_scale_template=None,
):
    pass
