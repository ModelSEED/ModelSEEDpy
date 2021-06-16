import logging
from modelseedpy.core.rast_client import RastClient
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msmodel import get_gpr_string, get_reaction_constraints_from_direction


import re
import copy
import cobra
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core import FBAHelper
from modelseedpy.fbapkg import GapfillingPkg, KBaseMediaPkg

SBO_ANNOTATION = "sbo"

#from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)

### temp stuff ###
core_biomass = {
    'cpd00032_c0': -1.7867,
    'cpd00005_c0': -1.8225,
    'cpd00169_c0': -1.496,
    'cpd11416_c0': 1,
    'cpd00003_c0': -3.547,
    'cpd00008_c0': 41.257,
    'cpd00024_c0': -1.0789,
    'cpd00009_c0': 41.257,
    'cpd00102_c0': -0.129,
    'cpd00101_c0': -0.8977,
    'cpd00236_c0': -0.8977,
    'cpd00002_c0': -41.257,
    'cpd00022_c0': -3.7478,
    'cpd00020_c0': -2.8328,
    'cpd00006_c0': 1.8225,
    'cpd00001_c0': -41.257,
    'cpd00072_c0': -0.0709,
    'cpd00010_c0': 3.7478,
    'cpd00004_c0': 3.547,
    'cpd00061_c0': -0.5191,
    'cpd00067_c0': 46.6265,
    'cpd00079_c0': -0.205
}

core_atp = {
    'cpd00067_c0': 46.6265,
    'cpd00002_c0': -41.257,
    'cpd00008_c0': 41.257,
    'cpd00001_c0': -41.257,
    'cpd00009_c0': 41.257,
}

gramneg = {
    'cpd00166_c0': -0.00280615915959131,
    'cpd00087_c0': -0.00280615915959131,
    'cpd15560_c0': -0.00280615915959131,
    'cpd00028_c0': -0.00280615915959131,
    'cpd10515_c0': -0.00280615915959131,
    'cpd15665_c0': -0.0250105977108944,
    'cpd12370_c0': 0.00280615915959131,
    'cpd15500_c0': -0.00280615915959131,
    'cpd00220_c0': -0.00280615915959131,
    'cpd00003_c0': -0.00280615915959131,
    'cpd00557_c0': -0.00280615915959131,
    'cpd00002_c0': -40.1101757365074,
    'cpd00023_c0': -0.219088153012743,
    'cpd00062_c0': -0.0908319049068452,
    'cpd00050_c0': -0.00280615915959131,
    'cpd00008_c0': 40,
    'cpd00264_c0': -0.00280615915959131,
    'cpd00010_c0': -0.00280615915959131,
    'cpd15533_c0': -0.0311453449430676,
    'cpd11416_c0': 1,
    'cpd15540_c0': -0.0311453449430676,
    'cpd00048_c0': -0.00280615915959131,
    'cpd00035_c0': -0.427934380173264,
    'cpd17042_c0': -1,
    'cpd00030_c0': -0.00280615915959131,
    'cpd00034_c0': -0.00280615915959131,
    'cpd00161_c0': -0.211072732780569,
    'cpd00201_c0': -0.00280615915959131,
    'cpd00016_c0': -0.00280615915959131,
    'cpd00104_c0': -0.00280615915959131,
    'cpd00067_c0': 40,
    'cpd11493_c0': -0.00280615915959131,
    'cpd00051_c0': -0.246696822701341,
    'cpd00017_c0': -0.00280615915959131,
    'cpd00357_c0': -0.0157642107352084,
    'cpd17041_c0': -1,
    'cpd00038_c0': -0.135406821203723,
    'cpd00107_c0': -0.375388847540127,
    'cpd00042_c0': -0.00280615915959131,
    'cpd00149_c0': -0.00280615915959131,
    'cpd00058_c0': -0.00280615915959131,
    'cpd00041_c0': -0.200830806928348,
    'cpd00129_c0': -0.184354665339991,
    'cpd15432_c0': -0.0250105977108944,
    'cpd00052_c0': -0.0841036156544863,
    'cpd00012_c0': 0.484600235732628,
    'cpd15352_c0': -0.00280615915959131,
    'cpd00322_c0': -0.241798510337235,
    'cpd00053_c0': -0.219088153012743,
    'cpd00006_c0': -0.00280615915959131,
    'cpd00345_c0': -0.00280615915959131,
    'cpd00063_c0': -0.00280615915959131,
    'cpd00033_c0': -0.509869786991038,
    'cpd00066_c0': -0.154519490031345,
    'cpd17043_c0': -1,
    'cpd00118_c0': -0.00280615915959131,
    'cpd00009_c0': 39.9971938408404,
    'cpd15793_c0': -0.0311453449430676,
    'cpd00356_c0': -0.01627686799489,
    'cpd01997_c0': 0.00280615915959131,
    'cpd00132_c0': -0.200830806928348,
    'cpd00060_c0': -0.127801422590767,
    'cpd00037_c0': -0.00280615915959131,
    'cpd00115_c0': -0.0157642107352084,
    'cpd00099_c0': -0.00280615915959131,
    'cpd00156_c0': -0.352233189091625,
    'cpd02229_c0': -0.0250105977108944,
    'cpd00069_c0': -0.120676604606612,
    'cpd00065_c0': -0.0472019191450218,
    'cpd00241_c0': -0.01627686799489,
    'cpd15666_c0': 0.0250105977108944,
    'cpd10516_c0': -0.00280615915959131,
    'cpd00084_c0': -0.0761464922056484,
    'cpd00056_c0': -0.00280615915959131,
    'cpd00119_c0': -0.0792636000737159,
    'cpd00001_c0': -35.5403092430435,
    'cpd03422_c0': 0.00280615915959131,
    'cpd00015_c0': -0.00280615915959131,
    'cpd00054_c0': -0.179456352975885,
    'cpd00205_c0': -0.00280615915959131,
    'cpd00039_c0': -0.285438020490179,
    'cpd00254_c0': -0.00280615915959131
}

grampos = {
    'cpd00241_c0': -0.0116907079028565,
    'cpd00017_c0': -0.00719527989638797,
    'cpd00033_c0': -0.409331301687739,
    'cpd00066_c0': -0.176188648374102,
    'cpd17043_c0': -1,
    'cpd03422_c0': 0.00719527989638797,
    'cpd17041_c0': -1,
    'cpd00557_c0': -0.00719527989638797,
    'cpd00129_c0': -0.161028229793075,
    'cpd00166_c0': -0.00719527989638797,
    'cpd00030_c0': -0.00719527989638797,
    'cpd00087_c0': -0.00719527989638797,
    'cpd00015_c0': -0.00719527989638797,
    'cpd00065_c0': -0.0544955586831525,
    'cpd00357_c0': -0.0151844826784228,
    'cpd00009_c0': 41.2498047201036,
    'cpd00038_c0': -0.0424026391792249,
    'cpd15667_c0': -0.00309563020839783,
    'cpd00069_c0': -0.111039822579957,
    'cpd15540_c0': -0.0251172136637642,
    'cpd00161_c0': -0.186841915485094,
    'cpd15748_c0': -0.00309563020839783,
    'cpd00035_c0': -0.267560900902997,
    'cpd00048_c0': -0.00719527989638797,
    'cpd12370_c0': 0.00719527989638797,
    'cpd00052_c0': -0.0261242266150642,
    'cpd15757_c0': -0.00309563020839783,
    'cpd00053_c0': -0.261005044219309,
    'cpd15533_c0': -0.0251172136637642,
    'cpd00002_c0': -41.2913947104178,
    'cpd00006_c0': -0.00719527989638797,
    'cpd00084_c0': -0.0569540049395353,
    'cpd10515_c0': -0.00719527989638797,
    'cpd00104_c0': -0.00719527989638797,
    'cpd00051_c0': -0.193397772168782,
    'cpd00028_c0': -0.00719527989638797,
    'cpd00118_c0': -0.00719527989638797,
    'cpd00107_c0': -0.347460404235438,
    'cpd00037_c0': -0.00719527989638797,
    'cpd15793_c0': -0.0251172136637642,
    'cpd00010_c0': -0.00719527989638797,
    'cpd11493_c0': -0.00719527989638797,
    'cpd00264_c0': -0.00719527989638797,
    'cpd15766_c0': -0.00309563020839783,
    'cpd00041_c0': -0.14832625746843,
    'cpd00056_c0': -0.00719527989638797,
    'cpd01997_c0': 0.00719527989638797,
    'cpd15668_c0': -0.00309563020839783,
    'cpd00254_c0': -0.00719527989638797,
    'cpd11416_c0': 1,
    'cpd02229_c0': -0.00309563020839783,
    'cpd00003_c0': -0.00719527989638797,
    'cpd00008_c0': 41.257,
    'cpd17042_c0': -1,
    'cpd00023_c0': -0.261005044219309,
    'cpd15665_c0': -0.00309563020839783,
    'cpd11459_c0': -0.00309563020839783,
    'cpd15666_c0': 0.0123825208335913,
    'cpd00115_c0': -0.0151844826784228,
    'cpd00050_c0': -0.00719527989638797,
    'cpd00063_c0': -0.00719527989638797,
    'cpd00205_c0': -0.00719527989638797,
    'cpd00054_c0': -0.216753011604418,
    'cpd00042_c0': -0.00719527989638797,
    'cpd00034_c0': -0.00719527989638797,
    'cpd15500_c0': -0.00719527989638797,
    'cpd00156_c0': -0.307715523090583,
    'cpd00132_c0': -0.14832625746843,
    'cpd00067_c0': -41.257,
    'cpd15775_c0': -0.00309563020839783,
    'cpd00119_c0': -0.0819482085460939,
    'cpd00060_c0': -0.11349826883634,
    'cpd00001_c0': 45.354000686262,
    'cpd00099_c0': -0.00719527989638797,
    'cpd00356_c0': -0.0116907079028565,
    'cpd00220_c0': -0.00719527989638797,
    'cpd00322_c0': -0.27042908820211,
    'cpd00062_c0': -0.0282246669459237,
    'cpd00345_c0': -0.00719527989638797,
    'cpd00012_c0': 0.184896624320595,
    'cpd10516_c0': -0.00719527989638797,
    'cpd00039_c0': -0.323695423757071,
    'cpd00201_c0': -0.00719527989638797,
    'cpd15669_c0': -0.00309563020839783,
    'cpd15560_c0': -0.00719527989638797,
    'cpd00149_c0': -0.00719527989638797,
    'cpd00058_c0': -0.00719527989638797,
    'cpd00016_c0': -0.00719527989638797,
    'cpd15352_c0': -0.00719527989638797
}


def build_biomass(rxn_id, cobra_model, template, biomass_compounds):
    bio_rxn = Reaction(rxn_id, 'biomass', '', 0, 1000)
    metabolites = {}
    for cpd_id in biomass_compounds:
        if cpd_id in cobra_model.metabolites:
            cpd = cobra_model.metabolites.get_by_id(cpd_id)
            metabolites[cpd] = biomass_compounds[cpd_id]
        else:
            ccpd = template.compcompounds.get_by_id(cpd_id[:-1])
            template_cpd = template.compounds.get_by_id(cpd_id[:-3])
            compartment = ccpd.templatecompartment_ref.split('/')[-1]
            cpd = Metabolite(cpd_id, template_cpd.formula, template_cpd.name, ccpd.charge, compartment + str('0'))
            metabolites[cpd] = biomass_compounds[cpd_id]
    bio_rxn.add_metabolites(metabolites)
    bio_rxn.annotation[SBO_ANNOTATION] = "SBO:0000629"
    return bio_rxn


def aaaa(genome):
    search_name_to_genes = {}
    search_name_to_orginal = {}
    for f in genome.features:
        if 'RAST' in f.ontology_terms:
            functions = f.ontology_terms['RAST']
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
    roles = dict(map(lambda x: (x['id'], x), template.roles))
    for r in template.reactions:
        rxn_roles[r.id] = set()
        complex_roles = r.get_complex_roles()
        if len(complex_roles) > 0:
            for cpx_id in complex_roles:
                for role_id in complex_roles[cpx_id]:
                    rxn_roles[r.id].add(normalize_role(roles[role_id]['name']))
                    # print(role_id, normalize_role(roles[role_id]['name']))
    return rxn_roles


def build_reaction(reaction_id, gpr, template, index='0'):
    template_reaction = template.reactions.get_by_id(reaction_id)
    reaction_compartment = 'c'
    metabolites = {}
    for o in template_reaction.templateReactionReagents:
        comp_compound = template_reaction.template.compcompounds.get_by_id(o['templatecompcompound_ref'].split('/')[-1])
        compound = template_reaction.template.compounds.get_by_id(comp_compound['templatecompound_ref'].split('/')[-1])
        comparment = comp_compound.templatecompartment_ref.split('/')[-1] + str(index)
        cpd = Metabolite(comp_compound.id + str(index), compound.formula, compound.name, comp_compound.charge,
                         comparment)
        metabolites[cpd] = o['coefficient']
    lower_bound, upper_bound = get_reaction_constraints_from_direction(template_reaction.direction)
    reaction = Reaction(
        "{}{}".format(template_reaction.id, index),
        "{}_{}{}".format(template_reaction.name, reaction_compartment, index),
        '',
        lower_bound, upper_bound
    )
    gpr_ll = []
    for complex_id in gpr:
        complex_set = set()
        for gene_id in gpr[complex_id]:
            complex_set.add(gene_id)
        gpr_ll.append(list(complex_set))

    reaction.add_metabolites(metabolites)
    reaction.gene_reaction_rule = get_gpr_string(gpr_ll)
    return reaction


class MSBuilder:

    def __init__(self, genome, template):
        """
        for future methods with better customization
        """
        self.genome = genome
        self.template = template

    @staticmethod
    def build_metabolic_model(model_id, genome, template, index='0',
                              allow_all_non_grp_reactions=False, annotate_with_rast=True):

        if annotate_with_rast:
            rast = RastClient()
            res = rast.annotate_genome(genome)

        search_name_to_genes, search_name_to_orginal = aaaa(genome)
        rxn_roles = aux_template(template)  # needs to be fixed to actually reflect template GPR rules

        metabolic_reactions = {}
        genome_search_names = set(search_name_to_genes)
        for rxn_id in rxn_roles:
            sn_set = set(genome_search_names & rxn_roles[rxn_id])
            if len(sn_set) > 0:
                genes = set()
                # print(rxn_id, sn_set)
                for sn in sn_set:
                    if sn in search_name_to_genes:
                        genes |= search_name_to_genes[sn]
                metabolic_reactions[rxn_id] = genes

        # temporary hack until proper complexes from template
        metabolic_reactions_2 = {}
        cpx_random = 0
        for rxn_id in metabolic_reactions:
            metabolic_reactions_2[rxn_id] = {}
            for gene_id in metabolic_reactions[rxn_id]:
                metabolic_reactions_2[rxn_id]['complex' + str(cpx_random)] = {gene_id}
                cpx_random += 1

        reactions = []
        for rxn_id in metabolic_reactions_2:
            reaction = build_reaction(rxn_id, metabolic_reactions_2[rxn_id], template, index)
            reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"
            reaction.annotation["seed.reaction"] = rxn_id
            reactions.append(reaction)
        cobra_model = Model(model_id)
        cobra_model.add_reactions(reactions)

        reactions_no_gpr = []
        reactions_in_model = set(map(lambda x: x.id, cobra_model.reactions))
        metabolites_in_model = set(map(lambda x: x.id, cobra_model.metabolites))

        for rxn in template.reactions:
            if rxn.data['type'] == 'universal' or rxn.data['type'] == 'spontaneous':
                reaction = build_reaction(rxn.id, {}, template, index)
                reaction_metabolite_ids = set(map(lambda x: x.id, set(reaction.metabolites)))
                if (len(metabolites_in_model & reaction_metabolite_ids) > 0 or allow_all_non_grp_reactions) and \
                        reaction.id not in reactions_in_model:
                    reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"
                    reaction.annotation["seed.reaction"] = rxn_id
                    reactions_no_gpr.append(reaction)
        cobra_model.add_reactions(reactions_no_gpr)

        reactions_exchanges = []
        for m in cobra_model.metabolites:
            if m.compartment == 'e0':
                rxn_exchange = Reaction('EX_' + m.id, 'Exchange for ' + m.name, 'exchanges', -1000, 1000)
                rxn_exchange.add_metabolites({m: -1})
                rxn_exchange.annotation[SBO_ANNOTATION] = "SBO:0000627"
                reactions_exchanges.append(rxn_exchange)
        cobra_model.add_reactions(reactions_exchanges)

        if template.name.startswith('CoreModel'):
            bio_rxn1 = build_biomass('bio1', cobra_model, template, core_biomass)
            bio_rxn2 = build_biomass('bio2', cobra_model, template, core_atp)
            cobra_model.add_reactions([bio_rxn1, bio_rxn2])
            cobra_model.objective = 'bio1'
        if template.name.startswith('GramNeg'):
            bio_rxn1 = build_biomass('bio1', cobra_model, template, gramneg)
            cobra_model.add_reactions([bio_rxn1])
            cobra_model.objective = 'bio1'
        if template.name.startswith('GramPos'):
            bio_rxn1 = build_biomass('bio1', cobra_model, template, grampos)
            cobra_model.add_reactions([bio_rxn1])
            cobra_model.objective = 'bio1'

        reactions_sinks = []
        for cpd_id in ['cpd02701_c0', 'cpd11416_c0', 'cpd15302_c0']:
            if cpd_id in cobra_model.metabolites:
                m = cobra_model.metabolites.get_by_id(cpd_id)
                rxn_exchange = Reaction('SK_' + m.id, 'Sink for ' + m.name, 'exchanges', 0, 1000)
                rxn_exchange.add_metabolites({m: -1})
                rxn_exchange.annotation[SBO_ANNOTATION] = "SBO:0000627"
                reactions_sinks.append(rxn_exchange)
        cobra_model.add_reactions(reactions_sinks)

        return cobra_model

    @staticmethod
    def gapfill_model(original_mdl, target_reaction, template, media):
        FBAHelper.set_objective_from_target_reaction(original_mdl, target_reaction)
        model = cobra.io.json.from_json(cobra.io.json.to_json(original_mdl))
        gfp = GapfillingPkg(model)
        gfp.build_package({
            "default_gapfill_templates":[template],
            "gapfill_all_indecies_with_default_templates":1,
            "minimum_obj":0.01,
            "set_objective":1
        })
        kmp = KBaseMediaPkg(model)
        kmp.build_package(media)
        sol=model.optimize()
        gfresults = gfp.compute_gapfilled_solution()
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
