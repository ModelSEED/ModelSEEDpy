import logging
from modelseedpy.core.rast_client import RastClient
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msmodel import get_gpr_string, get_reaction_constraints_from_direction


import re
import copy
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)

### temp stuff ###
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
                reactions_no_gpr.append(reaction)
    cobra_model.add_reactions(reactions_no_gpr)

    return cobra_model