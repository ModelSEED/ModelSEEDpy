import math
import copy
import logging
logger = logging.getLogger(__name__)


cobra_metabolite = {
    'id': '13dpg_c',
    'name': 'compound',
    'compartment': 'c',
    'charge': 0,
    'formula': '',
    'annotation': {}
}

cobra_reaction = {
    'id': 'ACALD',
    'name': 'Acetaldehyde dehydrogenase (acetylating)',
    'metabolites':
        {
            'acald_c': -1.0,
            'accoa_c': 1.0,
            'coa_c': -1.0,
            'h_c': 1.0,
            'nad_c': -1.0,
            'nadh_c': 1.0
        },
    'lower_bound': -1000.0,
    'upper_bound': 1000.0,
    'gene_reaction_rule': '',
    'annotation': {}
}


def modelseed_to_cobra_reaction(rxn, compartment=None, rxn_cmp_suffix='c',
                                lb=None, ub=None, rev='direction',
                                cs_override=None):
    if compartment is None:
        compartment = {'0': None}
    crxn = copy.deepcopy(cobra_reaction)
    if rxn_cmp_suffix is not None and not len(rxn_cmp_suffix.strip()) == 0:
        crxn['id'] = rxn.id + '_' + rxn_cmp_suffix
    else:
        crxn['id'] = rxn.id
    crxn['name'] = "{} [{}]".format(rxn.id, rxn_cmp_suffix)

    metabolites = {}
    cpd_cmp = {}
    # cs = get_cstoichiometry_plant(rxn.data)
    cs = rxn.cstoichiometry
    if not cs_override == None:
        logger.debug("cs_override %s -> %s", cs, cs_override)
        cs = cs_override

    logger.debug("%s: %s", rxn.id, cs)
    logger.debug("%s", compartment)
    for p in cs:
        value = cs[p]
        met_id = p[0]
        if not p[1] in compartment:
            print(rxn.data)
        cmp = compartment[p[1]]
        if cmp is not None and not len(cmp.strip()) == 0:
            met_id += '_' + cmp.strip()
        else:
            cmp = 'z'
        if met_id not in metabolites:
            metabolites[met_id] = value
            if not p[0] in cpd_cmp:
                cpd_cmp[p[0]] = set()
            cpd_cmp[p[0]].add(cmp)
        else:
            print('!', rxn.id)

    rev = rxn.data[rev]  # direction
    if lb is None and ub is None:

        lb = -1000
        ub = 1000
        # print('reversibility', rxn.data['reversibility'], 'direction', rxn.data['direction'])
        if rev == '>':
            lb = 0
        elif rev == '<':
            ub = 0
    crxn['lower_bound'] = lb
    crxn['upper_bound'] = ub
    crxn['metabolites'] = metabolites

    logger.debug("%s: %s", rxn.id, metabolites)

    return crxn, cpd_cmp
