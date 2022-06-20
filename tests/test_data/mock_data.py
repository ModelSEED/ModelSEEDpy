def mock_atp_biomass():
    return {
        'cellwall': 0,
        'cofactor': 0,
        'dna': 0,
        'energy': 1,
        'id': 'bio1',
        'lipid': 0,
        'name': 'EnergyProduction',
        'other': 0,
        'protein': 0,
        'rna': 0,
        'templateBiomassComponents': [
            {
                'class': 'energy',
                'coefficient': -1,
                'coefficient_type': 'MULTIPLIER',
                'link_coefficients': [1, -1, -1, -1],
                'linked_compound_refs': [
                    '~/compcompounds/id/cpd00001_c',
                    '~/compcompounds/id/cpd00008_c',
                    '~/compcompounds/id/cpd00009_c',
                    '~/compcompounds/id/cpd00067_c'
                ],
                'templatecompcompound_ref': '~/compcompounds/id/cpd00002_c'
            }
        ],
        'type': 'energy'
    }

def mock_biomass():
    return {
        'cellwall': 0,
        'cofactor': 0,
        'dna': 0,
        'energy': 41.257,
        'id': 'bio2',
        'lipid': 0,
        'name': 'CoreBiomass',
        'other': 1,
        'protein': 0,
        'rna': 0,
        'templateBiomassComponents': [
            {
                'class': 'other',
                'coefficient': -0.0709,
                'coefficient_type': 'EXACT',
                'link_coefficients': [],
                'linked_compound_refs': [],
                'templatecompcompound_ref': '~/compcompounds/id/cpd00072_c'
            },
            {
                'class': 'other',
                'coefficient': -0.8977,
                'coefficient_type': 'EXACT',
                'link_coefficients': [],
                'linked_compound_refs': [],
                'templatecompcompound_ref': '~/compcompounds/id/cpd00101_c'
            },
            {'class': 'other',
              'coefficient': -0.8977,
              'coefficient_type': 'EXACT',
              'link_coefficients': [],
              'linked_compound_refs': [],
              'templatecompcompound_ref': '~/compcompounds/id/cpd00236_c'},
             {'class': 'other',
              'coefficient': -0.129,
              'coefficient_type': 'EXACT',
              'link_coefficients': [],
              'linked_compound_refs': [],
              'templatecompcompound_ref': '~/compcompounds/id/cpd00102_c'},
             {'class': 'other',
              'coefficient': -1.496,
              'coefficient_type': 'EXACT',
              'link_coefficients': [],
              'linked_compound_refs': [],
              'templatecompcompound_ref': '~/compcompounds/id/cpd00169_c'},
             {'class': 'other',
              'coefficient': -0.5191,
              'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00061_c'},
                                                     {'class': 'other',
                                                      'coefficient': -2.8328,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00020_c'},
                                                     {'class': 'other',
                                                      'coefficient': -3.7478,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00022_c'},
                                                     {'class': 'other',
                                                      'coefficient': -1.7867,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00032_c'},
                                                     {'class': 'other',
                                                      'coefficient': -1.0789,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00024_c'},
                                                     {'class': 'other',
                                                      'coefficient': -0.205,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00079_c'},
                                                     {'class': 'other',
                                                      'coefficient': -1.8225,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [-1],
                                                      'linked_compound_refs': ['~/compcompounds/id/cpd00067_c'],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00005_c'},
                                                     {'class': 'other',
                                                      'coefficient': -3.547,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [-1],
                                                      'linked_compound_refs': ['~/compcompounds/id/cpd00067_c'],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00003_c'},
                                                     {'class': 'other',
                                                      'coefficient': 3.7478,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00010_c'},
                                                     {'class': 'other',
                                                      'coefficient': 1,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd11416_c'},
                                                     {'class': 'other',
                                                      'coefficient': 1.8225,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00006_c'},
                                                     {'class': 'other',
                                                      'coefficient': 3.547,
                                                      'coefficient_type': 'EXACT',
                                                      'link_coefficients': [],
                                                      'linked_compound_refs': [],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00004_c'},
                                                     {'class': 'energy',
                                                      'coefficient': -1,
                                                      'coefficient_type': 'MULTIPLIER',
                                                      'link_coefficients': [1, -1, -1, -1],
                                                      'linked_compound_refs': ['~/compcompounds/id/cpd00001_c',
                                                                               '~/compcompounds/id/cpd00008_c',
                                                                               '~/compcompounds/id/cpd00009_c',
                                                                               '~/compcompounds/id/cpd00067_c'],
                                                      'templatecompcompound_ref': '~/compcompounds/id/cpd00002_c'}],
                       'type': 'growth'}


def mock_template():
    template_data = {
        '__VERSION__': 1,
        'id': 'CoreModelTemplateTest',
        'biochemistry_ref': '1/2/3',
        'type': 'GenomeScale',
        'domain': 'Bacteria',
        'biomasses': [mock_biomass(), mock_atp_biomass()],
        'compartments': [
            {
                'aliases': [],
                          'hierarchy': 3,
                          'id': 'c',
                          'index': '0',
                          'name': 'Cytosol',
                          'pH': 7},
            {'aliases': [],
                          'hierarchy': 0,
                          'id': 'e',
                          'index': '1',
                          'name': 'Extracellular',
             }],

    }
    return template_data


def mock_genome_rast():
    pass


def mock_model():
    pass


def remap(model, bigg_to_seed_cpd, bigg_to_seed_rxn, index='0'):
    from cobra.core import Model, Metabolite, Reaction
    from modelseedpy.core.msmodel import get_cmp_token
    model_result = Model(model.id)
    metabolites = {}
    for m in model.metabolites:
        metabolite_id = m.id
        if m.id[:-2] in bigg_to_seed_cpd:
            metabolite_id = f'{bigg_to_seed_cpd[m.id[:-2]]}_{m.compartment}{index}'
        m_result = Metabolite(metabolite_id, m.formula, m.name, m.charge, m.compartment)
        metabolites[m.id] = m_result

    model_result.add_metabolites(metabolites.values())
    reactions = {}
    for r in model.reactions:
        reaction_id = r.id
        if r.id in bigg_to_seed_rxn:
            compartment = get_cmp_token(r.compartments)
            reaction_id = f'{bigg_to_seed_rxn[r.id]}_{compartment}{index}'
        r_result = Reaction(reaction_id, r.name, r.subsystem, r.lower_bound, r.upper_bound)
        s = {}
        for m, value in r.metabolites.items():
            if m.id in metabolites:
                s[model_result.metabolites.get_by_id(metabolites[m.id].id)] = value
            else:
                s[model_result.metabolites.get_by_id(m.id)] = value
        r_result.add_metabolites(s)
        reactions[r.id] = r_result
    model_result.add_reactions(reactions.values())
    model_result.objective = model.objective
    return model_result


def mock_model_ecoli_core(seed=True):
    from cobra.io import load_json_model
    from os import path
    model = load_json_model(path.join(path.dirname(__file__),'e_coli_core.json'))
    if not seed:
        return model
    bigg_to_seed_cpd = {'pyr': 'cpd00020',
     '3pg': 'cpd00169',
     'q8h2': 'cpd15561',
     'pi': 'cpd00009',
     'xu5p__D': 'cpd00198',
     '2pg': 'cpd00482',
     'glu__L': 'cpd00023',
     'glc__D': 'cpd00027',
     'g3p': 'cpd00102',
     'e4p': 'cpd00236',
     '6pgl': 'cpd00911',
     'lac__D': 'cpd00221',
     'acon_C': 'cpd00331',
     'oaa': 'cpd00032',
     'coa': 'cpd00010',
     'ru5p__D': 'cpd00171',
     'atp': 'cpd00002',
     'nadh': 'cpd00004',
     'o2': 'cpd00007',
     'gln__L': 'cpd00053',
     'nadp': 'cpd00006',
     'acald': 'cpd00071',
     'h': 'cpd00067',
     'nh4': 'cpd00013',
     'nadph': 'cpd00005',
     'h2o': 'cpd00001',
     'co2': 'cpd00011',
     'accoa': 'cpd00022',
     'adp': 'cpd00008',
     'dhap': 'cpd00095',
     'cit': 'cpd00137',
     'fum': 'cpd00106',
     'f6p': 'cpd00072',
     'r5p': 'cpd00101',
     'etoh': 'cpd00363',
     '6pgc': 'cpd00284',
     'akg': 'cpd00024',
     'fdp': 'cpd00290',
     'amp': 'cpd00018',
     'for': 'cpd00047',
     'nad': 'cpd00003',
     '13dpg': 'cpd00203',
     'succoa': 'cpd00078',
     'succ': 'cpd00036',
     'glx': 'cpd00040',
     'icit': 'cpd00260',
     'mal__L': 'cpd00130',
     'q8': 'cpd15560',
     'actp': 'cpd00196',
     'fru': 'cpd00082',
     'g6p': 'cpd00079',
     'ac': 'cpd00029',
     's7p': 'cpd00238',
     'pep': 'cpd00061'}
    bigg_to_seed_rxn = {
        'FBP': 'rxn00549', 'PIt2r': 'rxn05312', 'PPCK': 'rxn00247', 'FBA': 'rxn00786', 'MALt2_2': 'rxn08865',
        'TKT1': 'rxn01200', 'FORt': 'rxn08525', 'FRD7': 'rxn09272_copy1', 'PFK': 'rxn00545', 'SUCCt3': 'rxn09270',
        'PTAr': 'rxn00173', 'CYTBD': 'rxn08288', 'ENO': 'rxn00459', 'NADTRHD': 'rxn00083', 'ADK1': 'rxn00097',
        'PPC': 'rxn32421', 'PDH': 'rxn00154', 'CS': 'rxn19316', 'GLNS': 'rxn00187', 'PGI': 'rxn00558',
        'H2Ot': 'rxn05319', 'TPI': 'rxn00747', 'O2t': 'rxn05468', 'PYK': 'rxn35242', 'FRUpts2': 'rxn08535',
        'GND': 'rxn01115', 'ACKr': 'rxn00225',
     'LDH_D': 'rxn00500',
     'D_LACt2': 'rxn08350',
     'MALS': 'rxn20162',
     'TALA': 'rxn01333',
     'FUM': 'rxn30674',
     'PGK': 'rxn01100',
     'NH4t': 'rxn39085',
     'ATPS4r': 'rxn08173',
     'SUCOAS': 'rxn00285',
     'ME2': 'rxn00161',
     'ALCD2x': 'rxn00543',
     'GLUt2r': 'rxn05297',
     'CO2t': 'rxn05467',
     'ACt2r': 'rxn05488',
     'GLUSy': 'rxn00085',
     'ATPM': 'rxn00062',
     'PFL': 'rxn34258',
     'ACALD': 'rxn00171',
     'FORt2': 'rxn05559',
     'SUCDi': 'rxn09272',
     'G6PDH2r': 'rxn00604',
     'ME1': 'rxn00159',
     'NADH16': 'rxn08972',
     'FUMt2_2': 'rxn08542',
     'PGM': 'rxn01106',
     'ICDHyr': 'rxn00198',
     'PGL': 'rxn01476',
     'GAPD': 'rxn00781',
     'THD2': 'rxn09295',
     'TKT2': 'rxn00785',
     'ACALDt': 'rxn08032',
     'GLUN': 'rxn00189',
     'ICL': 'rxn00336',
     'PPS': 'rxn00147',
     'SUCCt2_2': 'rxn09269',
     'MDH': 'rxn00248',
     'ACONTa': 'rxn00974',
     'AKGt2r': 'rxn05493',
     'RPI': 'rxn00777',
     'AKGDH': 'rxn08094',
     'GLCpts': 'rxn05226',
     'GLNabc': 'rxn05155',
     'ETOHt2r': 'rxn08427',
     'GLUDy': 'rxn00184',
     'PYRt2': 'rxn05469',
     'ACONTb': 'rxn32587',
     'RPE': 'rxn01116'}

    return remap(model, bigg_to_seed_cpd, bigg_to_seed_rxn)
