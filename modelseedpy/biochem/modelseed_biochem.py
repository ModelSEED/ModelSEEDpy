import logging
import pandas as pd
from cobra.core.dictlist import DictList
from modelseedpy.biochem.modelseed_compound import ModelSEEDCompound, ModelSEEDCompound2
from modelseedpy.biochem.modelseed_reaction import ModelSEEDReaction, ModelSEEDReaction2

logger = logging.getLogger(__name__)

ALIAS_CPD_IDENTIFIERS_ORG = {
    'BiGG': 'bigg.metabolite',
    'KEGG': 'kegg.compound',
    'MetaCyc': 'metacyc.compound',
    'metanetx.chemical': 'metanetx.chemical'
}
ALIAS_RXN_IDENTIFIERS_ORG = {
    'KEGG': 'kegg.reaction',
    'MetaCyc': 'metacyc.reaction',
    'metanetx.reaction': 'metanetx.reaction',
    'BiGG': 'bigg.reaction',
    'rhea': 'rhea'
}
ALIAS_MODELS = {
    'iAF1260', 'iAF692', 'iAG612', 'iAO358', 'iAbaylyiv4', 'iGT196', 'iIN800',
    'iIT341', 'iJN746', 'iJR904', 'iMA945', 'iMEO21', 'iMM904', 'iMO1053-PAO1',
    'iMO1056', 'iND750', 'iNJ661', 'iPS189', 'iRR1083', 'iRS1563', 'iRS1597',
    'iSB619', 'iSO783', 'iYO844', 'Maize_C4GEM', 'AlgaGEM', 'AraGEM',
    'JP_Creinhardtii_MSB', 'JP_Creinhardtii_NMeth', 'JM_Creinhardtii', 'TS_Athaliana'
}


def get_low(ids):
    low = None
    ret = None
    for id in ids:
        if low is None or int(id[3:]) < low:
            low = int(id[3:])
            ret = id
    return ret


def make_alias_dict(compound_aliases):
    compound_alias = {}
    for row_id, d in compound_aliases.iterrows():
        alias = d['External ID']
        seed_id = d['ModelSEED ID']
        if alias not in compound_alias:
            compound_alias[alias] = {'seed_id': seed_id}
    return compound_alias


def process_aliases(aliases, alias_mapper, model_ids):
    annotation = {}
    genome_scale_models = {}
    others = {}
    for alias_key in aliases:
        if alias_key in alias_mapper:
            annotation[alias_mapper[alias_key]] = aliases[alias_key]
        elif alias_key in model_ids:
            genome_scale_models[alias_key] = aliases[alias_key]
        else:
            others[alias_key] = aliases[alias_key]
    return annotation, genome_scale_models, others


def load_metabolites_from_df(df: pd.DataFrame, names: dict, aliases: dict, structures: dict):
    compounds = []
    for t in df.itertuples():
        try:
            cpd_id = t[1]
            abbr = t[2]
            name = t[3]
            formula = t[4]
            mass = None if t[5] == 'None' or pd.isna(t[5]) else float(t[5])
            source = None if pd.isna(t[6]) else t[6]
            inchi_key = t[7]
            charge = t[8]
            is_core = True if t[9] else False
            is_obsolete = True if t[10] else False
            linked_compound = t[11]
            is_cofactor = t[12]
            delta_g = t[13]
            delta_g_err = t[14]
            pka = t[15]
            pkb = t[16]
            is_abstract = t[17]
            if pd.isna(is_abstract):
                is_abstract = False
            else:
                is_abstract = True if is_abstract else False
            comprised_of = t[18]
            smiles = t[20]
            flags = t[21]
            inchi = None
            if flags and not pd.isna(flags):
                flags = flags.split('|')
            elif pd.isna(flags):
                flags = set()
            aliases_annotation = {}
            aliases_models = {}
            aliases_others = {}
            if cpd_id in aliases:
                aliases_annotation, aliases_models, aliases_others = process_aliases(aliases[cpd_id],
                                                                                     ALIAS_CPD_IDENTIFIERS_ORG,
                                                                                     ALIAS_MODELS)
            if cpd_id in structures:
                if 'SMILE' in structures[cpd_id]:
                    smiles = structures[cpd_id]['SMILE']
                if 'InChI' in structures[cpd_id]:
                    inchi = structures[cpd_id]['InChI']
                if 'InChIKey' in structures[cpd_id]:
                    inchi_key = structures[cpd_id]['InChIKey']
            inchi_key = None if pd.isna(inchi_key) or len(inchi_key) == 0 else inchi_key
            other_names = set()
            if cpd_id in names:
                other_names = set(names[cpd_id])
            cpd = ModelSEEDCompound2(cpd_id, formula, name, charge, 'z',
                                     abbr, other_names,
                                     mass, delta_g, delta_g_err,
                                     smiles, inchi_key, inchi,
                                     is_core, is_obsolete, is_cofactor, is_abstract,
                                     pka, pkb, source, flags)
            cpd.annotation.update(aliases_annotation)
            cpd.notes['models'] = aliases_models
            cpd.notes['other_aliases'] = aliases_others
            compounds.append(cpd)
        except Exception as e:
            logger.error('failed to read compound at Index: %s. %s', t[0], e)
    return compounds


def load_reactions_from_df(df: pd.DataFrame, database_metabolites: dict, names: dict, aliases: dict, reaction_ecs: dict):
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction
    from modelseedpy.biochem.modelseed_reaction import get_cstoichiometry
    reactions = []
    metabolites_indexed = {}
    for t in df.itertuples():
        try:
            rxn_id = t[1]
            abbr = t[2]
            name = t[3]
            stoichiometry = t[5]
            reversibility = t[9]
            direction = t[10]
            is_abstract = t[11]
            delta_g = t[15]
            delta_g_err = t[16]
            status = t[18]
            is_obsolete = t[19]
            linked_reactions = t[20]
            flags = t[21]
            source = t[22]
            lower_bound, upper_bound = get_reaction_constraints_from_direction(reversibility)
            stoichiometry = None if pd.isna(stoichiometry) else stoichiometry
            is_abstract = False if pd.isna(is_abstract) else True
            is_abstract = True if is_abstract else False

            c_stoichiometry = get_cstoichiometry(stoichiometry)
            metabolites = {}
            for (cpd_id, cmp_token), value in c_stoichiometry.items():
                cpd = database_metabolites[cpd_id]
                cpd_index_id = f'{cpd.id}_{cmp_token}'
                if cpd_index_id not in metabolites_indexed:
                    cpd_token = cpd.copy()
                    cpd_token.id = f'{cpd.id}_{cmp_token}'
                    cpd_token.base_id = cpd.id
                    cpd_token.compartment = cmp_token
                    metabolites_indexed[cpd_index_id] = cpd_token
                metabolites[metabolites_indexed[cpd_index_id]] = value

            if flags and not pd.isna(flags):
                flags = flags.split('|')
            elif pd.isna(flags):
                flags = set()

            other_names = set()
            if rxn_id in names:
                other_names = set(names[rxn_id])
            rxn = ModelSEEDReaction2(rxn_id, name, '', lower_bound, upper_bound, abbr, other_names, delta_g,
                                     delta_g_err, is_obsolete, is_abstract, status, source, flags)

            annotation = {}
            if rxn_id in reaction_ecs and 'Enzyme Class' in reaction_ecs[rxn_id]:
                annotation['ec-code'] = reaction_ecs[rxn_id]['Enzyme Class']
            if rxn_id in aliases:
                aliases_annotation, aliases_models, aliases_others = process_aliases(aliases[rxn_id],
                                                                                     ALIAS_RXN_IDENTIFIERS_ORG,
                                                                                     ALIAS_MODELS)
                annotation.update(aliases_annotation)
                rxn.notes['models'] = aliases_models
                rxn.notes['other_aliases'] = aliases_others

            rxn.annotation.update(annotation)
            rxn.add_metabolites(metabolites)
            reactions.append(rxn)
        except Exception as e:
            logger.error('failed to read reaction at Index: %s. %s', t[0], e)
    return reactions, list(metabolites_indexed.values())


class ModelSEEDDatabase:
    """
    ModelSEED database instance.
    """

    def __init__(self, compounds, reactions, compound_tokens):
        self.compounds = DictList()
        self.compound_tokens = DictList()
        self.reactions = DictList()
        self.compounds += compounds
        self.reactions += reactions
        self.reactions += compound_tokens
        self.inchi_key_lookup = {}
        self.metabolite_reactions = {}

    def compounds_by_alias(self, alias, value):
        pass

    def reactions_by_alias(self, alias, value):
        pass

    def find_compounds_by_inchi_key(self, inchi_key, exact=True):
        pass

    def find_reactions_by_compounds(self, compounds):
        pass

    def add_compound(self, cpd):
        if cpd.inchi_key:
            a, b, p = cpd.inchi_key.split('-')
            if a not in self.inchi_key_lookup:
                self.inchi_key_lookup[a] = {}
            if b not in self.inchi_key_lookup[a]:
                self.inchi_key_lookup[a][b] = set()
            self.inchi_key_lookup[a][b].add((cpd.id, p))

    def add_reaction(self, rxn):
        for m in rxn.metabolites:
            if m.seed_id not in self.metabolite_reactions:
                self.metabolite_reactions[m.seed_id] = set()
            self.metabolite_reactions[m.seed_id].add(rxn.id)


class ModelSEEDBiochem:
    
    def __init__(self, compounds, reactions,
                 compound_aliases=None,
                 reaction_aliases=None,
                 compound_structures=None,
                 reaction_ecs=None):
        if reaction_ecs is None:
            reaction_ecs = {}
        if compound_structures is None:
            compound_structures = {}
        if reaction_aliases is None:
            reaction_aliases = {}
        if compound_aliases is None:
            compound_aliases = {}

        self.compounds = compounds
        self.reactions = reactions
        self.compound_aliases = compound_aliases
        self.reaction_aliases = reaction_aliases
        self.compound_structures = compound_structures
        self.reaction_ecs = reaction_ecs

    def summary(self):
        print("cpds:", len(self.compounds), "rxns:", len(self.reactions))
        print('structures', len(self.compound_structures))
        
    def get_formula(self, seed_id):
        return self.get_attribute('formula', seed_id)
    
    def get_name(self, seed_id):
        return self.get_attribute('name', seed_id)

    def get_attribute(self, attr, seed_id):
        seed = self.get_seed_compound(seed_id)
        if seed is None:
            return None
        if attr in seed.data:
            value = seed.data[attr]
            if not value == 'null':
                return seed.data[attr]
        return None
    
    def get_non_obsolete(self, seed_reaction):
        if seed_reaction.is_obsolete and 'linked_reaction' in seed_reaction.data:
            for id in seed_reaction.data['linked_reaction'].split(';'):
                other = self.get_seed_reaction(id)
                if not other == None:
                    if not other.is_obsolete:
                        return other.id

        return seed_reaction.id

    def get_seed_compound(self, seed_id):
        if seed_id in self.compounds:
            return ModelSEEDCompound(self.compounds[seed_id], self)
        return None
    
    def get_seed_reaction(self, seed_id):
        if seed_id in self.reactions:
            return ModelSEEDReaction(self.reactions[seed_id], self)
        return None
    
    def get_seed_compound_by_alias(self, database, cpd_id):
        o = None
        for o_id in self.compounds:
            aliases_str = self.compounds[o_id]['aliases']
            alias = dict((a.split(':')[0], a.split(':')[1]) for a in aliases_str.split(';'))
            if database in alias and cpd_id in alias[database].split('|'):
                o = self.compounds[o_id]
                break
        return o
    
    def get_seed_reaction_by_alias(self, id):
        o = None
        for o_id in self.reactions:
            if id in self.reactions[o_id]['aliases']:
                o = self.reactions[o_id]
                break
        return o

    def get_mapping_external_to_seed(self, db):
        cobra_to_seed = {}
        # db = 'iJR904'
        for seed_id in self.compound_aliases:
            mapping = self.compound_aliases[seed_id]
            if db in mapping:
                for other_id in mapping[db]:
                    if other_id not in cobra_to_seed:
                        cobra_to_seed[other_id] = set()
                    cobra_to_seed[other_id].add(seed_id)
                    # for other_id in cobra_to_seed:
                    #    if len(cobra_to_seed[other_id]) > 1:
                    #        if len(cobra_to_seed[other_id] & kbase_rxns) == 1:
                    #            cobra_to_seed[other_id] = cobra_to_seed[other_id] & kbase_rxns
                    #        else:
                    1
        return cobra_to_seed


def get_aliases_from_df(df: pd.DataFrame):
    aliases = {}
    for t in df.itertuples():
        seed_id = t[1]
        database_id = t[2]
        database = t[3]
        if seed_id not in aliases:
            aliases[seed_id] = {}
        if database not in aliases[seed_id]:
            aliases[seed_id][database] = set()
        aliases[seed_id][database].add(database_id)
    return aliases


def from_github(commit, database_repo='https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase'):
    reactions_url = database_repo + '/%s/Biochemistry/reactions.tsv' % commit
    compounds_url = database_repo + '/%s/Biochemistry/compounds.tsv' % commit
    reactions_aliases_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt' % commit
    compounds_aliases_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt' % commit
    compounds_structures_url = database_repo + '/%s/Biochemistry/Structures/Unique_ModelSEED_Structures.txt' % commit
    reactions_ec_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt' % commit
    reactions_names_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Names.txt' % commit
    compounds_names_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Compound_Names.txt' % commit

    return load_database(compounds_url, reactions_url,
                         compounds_names_url, compounds_aliases_url, compounds_structures_url,
                         reactions_names_url, reactions_aliases_url, reactions_ec_url)


def from_local2(path):
    reactions_url = path + '/Biochemistry/reactions.tsv'
    compounds_url = path + '/Biochemistry/compounds.tsv'
    reactions_aliases_url = path + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt'
    compounds_aliases_url = path + '/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt'
    reactions_names_url = path + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Names.txt'
    compounds_names_url = path + '/Biochemistry/Aliases/Unique_ModelSEED_Compound_Names.txt'
    compounds_structures_url = path + '/Biochemistry/Structures/Unique_ModelSEED_Structures.txt'
    reactions_ec_url = path + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt'
    return load_database(compounds_url, reactions_url,
                         compounds_names_url, compounds_aliases_url, compounds_structures_url,
                         reactions_names_url, reactions_aliases_url, reactions_ec_url)


def load_database(compounds_url, reactions_url,
                  compounds_names_url, compounds_aliases_url, compounds_structures_url,
                  reactions_names_url, reactions_aliases_url, reactions_ec_url):
    logger.info("load: %s", compounds_structures_url)
    compound_structures = get_structures_from_df(pd.read_csv(compounds_structures_url, sep='\t'))
    logger.info("load: %s", compounds_names_url)
    compound_names = get_names_from_df(pd.read_csv(compounds_names_url, sep='\t'))
    logger.info("load: %s", compounds_aliases_url)
    compound_aliases = get_aliases_from_df(pd.read_csv(compounds_aliases_url, sep='\t'))

    logger.info("load: %s", compounds_url)
    df_compounds = pd.read_csv(compounds_url, sep='\t', low_memory=False)  # Columns (10,14,15) have mixed types.
    compounds = load_metabolites_from_df(df_compounds, compound_names, compound_aliases, compound_structures)

    logger.info("load: %s", reactions_names_url)
    reaction_names = get_names_from_df(pd.read_csv(reactions_names_url, sep='\t'))
    logger.info("load: %s", reactions_aliases_url)
    reaction_aliases = get_aliases_from_df(pd.read_csv(reactions_aliases_url, sep='\t'))
    logger.info("load: %s", reactions_ec_url)
    reaction_ecs = get_aliases_from_df(pd.read_csv(reactions_ec_url, sep='\t'))

    df_reactions = pd.read_csv(reactions_url, sep='\t', low_memory=False)
    reactions, metabolites_indexed = load_reactions_from_df(
        df_reactions, dict(map(lambda x: (x.id, x), compounds)), reaction_names,
        reaction_aliases, reaction_ecs)

    database = ModelSEEDDatabase(compounds, reactions, metabolites_indexed)
    return database


def from_local(path):
    database_repo = path
    reactions_url = database_repo + '/Biochemistry/reactions.tsv'
    compounds_url = database_repo + '/Biochemistry/compounds.tsv'
    reactions_aliases_url = database_repo + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt' 
    compounds_aliases_url = database_repo + '/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt' 
    reactions_names_url = database_repo + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Names.txt'
    compounds_names_url = database_repo + '/Biochemistry/Aliases/Unique_ModelSEED_Compound_Names.txt'
    compounds_structures_url = database_repo + '/Biochemistry/Structures/Unique_ModelSEED_Structures.txt'
    reactions_ec_url = database_repo + '/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt'
    
    compounds = {}
    reactions = {}

    logger.info("load: %s", reactions_url)
    seed_reactions = pd.read_csv(reactions_url, sep='\t', low_memory=False)
    logger.info("load: %s", compounds_url)
    seed_compounds = pd.read_csv(compounds_url, sep='\t', low_memory=False)  # Columns (10,14,15) have mixed types.
    
    for row_id, d in seed_reactions.iterrows():
        seed_reaction = build_reaction(d)
        reactions[seed_reaction['id']] = seed_reaction

    for row_id, d in seed_compounds.iterrows():
        seed_compound = build_compound(d)
        compounds[seed_compound['id']] = seed_compound

    logger.info("load: %s", compounds_structures_url)
    df_structures = pd.read_csv(compounds_structures_url, sep='\t')
    compound_structures = get_structures(df_structures)

    logger.info("load: %s", compounds_aliases_url)
    df_compound_aliases = pd.read_csv(compounds_aliases_url, sep='\t')
    logger.info("load: %s", reactions_aliases_url)
    df_reaction_aliases = pd.read_csv(reactions_aliases_url, sep='\t')
    logger.info("load: %s", reactions_ec_url)
    df_reaction_ecs = pd.read_csv(reactions_ec_url, sep='\t')
    
    compound_aliases = get_aliases_from_df(df_compound_aliases)
    reaction_aliases = get_aliases_from_df(df_reaction_aliases)
    reaction_ecs = get_aliases_from_df(df_reaction_ecs)
    
    modelseed = ModelSEEDBiochem(compounds, reactions, compound_aliases, reaction_aliases, compound_structures, reaction_ecs)
    
    return modelseed


def get_names_from_df(df):
    names = {}
    for t in df.itertuples():
        seed_id = t[1]
        v = t[2]
        if seed_id not in names:
            names[seed_id] = set()
        names[seed_id].add(v)
    return names


def get_structures_from_df(df: pd.DataFrame):
    compound_structures = {}

    for t in df.itertuples():
        seed_id = t[1]
        structure_type = t[2]
        value = t[6]
        if seed_id not in compound_structures:
            compound_structures[seed_id] = {}
        if structure_type in compound_structures[seed_id]:
            logger.warning('warning duplicate structure: %s (%s)', structure_type, seed_id)
        compound_structures[seed_id][structure_type] = value

    return compound_structures


def get_structures(seed_structures):
    compound_structures = {}
    
    for row_id, d in seed_structures.iterrows():
        seed_id = d['ID']
        t = d['Type']
        value = d['Structure']
        if seed_id not in compound_structures:
            compound_structures[seed_id] = {}
        if t in compound_structures[seed_id]:
            print('warning duplicate structure:', t, seed_id)
        compound_structures[seed_id][t] = value
    
    return compound_structures


def build_compound(d):
    seed_compound = {
        'id': d['id'],
        'abbreviation': d['abbreviation'],
        'name': d['name'],
        'formula': d['formula'],
        'mass': d['mass'],
        'source': d['source'],
        'inchikey': d['inchikey'],
        'charge': d['charge'],
        'is_core': d['is_core'],
        'is_obsolete': d['is_obsolete'],
        'linked_compound': d['linked_compound'],
        'is_cofactor': d['is_cofactor'],
        'deltag': d['deltag'],
        'deltagerr': d['deltagerr'],
        'pka': d['pka'],
        'pkb': d['pkb'],
        'abstract_compound': d['abstract_compound'],
        'comprised_of': d['comprised_of'],
        'aliases': d['aliases'],
        'smiles': d['smiles']
    }
    return seed_compound


def build_reaction(d):
    seed_reaction = {
            'id': d['id'],
            'abbreviation': d['abbreviation'],
            'name': d['name'],
            'code': d['code'],
            'stoichiometry': d['stoichiometry'],
            'is_transport': d['is_transport'],
            'equation': d['equation'],
            'definition': d['definition'],
            'reversibility': d['reversibility'],
            'direction': d['direction'],
            'abstract_reaction': d['abstract_reaction'],
            'pathways': d['pathways'],
            'aliases': d['aliases'],
            'ec_numbers': d['ec_numbers'],
            'deltag': d['deltag'],
            'deltagerr': d['deltagerr'],
            'compound_ids': d['compound_ids'],
            'status': d['status'],
            'is_obsolete': d['is_obsolete'],
            'linked_reaction': d['linked_reaction'],
            'notes': d['notes'],
            'source': d['source']
    }
    return seed_reaction
