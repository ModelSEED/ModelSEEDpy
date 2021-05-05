import logging
import pandas as pd
from modelseedpy.biochem.modelseed_compound import ModelSEEDCompound
from modelseedpy.biochem.modelseed_reaction import ModelSEEDReaction

logger = logging.getLogger(__name__)


def get_low(ids):
    low = None
    ret = None
    for id in ids:
        if low == None or int(id[3:]) < low:
            low = int(id[3:])
            ret = id
    return ret


def make_alias_dict(compound_aliases):
    compound_alias = {}
    for row_id, d in compound_aliases.iterrows():
        #print(row_id, d)
        alias = d['External ID']
        seed_id = d['ModelSEED ID']
        #print()
        if not alias in compound_alias:
            compound_alias[alias] = {'seed_id' : seed_id}
    return compound_alias


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


def get_aliases_from_df(df):
    aliases = {}
    for t in df.itertuples():
        seed_id = t[1]
        database_id = t[2]
        database = t[3]
        if not seed_id in aliases:
            aliases[seed_id] = {}
        if not database in aliases[seed_id]:
            aliases[seed_id][database] = set()
        aliases[seed_id][database].add(database_id)
    return aliases


def from_github(commit):
    print(commit)
    #https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/566a42f8c04a0dab4536b59377074ad23bcc08be/Biochemistry/reactions.tsv
    database_repo = 'https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase'
    reactions_url = database_repo + '/%s/Biochemistry/reactions.tsv' % commit
    compounds_url = database_repo + '/%s/Biochemistry/compounds.tsv' % commit
    reactions_aliases_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt' % commit
    compounds_aliases_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt' % commit
    compounds_structures_url = database_repo + '/%s/Biochemistry/Structures/Unique_ModelSEED_Structures.txt' % commit
    reactions_ec_url = database_repo + '/%s/Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt' % commit
    compounds = {}
    reactions = {}

    print('load:', reactions_url)
    seed_reactions = pd.read_csv(reactions_url, sep='\t', low_memory=False)
    
    for row_id, d in seed_reactions.iterrows():
        seed_reaction = build_reaction(d)
        reactions[seed_reaction['id']] = seed_reaction
        
    print('load:', compounds_url)
    seed_compounds = pd.read_csv(compounds_url, sep='\t', low_memory=False) #Columns (10,14,15) have mixed types.
    
    for row_id, d in seed_compounds.iterrows():
        seed_compound = build_compound(d)
        compounds[seed_compound['id']] = seed_compound
        
        
    compound_aliases = pd.read_csv(compounds_aliases_url, sep='\t')
    reaction_aliases = pd.read_csv(reactions_aliases_url, sep='\t')
    
    alias_to_database_to_cpd = get_aliases_from_df(compound_aliases)
    alias_to_database_to_rxn = get_aliases_from_df(reaction_aliases)
    print('load:', reactions_ec_url)
    df_reaction_ecs = pd.read_csv(reactions_ec_url, sep='\t')
    print('load:', compounds_structures_url)
    df_structures = pd.read_csv(compounds_structures_url, sep='\t')

    compound_structures = get_structures(df_structures)
    reaction_ecs = get_aliases_from_df(df_reaction_ecs)
    
    modelseed = ModelSEEDBiochem(compounds, reactions, alias_to_database_to_cpd, alias_to_database_to_rxn, compound_structures, reaction_ecs)
    return modelseed


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
    seed_compounds = pd.read_csv(compounds_url, sep='\t', low_memory=False) #Columns (10,14,15) have mixed types.
    
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


def get_structures(seed_structures):
    compound_structures = {}
    
    for row_id, d in seed_structures.iterrows():
        seed_id = d['ID']
        t = d['Type']
        value = d['Structure']
        if not seed_id in compound_structures:
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
