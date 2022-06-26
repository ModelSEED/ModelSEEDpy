from chemicals import periodic_table
from warnings import warn
from chemw import ChemMW


def get_mapping(df, valid_databases, pd):
    mapping = {}
    for id, row in df.iterrows():
        mapping[id] = {}
        #print(dir(row))
        for db in df.keys():
            #print(db)
            if db in valid_databases:
                value = row[db]
                if not pd.isna(value):
                    mapping[id][db] = value
    return mapping


def atom_count(formula):
    chem_mw = ChemMW()
    try:
        chem_mw.mass(formula)
        return chem_mw.proportions
    except ValueError:
        warn(f'The {formula} formula is invalid')
        return None


def is_valid_formula(f, pt):
    warn('The "utils.is_valid_formula" is deprecated in favor of "utils.molecular_weight(formula)"')
    return molecular_weight(f)


def molecular_weight(formula):
    chem_mw = ChemMW()
    try:
        return chem_mw.mass(formula)
    except ValueError:
        warn(f'The {formula} formula is invalid')
        return None


class PeriodicTable:
    
    def __init__(self):
        self.elements_v = {ele.symbol:ele.protons for ele in periodic_table} 
        self.elements = {ele.symbol:ele.name for ele in periodic_table}

    def get_element_name(self, e):
        if e in self.elements:
            return self.elements[e]
        return None
