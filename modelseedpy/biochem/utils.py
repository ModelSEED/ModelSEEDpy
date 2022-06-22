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
        return len(chem_mw.proportions)
    except:
        warn(f'The {formula} formula is invalid')
        return False
    
def molecular_weight(formula):
    chem_mw = ChemMW()
    try:
        return chem_mw.mass(formula)
    except:
        warn(f'The {formula} formula is invalid')

class PeriodicTable:
    
    def __init__(self):
        self.elements_v = self.elements = {}
        for element in periodic_table:
            self.elements_v[element.symbol] = element.protons
            self.elements[element.symbol] = element.name

    def get_element_name(self, e):
        if e in self.elements:
            return self.elements[e]
        return None