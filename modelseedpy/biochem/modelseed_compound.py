from modelseedpy.biochem.seed_object import ModelSEEDObject
import pandas as pd


class ModelSEEDCompound(ModelSEEDObject):

    @property
    def formula(self):
        return self.data['formula']

    @property
    def database(self):
        return "seed.compound"

    @property
    def inchi(self):
        if self.api is not None:
            if self.id in self.api.compound_structures and 'InChI' in self.api.compound_structures[self.id]:
                return self.api.compound_structures[self.id]['InChI']
        return None

    @property
    def inchikey(self):
        if self.api is not None:
            if self.id in self.api.compound_structures and 'InChIKey' in self.api.compound_structures[self.id]:
                return self.api.compound_structures[self.id]['InChIKey']
        if pd.isna(self.data['inchikey']):
            return None
        return self.data['inchikey']

    @property
    def smiles(self):
        if self.api is not None:
            if self.id in self.api.compound_structures and 'SMILE' in self.api.compound_structures[self.id]:
                return self.api.compound_structures[self.id]['SMILE']
        if pd.isna(self.data['smiles']):
            return None
        return self.data['smiles']

    @property
    def aliases(self):
        if self.api is not None:
            if self.id in self.api.compound_aliases:
                return self.api.compound_aliases[self.id]
        return {}

    @property
    def deltag(self):
        return self.data['deltag']

    @property
    def is_obsolete(self):
        if 'is_obsolete' in self.data:
            is_obsolete = self.data['is_obsolete']
            if is_obsolete == 0:
                return False
            else:
                return True
        return False
