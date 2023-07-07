# -*- coding: utf-8 -*-
from modelseedpy.biochem.seed_object import ModelSEEDObject
from modelseedpy.core.mstemplate import MSTemplateSpecies, MSTemplateMetabolite
from cobra.core import Metabolite
import pandas as pd

_SMILE_ALIAS = "SMILE"
_INCHI_ALIAS = "InChI"
_INCHI_KEY_ALIAS = "InChIKey"


class ModelSEEDCompound2(Metabolite):
    def __init__(
        self,
        cpd_id=None,
        formula=None,
        name="",
        charge=None,
        compartment=None,
        abbr=None,
        names=None,
        mass=None,
        delta_g=None,
        delta_g_error=None,
        is_core=False,
        is_obsolete=False,
        is_cofactor=False,
        is_abstract=False,
        pka=None,
        pkb=None,
        source=None,
        flags=None,
    ):

        super().__init__(cpd_id, formula, name, charge, compartment)
        self.seed_id = cpd_id
        self.abbr = abbr
        self.names = set()
        if names:
            self.names |= set(names)
        self.mass = mass
        self.source = source

        self.is_core = is_core
        self.is_obsolete = is_obsolete
        self.is_cofactor = is_cofactor
        self.is_abstract = is_abstract

        self.delta_g = delta_g
        self.delta_g_error = delta_g_error

        self.linked_compound = None
        self.pka = pka
        self.pkb = pkb
        self.flags = set()
        if flags:
            self.flags |= set(flags)

    def to_template_compartment_compound(self, compartment):
        cpd_id = f"{self.seed_id}_{compartment}"
        # build Template Compound
        metabolite = MSTemplateMetabolite(
            self.seed_id,
            self.formula,
            self.name,
            self.charge,
            self.mass,
            self.delta_g,
            self.delta_g_error,
            self.is_cofactor,
            self.abbr,
        )
        # build Template Compartment Compound
        res = MSTemplateSpecies(cpd_id, self.charge, compartment, metabolite.id)

        # assign Compound to Compartment Compound
        res._template_compound = metabolite
        res.annotation.update(self.annotation)
        return res

    @property
    def smiles(self):
        return (
            None
            if _SMILE_ALIAS not in self.annotation
            else self.annotation[_SMILE_ALIAS]
        )

    @property
    def inchi_key(self):
        return (
            None
            if _INCHI_KEY_ALIAS not in self.annotation
            else self.annotation[_INCHI_KEY_ALIAS]
        )

    @property
    def inchi(self):
        return (
            None
            if _INCHI_ALIAS not in self.annotation
            else self.annotation[_INCHI_ALIAS]
        )


class ModelSEEDCompound(ModelSEEDObject):
    @property
    def formula(self):
        return self.data["formula"]

    @property
    def database(self):
        return "seed.compound"

    @property
    def inchi(self):
        if self.api is not None:
            if (
                self.id in self.api.compound_structures
                and "InChI" in self.api.compound_structures[self.id]
            ):
                return self.api.compound_structures[self.id]["InChI"]
        return None

    @property
    def inchikey(self):
        if self.api is not None:
            if (
                self.id in self.api.compound_structures
                and "InChIKey" in self.api.compound_structures[self.id]
            ):
                return self.api.compound_structures[self.id]["InChIKey"]
        if pd.isna(self.data["inchikey"]):
            return None
        return self.data["inchikey"]

    @property
    def smiles(self):
        if self.api is not None:
            if (
                self.id in self.api.compound_structures
                and "SMILE" in self.api.compound_structures[self.id]
            ):
                return self.api.compound_structures[self.id]["SMILE"]
        if pd.isna(self.data["smiles"]):
            return None
        return self.data["smiles"]

    @property
    def aliases(self):
        if self.api is not None:
            if self.id in self.api.compound_aliases:
                return self.api.compound_aliases[self.id]
        return {}

    @property
    def deltag(self):
        return self.data["deltag"]

    @property
    def is_obsolete(self):
        if "is_obsolete" in self.data:
            is_obsolete = self.data["is_obsolete"]
            if is_obsolete == 0:
                return False
            else:
                return True
        return False
