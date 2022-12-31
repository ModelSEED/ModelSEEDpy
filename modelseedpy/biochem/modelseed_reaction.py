
# -*- coding: utf-8 -*-
import math
from modelseedpy.biochem.seed_object import ModelSEEDObject
from cobra.core import Reaction
from modelseedpy.core.mstemplate import MSTemplateReaction


def to_str2(rxn, cmp_replace=None, cpd_replace={}):
    direction = rxn.data["direction"]
    op = "<?>"
    if direction == "=":
        op = "<=>"
    elif direction == ">":
        op = "-->"
    elif direction == "<":
        op = "<--"
    else:
        op = "<?>"
    cstoichiometry = rxn.cstoichiometry
    l = []
    r = []
    for o in cstoichiometry:
        if cstoichiometry[o] > 0:
            r.append((o[0], o[1], math.fabs(cstoichiometry[o])))
        else:
            l.append((o[0], o[1], math.fabs(cstoichiometry[o])))
    l_str = print_stoich_block(l, cmp_replace, cpd_replace)
    r_str = print_stoich_block(r, cmp_replace, cpd_replace)
    return "{} {} {}".format(l_str, op, r_str)


def print_stoich_block(b, cmp_replace=None, cpd_replace={}):
    if cmp_replace is None:
        cmp_replace = {}
    l_text = []
    for o in b:
        # print(o)
        v = o[2]
        v_text = ""
        if not 1 == o[2]:
            v_text = str(o[2]) + " "
        cpd_text = cpd_replace[o[0]] if o[0] in cpd_replace else o[0]
        cmp_text = cmp_replace[o[1]] if o[1] in cmp_replace else o[1]
        v_text += "{} [{}]".format(cpd_text, cmp_text)
        l_text.append(v_text)

    # print(l_text, r_text)
    return " + ".join(l_text)


def get_lhs_rhs(seed_reaction, eval_func):
    result = {}
    result_no_zero = {}
    stoichiometry = seed_reaction["stoichiometry"]
    if type(stoichiometry) == str:
        for reagent_str in stoichiometry.split(";"):
            reagent_data = reagent_str.split(":")
            cpd_id = reagent_data[1]
            value = float(reagent_data[0])
            if eval_func(value):
                if not cpd_id in result:
                    result[cpd_id] = 0
                result[reagent_data[1]] += value
            elif value == 0:
                raise Exception("zero value stoich: " + stoichiometry)

    for cpd_id in result:
        if result[cpd_id] == 0:
            print(f'{seed_reaction["id"]} contains compound bein transported {cpd_id}')
        else:
            result_no_zero[cpd_id] = result[cpd_id]
    return result_no_zero


def get_stoichiometry(seed_reaction):
    result = {}
    result_no_zero = {}
    stoichiometry = seed_reaction["stoichiometry"]
    if type(stoichiometry) == str:
        for reagent_str in stoichiometry.split(";"):
            reagent_data = reagent_str.split(":")
            cpd_id = reagent_data[1]
            if not cpd_id in result:
                result[cpd_id] = 0
            result[reagent_data[1]] += float(reagent_data[0])

    for cpd_id in result:
        if result[cpd_id] == 0:
            print(f'{seed_reaction["id"]} contains compound bein transported {cpd_id}')
        else:
            result_no_zero[cpd_id] = result[cpd_id]
    return result_no_zero


def get_cstoichiometry(stoichiometry):
    if stoichiometry is None or len(stoichiometry) == 0:
        return {}
    result = {}
    result_no_zero = {}
    # print(stoichiometry)
    for reagent_str in stoichiometry.split(";"):
        reagent_data = reagent_str.split(":")
        value = reagent_data[0]
        cpd_id = reagent_data[1]
        cmp_index = reagent_data[2]
        if not (cpd_id, cmp_index) in result:
            result[(cpd_id, cmp_index)] = 0
        result[(cpd_id, cmp_index)] += float(value)

    for p in result:
        result_no_zero[p] = result[p]

    return result_no_zero


class ModelSEEDReaction2(Reaction):
    """
    Reaction instance from ModelSEED database.
    """

    def __init__(
        self,
        rxn_id,
        name="",
        subsystem="",
        lower_bound=0.0,
        upper_bound=None,
        abbr=None,
        names=None,
        delta_g=None,
        delta_g_error=None,
        is_obsolete=False,
        is_abstract=False,
        status=None,
        source=None,
        flags=None,
    ):

        super().__init__(rxn_id, name, subsystem, lower_bound, upper_bound)
        self.abbr = abbr
        self.names = set()
        if names:
            self.names |= set(names)

        self.source = source
        self.status = status

        self.is_obsolete = is_obsolete
        if self.is_obsolete:
            self.is_obsolete = True
        else:
            self.is_obsolete = False
        self.is_abstract = is_abstract

        self.delta_g = float(delta_g) if delta_g else None
        self.delta_g_error = float(delta_g_error) if delta_g_error else None

        # removing symbolic high values representing null/none
        if self.delta_g and self.delta_g > 10000:
            self.delta_g = None
        if self.delta_g_error and self.delta_g_error > 10000:
            self.delta_g_error = None

        self.flags = set()
        if flags:
            self.flags |= set(flags)

    @property
    def compound_ids(self):
        return None

    def to_template_reaction(self, compartment_setup=None):
        if compartment_setup is None:
            raise ValueError("invalid compartment setup")
        from modelseedpy.core.msmodel import get_cmp_token

        reaction_compartment = get_cmp_token(compartment_setup.values())
        rxn_id = f"{self.id}_{reaction_compartment}"
        name = f"{self.name}"
        metabolites = {}
        for m, v in self.metabolites.items():
            if m.compartment not in compartment_setup:
                raise ValueError(
                    f"invalid compartment setup missing key [{m.compartment}] in {compartment_setup}"
                )
            cpd = m.to_template_compartment_compound(compartment_setup[m.compartment])
            metabolites[cpd] = v
            # print(m.id, m.compartment, cpd, v)

        # if len(str(index)) > 0:
        #    name = f'{self.name} [{compartment}]'
        reaction = MSTemplateReaction(
            rxn_id, self.id, name, self.subsystem, self.lower_bound, self.upper_bound
        )
        reaction.add_metabolites(metabolites)
        reaction.annotation.update(self.annotation)
        return reaction

    @property
    def is_transport(self):
        pass

    @property
    def is_translocation(self):
        return len(self.compartments) > 1

    @property
    def ec_numbers(self):
        return self.annotation["ec-code"] if "ec-code" in self.annotation else []

    @property
    def aliases(self):
        return self.annotation


class ModelSEEDReaction(ModelSEEDObject):
    """
    def __init__(self, data, api=None):
        self.data = data
        self.api = api
        self.compounds = {}
    """

    @property
    def stoichiometry(self):
        return get_stoichiometry(self.data)

    @property
    def cstoichiometry(self):
        return get_cstoichiometry(self.data["stoichiometry"])

    @property
    def ec_numbers(self):
        if self.api is not None:
            if self.id in self.api.reaction_ecs:
                return self.api.reaction_ecs[self.id]
        return {}

    @property
    def database(self):
        return "seed.reaction"

    @property
    def lhs(self):
        return get_lhs_rhs(self.data, lambda x: x < 0)

    @property
    def rhs(self):
        return get_lhs_rhs(self.data, lambda x: x > 0)

    @property
    def is_transport(self):
        if "is_transport" in self.data:
            is_transport = self.data["is_transport"]
            if is_transport == 1:
                return True

        if not len(self.lhs.keys() & self.rhs.keys()) == 0:
            return True

        return False

    @property
    def aliases(self):
        if self.api is not None:
            if self.id in self.api.reaction_aliases:
                return self.api.reaction_aliases[self.id]
        return {}

    @property
    def is_obsolete(self):
        if "is_obsolete" in self.data:
            is_obsolete = self.data["is_obsolete"]
            if is_obsolete == 0:
                return False
            else:
                return True
        return False

    def build_reaction_string(
        self, use_metabolite_names=False, use_compartment_names=None
    ):
        cpd_name_replace = {}
        if use_metabolite_names:
            if self.api:
                for cpd_id in set(map(lambda x: x[0], self.cstoichiometry)):
                    cpd = self.api.get_seed_compound(cpd_id)
                    cpd_name_replace[cpd_id] = cpd.name
            else:
                return self.data["definition"]
        return to_str2(self, use_compartment_names, cpd_name_replace)

    def __str__(self):
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string()
        )
