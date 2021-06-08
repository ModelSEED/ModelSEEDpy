import math
from modelseedpy.biochem.seed_object import ModelSEEDObject


def to_str2(rxn, cmp_replace=None, cpd_replace={}):
    direction = rxn.data['direction']
    op =  '<?>'
    if direction == '=':
        op = '<=>'
    elif direction == '>':
        op = '-->'
    elif direction == '<':
        op = '<--'
    else:
        op =  '<?>'
    cstoichiometry = rxn.cstoichiometry
    l = []
    r = []
    for o in cstoichiometry:
        if cstoichiometry[o] > 0:
            r.append((o[0],o[1], math.fabs(cstoichiometry[o])))
        else:
            l.append((o[0],o[1], math.fabs(cstoichiometry[o])))
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
    return ' + '.join(l_text)


def get_lhs_rhs(seed_reaction, eval_func):
    result = {}
    result_no_zero = {}
    stoichiometry = seed_reaction['stoichiometry']
    if type(stoichiometry) == str:
        for reagent_str in stoichiometry.split(';'):
            reagent_data = reagent_str.split(':')
            cpd_id = reagent_data[1]
            value = float(reagent_data[0])
            if eval_func(value):
                if not cpd_id in result:
                    result[cpd_id] = 0
                result[reagent_data[1]] += value
            elif value == 0:
                raise Exception('zero value stoich: ' + stoichiometry)

    for cpd_id in result:
        if result[cpd_id] == 0:
            print(f'{seed_reaction["id"]} contains compound bein transported {cpd_id}')
        else:
            result_no_zero[cpd_id] = result[cpd_id]
    return result_no_zero


def get_stoichiometry(seed_reaction):
    result = {}
    result_no_zero = {}
    stoichiometry = seed_reaction['stoichiometry']
    if type(stoichiometry) == str:
        for reagent_str in stoichiometry.split(';'):
            reagent_data = reagent_str.split(':')
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


def get_cstoichiometry(seed_reaction):
    result = {}
    result_no_zero = {}
    stoichiometry = seed_reaction['stoichiometry']
    # print(stoichiometry)
    if type(stoichiometry) == str:
        for reagent_str in stoichiometry.split(';'):
            reagent_data = reagent_str.split(':')
            value = reagent_data[0]
            cpd_id = reagent_data[1]
            cmp_index = reagent_data[2]
            if not (cpd_id, cmp_index) in result:
                result[(cpd_id, cmp_index)] = 0
            result[(cpd_id, cmp_index)] += float(value)

    for p in result:
        if result[p] == 0:
            print(f'{seed_reaction["id"]} contains compound bein transported {p}')
        else:
            result_no_zero[p] = result[p]
    return result_no_zero


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
        return get_cstoichiometry(self.data)

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
        if 'is_transport' in self.data:
            is_transport = self.data['is_transport']
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
        if 'is_obsolete' in self.data:
            is_obsolete = self.data['is_obsolete']
            if is_obsolete == 0:
                return False
            else:
                return True
        return False

    def build_reaction_string(self, use_metabolite_names=False, use_compartment_names=None):
        cpd_name_replace = {}
        if use_metabolite_names:
            if self.api:
                for cpd_id in set(map(lambda x: x[0], self.cstoichiometry)):
                    cpd = self.api.get_seed_compound(cpd_id)
                    cpd_name_replace[cpd_id] = cpd.name
            else:
                return self.data['definition']
        return to_str2(self, use_compartment_names, cpd_name_replace)

    def __str__(self):
        return "{id}: {stoichiometry}".format(
            id=self.id, stoichiometry=self.build_reaction_string())
