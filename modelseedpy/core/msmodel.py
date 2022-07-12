import logging
import re
from cobra.core import Model
from pyeda.inter import expr  # wheels must be specially downloaded and installed for Windows https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda

logger = logging.getLogger(__name__)


def get_reaction_constraints_from_direction(direction: str) -> (float, float):
    """
    Converts direction symbols ( > or <) to lower and upper bound, any other value is returned as reversible bounds
    :param direction:
    :return:
    """
    if direction == '>':
        return 0, 1000
    elif direction == '<':
        return -1000, 0
    else:
        return -1000, 1000


def get_direction_from_constraints(lower, upper):
    if lower < 0 < upper:
        return '='
    elif upper > 0:
        return '>'
    elif lower < 0:
        return '<'
    logger.error(f'The [{lower}, {upper}] bounds are not amenable with a direction string.')
    return '?'

def get_gpr_string(gpr):
    ors = []
    for ands in gpr:
        a = []
        for g in ands:
            a.append(g)
        ors.append(" and ".join(a))
    gpr_string = "(" + (") or (".join(ors)) + ")"
    if gpr_string == "()":
        return ""

    return gpr_string


def split_compartment_from_index(cmp_str: str):
    """
    Splits index from compartment.
    Example: c0 returns ('c', '0'), cm10 returns ('cm', '10)
    :param cmp_str: a compartment string with or without index
    :return:
    """
    s = re.split(r'(\d+)', cmp_str)
    index_val = None
    cmp_val = None
    if len(s) == 3 and len(s[0]) > 0:
        cmp_val, index_val, empty = s
        if empty > 0:
            cmp_val = None  # set cmp_val to None to raise error if empty element is not empty
    elif len(cmp_str) == 1:
        cmp_val = s[0]
    if cmp_val is None:
        raise ValueError(f"""Bad value {cmp_str} Value of compartment string must start with letter(s)
         and ending (optional) with digits""")
    return cmp_val, index_val


def get_cmp_token(compartments):
    """

    :param compartments:
    :return:
    """
    if len(compartments) == 0:
        logger.warning('The compartments parameter is empty. The "c" parameter is assumed.')
        return 'c'
    if len(compartments) == 1:
        return list(compartments)[0]
    if len(compartments) == 2:
        if set(compartments) == {'e', 'p'}:
            return 'e'
        elif 'b' in compartments and 'e' in compartments:
            return 'b'
        elif 'e' in compartments and 'c' in compartments:
            return 'c'
        elif 'k' in compartments:
            return 'k'
        elif 'c' in compartments:
            return list(filter(lambda x: not x == 'c', compartments))[0]
    return None


def get_set_set(expr_str):   # !!! this currently returns dictionaries, not sets??
    if len(expr_str.strip()) == 0:
        return {}
    expr_str = expr_str.replace(' or ', ' | ')
    expr_str = expr_str.replace(' and ', ' & ')
    dnf = expr(expr_str).to_dnf()
    if len(dnf.inputs) == 1 or dnf.NAME == 'And':
        return {frozenset({str(x) for x in dnf.inputs})}
    else:
        return {frozenset({str(x) for x in o.inputs}) for o in dnf.xs}


class MSModel(Model):

    def __init__(self, id_or_model=None, genome=None, template=None):
        """
        Class representation for a ModelSEED model.
        """
        super().__init__(self, id_or_model)
        if genome:
            self.genome_object = genome
        if template:
            self.template_object = template

    @property
    def template(self):
        return self.template_object

    @template.setter
    def template(self, template):
        self.template_object = template

    @property
    def genome(self):
        return self.genome_object

    @genome.setter
    def genome(self, genome):
        self.genome_object = genome

    def _set_genome_to_model(self, genome):
        # TODO: implement genome assignment checks if features matches genes
        pass
