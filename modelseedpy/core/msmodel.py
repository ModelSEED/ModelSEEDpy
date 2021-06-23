from cobra.core import Model


def get_reaction_constraints_from_direction(direction):
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
