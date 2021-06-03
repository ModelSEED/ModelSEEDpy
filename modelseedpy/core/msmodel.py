def get_reaction_constraints_from_direction(direction):
    if direction == '>':
        return 0, 1000
    elif direction == '<':
        return -1000, 0
    else:
        return -1000, 1000


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
