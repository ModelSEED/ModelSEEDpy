# -*- coding: utf-8 -*-
from cobrakbase.core.kbasegenomesgenome import normalize_role
import pandas as pd


def is_transport_pts(rxn):
    pep = False
    pyr = False
    other = False
    cstoichiometry = rxn.cstoichiometry
    parts = {}
    for o in cstoichiometry:
        if not o[1] in parts:
            parts[o[1]] = {}
        parts[o[1]][o[0]] = cstoichiometry[o]
        if o[0] == "cpd00061":
            pep = True
        elif o[0] == "cpd00020":
            pyr = True
        else:
            other = True
    if not len(parts) == 2:
        return (False, None, None, None, None)
    # print(rxn.data['definition'])
    inside = None
    outside = None
    for o in parts:
        compounds = set(parts[o])
        # print(compounds)
        if "cpd00061" in compounds and "cpd00020" in compounds:
            if inside == None:
                inside = o
            else:
                print("multiple match")
        else:
            outside = o
    # print(inside, outside)
    if inside == None or outside == None:
        return (False, None, None, None, None)

    outside_compounds = set(filter(lambda x: not x == "cpd00067", set(parts[outside])))
    transport_compound = set(parts[outside]).pop()
    t = "import"  # if parts[inside][transport_compound] < 0 else 'import'
    # print(transport_compound)
    if len(outside_compounds) == 1:
        return True, inside, outside, t, transport_compound
    else:
        return (False, None, None, None, None)


def is_transport_abc(rxn):
    atp = False
    adp = False
    h2o = False
    pi = False
    cstoichiometry = rxn.cstoichiometry
    parts = {}

    for o in cstoichiometry:
        if not o[1] in parts:
            parts[o[1]] = {}
        parts[o[1]][o[0]] = cstoichiometry[o]
        if o[0] == "cpd00001":
            h2o = True
        elif o[0] == "cpd00002":
            atp = True
        elif o[0] == "cpd00008":
            adp = True
        elif o[0] == "cpd00009":
            pi = True
        else:
            pass
    # print(parts)
    if not len(parts) == 2:
        return (False, None, None, None, None)
    inside = None
    outside = None
    for o in parts:
        compounds = set(parts[o])
        # print(compounds)
        if (
            "cpd00002" in compounds
            and "cpd00008" in compounds
            and "cpd00009" in compounds
        ):
            if inside == None:
                inside = o
            else:
                print("multiple match")
        else:
            outside = o
    # print(inside, outside)
    if inside == None or outside == None:
        return (False, None, None, None, None)
    # print(parts)
    outside_compounds = set(filter(lambda x: not x == "cpd00067", set(parts[outside])))
    transport_compound = set(parts[outside]).pop()

    # print(transport_compound, parts[outside], parts[inside])
    if (
        transport_compound == None
        or not transport_compound in parts[outside]
        or not transport_compound in parts[inside]
    ):
        return (False, None, None, None, None)
    t = "export" if parts[inside][transport_compound] < 0 else "import"
    if len(outside_compounds) == 1:
        return atp and adp and h2o and pi, inside, outside, t, transport_compound
    else:
        return (False, None, None, None, None)


class TransportersPipeline:
    def __init__(self, modelseed_database, annotation_api):
        self.db = modelseed_database
        self.annotation_api = annotation_api
        self.t_permeases = {}
        self.t_pts = {}
        self.t_abc = {}
        self.t_others = {}
        self.t_abc, self.t_pts, self.t_permeases, self.t_others = self.get_t_class()

    def get_cpd_alias(self):
        cpd_alias = {}
        for cpd_id in self.db.compounds:
            cpd = self.db.get_seed_compound(cpd_id)
            cpd_alias[cpd.id] = cpd.name
        return cpd_alias

    def collect_roles(self, query):
        """
        WHERE n.key CONTAINS 'symporter'
        WHERE n.key CONTAINS 'antiporter'
        WHERE n.key CONTAINS 'ABC transporter'
        WHERE n.key CONTAINS 'PTS'
        """
        abc_roles = {}
        loop = True
        page = 0
        while loop:
            a, b = self.annotation_api.page_nodes2("Function", page, 10, query)
            if a is None:
                loop = False
            else:
                for o in a:
                    nfunction = self.annotation_api.get_function_by_uid(o["n"].id)
                    subfunctions = nfunction.sub_functions
                    if len(subfunctions) > 0:
                        for subfunction in subfunctions:
                            role_str = subfunction.value
                            role_sn = normalize_role(role_str)
                            if role_sn not in abc_roles:
                                abc_roles[role_sn] = set()
                            abc_roles[role_sn].add(role_str)
                    else:
                        role_str = o["n"]["key"]
                        role_sn = normalize_role(role_str)
                        if role_sn not in abc_roles:
                            abc_roles[role_sn] = set()
                        abc_roles[role_sn].add(role_str)
                page += 1
        return abc_roles

    def get_t_class(self):
        t_permeases = {}
        t_pts = {}
        t_abc = {}
        t_others = {}
        for seed_id in self.db.reactions:
            rxn = self.db.get_seed_reaction(seed_id)
            if not rxn.is_obsolete and rxn.is_transport:
                cstoichiometry = rxn.cstoichiometry
                if len(cstoichiometry) == 2:
                    t_permeases[seed_id] = rxn
                elif is_transport_abc(rxn)[0]:
                    t_abc[seed_id] = rxn
                elif is_transport_pts(rxn)[0]:
                    t_pts[seed_id] = rxn
                else:
                    t_others[seed_id] = rxn
        return t_abc, t_pts, t_permeases, t_others

    def get_compounds_to_reaction(self, t_abc, ban, is_transport_function):
        compounds_to_reaction = {}
        for rxn_id in t_abc:
            if rxn_id not in ban:
                rxn = t_abc[rxn_id]
                is_transport, inside, outside, t, cpd_id = is_transport_function(rxn)
                # print(rxn_id, b, inside, outside, t, cpd_id)
                if is_transport and t == "import":
                    if cpd_id not in compounds_to_reaction:
                        compounds_to_reaction[cpd_id] = (rxn, inside, outside)
                    else:
                        print("!", rxn_id, cpd_id, compounds_to_reaction[cpd_id][0].id)
        return compounds_to_reaction

    def f1(self, df_abc, t_abc, ban, is_transport_function):
        compounds_to_reaction = self.get_compounds_to_reaction(
            t_abc, ban, is_transport_function
        )
        missing_abc_reaction = set()
        reaction_roles = {}
        role_reactions = {}
        for row_id, d in df_abc.iterrows():
            if not pd.isna(d["compound"]):
                for cpd_id in d["compound"].split("/"):
                    if not d["roles"] in role_reactions:
                        role_reactions[d["roles"]] = {}
                    if cpd_id in compounds_to_reaction:
                        rxn, inside, outside = compounds_to_reaction[cpd_id]
                        if not rxn.id in reaction_roles:
                            reaction_roles[rxn.id] = set()
                        reaction_roles[rxn.id].add(d["roles"])
                        role_reactions[d["roles"]][cpd_id] = compounds_to_reaction[
                            cpd_id
                        ]
                    else:
                        role_reactions[d["roles"]][cpd_id] = None
                        missing_abc_reaction.add(cpd_id)
        return missing_abc_reaction, reaction_roles, role_reactions

    def report(self, role_reactions):
        for role_str in role_reactions:
            print(role_str)
            for cpd_id in role_reactions[role_str]:
                print("\t", cpd_id)
                if role_reactions[role_str][cpd_id] == None:
                    print("\t", "missing reaction")
                else:
                    rxn, inside, outside = role_reactions[role_str][cpd_id]
                    cmp_alias = {inside: "in", outside: "out"}
                    # print('\t', rxn.id, ':', to_str2(rxn, cmp_alias, cpd_alias))
                    print("\t", rxn.id, ":", rxn.build_reaction_string(True, cmp_alias))
