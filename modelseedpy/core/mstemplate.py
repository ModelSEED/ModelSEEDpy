# -*- coding: utf-8 -*-
import logging
import copy
import math
from enum import Enum
from cobra.core import Metabolite, Reaction
from cobra.core.dictlist import DictList
from cobra.util import format_long_string
from modelseedpy.core.msmodel import (
    get_direction_from_constraints,
    get_reaction_constraints_from_direction,
    get_cmp_token,
)
from cobra.core.dictlist import DictList

# from cobrakbase.kbase_object_info import KBaseObjectInfo

logger = logging.getLogger(__name__)


class AttrDict(dict):
    """
    Base object to use for subobjects in KBase objects
    """

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class TemplateReactionType(Enum):
    CONDITIONAL = "conditional"
    UNIVERSAL = "universal"
    SPONTANEOUS = "spontaneous"
    GAPFILLING = "gapfilling"


class TemplateBiomassCoefficientType(Enum):
    MOLFRACTION = "MOLFRACTION"
    MOLSPLIT = "MOLSPLIT"
    MULTIPLIER = "MULTIPLIER"
    EXACT = "EXACT"


class MSTemplateMetabolite:
    def __init__(
        self,
        cpd_id,
        formula=None,
        name="",
        default_charge=None,
        mass=None,
        delta_g=None,
        delta_g_error=None,
        is_cofactor=False,
        abbreviation="",
        aliases=None,
    ):
        self.id = cpd_id
        self.formula = formula
        self.name = name
        self.abbreviation = abbreviation
        self.default_charge = default_charge
        self.mass = mass
        self.delta_g = delta_g
        self.delta_g_error = delta_g_error
        self.is_cofactor = is_cofactor
        self.aliases = aliases or []
        self.species = set()
        self._template = None

    @staticmethod
    def from_dict(d):
        return MSTemplateMetabolite(
            d["id"],
            d["formula"],
            d["name"],
            d["defaultCharge"],
            d["mass"],
            d["deltaG"],
            d["deltaGErr"],
            d["isCofactor"] == 1,
            d["abbreviation"],
            d["aliases"],
        )

    def get_data(self):
        return {
            "id": self.id,
            "name": self.name,
            "abbreviation": self.abbreviation if self.abbreviation else "",
            "aliases": [],
            "defaultCharge": self.default_charge if self.default_charge else 0,
            "deltaG": self.delta_g if self.delta_g else 10000000,
            "deltaGErr": self.delta_g_error if self.delta_g_error else 10000000,
            "formula": self.formula if self.formula else "R",
            "isCofactor": 1 if self.is_cofactor else 0,
            "mass": self.mass if self.mass else 0,
        }

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def __str__(self):
        return "{}:{}".format(self.id, self.name)

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Compound identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{format_long_string(self.name)}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{'0x0%x' % id(self)}</td>
            </tr><tr>
                <td><strong>Formula</strong></td><td>{self.formula}</td>
            </tr><tr>
                <td><strong>In {len(self.species)} species</strong></td><td>
                    {format_long_string(', '.join(r.id for r in self.species), 200)}</td>
            </tr>
        </table>"""


class MSTemplateSpecies(Metabolite):
    def __init__(
        self,
        comp_cpd_id: str,
        charge: int,
        compartment: str,
        cpd_id,
        max_uptake=0,
        template=None,
    ):
        self._template_compound = None
        super().__init__(comp_cpd_id, "", "", charge, compartment)
        self._template = template
        self.cpd_id = cpd_id
        self.max_uptake = max_uptake
        if self._template:
            if self.cpd_id in self._template.compounds:
                self._template_compound = self._template.compounds.get_by_id(
                    self.cpd_id
                )

    def to_metabolite(self, index="0"):
        """
        Create cobra.core.Metabolite instance
        :param index: compartment index
        :return: cobra.core.Metabolite
        """
        if index is None:
            index = ""
        cpd_id = f"{self.id}{index}"
        compartment = f"{self.compartment}{index}"
        name = f"{self.name}"
        if len(str(index)) > 0:
            name = f"{self.name} [{compartment}]"
        metabolite = Metabolite(cpd_id, self.formula, name, self.charge, compartment)
        return metabolite

    @property
    def compound(self):
        return self._template_compound

    @property
    def name(self):
        if self._template_compound:
            return self._template_compound.name
        return ""

    @name.setter
    def name(self, value):
        if self._template_compound:
            self._template_compound.name = value

    @property
    def formula(self):
        if self._template_compound:
            return self._template_compound.formula
        return ""

    @formula.setter
    def formula(self, value):
        if self._template_compound:
            self._template_compound.formula = value

    @staticmethod
    def from_dict(d, template=None):
        return MSTemplateSpecies(
            d["id"],
            d["charge"],
            d["templatecompartment_ref"].split("/")[-1],
            d["templatecompound_ref"].split("/")[-1],
            d["maxuptake"],
            template,
        )

    def get_data(self):
        return {
            "charge": self.charge,
            "id": self.id,
            "maxuptake": self.max_uptake,
            "templatecompartment_ref": "~/compartments/id/" + self.compartment,
            "templatecompound_ref": "~/compounds/id/" + self.cpd_id,
        }


class MSTemplateReaction(Reaction):
    def __init__(
        self,
        rxn_id: str,
        reference_id: str,
        name="",
        subsystem="",
        lower_bound=0.0,
        upper_bound=None,
        reaction_type=TemplateReactionType.CONDITIONAL,
        gapfill_direction="=",
        base_cost=1000,
        reverse_penalty=1000,
        forward_penalty=1000,
        status="OK",
        reference_reaction_id=None,
    ):
        """

        :param rxn_id:
        :param reference_id:
        :param name:
        :param subsystem:
        :param lower_bound:
        :param upper_bound:
        :param reaction_type:
        :param gapfill_direction:
        :param base_cost:
        :param reverse_penalty:
        :param forward_penalty:
        :param status:
        :param reference_reaction_id: DO NOT USE THIS duplicate of reference_id
        :param template:
        """
        super().__init__(rxn_id, name, subsystem, lower_bound, upper_bound)
        self.reference_id = reference_id
        self.GapfillDirection = gapfill_direction
        self.base_cost = base_cost
        self.reverse_penalty = reverse_penalty
        self.forward_penalty = forward_penalty
        self.status = status
        self.type = reaction_type.value if isinstance(reaction_type, TemplateReactionType) else reaction_type
        self.complexes = DictList()
        self.templateReactionReagents = {}
        self._template = None

    @property
    def gene_reaction_rule(self):
        return " or ".join(map(lambda x: x.id, self.complexes))

    @gene_reaction_rule.setter
    def gene_reaction_rule(self, gpr):
        pass

    @property
    def compartment(self):
        """

        :return:
        """
        return get_cmp_token(self.compartments)

    def to_reaction(self, model=None, index="0"):
        if index is None:
            index = ""
        rxn_id = f"{self.id}{index}"
        compartment = f"{self.compartment}{index}"
        name = f"{self.name}"
        metabolites = {}
        for m, v in self.metabolites.items():
            if model and m.id in model.metabolites:
                metabolites[model.metabolites.get_by_id(m.id)] = v
            else:
                metabolites[m.to_metabolite(index)] = v

        if len(str(index)) > 0:
            name = f"{self.name} [{compartment}]"
        reaction = Reaction(
            rxn_id, name, self.subsystem, self.lower_bound, self.upper_bound
        )
        reaction.add_metabolites(metabolites)
        reaction.annotation["seed.reaction"] = self.reference_id
        return reaction

    @staticmethod
    def from_dict(d, template):
        metabolites = {}
        complexes = set()
        for o in d["templateReactionReagents"]:
            comp_compound = template.compcompounds.get_by_id(
                o["templatecompcompound_ref"].split("/")[-1]
            )
            metabolites[comp_compound] = o["coefficient"]
        for o in d["templatecomplex_refs"]:
            protein_complex = template.complexes.get_by_id(o.split("/")[-1])
            complexes.add(protein_complex)
        lower_bound, upper_bound = get_reaction_constraints_from_direction(
            d["direction"]
        )
        if "lower_bound" in d and "upper_bound" in d:
            lower_bound = d["lower_bound"]
            upper_bound = d["upper_bound"]
        reaction = MSTemplateReaction(
            d["id"],
            d["reaction_ref"].split("/")[-1],
            d["name"],
            "",
            lower_bound,
            upper_bound,
            d["type"],
            d["GapfillDirection"],
            d["base_cost"],
            d["reverse_penalty"],
            d["forward_penalty"],
            d["status"] if "status" in d else None,
            d["reaction_ref"].split("/")[-1],
        )
        reaction.add_metabolites(metabolites)
        reaction.add_complexes(complexes)
        return reaction

    def add_complexes(self, complex_list):
        self.complexes += complex_list

    @property
    def cstoichiometry(self):
        return {(met.id, met.compartment): coefficient for (met, coefficient) in self.metabolites.items()}

    def remove_role(self, role_id):
        pass

    def remove_complex(self, complex_id):
        pass

    def get_roles(self):
        """

        :return:
        """
        roles = set()
        for cpx in self.complexes:
            roles.update(cpx.roles)
        return roles

    def get_complexes(self):
        return self.complexes

    def get_complex_roles(self):
        res = {}
        for complexes in self.data["templatecomplex_refs"]:
            complex_id = complexes.split("/")[-1]
            res[complex_id] = set()
            if self._template:
                cpx = self._template.get_complex(complex_id)

                if cpx:
                    for complex_role in cpx["complexroles"]:
                        role_id = complex_role["templaterole_ref"].split("/")[-1]
                        res[complex_id].add(role_id)
                else:
                    print(f'The complex for ID {complex_id} does not exist.')
        return res

    def get_data(self):
        template_reaction_reagents = list(
            map(
                lambda x: {
                    "coefficient": x[1],
                    "templatecompcompound_ref": "~/compcompounds/id/" + x[0].id,
                },
                self.metabolites.items(),
            )
        )
        return {
            "id": self.id,
            "name": self.name,
            "GapfillDirection": self.GapfillDirection,
            "base_cost": self.base_cost,
            "reverse_penalty": self.reverse_penalty,
            "forward_penalty": self.forward_penalty,
            "upper_bound": self.upper_bound,
            "lower_bound": self.lower_bound,
            "direction": get_direction_from_constraints(
                self.lower_bound, self.upper_bound
            ),
            "maxforflux": self.upper_bound,
            "maxrevflux": 0 if self.lower_bound > 0 else math.fabs(self.lower_bound),
            "reaction_ref": "kbase/default/reactions/id/" + self.reference_id,
            "templateReactionReagents": template_reaction_reagents,
            "templatecompartment_ref": "~/compartments/id/" + self.compartment,
            "templatecomplex_refs": list(
                map(lambda x: "~/complexes/id/" + x.id, self.complexes)
            ),
            # 'status': self.status,
            "type": self.type,
        }

    # def build_reaction_string(self, use_metabolite_names=False, use_compartment_names=None):
    #    cpd_name_replace = {}
    #    if use_metabolite_names:
    #        if self._template:
    #            for cpd_id in set(map(lambda x: x[0], self.cstoichiometry)):
    #                name = cpd_id
    #                if cpd_id in self._template.compcompounds:
    #                    ccpd = self._template.compcompounds.get_by_id(cpd_id)
    #                    cpd = self._template.compounds.get_by_id(ccpd['templatecompound_ref'].split('/')[-1])
    #                    name = cpd.name
    #                cpd_name_replace[cpd_id] = name
    #        else:
    #            return self.data['definition']
    #    return to_str2(self, use_compartment_names, cpd_name_replace)

    # def __str__(self):
    #    return "{id}: {stoichiometry}".format(
    #        id=self.id, stoichiometry=self.build_reaction_string())


class MSTemplateBiomass:
    def __init__(self):
        pass

    @staticmethod
    def from_dict(d):
        pass

    def add_biomass_component(self):
        pass

    def to_reaction(self, model=None, index="0"):
        pass

    def get_data(self):
        pass


class NewModelTemplateRole:
    def __init__(self, role_id, name, features=None, source="", aliases=None):
        """

        :param role_id:
        :param name:
        :param features:
        :param source:
        :param aliases:
        """
        self.id = role_id
        self.name = name
        self.source = source
        self.features = [] if features is None else features
        self.aliases = [] if aliases is None else aliases
        self._complexes = set()
        self._template = None

    @staticmethod
    def from_dict(d):
        return NewModelTemplateRole(
            d["id"], d["name"], d["features"], d["source"], d["aliases"]
        )

    def get_data(self):
        return {
            "id": self.id,
            "name": self.name,
            "aliases": self.aliases,
            "features": self.features,
            "source": self.source,
        }

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def __str__(self):
        return "{}:{}".format(self.id, self.name)

    def _repr_html_(self):
        return """
        <table>
            <tr>
                <td><strong>Role identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Function</strong></td><td>{format_long_string(self.name)}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{'0x0%x' % id(self)}</td>
            </tr><tr>
                <td><strong>In {len(self._complexes)} complexes</strong></td><td>
                    {format_long_string(
                               ', '.join(r.id for r in self._complexes), 200)}</td>
            </tr>
        </table>"""


class NewModelTemplateComplex:
    def __init__(
        self, complex_id, name, source="", reference="", confidence=0, template=None
    ):
        """

        :param complex_id:
        :param name:
        :param source:
        :param reference:
        :param confidence:
        :param template:
        """
        self.id = complex_id
        self.name = name
        self.source = source
        self.reference = reference
        self.confidence = confidence
        self.roles = {}
        self._template = template

    @staticmethod
    def from_dict(d, template):
        protein_complex = NewModelTemplateComplex(
            d["id"], d["name"], d["source"], d["reference"], d["confidence"], template
        )
        for o in d["complexroles"]:
            role = template.roles.get_by_id(o["templaterole_ref"].split("/")[-1])
            protein_complex.add_role(
                role, o["triggering"] == 1, o["optional_role"] == 1
            )
        return protein_complex

    def add_role(self, role: NewModelTemplateRole, triggering=True, optional=False):
        """
        Add role (function) to the complex
        :param role:
        :param triggering:
        :param optional:
        :return:
        """
        self.roles[role] = (triggering, optional)

    def get_data(self):
        complex_roles = []
        for role in self.roles:
            triggering, optional = self.roles[role]
            complex_roles.append(
                {
                    "triggering": 1 if triggering else 0,
                    "optional_role": 1 if optional else 0,
                    "templaterole_ref": "~/roles/id/" + role.id,
                }
            )
        return {
            "id": self.id,
            "name": self.name,
            "reference": self.reference,
            "confidence": self.confidence,
            "source": self.source,
            "complexroles": complex_roles,
        }

    def __str__(self):
        return " and ".join(
            ["{}{}{}".format(role[0].id, ":trig" if role[1][0] else "", ":optional" if role[1][1] else "") for role in self.roles.items()])

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def _repr_html_(self):
        complexes=format_long_string(', '.join("{}:{}:{}:{}".format(r[0].id, r[0].name, r[1][0], r[1][1]) for r in self.roles.items()), 200)
        return f"""
        <table>
            <tr>
                <td><strong>Complex identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{format_long_string(self.name)}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{'0x0%x' % id(self)}</td>
            </tr><tr>
                <td><strong>Contains {len(self.roles)} role(s)</strong></td><td>
                    {complexes}</td>
            </tr>
        </table>"""
        

class MSTemplateCompartment:
    def __init__(
        self, compartment_id: str, name: str, ph: float, hierarchy=0, aliases=None
    ):
        self.id = compartment_id
        self.name = name
        self.ph = ph
        self.hierarchy = hierarchy
        self.aliases = [] if aliases is None else list(aliases)
        self._template = None

    @staticmethod
    def from_dict(d):
        return MSTemplateCompartment(
            d["id"], d["name"], d["pH"], d["hierarchy"], d["aliases"]
        )

    def get_data(self):
        return {
            "id": self.id,
            "name": self.name,
            "pH": self.ph,
            "aliases": self.aliases,
            "hierarchy": self.hierarchy,
        }


class MSTemplate:
    def __init__(
        self,
        template_id,
        name="",
        domain="",
        template_type="",
        version=1,
        info=None,
        args=None,
    ):
        self.id = template_id
        self.name = name
        self.domain = domain
        self.template_type = template_type
        self.__VERSION__ = version
        self.biochemistry_ref = ""
        self.compartments, self.biomasses, self.reactions = DictList(), DictList(), DictList()
        self.compounds, self.pathways, self.subsystems = DictList(), DictList(), DictList()
        self.roles, self.complexes, self.compcompounds = DictList(), DictList(), DictList()

    def add_compartments(self, compartments: list):
        """

        :param compartments:
        :return:
        """
        duplicates = list(set(self.compartments).intersection(compartments))
        if len(duplicates) > 0:
            logger.error(f"The duplicate compartments {duplicates} cannot be added to the template")
            return None  #!!! Should the non-duplicate compartments still be added?

        for x in compartments:
            x._template = self
        self.compartments += compartments

    def add_roles(self, roles: list):
        """

        :param roles:
        :return:
        """
        duplicates = list(set(self.roles).intersection(roles))
        if len(duplicates) > 0:
            logger.error(f"The duplicate roles {duplicates} cannot be added to the template")
            return None  # !!! Should the non-duplicate compartments still be added?

        for x in roles:
            x._template = self
        self.roles += roles

    def add_complexes(self, complexes: list):
        """

        :param complexes:
        :return:
        """
        duplicates = list(set(self.complexes).intersection(complexes))
        if len(duplicates) > 0:
            logger.error(f"The duplicate complexes {duplicates} cannot be added to the template")
            return None  #!!! Should the non-duplicate compartments still be added?

        roles_to_add = []
        for complx in complexes:
            complx._template = self
            roles_rep = {}
            for role in complx.roles:
                r = role
                if role.id not in self.roles:
                    roles_to_add.append(role)
                else:
                    r = self.roles.get_by_id(role.id)
                roles_rep[r] = complx.roles[role]
                r._complexes.add(complx)
            complx.roles = roles_rep

        self.roles += roles_to_add
        self.complexes += complexes

    def add_compounds(self, compounds: list):
        """

        :param compounds:
        :return:
        """
        duplicates = list(set(self.compounds).intersection(compounds))
        if len(duplicates) > 0:
            logger.error(f"The duplicate compounds {duplicates} cannot be added to the template")
            return None  #!!! Should the non-duplicate compartments still be added?

        for cpd in compounds:
            cpd._template = self
        self.compounds += compounds

    def add_comp_compounds(self, comp_compounds: list):
        """
        Add a compartment compounds (i.e., species) to the template
        :param comp_compounds:
        :return:
        """
        duplicates = list(set(self.compcompounds).intersection(comp_compounds))
        if len(duplicates) > 0:
            logger.error(f"The duplicate comp compounds {duplicates} cannot be added to the template")
            return None  #!!! Should the non-duplicate compartments still be added?

        for comp_cpd in comp_compounds:
            comp_cpd._template = self
            if comp_cpd.cpd_id in self.compounds:
                comp_cpd._template_compound = self.compounds.get_by_id(comp_cpd.cpd_id)
                comp_cpd._template_compound.species.add(comp_cpd)
        self.compcompounds += comp_compounds

    def add_reactions(self, reaction_list: list):
        """

        :param reaction_list:
        :return:
        """
        duplicates = list(set(self.reactions).intersection(reaction_list))
        if len(duplicates) > 0:
            logger.error("unable to add reactions [%s] already present in the template", duplicates)
            return None  # !!! Should the non-duplicate compartments still be added?

        for rxn in reaction_list:
            metabolites_replace = {}
            complex_replace = set()
            rxn._template = self
            for comp_cpd, coefficient in rxn.metabolites.items():
                if comp_cpd.id not in self.compcompounds:
                    self.add_comp_compounds([comp_cpd])
                metabolites_replace[self.compcompounds.get_by_id(comp_cpd.id)] = coefficient
            for cpx in rxn.complexes:
                if cpx.id not in self.complexes:
                    self.add_complexes([cpx])
                complex_replace.add(self.complexes.get_by_id(cpx.id))
            rxn._metabolites = metabolites_replace
            rxn.complexes = complex_replace
        self.reactions += reaction_list

    def get_role_sources(self):
        pass

    def get_complex_sources(self):
        pass

    def get_complex_from_role(self, roles):
        cpx_role_str = ";".join(sorted(roles))
        if cpx_role_str in self.role_set_to_cpx:
            return self.role_set_to_cpx[cpx_role_str]
        return None

    @staticmethod
    def get_last_id_value(object_list, prefix):
        last_id = 0
        for obj in object_list:
            if obj.id.startswith(prefix):
                number_part = id[len(prefix):]
                if len(number_part) == 5:
                    last_id = max(last_id, int(number_part))  
        return last_id

    def get_complex(self, obj_id):
        return self.complexes.get_by_id(obj_id)

    def get_reaction(self, obj_id):
        return self.reactions.get_by_id(obj_id)

    def get_role(self, obj_id):
        return self.roles.get_by_id(obj_id)

    # def _to_object(self, key, data):
    #    if key == 'compounds':
    #        return NewModelTemplateCompound.from_dict(data, self)
    #    if key == 'compcompounds':
    #        return NewModelTemplateCompCompound.from_dict(data, self)
    #    #if key == 'reactions':
    #    #    return NewModelTemplateReaction.from_dict(data, self)
    #    if key == 'roles':
    #        return NewModelTemplateRole.from_dict(data, self)
    #    if key == 'subsystems':
    #        return NewModelTemplateComplex.from_dict(data, self)
    #    return super()._to_object(key, data)

    def get_data(self):
        """
        typedef structure {
            modeltemplate_id id;
            string name;
            string type;
            string domain;
            Biochemistry_ref biochemistry_ref;
            list < TemplateRole > roles;
            list < TemplateComplex > complexes;
            list < TemplateCompound > compounds;
            list < TemplateCompCompound > compcompounds;
            list < TemplateCompartment > compartments;
            list < NewTemplateReaction > reactions;
            list < NewTemplateBiomass > biomasses;
            list < TemplatePathway > pathways;
        } NewModelTemplate;
        """

        d = {
            "__VERSION__": self.__VERSION__,
            "id": self.id,
            "name": self.name,
            "domain": self.domain,
            "biochemistry_ref": self.biochemistry_ref,
            "type": "Test",
            "compartments": list(x.get_data() for x in self.compartments),
            "compcompounds": list(x.get_data() for x in self.compcompounds),
            "compounds": list(x.get_data() for x in self.compounds),
            "roles": list(x.get_data() for x in self.roles),
            "complexes": list(x.get_data() for x in self.complexes),
            "reactions": list(x.get_data() for x in self.reactions),
            "biomasses": list(self.biomasses),
            "pathways": [],
            "subsystems": [],
        }

        if self.drains is not None:
            d["drain_list"] = {c.id: t for c, t in self.drains.items()}

        return d

    def _repr_html_(self):
        """
        taken from cobra.core.Model :)
        :return:
        """
        return """
        <table>
            <tr>
                <td><strong>ID</strong></td>
                <td>{id}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{address}</td>
            </tr><tr>
                <td><strong>Number of metabolites</strong></td>
                <td>{num_metabolites}</td>
            </tr><tr>
                <td><strong>Number of species</strong></td>
                <td>{num_species}</td>
            </tr><tr>
                <td><strong>Number of reactions</strong></td>
                <td>{num_reactions}</td>
            </tr><tr>
                <td><strong>Number of biomasses</strong></td>
                <td>{num_bio}</td>
            </tr><tr>
                <td><strong>Number of roles</strong></td>
                <td>{num_roles}</td>
            </tr><tr>
                <td><strong>Number of complexes</strong></td>
                <td>{num_complexes}</td>
            </tr>
          </table>""".format(
            id=self.id,
            address="0x0%x" % id(self),
            num_metabolites=len(self.compounds),
            num_species=len(self.compcompounds),
            num_reactions=len(self.reactions),
            num_bio=len(self.biomasses),
            num_roles=len(self.roles),
            num_complexes=len(self.complexes),
        )


class MSTemplateBuilder:
    def __init__(
        self,
        template_id,
        name="",
        domain="",
        template_type="",
        version=1,
        info=None,
        biochemistry=None,
        biomasses=None,
        pathways=None,
        subsystems=None,
    ):
        self.id = template_id
        self.version = version
        self.name = name
        self.domain = domain
        self.template_type = template_type
        self.compartments, self.biomasses, self.roles, self.complexes = [], [], [], []
        self.compounds, self.compartment_compounds, self.reactions = [], [], []
        self.info = info
        self.biochemistry_ref = None
        self.drains = {}

    @staticmethod
    def from_dict(d, info=None, args=None):
        """

        :param d:
        :param info:
        :param args:
        :return:
        """
        builder = MSTemplateBuilder(
            d["id"], d["name"], d["domain"], d["type"], d["__VERSION__"], None
        )
        builder.compartments = d["compartments"]
        builder.roles = d["roles"]
        builder.complexes = d["complexes"]
        builder.compounds = d["compounds"]
        builder.compartment_compounds = d["compcompounds"]
        builder.reactions = d["reactions"]
        builder.biochemistry_ref = d["biochemistry_ref"]
        builder.biomasses = d["biomasses"]

        return builder

    @staticmethod
    def from_template(template):
        builder = MSTemplateBuilder()
        for compartment in template.compartments:
            builder.compartments.append(copy.deepcopy(compartment))

        return builder

    def with_compound_modelseed(self, seed_id, modelseed):
        pass

    def with_role(self, template_rxn, role_ids, auto_complex=False):
        # TODO: copy from template curation
        complex_roles = template_rxn.get_complex_roles()
        role_match = {role_id:False for role_id in role_ids}
        for complex_id in complex_roles:
            for role in role_match:
                if role in complex_roles[complex_id]:
                    role_match[role] = True
        all_roles_present = True
        for role in role_match:
            all_roles_present &= role_match[role]
        if all_roles_present:
            logger.debug(f'At least one complex does not express all one role of {role_ids}.')
            return None
        complex_id = self.template.get_complex_from_roles(role_ids)
        if complex_id is None:
            logger.warning(f'A corresponding complex for the roles {role_ids} cannot be found.')
            if auto_complex:
                role_names = set([self.template.get_role(role_id)['name'] for role_id in role_ids])
                logger.warning(f'The complex {role_names} will be added to the template.')
                complex_id = self.template.add_complex_from_role_names(role_names)
            else:
                return None
        complex_ref = '~/complexes/id/' + complex_id
        if complex_ref in template_rxn.data['templatecomplex_refs']:
            logger.debug(f'The template already contains a complex reference {complex_ref} for complex {complex_id}.')
            return None
        return complex_ref

    def with_compound(self):
        pass

    def with_compound_compartment(self):
        pass

    def with_compartment(self, cmp_id, name, ph=7, index='0'):
        res = list(x for x in self.compartments if x['id'] == cmp_id)
        if len(res) > 0:
            return res[0]

        self.compartments.append(
            {
                "id": cmp_id,
                "name": name,
                "aliases": [],
                "hierarchy": 3,  # TODO: what is this?
                "index": index,
                "pH": ph,
            }
        )

        return self

    def build(self):
        template = MSTemplate(self.id, self.name, self.domain, self.template_type, self.version)
        template.add_compartments([MSTemplateCompartment.from_dict(x) for x in self.compartments])
        template.add_compounds([MSTemplateMetabolite.from_dict(x) for x in self.compounds])
        template.add_comp_compounds([MSTemplateSpecies.from_dict(x) for x in self.compartment_compounds])
        template.add_roles([NewModelTemplateRole.from_dict(x) for x in self.roles])
        template.add_complexes([NewModelTemplateComplex.from_dict(x, template) for x in self.complexes])
        template.add_reactions([MSTemplateReaction.from_dict(x, template) for x in self.reactions])
        template.biomasses += [AttrDict(x) for x in self.biomasses]  # TODO: biomass object
        return template
