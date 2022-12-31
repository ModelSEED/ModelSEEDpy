mstemplate
--------------------------------------

+++++++++++++++++++++
TemplateReactionType
+++++++++++++++++++++

A class that defines reaction types, where the attributes ``CONDITIONAL``, ``UNIVERSAL``, ``SPONTANEOUS``, and ``GAPFILLING`` correspond to the strings ``"conditional"``, ``"universal"``, ``"spontaneous"``, and ``"gapfilling"``, respectively:

+++++++++++++++++++++
MSTemplateMetabolite
+++++++++++++++++++++

A function that assembles a unique list of features for the specified genome:

.. code-block:: python

 met_template = MSTemplateMetabolite(cpd_id, formula=None, name='', default_charge=None, mass=None,
                   delta_g=None, delta_g_error=None, is_cofactor=False, abbreviation='', aliases=None)

- *cpd_id*, *name*, *formula*, *abbreviation*, & *aliases* ``str``: The ID, name, formula, abbreviation, and alias of the ModelSEED compound that will be constructed into a template.
- *default_charge*, *mass*, *delta_g*, & *delta_g_error* ``float``: The charge, mass, :math:`\Delta` g, and error in the :math:`\Delta` g that describe the metabolite that will be constructed into a template.
- *is_cofactor* ``bool``: A description of whether the metabolite is a cofactor in metabolic reactions.


------------------
from_dict()
------------------

Returns a ``MSTemplateMetabolite`` object from a metabolite dictionary:

.. code-block:: python

 metabolite_template = met_template.from_dict(met_dict)

- *met_dict* ``dict``: A dictionary description of the ModelSEED compound that possesses the following keys ``"id"``, ``"formula"``, ``"name"``, ``"defaultCharge"``, ``"mass"``, ``"deltaG"``, ``"deltaGErr"``, ``"isCofactor"``, ``"abbreviation"``, & ``"aliases"``.

**Returns** *metabolite_template* ``modelseedpy.core.mstemplate.MSTemplateMetabolite``: The metabolite object that embodies the content from the dictionary.

------------------
get_data()
------------------

Returns the metabolite template content.

.. code-block:: python

 met_template_dict = met_template.get_data()

**Returns** *met_template_dict* ``dict``: A dictionary of the ``MSTemplateMetabolite`` content.

---------------------------------------
__repr__() & __str__() & _repr_html_()
---------------------------------------

Constructs strings of the metabolite template content.

.. code-block:: python

 met_template_str = met_template.__repr__()
 met_template_str = met_template.__str__()
 met_template_str = met_template.__repr_html__()

**Returns** *met_template_str* ``str``: formulations of the ``MSTemplateMetabolite`` content.


+++++++++++++++++++++
MSTemplateSpecies
+++++++++++++++++++++

A class that defines a metabolite of a template:

.. code-block:: python

 met_template = MSTemplateSpecies(comp_cpd_id: str, charge: int, compartment: str, cpd_id, max_uptake=0, template=None)

- *cobra_cpd_id*, *compartment*, & *cpd_id* ``str``: The COBRA ID, compartment, and ModelSEED ID for the respective compound.
- *charge* & *max_uptake* ``int``: The compound charge and max uptake of the compound by the examined model.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The model template in which the compound will be searched.

----------------
to_metabolite()
----------------

Creates a COBRA object for a ModelSEEDpy metabolite:

.. code-block:: python

 met = met_template.to_metabolite(index='0')

- *index* ``string``: The compartment index of the respective metabolite.

**returns** *met* ``cobra.core.metabolite.Metabolite``: The COBRA metabolite object of the respective metabolite.

--------------------------------
compound(), name(), & formula()
--------------------------------

Property methods that return the object, name, or formula of the template compound, respectively:

.. code-block:: python

 template_compound = met_template.compound()
 compound_name = met_template.name()
 compound_formula = met_template.formula()

--------------------------------
name() & formula()
--------------------------------

Property setter methods that set the name or formula of the template compound, respectively:

.. code-block:: python

 met_template.name(name)
 met_template.formula(formula)

- *name* & *formula* ``string``: The name and formula that will be assigned to the respective metabolite.

--------------------------------
from_dict()
--------------------------------

Returns methods that return the template compound object:

.. code-block:: python

 template_compound = met_template.from_dict(met_dict, template)

- *met_dict* ``dict``: A dictionary description of the ModelSEED compound that possesses the following keys ``"id"``, ``"charge"``, ``"templatecompartment_ref"``, & ``"templatecompound_ref"``.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The model template in which the compound will be searched.

------------------
get_data()
------------------

Returns the metabolite template content.

.. code-block:: python

 met_template_dict = met_template.get_data()

**Returns** *met_template_dict* ``dict``: A dictionary of the ``MSTemplateMetabolite`` content.


+++++++++++++++++++++
MSTemplateReaction
+++++++++++++++++++++

A class that defines a metabolite of a template:

.. code-block:: python

 rxn_template = MSTemplateSpecies(rxn_id: str, reference_id: str, name='', subsystem='', lower_bound=0.0,
                 upper_bound=None,reaction_type=TemplateReactionType.CONDITIONAL, gapfill_direction='=',
                 base_cost=1000, reverse_penalty=1000, forward_penalty=1000, status='OK')

- *rxn_id*, *reference_id*, *name*, *subsystem*, & *gapfill_direction* ``str``: The COBRA ID, KBase reference ID, name, subsystem, and gapfilling direction of the respective reaction.
- *lower_bound* & *upper_bound* ``int``: The reaction flux limitations.
- *reaction_type* ``modelseedpy.core.mstemplate.TemplateReactionType``: A description of the reaction type from the set of four options that are offered in the ``TemplateReactionType`` class.
- *base_cost*, *reverse_penalty*, & *forward_penalty* ``float``: defines the minimal flux cost and the specific costs of reverse and forward fluxes, respectively.
- *status* ``str``: specifies the gapfilling status.

----------------------
gene_reaction_rule()
----------------------

Property methods that return the gene complexes for the reaction:

.. code-block:: python

 gene_rules = rxn_template.gene_reaction_rule()

**returns** *gene_rules* ``str``: The set of gene complexes, delimited by ``" or "``.

--------------------------
compartment()
--------------------------

Property methods that return the interesting compartment of the respective reaction:

 comptment = rxn_template.compartment()

**returns** *comptment* ``str``: The interesting compartment character from the reaction.

--------------
to_reaction()
--------------

Creates a COBRA object for a ModelSEEDpy reaction:

.. code-block:: python

 reaction = rxn_template.to_reaction(model=None, index='0')

- *model* ``cobra.core.model.Model``: The CobraKBase model in which the examined reaction exists.
- *index* ``string``: The compartment within which the reaction executes.

**returns** *reaction* ``cobra.core.reaction.Reaction``: The COBRA metabolite object of the respective metabolite.

--------------------------------
from_dict()
--------------------------------

Returns methods that return the template compound object and the name of the template compound, respectively:

.. code-block:: python

 reaction = rxn_template.from_dict(rxn_dict, template)

- *rxn_dict* ``dict``: A dictionary description of the ModelSEED reaction that possesses the following keys ``"id"``, ``"reaction_ref"``, ``"name"``, ``"type"``, ``"GapfillDirection"``, ``"base_cost"``, ``"reverse_penalty"``, ``"forward_penalty"``, & ``"status"``.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The model template in which the reaction will be searched.

------------------------------------
add_complexes() & get_complexes()
------------------------------------

Concatenates a list of complexes to the existing list of complexes, and returns the list of complexes, respectively.

.. code-block:: python

 rxn_template.add_complexes(complex_list)
 complexes = rxn_template.get_complexes()

- *complex_list* ``list``: The list of complexes that will be extended to the existing list of complexes.

**returns** *complexes* ``list``: The collection of complexes in the ``MSTemplateReaction`` object.

------------------
cstoichiometry()
------------------

Property methods that return a dictionary of stoichiometric coefficients for each metabolite in the reaction:

.. code-block:: python

 rxn_stoichiometry = rxn_template.cstoichiometry()

**returns** *rxn_stoichiometry* ``dict``: The stoichiometry of each metabolite in the reaction (``value``) for each metabolite ID and compartment as a tuple (``key``).

--------------
get_roles()
--------------

The set of all roles in all complexes are returned:

.. code-block:: python

 roles = rxn_template.get_roles()

**returns** *roles* ``set``: The set of all roles in the complexes of the ``MSTemplateReaction`` object.

----------------------
get_complex_roles()
----------------------

The creates a dictionary of the roles for each complex:

.. code-block:: python

 roles = rxn_template.get_complex_roles()

**returns** *roles* ``dict``: The set of all roles (``keys``) for each complex in the ``MSTemplateReaction`` object.

------------------
get_data()
------------------

Returns the reaction template content.

.. code-block:: python

 rxn_template_dict = rxn_template.get_data()

**Returns** *rxn_template_dict* ``dict``: A dictionary of the ``MSTemplateReaction`` content.

+++++++++++++++++++++
NewModelTemplateRole
+++++++++++++++++++++

A class that defines a template role for a model:

.. code-block:: python

 new_model_tmp = NewModelTemplateRole(role_id, name, features=None, source='', aliases=None)

- *role_id* & *name* ``str``: The ID and name of the role that will be refined into a template.
- *features* & *aliases* ``list``: The collections of features and aliases of the role that will be translated into a template.
- *source* ``str``: The source of the role.

----------------
from_dict()
----------------

Returns a role template object that is constructed from a dictionary:

.. code-block:: python

 role_template = new_model_tmp.from_dict(role_dict)

- *role_dict* ``dict``: A dictionary description of the ModelSEED compound that possesses the following keys ``"id"``, ``"name"``, ``"features"``, ``"source"``, & ``"aliases"``.

**Returns** *role_template* ``modelseedpy.core.mstemplate.NewModelTemplateRole``: A template role object.

------------------
get_data()
------------------

Returns the reaction template content.

.. code-block:: python

 role_template_dict = new_model_tmp.get_data()

**Returns** *role_template_dict* ``dict``: A dictionary of the ``NewModelTemplateRole`` content.

---------------------------------------
__repr__() & __str__() & _repr_html_()
---------------------------------------

Constructs strings of the template content.

.. code-block:: python

 role_template_str = new_model_tmp.__repr__()
 role_template_str = new_model_tmp.__str__()
 role_template_str = new_model_tmp.__repr_html__()

**Returns** *role_template_str* ``str``: formulations of the ``NewModelTemplateRole`` content.

+++++++++++++++++++++++++++++
NewModelTemplateComplex
+++++++++++++++++++++++++++++

A class that defines a template for a protein complex:

.. code-block:: python

 complex_template = NewModelTemplateComplex(complex_id, name, source='', reference='', confidence=0, template=None)

- *complex_id*, *name*, *source*, *reference* ``str``: The ID, name, source, and reference of the complex that will be refined into a template.
- *confidence* ``int``: A confidence rating of the
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template upon which the complex will be added.

----------------
from_dict()
----------------

Returns a complex template object that is constructed from a dictionary:

.. code-block:: python

 complex = complex_template.from_dict(complex_dict, template)

- *complex_dict* ``dict``: A dictionary description of the ModelSEED compound that possesses the following keys ``"id"``, ``"name"``, ``"source"``, ``"reference"``, & ``"confidence"``.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template upon which the complex will be added.

**Returns** *complex* ``modelseedpy.core.mstemplate.NewModelTemplateComplex``: A complex template object.

-------------
add_role()
-------------

Adds triggering and optional functions of a role to the dictionary of roles for the respective complex:

.. code-block:: python

 complex_template.add_role(role, triggering=True, optional=False)

- *role* ``modelseedpy.core.mstemplate.NewModelTemplateRole``: The role that will be added to the collection of roles for the complex.
- *triggering* & *optional* ``bool``: Descriptions of the role that will be added.

-------------
get_data()
-------------

Returns the complex template information:

.. code-block:: python

 complex_data = complex_template.get_data()

**Returns** *complex_data* ``dict``: The complex template information.

---------------------------------------
__repr__() & __str__() & _repr_html_()
---------------------------------------

Constructs strings of the template content.

.. code-block:: python

 complex_template_str = complex_template.__repr__()
 complex_template_str = complex_template.__str__()
 complex_template_str = complex_template.__repr_html__()

**Returns** *complex_template_str* ``str``: formulations of the ``NewModelTemplateComplex`` content.


+++++++++++++++++++++++++++++
MSTemplateCompartment
+++++++++++++++++++++++++++++

A class that defines template compartments:

.. code-block:: python

 complex_template = NewModelTemplateComplex(compartment_id: str, name: str, ph: float, hierarchy=0, aliases=None)

- *compartment_id* & *name* ``str``: The ID and name of the compartment that will be refined into a template.
- *ph* ``float``: The pH of the compartment.
- *hierarchy* ``float``: The pH of the compartment.
- *aliases* ``list``: The collection of alternative identifications for the compartment.

----------------
from_dict()
----------------

Returns a compartment template object that is constructed from a dictionary:

.. code-block:: python

 compartment = complex_template.from_dict(compartment_dict)

- *compartment_dict* ``dict``: A dictionary description of the ModelSEED compound that possesses the following keys ``"id"``, ``"name"``, ``"pH"``, ``"hierarchy"``, & ``"aliases"``.

**Returns** *compartment* ``modelseedpy.core.mstemplate.MSTemplateCompartment``: A compartment template object.

-------------
get_data()
-------------

Returns the compartment template information:

.. code-block:: python

 complex_data = complex_template.get_data()

**Returns** *complex_data* ``dict``: The complex template information.


+++++++++++++++++++++++++++++
MSTemplate
+++++++++++++++++++++++++++++

A class that defines model templates, while leveraging the aforementioned classes:

.. code-block:: python

 template = MSTemplate(template_id, name='', domain='', template_type='', version=1, info=None, args=None)

- *template_id*, *name*, *domain*, & *template_type* ``str``: The ID, name, domain, and type of the template that will be constructed.
- *version* ``int``: The version of the template.

-----------------------------------------------------------------------------------------------------------
add_compartments(), add_roles(), add_complexes(), add_compounds(), add_comp_compounds(), & add_reactions()
-----------------------------------------------------------------------------------------------------------

Functions that add compartments, roles, complexes, compartment compounds, and reactions, respectively, to the developing template. These functions will only add the provided values to the template when they are all unique:

.. code-block:: python

 template.add_compartments(compartments)
 template.add_roles(roles)
 template.add_complexes(complexes)
 template.add_compounds(compounds)
 template.add_comp_compounds(comp_compounds)
 template.add_reactions(reactions)

- *compartments*, *roles*, *complexes*, *compounds*, *comp_compounds*, & *reactions* ``list``: The collections of compartments, roles, complexes, compounds, comp_compounds, and reactions that will be added to the template, provided that all list elements are not extant in the model.

------------------------------
get_complex_from_roles()
------------------------------

A function that yields a complex based upon a descriptive set of complex roles:

.. code-block:: python

 complex = template.get_complex_from_roles(roles)

- *roles* ``list``: The collection of complex roles that will be used to discern the associated complex.

**Returns** *complex* ``modelseedpy.core.mstemplate.NewModelTemplateComplex``: The complex that is discerned from the collection of roles.

------------------------------
get_last_id_value()
------------------------------

A function that yields the largest id from a collection of COBRA objects:

.. code-block:: python

 last_id = template.get_complex_from_roles(objects)

- *objects* ``list``: The collection of COBRA objects whose IDs will be examined.

**Returns** *last_id* ``int``: The largest ID from the collection of COBRA objects.

------------------------------------------------
get_complex(), get_reaction(), & get_role()
------------------------------------------------

A function that yields the largest id from a collection of COBRA objects:

.. code-block:: python

 complex = template.get_complex(obj_id)
 reaction = template.get_reaction(obj_id)
 role = template.get_role(obj_id)

- *obj_id* ``str``: The COBRA ID whose associated complex, reaction, and role will be examined.

**Returns** *complex* ``modelseedpy.core.mstemplate.NewModelTemplateComplex``: The complex that matches the COBRA object ID.
**Returns** *reaction* ``cobra.core.reaction.Reaction``: The COBRA reaction that matches the ID.
**Returns** *role* ``modelseedpy.core.mstemplate.NewModelTemplateRole``: The role that matches the COBRA objects.

------------
get_data()
------------

A function that returns the template data:

.. code-block:: python

 template_data = template.get_data()

**Returns** *template_data* ``dict``: The template data organized into a dictionary structure.

-----------------
_repr_html_()
-----------------

Constructs and returns strings of the template content:

.. code-block:: python

 template_html = template.__repr_html__()

**Returns** *template_html* ``str``: A str of the template data organized into HTML.


+++++++++++++++++++++++++++++
MSTemplateBuilder
+++++++++++++++++++++++++++++

A class that defines model templates, while leveraging the aforementioned classes:

.. code-block:: python

 template = MSTemplateBuilder(template_id, name='', domain='', template_type='', version=1, info=None,
                 biochemistry=None, biomasses=None, pathways=None, subsystems=None)

- *template_id*, *name*, *domain*, & *template_type* ``str``: The ID, name, domain, and type of the template that will be constructed.
- *version* ``int``: The version of the template.
- *info* ``str``: A description of the template that will be stored with the constructed template.

----------------
from_dict()
----------------

Returns a template builder object that is constructed from a dictionary:

.. code-block:: python

 builder = complex_template.from_dict(template_dict)

- *template_dict* ``dict``: A dictionary description of the template, which possesses keys of ``"id"``, ``"name"``, ``"domain"``, ``"type"``, ``"__VERSION__"``, ``"compartments"``, ``"roles"``, ``"complexes"``, ``"compounds"``, ``"compcompounds"``, ``"reactions"``, ``"biochemistry_ref"``, & ``"biomasses"``.

**Returns** *builder* ``modelseedpy.core.mstemplate.MSTemplateBuilder``: The template builder object that was constructed from the dictionary.

-----------------
from_template()
-----------------

Returns a template builder object whose compartments are copied from an existing template:

.. code-block:: python

 builder = complex_template.from_dict(template)

- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template upon which the complex will be added.

**Returns** *builder* ``modelseedpy.core.mstemplate.MSTemplateBuilder``: The template builder object that was partly copied from the existing template.

-----------------
with_role()
-----------------

Returns the complex reference for the given reaction and role IDs:

.. code-block:: python

 complex_ref = complex_template.with_role(template_rxn, role_ids, auto_complex=False)

- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template upon which the complex will be added.
- *role_ids* ``list``: The collection of role IDs for the complex that will identify the corresponding complex ID.
- *auto_complex* ``bool``: specifies whether a complex will be added from the roles if roles are not identified with an associated complex.

**Returns** *complex_ref* ``str``: The complex reference path with the determined complex ID.

----------------------
with_compartment()
----------------------

Returns a matched compartment with the provided ID, otherwise the compartment is added to the MSTemplateBuilder object and the MSBuilder object is returned:

.. code-block:: python

 compartment = complex_template.with_compartment(cmp_id, name, ph=7, index='0')

- *cmp_id*, *name*, & *index* ``str``: The ID, name, and index of the compartment that will be returned or added to the template.
- *ph* ``float``: The pH of the corresponding compartment.

**Returns** *compartment* ``str``: The compartment, or the first of numerous compartments, that matches the provided ID.

-----------
build()
-----------

The function that amalgamates the content of the MSTemplateBuilder object into a MSTemplate object:

.. code-block:: python

 template = complex_template.build()

**Returns** *template* ``modelseedpy.core.mstemplate.MSTemplate``: The MSTemplate object that is constructed from the content of the MSTemplateBuilder object.
