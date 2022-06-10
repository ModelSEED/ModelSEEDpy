Template packages
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

**metabolite_template** ``modelseedpy.core.mstemplate.MSTemplateMetabolite``: The metabolite object that embodies the content from the dictionary.

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

**Returns** *met_template_str* ``str``: formulations of the ``MSTemplateMetabolite`` content.