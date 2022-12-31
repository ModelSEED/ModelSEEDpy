fbahelper
-------------------

+++++++++++++++++++++
FBAHelper()
+++++++++++++++++++++

This class offers a suite of static method functions that assist users in editing and expanding COBRA models:

.. code-block:: python

 from modelseedpy.core import FBAHelper

---------------------------------------------------
add_autodrain_reactions_to_community_model()
---------------------------------------------------

Drain reactions are added to a model:

.. code-block:: python

 fbahelper.add_autodrain_reactions_to_community_model(model,auto_sink = ["cpd02701", "cpd15302"])

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *auto_sink* ``list``: The collection of metabolites for which drain reactions will be created and added to the model.

--------------------------------
add_drain_from_metabolite_id()
--------------------------------

A drain reaction is constructed for a specified metabolite:

.. code-block:: python

 drain_reaction = fbahelper.add_drain_from_metabolite_id(model, cpd_id, uptake, excretion, prefix='EX_', prefix_name='Exchange for ')

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *cpd_id* ``str``: The ID of the compound for which a drain reaction will be constructed.
- *uptake* & *excretion* ``float``: The magnitudes of metabolite uptake and excretion that define the lower and upper bounds of the drain reaction.
- *prefix* & *prefix_name* ``str``: The prefixes for the drain reaction ID and name, respectively.

**returns** *drain_reaction* ``cobra.core.model.Reaction``: The drain reaction of the respective metabolite.

--------------------------------
test_condition_list()
--------------------------------

A collection of simulation conditions are examined for the specified model:

.. code-block:: python

 test_result = fbahelper.test_condition_list(model, condition_list, pkgmgr)

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *condition_list* ``list``: A list of simulation conditions that will be each examined during the simulation.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.

**returns** *test_result* ``bool``: specifies whether all of the conditions yielded feasible solutions and objective values that surpassed the provided threshold with each condition.

--------------------------------
reaction_expansion_test()
--------------------------------

The reactions that cause simulations to not pass all conditions is collected and returned:

.. code-block:: python

 filtered_list = fbahelper.reaction_expansion_test(model, reaction_list, condition_list, pkgmgr)

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *reaction_list* ``list``: The collections of COBRA reactions that will be knocked-out and then are iteratively restored while assessing whether all of the conditions pass the reaction.
- *condition_list* ``list``: and simulation conditions that will be simulated.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.

**returns** *filtered_list* ``list``: The collection of reactions from the ``reaction_list`` for which at least one of the simulation conditions failed.

-------------------------------------
set_reaction_bounds_from_direction()
-------------------------------------

The reaction bounds are set based upon the reaction direction:

.. code-block:: python

 fbahelper.set_reaction_bounds_from_direction(reaction, direction, add=False)

- *reaction* ``cobra.core.model.Reaction``: A COBRA reaction whose flux bounds will be adjusted based upon the direction of the reaction.
- *direction* ``str``: The ``<`` or ``>`` designation of the reaction direction.
- *add* ``bool``: specifies whether the upper bound for ``<`` directions or the lower bound for ``>`` directions will be assigned to zero.

-------------------------------------
set_objective_from_target_reaction()
-------------------------------------

The FBA reaction objective is defined:

.. code-block:: python

 target_reaction = fbahelper.set_objective_from_target_reaction(model,target_reaction,minimize = False)

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *target_reaction* ``str``: The ID of the COBRA reaction in the parameterized model that will be set as the model objective.
- *minimize* ``bool``: specifies whether the simulation will minimize the objective.

**returns** *target_reaction* ``cobra.core.model.Reaction``: The reaction that is specified to be the simulation objective.

-------------------------------------
compute_flux_values_from_variables()
-------------------------------------

Defines the reaction fluxes for all model reactions:

.. code-block:: python

 flux_values = fbahelper.compute_flux_values_from_variables(model)

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which will be interpreted for the fluxes.

**returns** *flux_values* ``cobra.core.model.Reaction``: The reaction that is specified to be the simulation objective.

-------------------------------------
modelseed_id_from_cobra_metabolite()
-------------------------------------

A ModelSEED compound ID is determined from a COBRA metabolite ID:

.. code-block:: python

 msid = fbahelper.modelseed_id_from_cobra_metabolite(metabolite)

- *model* ``cobra.core.model.Metabolite``: The COBRA metabolite whose ModelSEED ID will be returned from a COBRA ID.

**returns** *msid* ``str``: The ModelSEED metabolite ID that is parsed from the COBRA ID.

-------------------------------------
modelseed_id_from_cobra_reaction()
-------------------------------------

A ModelSEED reaction ID is determined from a COBRA reaction ID:

.. code-block:: python

 msid = fbahelper.modelseed_id_from_cobra_reaction(reaction)

- *model* ``cobra.core.model.Reaction``: The COBRA reaction whose ModelSEED ID will be returned from a COBRA ID.

**returns** *msid* ``str``: The ModelSEED reaction ID that is parsed from the COBRA reaction ID.

-------------------------------------
metabolite_mw()
-------------------------------------

The molecular weight of a metabolite is calculated from its elemental composition, or its chemical formula when the elements are unavailable:

.. code-block:: python

 mw = fbahelper.metabolite_mw(metabolite)

- *model* ``cobra.core.model.Metabolite``: The COBRA metabolite whose molecular weight will be calculated.

**returns** *mw* ``float``: The molecular weight of the parameterized metabolite.

-------------------------------------
elemental_mass()
-------------------------------------

**returns** *elementmasses* ``dict``: A dictionary of all elemental masses (``values``) for all chemical symbols (``keys``).

-------------------------------------
get_modelseed_db_api()
-------------------------------------

**returns** *modelseed_api* ``ModelSEED``: The ModelSEED Database that can be used for ModelSEEDpy operations.

-------------------------------------
is_ex() & is_biomass()
-------------------------------------

Functions that determine whether parameterized reactions are exchange or biomass reactions, respectively:

.. code-block:: python

 result = fbahelper.is_ex(reaction)
 result = fbahelper.is_biomass(reaction)

- *reaction* ``cobra.core.model.Reaction``: The COBRA reaction that will be examined as being either an exchange or biomass reaction, respectively.

**returns** *result* ``bool``: specifies whether the parameterized reaction is one of the two reaction types.

-------------------------------------
find_reaction()
-------------------------------------

Identifies a reaction in a model based upon reaction stoichiometry:

.. code-block:: python

 reaction = fbahelper.find_reaction(model, stoichiometry)

- *model* ``cobra.core.model.Model``: The COBRA model in which the reaction will be searched.
- *stoichiometry* ``dict``: The stoichiometry of the reaction that will be searched in the model.

**returns** *reaction* ``cobra.core.model.Model``: The located COBRA reaction, if its reaction string is identified in the model, or ``None`` otherwise.

-------------------------------------
msid_hash()
-------------------------------------

Assembles all of the COBRA metabolites that correspond to the same ModelSEED metabolite ID:

.. code-block:: python

 reaction = fbahelper.msid_hash(model)

- *model* ``cobra.core.model.Model``: The COBRA model where all of the metabolites will be searched.

**returns** *metabolites* ``dict``: Lists of all COBRA metabolites (``values``) that are represented by a ModelSEED compound ID.

-------------------------------------
rxn_hash()
-------------------------------------

Pairs all reaction strings, in both directions, with their corresponding COBRA reaction object:

.. code-block:: python

 reaction = fbahelper.rxn_hash(model)

- *model* ``cobra.core.model.Model``: The COBRA model where all of the metabolites will be searched.

**returns** *reactions* ``dict``: Lists of all COBRA reaction objects with a designation of their directionality (``values``) according to their reaction strings (``keys``).

-------------------------------------
rxn_compartment()
-------------------------------------

Determines the non-extracellular compartment of the parameterized reaction:

.. code-block:: python

 compartment = fbahelper.rxn_compartment(reaction)

- *model* ``cobra.core.model.Reaction``: The COBRA reaction whose non-extracellular compartment will be provided.

**returns** *compartment* ``str``: The non-extracellular reaction compartment.

-------------------------------------
add_atp_hydrolysis()
-------------------------------------

Adds an ATP Hydrolysis reaction to the parameterized model in the specified compartment:

.. code-block:: python

 reaction_dict = fbahelper.add_atp_hydrolysis(model,compartment)

- *model* ``cobra.core.model.Model``: The COBRA model into which an ATP hydrolysis reaction will be added.
- *compartment* ``str``: The compartment of the ATP hydrolysis reaction, which is sensitive to community models where each species is represented with a unique compartment.

**returns** *reaction_dict* ``dict``: The reaction string, direction, and newness of the reaction (``values``) are specified for the constructed ATP hydrolysis reaction (``keys``).

-------------------------------------
parse_id()
-------------------------------------

Determines the non-extracellular compartment of the parameterized reaction:

.. code-block:: python

 ID_components = fbahelper.parse_id(cobra_obj)

- *cobra_obj* ``str``: The COBRA object ID that will be parsed.

**returns** *ID_components* ``tuple``: The basename, compartment, and index of the COBRA object in a single tuple, respectively.

-------------------------------------
medianame()
-------------------------------------

**returns** *media_id* ``str``: The ID of a media, where "Complete" is provided by default.

-------------------------------------
validate_dictionary()
-------------------------------------

Validates a dictionary based upon the requirements and optional default values are added to the dictionary:

.. code-block:: python

 validated_dictionary = fbahelper.validate_dictionary(dictionary, required_keys, defaults)

- *dictionary* ``dict``: The dictionary that will be validated.
- *required_keys* ``list``: The collection of keys that must be provided in the dictionary.
- *defaults* ``dict``: The default entries that will be added to the dictionary.

**returns** *validated_dictionary* ``dict``: The dictionary that has been validated through the function.
