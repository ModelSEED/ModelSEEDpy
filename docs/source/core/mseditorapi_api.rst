mseditor
------------

+++++++++++++++++++++
MSEditorAPI
+++++++++++++++++++++

A class of static methods that offer various editing features

------------------------------------
remove_reactions()
------------------------------------

``staticMethod`` Removes a list of reactions from a model:

.. code-block:: python

 MSEditorAPI.remove_reactions(model, rxn_id_list = [])

- *model* ``core.core.model.Model``: The model that will be edited by removing a reaction.
- *rxn_id_list* ``list``: The IDs of all reaction that will be removed from the model.

----------------------
edit_reaction()
----------------------

``staticMethod`` Defines genome features from a FASTA file:

.. code-block:: python

 MSEditorAPI.edit_reaction(model, rxn_id, direction=None, gpr=None)

- *model* ``core.core.model.Model``: The model that will be edited by removing a reaction.
- *rxn_id* ``str``: The ID of the reaction that will be edited.
- *direction* ``str``: The "=>", "<=", or "<=>" reaction direction that represents thermodynamic favorability.
- *gpr* ``str``: The reaction GPR that will be set in the COBRA model.

------------------------------------
edit_biomass_compound()
------------------------------------

``staticMethod`` Adds a compound to the biomass reaction, or creates a biomass reaction with the compound where a biomass reaction is not defined:

.. code-block:: python

 MSEditorAPI.edit_biomass_compound(model,biomass_id,cpd_id,new_coef)

- *model* ``core.core.model.Model``: The model whose biomass reaction will be edited.
- *biomass_id* & *cpd_id* ``str``: The IDs of the biomass reaction and the compound that will be added to the biomass reaction.
- *new_coef* ``float``: The stoichiometric coefficient of the added compound to the biomass reaction.

---------------------------
compute_molecular_weight()
---------------------------

``staticMethod`` Returns the molecular weight of a model metabolite:

.. code-block:: python

 met_mw = MSEditorAPI.compute_molecular_weight(model, metabolite_id)

- *model* ``core.core.model.Model``: The model that contains the metabolite and thus the molecular formula.
- *metabolite_id* ``str``: The ID of the metabolite whose MW will be calculated.

**Returns** *met_mw* ``float``: The MW of the metabolite.

------------------------
add_custom_reaction()
------------------------

``staticMethod`` Adds a reaction to the model with the parameterized characteristics:

.. code-block:: python

 MSEditorAPI.add_custom_reaction(model,rxn_id,MSEquation,gpr = None)

- *model* ``core.core.model.Model``: The model that will be expanded with the reaction.
- *rxn_id* ``str``: The ID of the reaction that will be added to the model.
- *MSEquation* ``modelseedpy.core.mseditorapi.msequation``: The ModelSEED reaction object that contains both the stoichiometry and direction.
- *gpr* ``str``: The reaction GPR that will be set in the COBRA model.

------------------------
add_ms_reaction()
------------------------

``staticMethod`` Adds a reaction with ModelSEED parameters to a model:

.. code-block:: python

 MSEditorAPI.add_ms_reaction(model, rxn_id, modelseed, compartment_equivalents = {'0':'c0', '1':'e0'}, direction = '>')

- *model* ``core.core.model.Model``: The model that will be expanded with the reaction.
- *rxn_id* ``str``: The ID of the reaction that will be added to the model.
- *modelseed* ``ModelSEED Database``: The ModelSEED Database object that will be used to acquire reaction information.
- *compartment_equivalents* ``dict``: The compartments and their indicies that are used in the reaction.
- *direction* ``str``: The "<", "=", or ">" direction of the reaction.

------------------------
copy_model_reactions()
------------------------

``staticMethod`` Adds specified reactions from a source model to a second model:

.. code-block:: python

 MSEditorAPI.copy_model_reactions(model,source_model,rxn_id_list = [])

- *model* ``core.core.model.Model``: The model that will be expanded with additional reactions.
- *source_model* ``core.core.model.Model``: The model whose reactions will be added to the ``model``.
- *rxn_id_list* ``list``: The list of reactions that may be potentially added, provided that they are in the ``source_model``.

-----------------------------
copy_all_model_reactions()
-----------------------------

``staticMethod`` Adds all new reactions from a source model to a second model:

.. code-block:: python

 MSEditorAPI.copy_model_reactions(model,source_model)

- *model* ``core.core.model.Model``: The model that will be expanded with all new reactions.
- *source_model* ``core.core.model.Model``: The model whose reactions will be added to the ``model``.

+++++++++++++++++++++
MSEquation
+++++++++++++++++++++

A class that constructs and organizes reaction information:

.. code-block:: python

 mse = MSEquation(stoichiometry, direction)

- *stoichiometry* ``dict``: The collection of stoichiometry (``values``) for all metabolites in the reaction (``keys``).
- *direction* ``str``: The "<", "=", or ">" reaction directionality.

------------------------------------
build_from_palsson_string()
------------------------------------

``staticMethod`` Parses a BiGG-formatted reaction string into an amenable form for ModelSEEDpy:

.. code-block:: python

 reaction_object = MSEquation.build_from_palsson_string(equation_string, default_group='c')

- *equation_string* ``str``: The BiGG reaction string.
- *default_group* ``str``: The reactant identifier.

**Returns** *reaction_object* ``modelseedpy.core.mseditorapi.MSEquation``: The ``MSEquation`` object version of a BiGG reaction.
