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

``staticMethod`` Returns the sought gene based upon a query term of features or aliases:

.. code-block:: python

 MSEditorAPI.add_custom_reaction(query)
 
- *query* ``str``: The search term of a feature ID or gene alias.

**Returns** *gene* ``modelseedpy.core.msgenome.MSGenome``: The gene that matches the search term, where ``None`` signifies that no match was discerned.