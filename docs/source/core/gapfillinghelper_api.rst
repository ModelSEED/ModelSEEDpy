gapfillinghelper
------------------

+++++++++++++++++++++
GapfillingHelper
+++++++++++++++++++++

A class of functions that assist modifying models:

.. code-block:: python

 gfhelper = GapfillingHelper(blacklist = [], auto_sink = ["cpd02701_c", "cpd11416_c0", "cpd15302_c"])

- *blacklist* ``list``: The collection of IDs for reactions that will not be examined during gapfilling.
- *auto_sink* ``list``: The collection of IDs for reactions that represent biomass growth.

------------------------------------------
test_reaction_additions_againt_limits()
------------------------------------------

Returns the collection of genes for all roles of all complexes from the template reaction:

.. code-block:: python

 filtered = gfhelper.test_reaction_additions_againt_limits(model,reactions,tests)

- *model* ``cobra.core.model.Model``: The COBRA model whose reactions will be tested.
- *reactions* ``list``: The collection of COBRA reactions in the ``model`` that will be tested.
- *tests* ``list``: The collection of dictionaries for the tests, with keys of ``"media"``, ``"default_uptake"``, ``"default_excretion"``, ``"target"``, & ``"maximize"``.

**Returns** *filtered* ``cobra.core.dictlist.DictList``: The collection of tests that contained an objective maximum of 1 and greater than the test limit.

--------------------------------------------
build_model_extended_for_gapfilling()
--------------------------------------------

Extends a model with reactions and metabolites from collections of source models and possibly templates:

.. code-block:: python

 gapfilling_penalties = gfhelper.build_model_extended_for_gapfilling(extend_with_template = True, source_models = [], input_templates = [], model_penalty = 1, reaction_scores = {})

- *extend_with_template* ``bool``: specifies whether the gapfilling penalties will be extended with a template.
- *source_models* & *input_templates* ``list``: The collections of models and templates whose reactions and metabolites be extend the object model before its gapfilling.
- *model_penalty* ``int``: The gapfilling penalty for reaction flux, with an equal weighting in both directions.
- *reaction_scores* ``dict``: The gapfilling reaction scores (``values``) for each gene of each reaction (``keys``).

**Returns** *gapfilling_penalties* ``dict``: The gapfilling penalties of the extended model.

---------------------------------------------------------
convert_modelreaction() & convert_modelcompound()
---------------------------------------------------------

Formats COBRA reactions and metabolites for ModelSEED operations, respectively:

.. code-block:: python

 cobra_rxn = gfhelper.convert_modelreaction(reaction, bigg=False)
 cobra_met = gfhelper.convert_modelreaction(metabolite, bigg=False)

- *reaction* ``cobra.core.reaction.Reaction``: The COBRA reaction that will be reformatted.
- *metabolite* ``cobra.core.metabolite.Metabolite``: The COBRA metabolite that will be reformatted.
- *bigg* ``bool``: specifies whether the COBRA object originates from a BiGG model, which requires an additional reformulation.

**Returns** *cobra_rxn* ``cobra.core.reaction.Reaction``: The reaction that is generated from the ModelSEED reaction.
**Returns** *cobra_met* ``cobra.core.metabolite.Metabolite``: The metabolite that is generated from the ModelSEED metabolite.

--------------------------------------
binary_check_gapfilling_solution()
--------------------------------------

Constructs binary variables for the direction of all model reactions, the sum of which are minimized and the resulting fluxes are returned:

.. code-block:: python

 flux_values = gfhelper.binary_check_gapfilling_solution(gapfilling_penalties,add_solution_exclusion_constraint)

- *gapfilling_penalties* ``dict``: The collection of gapfilling penalties (``values``) for each direction of all reaction IDs (``keys``).
- *add_solution_exclusion_constraint* ``bool``: specifies whether a binary exclusion constraint will be added based upon the primal flux values, which renders a gapfilled solution infeasible and thus permits the determination of a new solution.

**Returns** *flux_values* ``dict``: The collection of all primal flux values (``values``) for each direction of all reaction IDs (``keys``).

--------------------------------------
create_minimal_reaction_objective()
--------------------------------------

Constructs an objective function that minimizes the flux of gapfilled reactions:

.. code-block:: python

 gene = gfhelper.create_minimal_reaction_objective(penalty_hash, default_penalty = 0)

- *penalty_hash* ``dict``: The collection of gapfilling penalties (``values``) for each direction of all reaction IDs (``keys``), which will be minimized through this function.
- *default_penalty* ``str``: The default gapfill penalty and is the default flux coefficient in the objective function for all reactions.

-----------------------------------------
convert_cobra_compound_to_kbcompound()
-----------------------------------------

Constructs a metadata dictionary of a COBRA Metabolite object that is returned and can be added to a KBase model:

.. code-block:: python

 cpd_data = gfhelper.convert_cobra_compound_to_kbcompound(cpd, kbmodel=None)

- *cpd* ``cobra.core.metabolite.Metabolite``: The COBRA Metabolite that will be converted into a KBase Metabolite.
- *kbmodel* ``cobrakbase model``: The KBase model that will be expanded with ``cpd`` metadata, where ``None`` specifies that the compound will not be added.

**Returns** *cpd_data* ``dict``: The collection of ``cpd`` attributes in key-value pairs.

-----------------------------------------
convert_cobra_reaction_to_kbreaction()
-----------------------------------------

Constructs a metadata dictionary of a COBRA Reaction object that is returned and can be added to a KBase model:

.. code-block:: python

 rxn_data = gfhelper.convert_cobra_reaction_to_kbreaction(rxn, kbmodel, direction="=", add_to_model=True)

- *rxn* ``cobra.core.metabolite.Metabolite``: The COBRA Metabolite that will be converted into a KBase Reaction.
- *kbmodel* ``cobrakbase model``: The KBase model that contains ``rxn``.
- *direction* ``str``: The "<", "=", or ">" direction of ``rxn``.
- *add_to_model* ``bool``: specifies whether the reaction metadata will be added to the KBase model.

**Returns** *rxn_data* ``dict``: The collection of ``rxn`` attributes in key-value pairs.

-----------------------------------------
convert_objective_to_constraint()
-----------------------------------------

Coverts an old objective function into a variable and constructs a new constraint that the new objective must equate the old object. The variable and constraint are added to the ``cobramodel`` in the extant object:

.. code-block:: python

 gfhelper.convert_objective_to_constraint(lower_bound, upper_bound)

- *lower_bound* & *upper_bound* ``float``: The bounds that will contrain the objective function variable.

-----------------------------------------
compute_gapfilled_solution()
-----------------------------------------

Returns the direction for all gapfilled reactions in a model:

.. code-block:: python

 directions = gfhelper.compute_gapfilled_solution(penalty_hash, flux_values = None)

- *penalty_hash* ``dict``: The collection of gapfilling penalties (``values``) for each direction of all reaction IDs (``keys``), which will be minimized through this function.
- *flux_values* ``dict``: The collection of all primal flux values (``values``) for each direction of all reaction IDs (``keys``), where ``None`` constructs *flux_values* from the from ``cobramodel`` in class object.

**Returns** *directions* ``dict``: The collection of directions (``values``) for all reactions in a ``cobramodel`` that are stored in ``penalty_hash``.

-----------------------------------------
add_gapfilling_solution_to_kbase_model()
-----------------------------------------

The gapfilled reactions of a solution are added to a model:

.. code-block:: python

 gfhelper.add_gapfilling_solution_to_kbase_model(newmodel, penalty_hash, media_ref)

- *newmodel* ``cobrakbase Model``: The model to which the gapfilled content will be added.
- *penalty_hash* ``dict``: The collection of gapfilling penalties (``values``) for each direction of all reaction IDs (``keys``), which will be minimized through this function.
- *media_ref* ``str``: The reference of the media that was used to gapfill the model.

-----------------------------------------
compute_reaction_scores()
-----------------------------------------

Returns the gapfilling reaction scores for all events, with possible weighting:

.. code-block:: python

 reaction_genes = gfhelper.compute_reaction_scores(weights=None)

- *weights* ``dict``: The collection of gapfill-weightings (``values``) for each event, via ``"description"``, ``"event_id"``, or ``"id"`` attributes of the event (``keys``). An argument of ``None`` specifies that all events will be equally weighted.

**Returns** *reaction_genes* ``dict``: The collection of reaction scores (``values``) for each gene of all reactions over all ontological events in ``fbamodel``.

-----------------------------------------
replicate_model()
-----------------------------------------

Returns a new model that contains a parameterized number of duplicate content of the ``cobramodel`` in the class object:

.. code-block:: python

 newmodel = gfhelper.replicate_model(count)

- *count* ``int``: The number of copies of the ``cobramodel`` that are added to the new model.

**Returns** *newmodel* ``cobra.core.model.Model``: The duplicated COBRA model.

-----------------------------------------
test_reaction_additions_againt_limits()
-----------------------------------------

Returns a new model that contains a parameterized number of duplicate content of the ``cobramodel`` in the class object:

.. code-block:: python

 newmodel = gfhelper.replicate_model(reactions, directions, tests)

- *reactions* ``dict``: The "<" or ">" reaction directions (``values``) for all COBRA reactions that will be tested (``keys``).
- *tests* ``list``: The collection of tests that will be examined for the reactions in the ``cobrakbase`` model.

**Returns** *filtered_tests* ``dict``: The collection of reaction directions and reaction objects in key-value pairs.
