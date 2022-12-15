msgapfill
---------------------

+++++++++++++++++++++
MSGapfill
+++++++++++++++++++++

A class that gapfills a given model:

.. code-block:: python

 msgap = MSGapfill(model, default_gapfill_templates=[], default_gapfill_models=[], test_conditions=[], reaction_scores={}, blacklist=[])

- *model* ``cobra.core.model.Model``: The model whose ATP Hydrolysis reaction will be corrected.
- *default_gapfill_templates* & *default_gapfill_models* ``list``: The collection of ``modelseedpy.core.mstemplate.MSTemplate`` templates and ``cobra.core.model.Model`` models, respectively, whose gapfilling penalties will be determined and updated during gapfilling.
- *test_conditions* ``list``: The collection of conditions in which the model is examined for feasibility.
- *reaction_scores* ``dict``: The collection of reaction scores (``values``) for all genes of each core reaction ID (``keys``).
- *blacklist* ``list``: The collection of reaction IDs that will not be included while gapfilling.

------------------------------
run_gapfilling()
------------------------------

Executes gapfilling for the specified model:

.. code-block:: python

 check_solution = msgap.run_gapfilling(media=None, target=None, minimum_obj=0.01, binary_check=True, solver = 'optland-cplex')
 solution = msgap.run_gapfilling(media=None, target=None, minimum_obj=0.01, binary_check=False, solver = 'optland-cplex')

- *media* ``modelseedpy.core.msmedia.MSMedia``: The media in which gapfilling will occur.
- *target* ``str``: The ID of the reaction that will be set as the gapfilling objective.
- *minimum_obj* ``float``: The minimum tolerable objective value, which is constrained as the lower bound of the objective function.
- *binary_check* ``bool``: specifies whether the reaction directions for all gapfilled reactions will be returned.
- *solver* ``str``: The ID specification of the linear programming solver that is used during optimization.

**Returns** *check_solution* ``dict``: The collection of "<" or ">" directions for all reversed reactions in the model that are described with gapfilling penalties.

**Returns** *solution* ``cobra.core.solution.Solution``: The COBRA optimization solution.

------------------------------
integrate_gapfill_solution()
------------------------------

Embeds a gapfilling solution into a model:

.. code-block:: python

 new_model = msgap.run_gapfilling(solution)

- *solution* ``cobra.core.solution.Solution``: The optimization solution that will be embedded in the model within the ``MSGapfill`` object.

**Returns** *new_model* ``cobra.core.model.Model``: The COBRA model that is updated with the gapfilling optimization solution.


------------------------------
gapfill()
------------------------------

``staticMethod`` Executes gapfilling of the specified :

.. code-block:: python

 new_model = MSGapfill.gapfill(model, media=None, target_reaction="bio1", default_gapfill_templates=[], default_gapfill_models=[], test_conditions=[], reaction_scores={}, blacklist=[])

- *model* ``cobra.core.model.Model``: The model whose ATP Hydrolysis reaction will be corrected.
- *media* ``modelseedpy.core.msmedia.MSMedia``: The media in which gapfilling will occur.
- *target_reaction* ``str``: The ID of the reaction that will be set as the gapfilling objective.
- *default_gapfill_templates* & *default_gapfill_models* ``list``: The collection of ``modelseedpy.core.mstemplate.MSTemplate`` templates and ``cobra.core.model.Model`` models, respectively, whose gapfilling penalties will be determined and updated during gapfilling.
- *test_conditions* ``list``: The collection of conditions in which the model is examined for feasibility.
- *reaction_scores* ``dict``: The collection of reaction scores (``values``) for all genes of each core reaction ID (``keys``).
- *blacklist* ``list``: The collection of reaction IDs that will not be included while gapfilling.

**Returns** *new_model* ``cobra.core.model.Model``: The COBRA model that is updated with the gapfilling optimization solution.
