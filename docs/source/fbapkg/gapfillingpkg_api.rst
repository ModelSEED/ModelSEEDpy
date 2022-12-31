gapfillingpkg
--------------------------------------

+++++++++++++++++++++
GapfillingPkg()
+++++++++++++++++++++

This class constrains the sum progressions of sets of reactions:

.. code-block:: python

 from modelseedpy.fbapkg import ReactionUsePkg
 gapfill = GapfillingPkg(model)

- *model* ``cobra.core.model.Model``: The CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build()
----------------------

The model is gapfilled with default parameters

.. code-block:: json

 {
            "gapfill_all_indecies_with_default_templates": 1,
            "set_objective": 1
 }

and specifications of the template and the minimum objective value:

.. code-block:: python

 gapfill.build(template, minimum_objective=0.01)

- *template* ``list``: The collection of templates that will be used to gap-fill the model.
- *minimum_objective* ``float``: The minimum permissible objective value.

----------------------
build_package()
----------------------

The model is gapfilled:

.. code-block:: python

 gapfill.build_package(parameters)

- *parameters* ``dict``: The parameters that will supplant default values

.. code-block:: json

 {
            "auto_sink": ["cpd02701", "cpd11416", "cpd15302"],
            "extend_with_template":1,
            "model_penalty":1,
            "default_gapfill_models":[],
            "default_gapfill_templates":[],
            "gapfill_templates_by_index":{},
            "gapfill_models_by_index":{},
            "reaction_scores":{},
            "gapfill_all_indecies_with_default_templates":1,
            "gapfill_all_indecies_with_default_models":1,
            "default_excretion":100,
            "default_uptake":-100,
            "minimum_obj":0.01,
            "set_objective":1,
            "blacklist":"default_blacklist"
 }

where the ``default_blacklist`` value is a list of approximately 100 reaction IDs that will not be included while gapfilling.

----------------------------------------------
extend_model_with_model_for_gapfilling()
----------------------------------------------

The reactions and metabolites from a source COBRA model are introduced to the autosink and exchange reactions of the model that is initiated by this class:

.. code-block:: python

 gapfill.extend_model_with_model_for_gapfilling(source_model, index)

- *source_model* ``cobra.core.model.Model``: The COBRA model whose reactions and metabolites will be imposed in the initiated model of the class.
- *index* ``int``: The number that corresponds with the species, which is relevant for distinguishing species in a community model.

----------------------------------------------
extend_model_with_template_for_gapfilling()
----------------------------------------------

Adds new reactions and metabolites from a template to the exchange reactions of the model that is initiated by this class:

.. code-block:: python

 gapfill.extend_model_with_template_for_gapfilling(template, index)

- *template* ``modelseedpy.core.mstemplate.MSTemplateBuilder``: The templates that will be used to gap-fill the model.
- *index* ``int``: The number that corresponds with the species, which is relevant for distinguishing species in a community model.


----------------------------------------------
binary_check_gapfilling_solution()
----------------------------------------------

Redefining the objective to the minimum sum of the reaction fluxes that are in the parameterized COBRA solution:

.. code-block:: python

 check_solution = gapfill.binary_check_gapfilling_solution(solution=None, flux_values=None)

- *solution* ``cobra.core.solution.Solution``: The FBA solution from a simulation of the respective model.
- *flux_values* ``dict``: The forward and reverse fluxes (``values``) are stored within "forward" and "reverse" keys for the IDs of all reactions (``keys``).

**Returns** *check_solution* ``dict``: The collection of "<" or ">" directions for all reversed reactions in the model that are described with gapfilling penalties.

----------------------------------------------
run_test_conditions()
----------------------------------------------

Redefining the objective to the minimum sum of the reaction fluxes that are in the parameterized COBRA solution:

.. code-block:: python

 gapfill.run_test_conditions(condition_list, solution = None, max_iterations = 10)

- *condition_list* ``list``: A list of simulation conditions that will be each be examined during the simulation.
- *solution* ``cobra.core.solution.Solution``: The COBRA solution that contains the reactions that will be examined in simulations for all conditions.
- *max_iterations* ``int``: The number of iterations through which the solution and conditions will be examined.


----------------------------------------------
filter_database_based_on_tests()
----------------------------------------------

Silencing reactions that have associated gapfilling penalties:

.. code-block:: python

 gapfill.run_test_conditions(condition_list)

- *condition_list* ``list``: A list of simulation conditions that are examined during the simulation to acquire the list of reactions to be silenced.


----------------------------------------------
filter_database_based_on_tests()
----------------------------------------------

Silencing reactions that have associated gapfilling penalties:

.. code-block:: python

 gapfill.run_test_conditions(condition_list)

- *condition_list* ``list``: A list of simulation conditions that are examined during the simulation to acquire the list of reactions to be silenced.

----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *gapfilling_penalties* ``dict``: A dictionary
- *new_metabolites* & *new_reactions* ``dict``: Dictionaries of metabolite and reaction COBRA objects (``values``) for all metabolite and reaction IDs (``keys``) that were added to the model.
