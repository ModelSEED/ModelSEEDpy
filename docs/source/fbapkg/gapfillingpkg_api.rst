GapFilling Package API
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

The drain reactions are created through this function:

.. code-block:: python

 gapfill.build_package(parameters)

- *parameters* ``dict``: The parameters that will populate the 

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

where the ``default_blacklist`` value is a list of approximately 100 reaction IDs that will not be 
           
----------------------------------
build_exclusion_constraint()
----------------------------------

The drain reactions are created through this function:

.. code-block:: python

 rxnuse.build_exclusion_constraint(flux_values = None)

- *flux_values* ``dict``: A dictionary of the fluxes (``values``) for the reactions IDs (``keys``) that will be constrained, which determines whether the reaction proceeds forwards or backwards. The ``None`` value defaults to determining the fluxes for all reactions in the model.

----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.