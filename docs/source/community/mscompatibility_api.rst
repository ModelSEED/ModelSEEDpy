mscompatibility
--------------------------

+++++++++++++++++++++
MSCompatibility()
+++++++++++++++++++++

This class compatibilizes a collection of individual metabolic models to facilitate syntrophy and accurate analysis of herefrom community models. This standardization leverages the ModelSEED Database as the arbitrator of metabolite and reaction IDs in the metabolic models:

.. code-block:: python

 from modelseedpy.community import MSCompatibility
 ms_compat = MSCompatibility(modelseed_db_path, printing = True)

- *modelseed_db_path* ``str``: the path to the ModelSEED Database, which is only required for the FullThermo, where ``None`` does not apply these constraints.
- *printing* ``bool``: specifies whether results will be printed.

----------------------
standardize()
----------------------

**Staticmethod**

The IDs and names of the metabolites and reactions models are standardized to those of the ModelSEED Database:

.. code-block:: python

 new_models (, optionally unknown_met_ids) = MSCompatibility.standardize(models, metabolites=True, exchanges=True, conflicts_file_name=None, model_names=None, 
                       export_directory=None, view_unknown_mets=True, printing=True, unknown_met_ids=None, changed_metabolites=None, changed_reactions=None)

- *models* ``Iterable o cobra.core.model.Model``: the collection of CobraKBase models that will be standardized. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *metabolites* ``bool``: specifies whether metabolites (``True``) or reactions (``False``) will be compared between the models.
- *exchanges* ``bool``: specifies whether only the exchanges of the models will be standardized.
- *conflicts_file_name* ``str``: the filename to where metabolite conflicts will be exported, where ``None`` does not export.
- *model_names* ``list``: the collection of model names that correspond with the indices of the ``models`` parameter, which is used to distinguish the exported files of each model.
- *export_directory* ``str``: specifies the directory to which all of the content will be exported.
- *view_unknown_mets* & *printing* ``bool``: specifies whether the unknown metabolite IDs and other results of the alignment functionality, respectively, are printed to the console for the user to review.
- *unknown_met_ids*, *changed_metabolites*, & *changed_reactions* ``Iterable``: collections of unknown metabolite IDs and corrected metabolites and reactions. These are internal entities that are passed as argumented to ``standardize()`` by ``align_exchanges()`` when the latter is provided ``True`` through the ``standardize`` parameter.

**returns** the collection of standardized COBRA models, and possibly the collection of unknown metabolite IDs

-----------------------------
align_exchanges()
-----------------------------

**Staticmethod**

The exchange reactions and metabolites of metabolic models are aligned to facilitate community assemblage and cross-feeding interactions:

.. code-block:: python

 new_models (, optionally extras) = MSCompatibility.align_exchanges(models, standardize=False, conflicts_file_name=None, 
                                                    model_names=None, export_directory=None, printing=True, extras=False)
 
- *models* ``Iterable o cobra.core.model.Model``: the collection of CobraKBase models that will be standardized. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *standardize* ``bool``: specifies whether the models will be standardized through the ``standardize`` function.
- *conflicts_file_name* ``str``: the filename to where metabolite conflicts will be exported, where ``None`` does not export.
- *model_names* ``Iterable``: the collection of model names that correspond with the indices of the ``models`` parameter, which is used to distinguish the exported files of each model.
- *export_directory* ``str``: specifies the directory to which all of the content will be exported.
- *printing* ``bool``: specifies whether results of the alignment functionality, respectively, are printed to the console for the user to review.
- *extras* ``bool``: specifies whether the ``unique_mets``, ``unknown_met_ids``, ``changed_metabolites``, and ``changed_reactions`` collections of internal data are provided in a tuple, respectively, as the second returned object from the function. This information is not provided, and the function returns one object, when this parameter is ``False``.

**returns** the collection of aligned COBRA models, and possibly the collection of the aforementioned ``extras``
