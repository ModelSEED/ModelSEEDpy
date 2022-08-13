mscompatibility
--------------------------

+++++++++++++++++++++
MSCompatibility()
+++++++++++++++++++++

This class determines the compatibility of individual models for the construction of a community model, as well as the standardization of models to the metabolite and reaction IDs and names of the ModelSEED Database:

.. code-block:: python

 from modelseedpy.community import MSCompatibility
 ms_compat = MSCompatibility(modelseed_db_path, printing = True)

- *modelseed_db_path* ``str``: the path to the ModelSEED Database, which is only required for the FullThermo, where ``None`` does not apply these constraints. 
- *printing* ``bool``: specifies whether results will be printed.

----------------------
standardize()
----------------------

The IDs and names of the metabolites and reactions of a model are standardized to those of the ModelSEED Database:

.. code-block:: python

 ms_compat.standardize(models, metabolites:bool=True, exchanges:bool=True, conflicts_file_name:str=None, 
                       model_names: list = None, export_directory:str=None)

- *models* ``list o cobra.core.model.Model``: the collection of CobraKBase models that will be standardized. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *metabolites* ``bool``: specifies whether metabolites (``True``) or reactions (``False``) will be compared between the models.
- *exchanges* ``bool``: specifies whether only the exchanges of the models will be standardized.
- *conflicts_file_name* ``str``: the filename to where metabolite conflicts will be exported, where ``None`` does not export.
- *model_names* ``list``: the collection of model names that correspond with the indices of the ``models`` parameter, which is used to distinguish the exported files of each model.
- *export_directory* ``str``: specifies the directory to which all of the content will be exported.

**returns** the standardized COBRA model

-----------------------------
align_exchanges()
-----------------------------

Determines the consistency of reaction or metabolite IDs and names between two models:

.. code-block:: python

 ms_compat.compare_models(model_1, model_2, metabolites = True, standardize = False)

- *models* ``list o cobra.core.model.Model``: the collection of CobraKBase models that will be standardized. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *standardize* ``bool``: specifies whether the models will be standardized through the ``standardize`` function.
- *conflicts_file_name* ``str``: the filename to where metabolite conflicts will be exported, where ``None`` does not export.
- *model_names* ``list``: the collection of model names that correspond with the indices of the ``models`` parameter, which is used to distinguish the exported files of each model.
- *export_directory* ``str``: specifies the directory to which all of the content will be exported.

**returns** the standardized COBRA model


----------------------
Accessible content
----------------------

Several objects within the ``MSCompatibility`` class may be useful for subsequent post-processing or troubleshooting of the simulation results:

- *compounds* & *reactions* ``dict``: the complete content of compounds and reactions, respectively, in the ModelSEED Database.
- *compound_names* & *compounds_id_indexed* ``OrderedDict``: ordered dictionaries of ModelSEED compound names and IDs (``values``) according to the corresponding compound IDs and names (``keys``), respectively.
- *reaction_ids* ``OrderedDict``: ordered dictionary of ModelSEED reaction IDs (``values``) according to the corresponding reaction names (``keys``).