MScompatibility
--------------------------

+++++++++++++++++++++
MSCompatibility()
+++++++++++++++++++++

This class determines the compatibility of individual models for the construction of a community model, as well as the standardization of models to the metabolite and reaction IDs and names of the ModelSEED Database:

.. code-block:: python

 from modelseedpy.core import MSCompatibility
 ms_compat = MSCompatibility(modelseed_db_path, printing = True)

- *modelseed_db_path* ``str``: the path to the ModelSEED Database, which is only required for the FullThermo, where ``None`` does not apply these constraints. 
- *printing* ``bool``: specifies whether results will be printed.

----------------------
standardize_MSD()
----------------------

The IDs and names of the metabolites and reactions of a model are standardized to those of the ModelSEED Database:

.. code-block:: python

 ms_compat.standardize_MSD(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be standardized. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 

**returns** the standardized COBRA model

-----------------------------
compare_models()
-----------------------------

Determines the consistency of reaction or metabolite IDs and names between two models:

.. code-block:: python

 ms_compat.compare_models(model_1, model_2, metabolites = True, standardize = False)

- *model_1* & *model_2* ``cobra.core.model.Model``: CobraKBase models that will be compared. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *metabolites* ``bool``: specifies whether metabolites (``True``) or reactions (``False``) will be compared between the models.
- *standardize* ``bool``: specifies whether the models will be standardized through the ``standardize_MSD`` function.

**returns** misaligned, model_1, model_2

- *misaligned* ``list``: the collection of the discrepancies between the two models, where the differing entries for each model are both provided.
- *model_1* & *model_2* ``cobra.core.model.Model``: the models that were compared, which is relevant only if the models were also standardized via ``standardize = True``.

----------------------
exchanges()
----------------------

Model variabilities in the exchange fluxes -- such as non-standard metabolite IDs (e.g NH4) and different metabolite IDs for each isomer (e.g. L-alanine is cpd00035 while beta-alanine is cpd00085) -- are systematically corrected to facilitate a compatible community of these models.

.. code-block:: python

 ms_compat.exchanges(model)

- *model_1* & *model_2* ``cobra.core.model.Model``: CobraKBase models that will be compared. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 

**returns** model, unknown_met_ids

- *model* ``cobra.core.model.Model``: the corrected model that is now compatible with other models to assemble a community.
- *unknown_met_ids* ``list``: the collection of non-standard metabolite IDs that were not able to be mapped with a standard ModelSEED Database ID.


----------------------
Accessible content
----------------------

Several objects within the ``MSCompatibility`` class may be useful for subsequent post-processing or troubleshooting of the simulation results:

- *compounds* & *reactions* ``dict``: the complete content of compounds and reactions, respectively, in the ModelSEED Database.
- *compound_names* & *compounds_id_indexed* ``OrderedDict``: ordered dictionaries of ModelSEED compound names and IDs (``values``) according to the corresponding compound IDs and names (``keys``), respectively.
- *reaction_ids* ``OrderedDict``: ordered dictionary of ModelSEED reaction IDs (``values``) according to the corresponding reaction names (``keys``).