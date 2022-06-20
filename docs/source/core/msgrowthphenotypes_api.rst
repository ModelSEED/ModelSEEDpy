msgrowthphenotypes
---------------------

+++++++++++++++++++++
MSGrowthPhenotype
+++++++++++++++++++++

A class that defines a growth phenotype and constructs media for the phenotype:

.. code-block:: python

 mspheno = MSGrowthPhenotype(obj_id, media=None, growth=None, gene_ko=[], additional_compounds=[], parent=None, name=None)
 
- *obj_id* & *name* ``str``: The ID and name of the growth phenotype.
- *media* ``modelseedpy.core.msmedia.MSMedia``: The media in which the phenotype will be simulated.
- *growth* ``float``: The objective value of the growth phenotype.
- *gene_ko* ``list``: The collection of genes that are knocked-out in association with this phenotype.
- *additional_compounds* ``list``: A collection of compounds that will construct a media through the ``build_media`` function.
- *parent* ``list``: The source of base media information, through the ``base_media``, ``base_uptake``, and ``base_excretion`` attributes.

-----------------------
build_media()
-----------------------

Returns a media that is constructed as the amalgamation of extant ``media`` and ``parent`` media and the ``additional_compounds``:

.. code-block:: python

 media = mspheno.build_media()

**Returns** *media* ``cobra.core.msmedia.MSMedia``: The media that is constituted from the ``additional_compounds`` and existing media in the ``media`` and ``parent`` attributes of the ``MSGrowthPhenotype`` object.

---------------
simulate()
---------------

Simulates the growht phenotype of a model that is defined within ``MSModelUtils`` object:

.. code-block:: python

 results = mspheno.simulate(modelutl, growth_threshold=0.001, add_missing_exchanges=False, save_fluxes=False, pfba=False)

- *modelutl* ``modelseedpy.core.msmodelutl.MSModelUtil``: A ``MSModelUtils`` object which possesses the model that will be manipulated and simulated.
- *growth_threshold* ``float``: The objective value threshold for the growth phenotype that is examined by this class.
- *add_missing_exchanges* ``bool``: specifies whether the missing exchange reactions will be added to the model.
- *save_fluxes* ``bool``: specifies whether the solution fluxes will be stored in the results dictionary.

**Returns** *results* ``dict``: The organization of simulation results and intermediate values in key-value pairs, including whether the growth predictions were correct or false.

---------------------------------
gapfill_model_for_phenotype()
---------------------------------

Formats COBRA reactions and metabolites for ModelSEED operations, respectively:

.. code-block:: python

 gpf_model = gfhelper.convert_modelreaction(modelutl, default_gapfill_templates, test_conditions, default_gapfill_models=[], blacklist=[], growth_threshold=0.001, add_missing_exchanges=False)
 
- *modelutl* ``modelseedpy.core.msmodelutl.MSModelUtil``: A ``MSModelUtils`` object which possesses the model that will be gapfilled.
- *default_gapfill_templates* & *test_conditions* ``list``: A collection of gapfilling templates and test conditions that will be used to gapfill the model.
- *default_gapfill_models* ``list``: The collection of models that will extend ``modelutl.model`` for gapfilling.
- *blacklist* ``list``: The collection of reactions that will not be included during gapfilling.
- *growth_threshold* ``float``: The objective value threshold for the growth phenotype that is examined by this class.
- *add_missing_exchanges* ``bool``: specifies whether the missing exchange reactions will be added to the model.

**Returns** *gpf_model* ``cobra.core.model.Model``: The gapfilled model.

+++++++++++++++++++++
MSGrowthPhenotypes
+++++++++++++++++++++

A class that defines a growth phenotype and combines phenotypes that are defined through ``MSGrowthPhenotype``:

.. code-block:: python

 mspheno = MSGrowthPhenotype(base_media=None, base_uptake=0, base_excretion=1000)
 
- *base_media* ``modelseedpy.core.msmedia.MSMedia``: The media that is associated with the growth phenotype.
- *base_uptake* & *base_excretion* ``int``: The uptake and excretion fluxes for the examined phenotype.

-----------------------
from_compound_hash()
-----------------------

``staticMethod`` Returns a ``MSGrowthPhenotypes`` object that is constructed from a dictionary that describes the compounds of the phenotype:

.. code-block:: python

 growthpheno = MSGrowthPhenotypes.from_compound_hash(compounds, base_media, base_uptake=0, base_excretion=1000)

- *compounds* ``list``: The collection of compounds that will comprise the growth phenotype.
- *base_media* ``str``: The media that is associated with the growth phenotype.
- *base_uptake* & *base_excretion* ``int``: The uptake and excretion fluxes for the examined phenotype.

**Returns** *growthpheno* ``modelseedpy.core.msgrowthphenotypes.MSGrowthPhenotypes``: A growth phenotype that is constructed from a dictionary that describes a compound.

--------------------------------------
from_kbase_object()
--------------------------------------

``staticMethod`` Returns a ``MSGrowthPhenotypes`` object that is constructed from KBase through the kbase API object:

.. code-block:: python

 growthpheno = MSGrowthPhenotypes.from_compound_hash(data, kbase_api)
 
- *data* ``dict``: The collection of phenotypes that will be defined and examined (``values``), under the ``phenotypes`` key.
- *kbase_api* ``KBase API``: The KBase API object that can acquire media information from a KBase reference from each phenotype.

**Returns** *growthpheno* ``modelseedpy.core.msgrowthphenotypes.MSGrowthPhenotypes``: The collective of growth phenotypes that are defined from the ``data`` dictionary.

--------------------------------------
from_kbase_file()
--------------------------------------

``staticMethod`` Returns a ``MSGrowthPhenotypes`` object that is constructed from a KBase TSV file:

.. code-block:: python

 growthpheno = MSGrowthPhenotypes.from_kbase_file(filename, base_media, kbase_api)
 
- *filename* ``str``: The name of the TSV file -- with a header of "media    mediaws    growth    geneko    addtlCpd" -- that will be parsed into a ``MSGrowthPhenotypes`` object.
- *base_media* ``str``: The media that is associated with the growth phenotype.
- *kbase_api* ``KBase API``: The KBase API object that can acquire media information from a KBase reference from each phenotype.

**Returns** *growthpheno* ``modelseedpy.core.msgrowthphenotypes.MSGrowthPhenotypes``: The collective of growth phenotypes that are defined from the ``data`` dictionary.

----------------------
from_ms_file()
----------------------

``staticMethod`` Returns a ``MSGrowthPhenotypes`` object that is constructed from a ModelSEED CSV file:

.. code-block:: python

 growthpheno = MSGrowthPhenotypes.from_ms_file(filename, base_media, base_uptake=0, base_excretion=100)
 
- *filename* ``str``: The name of the CSV file -- with a header of "media    mediaws    growth    geneko    addtlCpd" -- that will be parsed into a ``MSGrowthPhenotypes`` object.
- *base_media* ``str``: The media that is associated with the growth phenotype.
- *base_uptake* & *base_excretion* ``int``: The uptake and excretion fluxes for the examined phenotype.

**Returns** *growthpheno* ``modelseedpy.core.msgrowthphenotypes.MSGrowthPhenotypes``: The collective of growth phenotypes that are defined from the ``data`` dictionary.

--------------------
add_phenotypes()
--------------------

Constructs a metadata dictionary of a COBRA Reaction object that is returned and can be added to a KBase model:

.. code-block:: python

 MSGrowthPhenotypes.add_phenotypes(new_phenotypes)
 
- *new_phenotypes* ``list``: The collection of phenotypes that will be added to the ``MSGrowthPhenotypes`` object list of phenotypes.

----------------------------
simulate_phenotypes()
----------------------------

Coverts an old objective function into a variable and constructs a new constraint that the new objective must equate the old object. The variable and constraint are added to the ``cobramodel`` in the extant object:

.. code-block:: python

 gfhelper.convert_objective_to_constraint(model, biomass, add_missing_exchanges=False, correct_false_negatives=False, template=None, growth_threshold=0.001)
  
- *model* ``cobra.core.model.Model``: The model wqhose phenotypes will be simulated.
- *biomass* ``cobra.core.reaction.Reaction``: The biomass reaction which is set as the model objective.
- *add_missing_exchanges* ``bool``: specifies whether the missing exchange reactions will be added to the model.
- *correct_false_negatives* ``bool``: specifies whether false negatives from each phenotype simulation will be corrected.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The model template that is used to gapfill the model.
- *growth_threshold* ``float``: The objective value threshold for the growth phenotypes.
