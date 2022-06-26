msbuilder
------------

+++++++++++++++++++++
MSBuilder
+++++++++++++++++++++

A class that constructs COBRA metabolic models from genomes and model templates:

.. code-block:: python

 new_model = MSBuilder(genome, template=None)

- *genome* ``modelseedpy.core.msgenome.MSGenome``: The genome that will form the COBRA model.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template that will construct the COBRA model.

------------------------------------
get_gpr_from_template_reaction()
------------------------------------

Returns the collection of genes for all roles of all complexes from the template reaction:

.. code-block:: python

 gpr_dict = new_model.get_gpr_from_template_reaction(template_reaction, allow_incomplete_complexes=True)

- *template_reaction* ``modelseedy.core.mstemplate.MSTemplateReaction``: The reaction whose GPR relationships will be discerned.
- *allow_incomplete_complexes* ``bool``: specifies whether the complexes will be built regardless of total complex determination.

**Returns** *gpr_dict* ``dict``: The collection of genes (``values``) for all roles of all complexes (``keys``) of the parameterized reaction.

------------------------------------
build_exchanges()
------------------------------------

Constructs exchange reactions from model reactions whose metabolites exist in the extracellular compartment:

.. code-block:: python

 reactions_exchanges = new_model.build_exchanges(model, extra_cell='e0')

- *model* ``cobra.core.model.Model``: The COBRA model that will be expanded with exchange reactions.
- *extra_cell* ``str``: The compartment of the excellular solution.

**Returns** *reactions_exchanges* ``list``: The collection of COBRA exchange reactions that were added to the model.

----------------------
build_biomasses()
----------------------

Constructs biomass reaction(s) for a model with a specified template:

.. code-block:: python

 biomass_reactions = new_model.build_biomasses(model, template, index)

- *model* ``cobra.core.model.Model``: The COBRA model that will be expanded with exchange reactions.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template within which the biomass reactions of the model will be constructed.
- *index* ``str``: The compartment index of the respective model.

**Returns** *biomass_reactions* ``list``: The collection of biomass reactions that are created and added to the model.

------------------------------------
auto_select_template()
------------------------------------

Returns the predicted class of a genome per the ``knn_ACNP_RAST_filter`` filter ID:

.. code-block:: python

 genome_class = new_model.auto_select_template()

**Returns** *genome_class* ``Pickle prediction``: The genome class that is predicted for the model based upon its genome and a ``knn_ACNP_RAST_filter`` filter ID.

------------------------------------
build_metabolic_reactions()
------------------------------------

Returns the collection of reactions that are constructed from the gpr set for each reaction in the template:

.. code-block:: python

 reactions = new_model.build_metabolic_reactions(index='0', allow_incomplete_complexes=True)
 
- *index* ``str``: The compartment index of the respective model.
- *allow_incomplete_complexes* ``bool``: specifies whether the complexes will be built regardless of total complex determination.

**Returns** *reactions* ``list``: The collection of formated template reactions that have associated gpr information.

------------------------------------
build_non_metabolite_reactions()
------------------------------------

Returns the collection of reactions that lack gpr information:

.. code-block:: python

 reactions_no_gpr = new_model.build_non_metabolite_reactions(cobra_model, index='0', allow_all_non_grp_reactions=False)
 
- *model* ``cobra.core.model.Model``: The COBRA model that will be expanded with exchange reactions.
- *index* ``str``: The compartment index of the respective model.

**Returns** *reactions_no_gpr* ``list``: The collection of formated template reactions that lack associated gpr information.

------------
build()
------------

Constructs a COBRA model based upon the genome in the MSBuilder class and the provided model ID:

.. code-block:: python

 model = new_model.build(model_id, index='0', annotate_with_rast=True)
 
- *model_id* ``str``: The ID of the model that will be constructed.
- *index* ``str``: The compartment index of the respective model.
- *allow_all_non_grp_reactions* ``bool``: specifies whether non-metabolite reactions will be added to the model.
- *annotate_with_rast* ``bool``: specifies whether the genome will be ontologically annotated via RAST.

**Returns** *cobra_model* ``cobra.core.model.Model``: The COBRA model that is generated from the provided genome and model ID.

------------------------------------
build_full_template_model()
------------------------------------

Constructs a COBRA model from a template:

.. code-block:: python

 model = new_model.build_full_template_model(template, model_id=None, index='0')
 
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template that will guide the model construction.
- *model_id* ``str``: The ID of the model that will be constructed.
- *index* ``str``: The compartment index of the respective model.

**Returns** *cobra_model* ``cobra.core.model.Model``: The COBRA model that is generated from the provided template and model ID.

------------------------------------
build_metabolic_model()
------------------------------------

A concise function that develops a COBRA metabolic model from a genome and various specifications:

.. code-block:: python

 model = new_model.build_full_template_model(model_id, genome, gapfill_media=None, template=None, index='0',
                              allow_all_non_grp_reactions=False, annotate_with_rast=True)
 
- *model_id* ``str``: The ID of the model that will be constructed.
- *genome* ``modelseedpy.core.msgenome.MSGenome``: The genome that will form the COBRA model.
- *gapfill_media* ``modelseedpy.core.msgapfill.MSGapfill``: The media that will be used to gapfill the model, where ``None`` specifies that the model will not be gapfilled.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template that will guide the model construction.
- *index* ``str``: The compartment index of the respective model.
- *allow_all_non_grp_reactions* ``bool``: specifies whether non-metabolite reactions will be added to the model.
- *annotate_with_rast* ``bool``: specifies whether the genome will be ontologically annotated via RAST.

**Returns** *cobra_model* ``cobra.core.model.Model``: The COBRA model that is generated from the provided template and model ID.

--------------------
gapfill_model()
--------------------

A model is gapfilled for a target reaction, extracellular media, and model template:

.. code-block:: python

 model = new_model.gapfill_model(original_mdl, target_reaction, template, media)
 
- *original_mdl* ``cobra.core.model.Model``: The model that will be gapfilled.
- *target_reaction* ``str``: The ID of the reaction that will be defined as the objective during gapfilling.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template that will guide the model gapfilling.
- *media* ``modelseedpy.core.msgapfill.MSGapfill``: The media that will be used to gapfill the model.

**Returns** *cobra_model* ``cobra.core.model.Model``: The gapfilled COBRA model.

