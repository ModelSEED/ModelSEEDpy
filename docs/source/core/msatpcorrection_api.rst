msatpcorrection
---------------------

+++++++++++++++++++++
MSATPCorrection
+++++++++++++++++++++

A class that corrects the ATP Hydrolysis reaction of a given model:

.. code-block:: python

 msatp = MSATPCorrection(model, core_template, atp_medias, compartment="c0", max_gapfilling=None, gapfilling_delta=0, atp_hydrolysis_id=None)

- *model* ``cobra.core.model.Model``: The model whose ATP Hydrolysis reaction will be corrected.
- *core_template* ``modelseedpy.core.mstemplate.MSTemplate``: The model template that is used to correct the model.
- *atp_medias* ``modelseedpy.core.msmedia.MSMedia``: The media that is associated with the growth phenotype.
- *compartment* ``str``: The compartment in which ATP will be added.
- *max_gapfilling* ``float``: The greatest extent of gapfilling for respective model.
- *gapfilling_delta* ``int``: The allowable margin above the best gapfilling score that permits selection of a growth media.
- *atp_hydrolysis_id* ``str``: The ID of the ATP Hydrolysis reaction, where ``None`` specifies a SEED reaction search.

------------------------------
atp_correction()
------------------------------

``staticMethod`` Restores the bounds on all noncore reactions:

.. code-block:: python

 msatpobj = msatp.run_atp_correction(model, coretemplate, atp_medias = None, max_gapfilling = None, gapfilling_delta = 0)

- *model* ``cobra.core.model.Model``: The model whose ATP Hydrolysis reaction will be corrected.
- *coretemplate* ``modelseedpy.core.mstemplate.MSTemplate``: The template that consists of the reactions that be used to assess which reactions in the model are categorized as core.
- *atp_medias* ``list``: The collection of media whose gapfilling requirements to produce ATP will be examined.
- *max_gapfilling* ``float``: The maximal gapfilling score, where ``None`` defaults to the best score from gapfilling the provided media.
- *gapfilling_delta* ``int``: The allowable margin above the best gapfilling score that permits selection of a growth media.

**Returns** *msatpobj* ``modelseedpy.core.msatpcorrection.MSATPCorrection``: The ``MSATPCorrection`` object that possesses the corrected model and all attributes of the correction process.
