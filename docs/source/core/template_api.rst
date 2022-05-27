Template Package
-------------------

+++++++++++++++++++++
Template()
+++++++++++++++++++++

This class reformats metabolites into compounds and reactions into an amenable format for templates:

.. code-block:: python

 from modelseedpy.core import Template
 bioplate = Template()

------------------------------
convert_template_compound()
------------------------------

A compound is reformatted for addition to a template:

.. code-block:: python

 met = bioplate.convert_template_compound(cpdid, index)

- *cpdid* ``cobra.core.model.Metabolite``: The compound that will be formatted for addition to a complate.
- *index* ``int``: The index of a respective compound, which will be the suffix of the reformatted compound ID.

**returns** *met* ``cobra.core.model.Metabolite``: The reformatted compound for addition to a model template.
           
--------------------------------
convert_template_reaction()
--------------------------------

A reaction is reformatted for addition to a template:

.. code-block:: python

 rxn = bioplate.convert_template_reaction(model, rxnid, index, for_gapfilling=True)

- *model* ``cobra.core.model.Model``: The COBRA model that contains the reactions which are to be reformatted for the template.
- *rxnid* ``str``: The ID for the reaction that will be reformatted.
- *index* ``int``: The index of a respective reactuin, which will be the suffix of the reformatted reaction ID and the constituent metabolite IDs.
- *for_gapfilling* ``bool``: specifies whether the reaction direction for gapfilling or not gapfilling is used in the formatted reaction.

**returns** *rxn* ``cobra.core.model.Reaction``: The reformatted reaction for addition to a model template.
           
----------------------
Accessible content
----------------------

The ``BilevelPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *compounds*, *compcompounds*, & *reactions* ``DictList``: The assemblies of base compounds, compartmentalized compounds, and reactions that were reformatted by the functions.