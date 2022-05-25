Flexible biomass 
--------------------------------------

+++++++++++++++++++++
FlexibleBiomassPkg()
+++++++++++++++++++++

This class constrains metabolites of the biomass reaction and adjusts the objective values:

.. code-block:: python

 from modelseedpy.fbapkg import FlexibleBiomassPkg
 flexbio = FlexibleBiomassPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
           
----------------------
build_package()
----------------------

The drain reactions are created through this function:

.. code-block:: python

 flexbio.build_package(parameters)

- *parameters* ``dict``: The parameters that govern how the constraints be created and implemented into the model. The only require ``key:value`` entry is specifying the ``bio_rxn_id`` of the biomass reaction. The default dictionary

.. code-block:: json

 {
            "flex_coefficient":0.75,
            "use_rna_class":[-0.75,0.75],
            "use_dna_class":[-0.75,0.75],
            "use_protein_class":[-0.75,0.75],
            "use_energy_class":[-0.1,0.1],
 }

can be supplanted in the ``parameters`` argument by listing the ``key`` to be changed with the new ``value``.

----------------------
Accessible content
----------------------

The ``FlexibleBiomassPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.