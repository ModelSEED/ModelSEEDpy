Metabolite FBA
--------------------------------------

+++++++++++++++++++++
MetaboFBAPkg()
+++++++++++++++++++++

This class constrains metabolites in the parameterized peaks to zero and adds simple thermodynamic constraints:

.. code-block:: python

 from modelseedpy.fbapkg import MetaboFBAPkg
 metFBA = MetaboFBAPkg(model)

- *model* ``cobra.core.model.Model``: The CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 

----------------------
build_package()
----------------------

The COBRA metabolite object is constrained:

.. code-block:: python

 metFBA.build_package(parameters)

- *parameters* ``dict``: The parameters that will govern how the model is constrained, with a default entry
 
.. code-block:: json
 
 {
            "set_objective":true,
 }

and a required key of ``peaks`` whose value is a list of paek data that will be used to identify the metabolites that are constrained.
       
----------------------
Accessible content
----------------------

The ``MetaboFBAPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the MetaboFBAPkg package.