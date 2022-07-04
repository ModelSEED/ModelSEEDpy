Proteome Fitting Package
--------------------------------------

+++++++++++++++++++++
ProteomeFittingPkg()
+++++++++++++++++++++

This class constrains the sum progressions of sets of reactions:

.. code-block:: python

 from modelseedpy.fbapkg import ProteomeFittingPkg
 proteofit = ProteomeFittingPkg(model)

- *model* ``cobra.core.model.Model``: The CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 

----------------------
build_package()
----------------------

The reactions that are defined in the protemoe are constrained and used to refine the objective:

.. code-block:: python

 proteofit.build_package(parameters)

- *parameters* ``dict``: The parameters that will supplant default values

.. code-block:: json

 {
            "flux_values":{},
            "kcat_values":{},
            "prot_coef" : 0.1,
            "totalflux" : 1,
            "kcat_coef" : 0.333,
            "obj_kfit":1,
            "obj_kvfit":1,
            "obj_vfit":1,
            "set_objective":1,
            "rescale_vfit_by_flux":true,
            "default_rescaling":0.1,
            "default_expression":10
 }

where keys of ``proteome`` and ``condition``, which lack default values, must be defined.
       
----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.