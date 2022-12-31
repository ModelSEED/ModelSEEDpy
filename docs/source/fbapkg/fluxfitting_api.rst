fluxfittingpkg
--------------------------------------

+++++++++++++++++++++
FluxFittingPkg()
+++++++++++++++++++++

This class constrains metabolites of the biomass reaction and adjusts the objective values:

.. code-block:: python

 from modelseedpy.fbapkg import FluxFittingPkg
 flexbio = FluxFittingPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

The drain reactions are created through this function:

.. code-block:: python

 flexbio.build_package(parameters)

- *parameters* ``dict``: The parameters that govern how the constraints will be created and implemented into the model. The default dictionary

.. code-block:: json

 {
            "target_flux":{},
            "totalflux":0,
            "set_objective":1,
            "default_rescaling":0.1,
            "rescale_vfit_by_flux":true
 }

can be supplanted in the ``parameters`` argument by listing the ``key`` to be changed with the new ``value``. The ``target_flux`` sub-dictionary specifies the fluxes (``values``) for each reaction ID (``keys``) that will be constrained through this package.

----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.
