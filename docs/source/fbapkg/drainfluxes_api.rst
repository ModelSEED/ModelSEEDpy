Drain fluxes 
-------------------

+++++++++++++++++++++
DrainFluxPkg()
+++++++++++++++++++++

This class adds drain reactions for each specified drain compound:

.. code-block:: python

 from modelseedpy.fbapkg import DrainFluxPkg
 drainflux = DrainFluxPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be expanded with drain reactions. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
           
----------------------
build_package()
----------------------

The drain reactions are created through this function:

.. code-block:: python

 drainflux.build_package(parameters)

- *parameters* ``dict``: The parameters that govern how the drain reactions will be created and implemented into the model. The default dictionary

.. code-block:: json

 {
            "add_all_intracellular_drains":false,
            "default_uptake":0,
            "default_excretion":100,
            "drain_compounds":{},
            "set_minimal_drain_objective":false,
            "update_drain_fluxes":false
 }

can be supplanted in the ``parameters`` argument by listing the ``key`` to be changed with the new ``value``.

----------------------
Accessible content
----------------------

The ``DrainFluxPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters, respectively.