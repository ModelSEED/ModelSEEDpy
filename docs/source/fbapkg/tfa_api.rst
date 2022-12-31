fullthermopkg
-----------------------------

+++++++++++++++++++++
FullThermoPkg()
+++++++++++++++++++++

This class applies and simulates rigorous thermodynamic constraints upon COBRA models:

.. code-block:: python

 from modelseedpy.fbapkg import FullThermoPkg
 tfa = FullThermoPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be simulated. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

A cobrakbase model is simulated with the parameterized kinetics data over the defined time and conditions:

.. code-block:: python

 tfa.build_package(parameters, verbose = True)

- *parameters* ``dict``: specifies simulations parameters in base metric units -- e.g. Molar, Kilojoules per mol, and Kelvins -- that are used to populate the FBA simulation. The default concentrations,

.. code-block:: json

 {
            "cpd00067_c0":[0.0000001,0.0000001],
            "cpd00007_c0":[1E-07,8.2E-06],
            "cpd00011_c0":[1E-08,0.0014],
            "cpd00067_e0":[3.16228E-07,3.16228E-07],
            "cpd00009_e0":[0.056,0.056],
            "cpd00048_e0":[0.0030,0.0030],
            "cpd00013_e0":[0.019,0.019],
            "cpd00971_e0":[0.16,0.16],
            "cpd00205_e0":[0.022,0.022],
            "cpd10515_e0":[0.062,0.062],
            "cpd00011_e0":[0.00010,0.00010],
            "cpd00007_e0":[8.2E-06,8.2E-06],
            "cpd00027_e0":[0.020,0.020]
 }

and compartment potentials (the extracellular environment is 0 by definition to facilitate community modeling),

.. code-block:: json

 {
            "e0":0,
            "c0":-160
 }

are supplanted by specified parameters

.. code-block:: json

 {
            "default_max_conc":0.02,
            "default_min_conc":0.000001,
            "default_max_error":5,
            "custom_concentrations":{},
            "custom_deltaG_error":{},
            "compartment_potential":{},
            "temperature":298,
            "filter":null,
            "infeasible_model": false,
            "dgbin":false
 }

that can be adjusted through the ``parameters`` argument. The only required key in ``parameters`` that must be provided by the user is the ``modelseed_db_path`` that enables the import of the ModelSEED Database.

- *verbose* ``bool``: specifies whether simulation details and calculations will be printed, which is valuable for troubleshooting.

----------------------
Accessible content
----------------------

Several objects within the ``FullThermo`` class may be useful for subsequent post-processing or troubleshooting of the simulation results:

- *model* ``cobra.core.model.Model``: The cobrakbase model, with the corresponding constraints, that is simulated.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.
