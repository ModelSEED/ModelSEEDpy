elementuptakepkg
-----------------------

+++++++++++++++++++++
ElementUptakePkg()
+++++++++++++++++++++

This class applies constraints of elemental consumption to the nutrient intake, which can influence cross-feeding interactions within a community:

.. code-block:: python

 from modelseedpy.community import ElementUptakePkg
 eleup = ElementUptakePkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

The uptake constraints are applied to the model:

.. code-block:: python

 eleup.build_package(element_limits)

- *element_limits* ``dict``: a dictionary of the uptake limits (``values``) for each elemental symbol (``keys``).

----------------------
Accessible content
----------------------

The ``ElementUptakePkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model, with the corresponding constraints, that is simulated.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters, respectively.
