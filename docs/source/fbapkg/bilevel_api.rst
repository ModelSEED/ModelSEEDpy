bilevel
-------------------

+++++++++++++++++++++
BilevelPkg()
+++++++++++++++++++++

This class applies constraints that consider bilevel interactions:

.. code-block:: python

 from modelseedpy.fbapkg import BilevelPkg
 bilevel = BilevelPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

The constraints are applied to the model:

.. code-block:: python

 bilevel.build_package(binary_variable_count = 0)

- *binary_variable_count* ``int``: The quantity of binary variables that will be defined in the model, where ``0`` signifies that no binary variables will be created.

----------------------
Accessible content
----------------------

The ``BilevelPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model, with the corresponding constraints, that is simulated.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters, respectively.
