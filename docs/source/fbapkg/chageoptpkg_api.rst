changeoptpkg
---------------------

+++++++++++++++++++++
ChangeOptPkg()
+++++++++++++++++++++

This class applies constraints that the objective function:

.. code-block:: python

 from modelseedpy.fbapkg import ChangeOptPkg
 bilevel = ChangeOptPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

The constraints are applied to the model:

.. code-block:: python

 bilevel.build_package(target_values = {}, build_objective = True)

- *target_values* ``dict``: the collection of objective coefficients (the ``values`` within a ``objcoef`` key) for a set of reactions (``keys``).
- *build_objective* ``bool``: specifies whether the redefined objective is set to the model.

----------------------
Accessible content
----------------------

The ``ChangeOptPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model, with the corresponding constraints, that is simulated.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters, respectively.
