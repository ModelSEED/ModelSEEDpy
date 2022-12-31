commkineticpkg 
-------------------

+++++++++++++++++++++
CommKineticPkg()
+++++++++++++++++++++

This class applies kinetic constraints to the individual growth rates of community members:

.. code-block:: python

 from modelseedpy.community import CommKineticPkg
 commkin = CommKineticPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.

----------------------
build_package()
----------------------

The constraints are applied to the model:

.. code-block:: python

 commkin.build_package(kinetic_coef, community_model = None)

- *kinetic_coef* ``float``: the kinetic coefficient that will constrain the cross-feeding interactions amongst the community species.
- *community_model* ``float``: the ``MSCommunity`` model object that will be constrained, where ``None`` specifies the model from the initiation of this package.

----------------------
Accessible content
----------------------

The ``CommKineticPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model, with the corresponding constraints, that is constrained.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and the simulation parameters, respectively.
