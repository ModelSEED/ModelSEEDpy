Reaction Use Package
--------------------------------------

+++++++++++++++++++++
ReactionUsePkg()
+++++++++++++++++++++

This class constrains the sum progressions of sets of reactions:

.. code-block:: python

 from modelseedpy.fbapkg import ReactionUsePkg
 rxnuse = ReactionUsePkg(model)

- *model* ``cobra.core.model.Model``: The CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
           
----------------------
build_package()
----------------------

The drain reactions are created through this function:

.. code-block:: python

 rxnuse.build_package(rxn_filter = None, reversibility = 0)

- *rxn_filter* ``dict``: The reaction directions (``values``) for all reaction IDs (``keys``) that will be constrained, where ``None`` signifies that all reactions will be constrained as equilibria.
- *reversibility* ``bool``: specifies whether the constrained reactions are reversible.
           
----------------------------------
build_exclusion_constraint()
----------------------------------

The drain reactions are created through this function:

.. code-block:: python

 rxnuse.build_exclusion_constraint(flux_values = None)

- *flux_values* ``dict``: A dictionary of the fluxes (``values``) for the reactions IDs (``keys``) that will be constrained, which determines whether the reaction proceeds forwards or backwards. The ``None`` value defaults to determining the fluxes for all reactions in the model.

----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.