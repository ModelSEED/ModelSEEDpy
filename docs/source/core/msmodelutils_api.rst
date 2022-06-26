Model Utilities
--------------------------------------

+++++++++++++++++++++
MSModelUtil()
+++++++++++++++++++++

This class offers a suite of utility functions that support editing and manipulating FBA models:

.. code-block:: python

 from modelseedpy.fbapkg import MSModelUtil
 msutil = MSModelUtil(model)

- *model* ``cobra.core.model.Model``: The CobraKBase model that will be constrained. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 

----------------------
find_met()
----------------------

The COBRA metabolite object is located based upon the metabolite name:

.. code-block:: python

 metabolite = msutil.find_met(name)

- *name* ``str``: The name of the metabolite that is will be returned by the function.

**returns** *metabolite* ``cobra.core.model.Metabolite``: The metabolite that was located by the function, or an empty list if the metabolite was not discovered.

----------------------
exchange_list()
----------------------

The returns a list of all COBRA reaction objects for the instantiated model whose reaction IDs possess an "EX_", which signifies an exchange reaction:

.. code-block:: python

 exchange_reactions = msutil.exchange_list()

**returns** *exchange_reactions* ``list``: The collection of exchange reactions for the instantiated model.

----------------------
exchange_hash()
----------------------

The returns a dictionary of COBRA metabolite (``key``) : reaction (``value``) pairs for all metabolites in the exchange reactions from ``msutil.exchange_list()``:

.. code-block:: python

 exchange_reactions = msutil.exchange_hash()

**returns** *exchange_reactions* ``dict``: The COBRA metabolite (``key``) : reaction (``value``) pairs for all model exchange reactions.

-----------------------------
add_missing_exchanges()
-----------------------------

The media compounds that are defined in the exchange reactions of the model are defined with drain reactions:

.. code-block:: python

 media_compounds = msutil.add_missing_exchanges()

**returns** *media_compounds* ``list``: The collection of media compounds that are represented in the exchange reactions of the model and for which drain reactions have been created.
       
-------------------------------------
add_exchanges_for_metabolites()
-------------------------------------

The media compounds that are defined in the exchange reactions of the model are defined with drain reactions:

.. code-block:: python

 drain_reactions = msutil.add_exchanges_for_metabolites(cpds, uptake=0, excretion=0, prefix='EX_', prefix_name='Exchange for ')

- *cpds* ``list``: The collection of media metabolites that represented in the model and will be added to the model as drain reactions.
- *uptake* & *excretion* ``int``: The magnitudes of the lower and upper bounds of the drain reaction, respectively.
- *prefix* & *prefix_name* ``str``: The prefixes of the compound ID and name, respectively.

**returns** *drain_reactions* ``list``: The drain reactions that were created for the parameterized list of compounds.
       
-------------------------------------------
convert_cobra_compound_to_kbcompound()
-------------------------------------------

The information of the parameterized compound will be organized into an amenable format for addition to the ``modelcompounds`` attribute of CobraKBase models:

.. code-block:: python

 cpd_data = msutil.convert_cobra_compound_to_kbcompound(cpd, kbmodel = None)

- *cpd* ``cobra.core.model.Metaoblite``: The COBRA metabolite whose information will be formatted as a KBase metabolite.
- *kbmodel* ``cobra.core.model.Model``: The CobraKBase model whose ``modelcompounds`` attribute will be appended with data from the COBRA metabolite, where ``None`` specifies that the defined dictionary of compound information will not be added to a model.

**returns** *cpd_data* ``dict``: The dictionary of compound information in the format of the ``modelcompounds`` attribute of CobraKBase models.
       
-------------------------------------------
convert_cobra_reaction_to_kbreaction()
-------------------------------------------

The information of the parameterized reaction will be organized into an amenable format for addition to the ``modelreactions`` attribute of CobraKBase models:

.. code-block:: python

 rxn_data = msutil.convert_cobra_reaction_to_kbreaction(rxn, kbmodel, cpd_hash, direction = "=", add_to_model = 1, reaction_genes = {})

- *rxn* ``cobra.core.model.Reaction``: The COBRA reaction whose information will be formatted as a KBase reaction.
- *kbmodel* ``cobra.core.model.Model``: The CobraKBase model whose ``modelreactions`` attribute will be appended with data from the COBRA reaction, where ``None`` specifies that the defined dictionary of compound information will not be added to a model.
- *direction* ``str``: Signification of the reversibility of the reaction as either "<", ">", or "=" as equilibrium.
- *reaction_genes* ``dict``: The collection of contribution (``values``) for each gene (``keys``) that contribute to each reaction (``keys``).

**returns** *rxn_data* ``dict``: The dictionary of reaction information in the format of the ``modelreactions`` attribute of CobraKBase models.
       
-------------------------------------------
add_gapfilling_solution_to_kbase_model()
-------------------------------------------

The parameterized CobraKBase model will be expanded with the content of a gapfilling solution:

.. code-block:: python

 rxn_data = msutil.add_gapfilling_solution_to_kbase_model(newmodel, gapfilled_reactions, gfid=None, media_ref = None, reaction_genes = None)

- *newmodel* ``cobra.core.model.Model``: The CobraKBase model whose information will be formatted as a KBase reaction.
- *gapfilled_reactions* ``dict``: The collection of COBRA reactions (``values``) in each "new" and "reversed" category of reactions (``keys``).
- *gfid* ``str``: The gapfilling ID, which defaults to "gf.#" where # is the smallest unused index.
- *media_ref* ``str``: The reference for the gapfilling media.
- *reaction_genes* ``dict``: The collection of contribution (``values``) for each gene (``keys``) that contribute to each reaction (``keys``).

**returns** *rxn_table* ``list``: A collection of dictionaries, one for each reaction that is added to the CobraKBase model.
       
----------------------
Accessible content
----------------------

The ``MSModelUtil`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the MSModelUtil package.
- *metabolite_hash* & *search_metabolite_hash* ``dict``: Lists of metabolite matches (``values``) for each metabolite name and refined search name, respectively.