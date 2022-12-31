mscommunity 
--------------------------

+++++++++++++++++++++++++
CommunityModelSpecies()
+++++++++++++++++++++++++

This class parses species in a community model based upon the composition of the model biomass reaction:

.. code-block:: python

 from modelseedpy.community import CommunityModelSpecies
 com_species = CommunityModelSpecies(community, biocpd, names)

- *community* ``cobra.core.model.Model``: the CobraKBase community model that will be simulated. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.
- *biocpd* ``cobrakbase.core.kbasefba.fbamodel_metabolite.ModelCompound``: the biomass metabolite that will be investigated in relation to the simulated community model.
- *names* ``list``: the list of species in the community model, which are indexed sequentially according to their community number.

----------------------
disable_species()
----------------------

A species can be practically muted from a community simulation by constraining all reaction fluxes of that species to be zero:

.. code-block:: python

 com_species.disable_species()

The code applies to the species and community that are loaded and parsed in the ``CommunityModel`` class.

-----------------------------
compute_max_biomass()
-----------------------------

The biomass production of the community model is determined when a reaction, in which the parameterized biomass metabolite is produced, is optimized.

.. code-block:: python

 com_species.compute_max_biomass()

The code applies to the species, community, and metabolite that are loaded and parsed in the ``CommunityModel`` class.

**returns** the optimization result

----------------------
compute_max_atp()
----------------------

The biomass production of the community model is determined when the ATP hydrolysis reaction, which is added to the model when it is not present, is optimized.

.. code-block:: python

 com_species.compute_max_atp()

The code applies to the species and community that are loaded and parsed in the ``CommunityModel`` class.

**returns** the optimization result

----------------------
Accessible content
----------------------

Several objects within the ``FullThermo`` class may be useful for subsequent post-processing or troubleshooting of the simulation results:

- *model* ``cobra.core.model.Model``: the cobrakbase model, with the corresponding constraints, that is simulated.
- *species_num* ``int``: the number that is assigned to the species under consideration in the community model.
- *id*  ``str``: the identification of the species under consideration in the community model.
- *abundance* ``int``: the abundance of the species under consideration in the community model.
- *biomasses* ``list``: the collection of community reactions, excluding transport reactions, that produce the investigated metabolite.
- *biomass_drain* ``cobra.core.model.Reaction``: the transport reaction that drains the investigated metabolite into the extracellular environment.

+++++++++++++++++++++
MSCommunity()
+++++++++++++++++++++

This class manipulates and simulates community models:

.. code-block:: python

 from modelseedpy.community import MSCommunity
 mscom = MSCommunity(model, names=[], abundances=None, pfba = True, lp_filename = None)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be simulated. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package.
- *names* ``list``: the list of species in the community model, which are indexed sequentially according to their community number.
- *abundances* ``dict``: the abundances (``values``) of the species in the community model (``keys``), in either absolute or relative terms.
- *pfba* ``bool``: signifies whether parsimonious FBA will be simulated.
- *lp_filename* ``str``: species the filename to which the LP file will be exported, where ``None`` does not export the LP file.

----------------------
set_abundance()
----------------------

The abundances of the community members are implemented in the model, and are normalized to relative abundances:

.. code-block:: python

 mscom.set_abundance(abundances)

- *abundances* ``dict``: the abundances (``values``) of the species in the community model (``keys``), in either absolute or relative terms.

----------------------
set_objective()
----------------------

The simulation objective for the community model is implemented:

.. code-block:: python

 mscom.set_objective(target = None, minimize = False)

- *target* ``str``: the ModelSEED id of the reaction for which the simulation will be optimized.
- *minimize* ``bool``: specifies whether the optimization will maximize or minimize the selected reaction, where ``False`` signifies maximization as the default.

----------------------
constrain()
----------------------

The simulation objective for the community model is implemented:

.. code-block:: python

 mscom.constrain(element_uptake_limit = None, kinetic_coeff = None, modelseed_db_path = None)

- *element_uptake_limit* ``dict``: the upper limits of consumption (``values``) for each element in the simulated system (where the element symbols are ``keys``), where ``None`` does not apply these constraints.
- *kinetic_coeff* ``float``: the kinetic coefficient of cross-feeding amongst members of the simulated community, where ``None`` does not apply this constraint.
- *modelseed_db_path* ``str``: the path to the ModelSEED Database, which is only required for the FullThermo, where ``None`` does not apply these constraints.

----------------------
print_lp()
----------------------

The Linear Programming file of the simulation is exported:

.. code-block:: python

 mscom.print_lp(filename= None)

- *filename* ``str``: the path to which the Linear Programming file of the simulation will be exported.

-----------------------------
compute_interactions()
-----------------------------

The cross-feeding interactions amongst all of the members of the community model are calculated:

.. code-block:: python

 cross_feeding = mscom.compute_interactions(solution = None, threshold=1)

- *solution* ``cobra.core.solution.Solution``: the simulation solution that will be parsed to calculate the cross-feeding interactions. The solution from the last simulation, which is stored within the class, is used when the argument is ``None``.
- *threshold* ``int``: the normalized flux threshold, above which the cross-feeding interactions will be considered.

**returns** *cross_feeding* ``pandas.core.frame.DataFrame`` A `Pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that provides the metabolite-level resolution of cross-feeding for each species in the community.


----------------------
gapfill()
----------------------

The community model will be gap-filled with specified media, templates, models, and conditions:

.. code-block:: python

 mscom.gapfill(media = None, target = None, minimize = False, default_gapfill_templates = [],
                     default_gapfill_models = [], test_conditions = [], reaction_scores = {}, blacklist = [])

- *media* ``str``: the media of the model that will be used for gap-filling, where ``None`` defaults to a complete media.
- *target* ``str``: the ModelSEED id of the reaction that will be optimized during the gap-filling.
- *default_gapfill_templates* & *default_gapfill_models* ``list``: collections of templates and models that will be used for gap-filling the community model.
- test_conditions ``list``: the collection of simulation conditions, including media and objective reactions and directions, that will be used to gap-fill the model.
- *reaction_scores* ``dict``: the highest score (``value``) of each gene (``key2``) for each reaction (``key1``), which rescales penalties via reaction scores and saving genes.
- *blacklist* ``list``: a collection of reaction ids that will not used for gap-filling.

**return** the gap-filled model


--------------------------------
test_individual_species()
--------------------------------

Examines the objective values of individual species in the simulated community:

.. code-block:: python

 mscom.test_individual_species(media = None, allow_interaction = True, run_atp = True, run_biomass = True)

- *media* ``str``: the media of the model that will be used for gap-filling, where ``None`` defaults to a complete media.
- *allow_cross_feeding* ``bool``: specifies whether cross-feeding is permitted.
- *run_atp* & *run_biomass* ``bool``: specify whether the species will be optimized for ATP and Biomass, respectively, and optimized.

**return** ``pandas.core.frame.DataFrame`` A `Pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that entails objective value for ATP and Biomass optimizations for each species in the community.

--------------------------------
atp_correction()
--------------------------------

The ATP hydrolysis reaction is defined in the model and the ``MSATPCorrection`` module is instantiated as ``mscom.atpcorrect`` for post-processing:

.. code-block:: python

 mscom.atp_correction(core_template, atp_medias, compartment="c0", max_gapfilling = None, gapfilling_delta = 0)

- *core_template* ``list``: the collection of templates that will be used to gap-fill the community model.
- *atp_medias* ``list``: the collection of media that will be used for gap-filling.
- *compartment* ``str``: specifies the model compartment to which the ATP hydrolysis reaction will be added.
- *max_gapfilling* & *gapfilling_delta* ``float``: specify the maximum graphfilling score and the acceptable variability from the best gapfilling score, below which a media will be selected for growth of the respective model.

--------------------------------
predict_abundances()
--------------------------------

The relative abundances of species members within a community are approximated from the biomass fluxes in the solution of the community objective:

.. code-block:: python

 mscom.predict_abundances(media = None, pfba = True, kinetic_coeff = None)

- *media* ``str``: the media of the model that will be used for gap-filling, where ``None`` defaults to a complete media.
- *pfba* ``bool``: signifies whether parsimonious FBA will be simulated.
- *kinetic_coeff* ``float``: the kinetic coefficient of cross-feeding amongst members of the simulated community. The combination of ``None`` for this argument and the absence of a defined ``kinetic_coeff`` in the ``MSCommunity`` class defaults to a value of 2000.

**return** ``pandas.core.frame.DataFrame`` A `Pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that provides the estimated abundance for each species in the community.

----------
run()
----------

The community model is simulated, with the :

.. code-block:: python

 solution = mscom.run(media = None, pfba = True)

- *media* ``str``: the media of the model that will be used for gap-filling, where ``None`` defaults to a complete media.
- *pfba* ``bool``: signifies whether parsimonious FBA will be simulated.

**return** *solution* ``cobra.core.solution.Solution`` The solution from simulation of the community model.


----------------------
Accessible content
----------------------

Several objects within the ``FullThermo`` class may be useful for subsequent post-processing or troubleshooting of the simulation results:

- *model* ``cobra.core.model.Model``: the cobrakbase model, with the corresponding constraints, that is simulated.
- *cross_feeding_df* ``pandas.core.frame.DataFrame``: the output DataFrame from the ``compute_interactions`` function that organizes metabolite-resolution of cross-feeding for each species in the community.
- *lp_filename* ``str``: the filename to which the Linear Programming problem is exported. This can alternatively be defined in the ``print_lp()`` function as an argument. The absence of a defined ``lp_filename`` prevents the LP problem from being exported.
- *gapfillings*  ``dict``: the collection of ``MSGapfil`` objects (``values``) for each combination of media and target objective that is parameterized in the function (``key``).
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.
- *solution* ``int``: the FBA solution from the most recent simulation of the model.
- *primary_biomass* & *biomass_drain* ``cobra.core.model.Reaction``: the COBRA model reactions that produce or excrete the biomass compound, respectively.
- *kinetic_coeff* ``float``: the kinetic coefficient that constrained cross-feeding amongst members of the simulated community.
- *element_uptake_limit* ``dict``: the upper limits of consumption (``values``) for each element in the simulated system (``keys``).
- *modelseed_db_path* ``str``: the path to the ModelSEED Database, if the FullThermo constraints were applied.
