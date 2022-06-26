Simulating a Model
________________________________________________________________________

Complete ``COBRAkbase`` models can be simulated through a multitude of FBA methods, and a subset of community methods. Individual and community ``COBRAkbase`` models can initially be constrained for various biochemical and physical processes.

++++++++++++
Constraints
++++++++++++

Elemental uptake
---------------------

The total uptake of specific elements can be constrained through the ``ElementUptakePkg`` package:

.. code-block:: python

 from modelseedpy.fbapkg import ElementUptakePkg
 eleup = ElementUptakePkg(model)
 eleup.build_package(element_limits)
 
The package applies uptake limits for any elemental symbol (``values`` and ``keys`` of the ``element_limits`` dictionary, respectively) to a given model. 

Reaction Thermodynamics
-------------------------

The thermodynamic free energy of each reaction, and variable formation energies of each metabolite, can be constrained through the ``FullThermoPkg`` package:

.. code-block:: python

 from modelseedpy.fbapkg import FullThermoPkg
 tfa = FullThermoPkg(model)
 tfa.build_package(parameters)
 
The package applies free energy constraints for all reactions in a given model, based upon the free energy data in the ModelSEED Database. The user must specify the ``"modelseed_db_path" `` in the argument (*parameters* ``dict``). The user is also able to redefine default values through *parameters*, which is detailed in the respective API documentation.

Community member growth rates
-------------------------------

The growth rate of community members can be constrained with a kinetic rate constant through the ``CommKineticPkg`` package:

.. code-block:: python

 from modelseedpy.community import CommKineticPkg
 commkin = CommKineticPkg(model)
 commkin.build_package(kinetic_coef)
 
The package applies the kinetic coefficient (*kinetic_coef* ``float``) to the biomass reaction of each species in a given community model. 



++++++++++++
FBA Methods
++++++++++++

``COBRAkbase`` models, either with or without additional constraints, can be simulated through a few FBA processes. 

Bilevel
---------------------



reactionuse
---------------------



dFBA
---------------------



metabofba
---------------------


+++++++++++++++++++
Community Methods
+++++++++++++++++++

``COBRAkbase`` community models, either with or without additional constraints, can be simulated through a few packages.

MSCommunity
---------------------

The ``MSCommunity`` package is an original package from ModelSEEDpy that resolves metabolite-level cross-feeding and provides a concise API for numerous constraints and FBA methods:

.. code-block:: python

 from modelseedpy.community import MSCommunity 
 mscom = MSCommunity(model, names=[], abundances=None, pfba = True, lp_filename = None)
 solution = mscom.run(media = None, pfba = True)

A community COBRA model is passed to the package with a list of the community members, which are indexed sequentially according to their community number. The abundances of the community members can be provided in relative or absolute terms (``values`` & ``keys``, respectively, in *abdundances* ``dict``). The community can then be simulated in an arbitrary KBase media, where ``None`` specifies a complete media.

Community cross-feeding is calculated through the ``compute_interactions`` interactions function:

.. code-block:: python

 cross_feeding_df = mscom.compute_interactions(solution = None, threshold=1)

An FBA solution, such as that from ``mscom.run``, is parsed to determine the cross-feeding interactions of the community that surpass a flux threshold (*threshold* ``int``). The function returns a `Pandas DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ , which conveniently permits user manipulation of the data.
