KBase Media Package 
--------------------------------------

+++++++++++++++++++++
KBaseMediaPkg()
+++++++++++++++++++++

This class constrains exchange reactions and media compounds:

.. code-block:: python

 from modelseedpy.fbapkg import KBaseMediaPkg
 kbmedia = KBaseMediaPkg(model)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be edited. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
           
----------------------
build_package()
----------------------

The drain reactions are created through this function:

.. code-block:: python

 kbmedia.build_package(media_or_parameters, default_uptake=None, default_excretion=None)

- *media_or_parameters* ``dict | cobrakbase.core.kbasebiochem.media.Media``: The parameters that govern flux bounds of the exchange reactions, or the media that will be simulated with the specified model. The default parameters 

.. code-block:: json

 {
                "default_uptake": 0,
                "default_excretion": 100,
                "media": null
 }

can be supplanted in the ``parameters`` argument by replacing the ``value`` of each ``key``.

----------------------
Accessible content
----------------------

The ``FluxFittingPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *model* ``cobra.core.model.Model``: The cobrakbase model that possesses the implemented drain reactions.
- *variables* & *parameters* ``dict``: Dictionaries of the linear programming variables and simulation parameters, respectively.
- *pkgmgr* ``modelseedpy.fbapkg.mspackagemanager.MSPackageManager``: The collection of associated classes that are used in the FullThermo package.
- *modelutl* ``modelseedpy.core.msmodelutl.MSModelUtil``: A utilities class that converts between COBRA and CobraKBase nomenclature and performs helpful tasks.