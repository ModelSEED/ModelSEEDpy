Core: editing models and assembling communities
________________________________________________________________________

|License|

.. |PyPI version| image:: https://img.shields.io/pypi/v/modelseedpy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/modelseedpy/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/modelseedpy/workflows/Test%20modelseedpy/badge.svg
   :target: https://github.com/freiburgermsu/modelseedpy/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/modelseedpy
   :target: https://pepy.tech/project/modelseedpy
   :alt: Downloads

The ModelSEEDpy packages that compatibilize the transport reactions of individual models and simulate cross-feeding amongst community members are organized in the ``community`` directory of the ModelSEEDpy library. 

----------------------
Import
----------------------

The core packages of ModelSEEDpy are locally imported::

 from modelseedpy.community import *   
   
----------------------
Contents
----------------------

The core set of ModelSEEDpy packages consist of the following:

.. toctree::

    mscommunity_api
    mscompatibility_api
    commkineticpkg_api
    elementuptakepkg_api