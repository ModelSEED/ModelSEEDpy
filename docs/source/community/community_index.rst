community
________________________________________________________________________

|PyPI version| |License| |Downloads|

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

The ModelSEEDpy packages that compatibilize the transport reactions of individual models and simulate cross-feeding amongst community members are organized in the ``community`` directory of the ModelSEEDpy library. These packages are imported via::

 from modelseedpy.community import *   

and include the following:

.. toctree::

    commkineticpkg_api
    dfbapkg_api
    mscommunity_api
    mscompatibility_api