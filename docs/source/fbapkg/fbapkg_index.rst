fbapkg
________________________________________________________________________

|docs|

.. |docs| image:: https://readthedocs.org/projects/modelseedpy/badge/?version=latest
   :target: https://modelseedpy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

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


The collection of packages that constrain and investigate, e.g. dFBA, metabolic models are assembled into the ``fbapkg`` directory of ModelSEEDpy. These packages are imported::

 from modelseedpy.fbapkg import *   
   
   
The core set of ModelSEEDpy packages consist of the following:

.. toctree::

   dfbapy_api
   tfa_api
   bilevel_api
   changeoptpkg_api
   flexiblebiomass_api
   fluxfitting_api
   