fbapkg
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


The collection of packages that constrain and investigate, e.g. dFBA, metabolic models are assembled into the ``fbapkg`` directory of ModelSEEDpy. These packages are imported::

 from modelseedpy.fbapkg import *   
   
   
and include the following:

.. toctree::

   bilevel_api
   changeoptpkg_api
   drainfluxes_api
   elementuptakepkg_api
   flexiblebiomass_api
   fluxfitting_api
   gapfillingpkg_api
   kbasemedia_api
   metabofba_api
   proteomefitting_api
   reactionuse_api
   tfa_api
   
