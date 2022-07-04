core
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

The ModelSEEDpy packages that parse and manipulate, e.g. gapfill, metabolic models are assembled into the ``core`` directory of ModelSEEDpy. These packages are imported via::

 from modelseedpy.core import *   
   
and include the following

.. toctree::

   biology_api
   fbahelper_api
   gapfillinghelper_api
   msatpcorrection_api
   msbuilder_api
   mseditorapi_api
   msgenome_api
   msgapfill_api
   msgenomeclassifier_api
   msgrowthphenotypes_api
   msmedia_api
   msmodel_api
   msmodelutils_api
   mstemplate_api
   rpcclient_api
   template_api
    
    
