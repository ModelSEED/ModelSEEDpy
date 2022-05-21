Metabolic modeling with the ModelSEED Database in Python
________________________________________________________________________

|Supported Python Versions| 

.. |Supported Python Versions| image:: https://img.shields.io/pypi/pyversions/modelseedpy)
   :target: https://pypi.org/project/modelseedpy/
   :alt: Python versions

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

Metabolic modeling is an pivotal method for computational research in synthetic biology and precision medicine. The metabolic models, such as the constrint-based flux balance analysis (FBA) algorithm, are improved with comprehensive datasets that capture more metabolic chemistry in the model and improve the accuracy of simulation predictions. We therefore developed ModelSEEDpy as a comprehensive suite of packages that bootstrap metabolic modeling with the ModelSEED Database (`Seaver et al., 2021 <https://academic.oup.com/nar/article/49/D1/D575/5912569?login=true>`_ ). These packages parse and manipulate (e.g. gapfill missing reactions or calculated chemical properties of metabolites), constrain (with kinetic, thermodynamics, and nutrient uptake), and simulate cobrakbase models (both individual models and communities). This is achieved by standardizing COBRA models through the   ``cobrakbase`` module into a form that is amenable with the KBase/ModelSEED ecosystem. These functionalities are exemplified in `Python Notebooks <https://github.com/ModelSEED/ModelSEEDpy/examples>`_ . Please submit errors, inquiries, or suggestions as `GitHub issues <https://github.com/ModelSEED/ModelSEEDpy/issues>`_ where they can be addressed by our developers.


----------------------
Installation
----------------------

ModelSEEDpy can be installed via ``pip`` through the ``PyPI`` channel and via ``conda`` through the ``conda-forge`` channel::

 pip install modelseedpy
 conda install modelseedpy
   
   
.. toctree::
    :hidden:

   core_index
   fbapkg_index
   contents
   
