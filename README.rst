.. image:: https://raw.githubusercontent.com/modelseed/modelseedpy/main/examples/ms-logo-horizontal.png?sanitize=true

Metabolic modeling with the ModelSEED Database
________________________________________________________________________

|PyPI version| |docs| |Downloads| |License|

.. |docs| image:: https://readthedocs.org/projects/modelseedpy/badge/?version=latest
   :target: https://modelseedpy.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |PyPI version| image:: https://img.shields.io/pypi/v/modelseedpy.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/modelseedpy/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/modelseed/modelseedpy/workflows/Test%20modelseedpy/badge.svg
   :target: https://github.com/modelseed/modelseedpy/actions
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

ModelSEEDpy will soon be installable via the ``PyPI`` channel::

 pip install modelseedpy
 
but, until then, the repository must cloned::

 git clone https://github.com/ModelSEED/ModelSEEDpy.git

and then locally installed with ``pip``::

 cd path/to/modelseedpy
 pip install .
   
The associated ModelSEED Database, which is required for a few packages, is simply downloaded by cloning the GitHub repository::

 git clone https://github.com/ModelSEED/ModelSEEDDatabase.git
   
and the path to this repository is passed as an argument to the corresponding packages. 
   
**Windows users** must separately install the ``pyeda`` module: 1) download the appropriate wheel for your Python version from `this website <https://www.lfd.uci.edu/~gohlke/pythonlibs/#pyeda>`_ ; and 2) install the wheel through the following commands in a command prompt/powershell console::

 cd path/to/pyeda/wheel
 pip install pyeda_wheel_name.whl
