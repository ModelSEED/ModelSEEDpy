ModelSEEDpy API
________________________________________________________________________

The detailed documentation of all user-operable classes and functions in the ModelSEEDpy library are provided.


core
---------

The ``core`` sub-library permits parsing and manipulating metabolic models::

 from modelseedpy.core import *   
 

.. toctree::
   :includehidden:
   
   core/biology_api
   core/fbahelper_api
   core/gapfillinghelper_api
   core/msatpcorrection_api
   core/msbuilder_api
   core/mseditorapi_api
   core/msgenome_api
   core/msgapfill_api
   core/msgenomeclassifier_api
   core/msgrowthphenotypes_api
   core/msmedia_api
   core/msmodel_api
   core/msmodelutils_api
   core/mstemplate_api
   core/rpcclient_api
   core/template_api
   
   
community
------------
   
The ``community`` sub-library compatibilizes transport reactions of individual models for a community model and simulates interactions amongst community members::

 from modelseedpy.community import *   


.. toctree::
   :includehidden:

   community/commkineticpkg_api
   community/dfbapkg_api
   community/mscommfitting_api
   community/mscommunity_api
   community/mscompatibility_api

   
fbapkg
------------
   
The ``fbapkg`` sub-library constrains and investigates metabolic models through Flux Balance Analysis methods::

 from modelseedpy.fbapkg import *   

.. toctree::
   :includehidden:

   fbapkg/bilevel_api
   fbapkg/changeoptpkg_api
   fbapkg/drainfluxes_api
   fbapkg/elementuptakepkg_api
   fbapkg/flexiblebiomass_api
   fbapkg/fluxfitting_api
   fbapkg/gapfillingpkg_api
   fbapkg/kbasemedia_api
   fbapkg/metabofba_api
   fbapkg/proteomefitting_api
   fbapkg/reactionuse_api
   fbapkg/tfa_api

   
ml
------------
      
The ``ml`` sub-library supports parsing features from genomes, with accessibility to KBase::

 from modelseedpy.ml import *   
   

.. toctree::
   :includehidden:
   
   ml/predict_phenotypes_api
