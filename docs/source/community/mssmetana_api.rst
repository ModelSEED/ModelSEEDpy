mssmetana
--------------------------

+++++++++++++++++++++
MSSmetana()
+++++++++++++++++++++

This class captures the scoring metrics of community interactions that were introduced by Zelezniak et al., 2015 (https://doi.org/10.1073/pnas.1421834112). Both class methods, which support multiple scores to be efficiently calculated for a single community of members, and static methods, which supports in-line calculation of single scores for a community ad hoc, are available for each score. The code will be expanded in the future to acommodate ``MSCommunity`` community objects:

.. code-block:: python

 from modelseedpy.community import MSSmetana
 smtna = MSSmetana(cobra_models, com_model, min_growth=0.1, n_solutions=100, 
                 abstol=1e-3, media_dict=None, printing=True, compatibilize=True, minimize_flux=True)

- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *com_model* ``cobra.core.model.Model``: the community model that combines the individual member models.
- *min_growth* ``float``: the minimal permissible community biomass objective value from simulation.
- *n_solutions* ``int``: the number of loops over which the MILP algorithms of the SMETANA scores search.
- *abstol* ``float``: the minimum flux above which the flux is considered to be non-zero.
- *media_dict* ``dict``: the minimal media of the community and its members, which adheres to the following structure:

.. code-block:: json

 {
     "community" : {
          "EX_cpd00149_e0": 0.0004692000000265238,
          "EX_cpd00205_e0": 0.019035799999983283,
          "EX_cpd01017_e0": 591.3730964939617,
          "EX_cpd00971_e0": 0.0004668000001402106,
          "EX_cpd00034_e0": 0.0004992000000265238,
          "EX_cpd00098_e0": 10.500407471304925,
          "EX_cpd00030_e0": 0.0005326000000265237,
          "EX_cpd00007_e0": 1000.0000223,
          "EX_cpd00254_e0": 0.0012921000000005733,
          "EX_cpd00063_e0": 0.0009620000000722939,
          "EX_cpd00166_e0": 0.0004891000000265238,
          "EX_cpd00244_e0": 3.07e-05,
          "EX_cpd03725_e0": 554.212955031386,
          "EX_cpd01570_e0": 138.92131569784812,
          "EX_cpd00099_e0": 0.0004952000000457701,
          "EX_cpd11574_e0": 2.1000000000000002e-06,
          "EX_cpd00058_e0": 6.739999999999998e-05
       },
     "members": {
         "Bacteroides_thetaiotaomicron_VPI-5482.fbamdl.23": {
             "media": {
                "EX_cpd00149_e0": 0.0004668000000265238,
                "EX_cpd00205_e0": 0.0004668000000265238,
                "EX_cpd01017_e0": 591.3730964939617,
                "EX_cpd00971_e0": 0.0004668000001402106,
                "EX_cpd00009_e0": 10.676809471303557,
                "EX_cpd00034_e0": 0.0004668000000265238,
                "EX_cpd00028_e0": 0.0004668000000265238,
                "EX_cpd00098_e0": 10.500407471304925,
                "EX_cpd00030_e0": 0.0004668000000265238,
                "EX_cpd00007_e0": 1000.0,
                "EX_cpd10515_e0": 0.0009336000000000001,
                "EX_cpd00254_e0": 0.0004668000000265238,
                "EX_cpd00063_e0": 0.0004668000000265238,
                "EX_cpd00166_e0": 0.0004668000000265238
             },
             "solution": {
                "12ETHDt_c0": -0.0018672000000000048,
                "12PPDRt_c0": 0.0,
                "12PPDt_c0": 0.0,
                "et cetera": 123,
                }
        },
        "iML1515": {
             "media": {
                "EX_cpd00244_e0": 3.07e-05,
                "EX_cpd00541_e0": 3.0000000000000004e-07,
                "EX_cpd00104_e0": 2e-07,
                "EX_cpd03725_e0": 554.212955031386,
                "EX_cpd01570_e0": 138.92131569784812,
                "EX_cpd00007_e0": 2.2300000000002873e-05,
                "EX_cpd00254_e0": 0.0008252999999740496,
                "EX_cpd00063_e0": 0.0004952000000457701,
                "EX_cpd00166_e0": 2.23e-05,
                "EX_cpd00099_e0": 0.0004952000000457701,
                "EX_cpd00205_e0": 0.01856899999995676,
                "EX_cpd00048_e0": 0.025150899998379828,
                "EX_cpd00149_e0": 2.4000000000000003e-06,
                "EX_cpd11574_e0": 2.1000000000000002e-06,
                "EX_cpd00034_e0": 3.24e-05,
                "EX_cpd00030_e0": 6.58e-05,
                "EX_cpd00058_e0": 6.739999999999998e-05
             },
             "solution": {
                "34dhpactex_e0": 0.0,
                "GUAtex_e0": 0.0,
                "rxn01256_c0": 0.03142369999992683,
                "rxn00411_c0": 3.64853774359535,
                "et cetera": 123,
    }
 }

The ``"community"`` key represents the community model that is being investigated, while the ``"members"`` key represents individual members of the community. The former contains only the minimal media of the community through key:value pairings of exchange reactions and their respective fluxes, with (+) denoting influx. The latter is organized by each member species, which contains both the minimal media and the full solution of the members in their minimal media, with (+) denoting influx and outflux in the minimal media and optimization solution, respectively.

- *printing* ``bool``: specifies whether debugging checkpoints are printed to the console.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.
- *minimize_flux* ``bool``: specifies whether the ``minimal_media`` will be determined via the minimal flux method, where ``False`` utilizes the ``minimal_components`` method.

----------------------------------------------------------------------------------------
mro_score(), mip_score(), mu_score(), mp_score(), sc_score(), smetana_score()
----------------------------------------------------------------------------------------

The individual SMETANA scores can be succinctly calculated in any order from the aforementioned class object, without the need for further parameters:

.. code-block:: python

 mro = smtna.mro_score()
 mip = smtna.mip_score()
 mu = smtna.mu_score()
 mp = smtna.mp_score()
 sc = smtna.sc_score()
 smetana = smtna.smetana_score()
 
 **returns** the respective score of the defined community system:

- *mro* & *mip* ``float``: The numerous scores from the MRO and MIP scores, respectively. 
- *mu*, *mp*, *sc*, & *smetana*  ``dict``: The collections of scores, organized by model IDs, for the MU, MP, SC, and SMETANA scores, respectively.
       
-----------------------------
Attributes
-----------------------------

The ``MSSmetana`` class object stores numerous attributes:
 
- *mro* & *mip* ``float``: The numerous scores from the MRO and MIP scores, respectively. 
- *mu*, *mp*, *sc*, & *smetana*  ``dict``: The collections of scores, organized by model IDs, for the MU, MP, SC, and SMETANA scores, respectively.
- *models* ``list o cobra.core.model.Model``: The collection of member models that comprise the community.
- *community* ``cobra.core.model.Model``: The community model.
- *media* ``dict``: The media object of the community.
- *printing* ``bool``: The setting for whether results of the alignment functionality, respectively, are printed to the console.

      
---------
mro()
---------

**Staticmethod**

The MRO SMETANA score can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.mro(cobra_models, min_growth=0.1, media_dict=None, compatibilize=True)

- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *min_growth* ``float``: the minimal permissible community biomass objective value from simulation.
- *media_dict* ``dict``: A dictionary of predetermined minimal media, per the above definition.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.


--------
mip()
--------

**Staticmethod**

The MIP SMETANA score can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.mip(com_model, cobra_models, min_growth=0.1, interacting_media_dict=None,
            noninteracting_media_dict=None, compatibilize=True)

- *com_model* ``cobra.core.model.Model``: the community model that combines the individual member models.
- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *min_growth* ``float``: the minimal permissible community biomass objective value from simulation.
- *interacting_media_dict* & *noninteracting_media_dict* ``dict``: Dictionaries of the predetermined minimal media that include and exclude cross-feeding (syntrophy), respectively. The MIP formulation essentially compares these two media, hence the calculation can be tremenedously expedited if both of these media objects are parameterized and need not be calculated in the logic.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.
       
       
---------
mu()
---------

**Staticmethod**

The MU SMETANA score can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.mu(cobra_models, n_solutions=100, abstol=1e-3, compatibilize=True)

- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *n_solutions* ``int``: the number of loops over which the MILP algorithms of the SMETANA scores search.
- *abstol* ``float``: the minimum flux above which the flux is considered to be non-zero.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.
       
---------
mp()
---------

**Staticmethod**

The MP SMETANA score can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.mp(cobra_models=None, com_model=None, abstol=1e-3, compatibilize=True)
       
- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *com_model* ``cobra.core.model.Model``: the community model that combines the individual member models.
- *abstol* ``float``: the minimum flux above which the flux is considered to be non-zero.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.
       
---------
sc()
---------

**Staticmethod**

The SC SMETANA score can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.sc(cobra_models=None, com_model=None, min_growth=0.1, 
              n_solutions=100, abstol=1e-3, compatibilize=True)
       
- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *com_model* ``cobra.core.model.Model``: the community model that combines the individual member models.
- *min_growth* ``float``: the minimal permissible community biomass objective value from simulation.
- *n_solutions* ``int``: the number of loops over which the MILP algorithms of the SMETANA scores search.
- *abstol* ``float``: the minimum flux above which the flux is considered to be non-zero.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.


-----------
smetana()
-----------

**Staticmethod**

The smetana SMETANA superscore can be specifically calculated without constructing a class object:

.. code-block:: python

 MSSmetana.smetana(cobra_models, com_model=None, min_growth=0.1, n_solutions=100, abstol=1e-6,
                prior_values=None, compatibilize=False, sc_coupling=False)

- *cobra_models* ``list o cobra.core.model.Model``: the collection of member models that comprise the examined community. 
- *com_model* ``cobra.core.model.Model``: the community model that combines the individual member models.
- *min_growth* ``float``: the minimal permissible community biomass objective value from simulation.
- *n_solutions* ``int``: the number of loops over which the MILP algorithms of the SMETANA scores search.
- *abstol* ``float``: the minimum flux above which the flux is considered to be non-zero.
- *prior_values* ``Iterable``: The collection of ``SC``, ``MU``, and ``MP`` score results that were previously calculated for the studied system, and thus do not need to be recalculated.
- *compatibilize* ``bool``: specifies whether the member models will be standardized to the ModelSEED Database conventions.
- *sc_coupling* ``bool``: specifies whether the SC score contributes to the calculation of the smetana score.
