mscommfitting
--------------------------

+++++++++++++++++++++++++
MSCommFitting()
+++++++++++++++++++++++++

This class contains the functions that load and parse experimental data, define a linear problem that represents the examined system, simulates the problem, graphs the results, and can optimize parameters to best fit the experimental system:

.. code-block:: python

 from modelseedpy.community import CommunityModelSpecies
 mscommfit = MSCommFitting()

----------------------
load_data()
----------------------

The experimental data, through multiple different channels, are imported and parsed into amenable forms for the fitting model:

.. code-block:: python

 mscommfit.load_data(community_members={}, kbase_token=None, solver='glpk', signal_tsv_paths={}, phenotypes_csv_path=None, media_conc_path=None,
                     species_abundance_path=None, carbon_conc_series={}, ignore_trials={}, ignore_timesteps=[], significant_deviation=2, zip_path=None)

- *community_members* ``dict``: the model of the community that was experimentally investigated and will be examined via fitting, which includes the permanent KBase ID of the media (e.g. 93465/3/1) that describe each respective community model.
- *kbase_token* ``str``: the KBase user token that must be provided to access permanent_media_id.
- *solver*  ``str``: the LP solver that will optimize the community model in the given media.
- *signal_tsv_paths* ``dict``: the dictionary of index names for each paths to signal TSV data that will be fitted.
- *phenotypes_csv_path* ``str``: a custom CSV of media phenotypic data.
- *media_conc_path* ``str``: a CSV of the media concentrations.
- *species_abundance_path* ``str``: a CSV over the series of species abundances for a range of trials.
- *carbon_conc_series* ``dict``: the concentrations (``values``) of carbon sources as ModelSEED IDs (``keys``) in the media.
- *ignore_trials* ``dict``: the trials (row+column coordinates) that will be ignored by the model.
- *ignore_timesteps* ``list``: timesteps that will be ignored.
- *significant_deviation* ``float``: the lowest multiple of a trial mean from its initial value that will permit its inclusion in the model fit.
- *zip_path* ``str``: specifies whether the input contents are in a zipped file.

-----------------------------
define_problem()
-----------------------------

The parsed experimental data is used to define and constrain a Global Linear Problem of the community system:

.. code-block:: python

 mscommfit.define_problem(parameters={}, zip_name=None, export_parameters=True, export_lp=True)

- *parameters* ``dict``: the parameters that will overwrite the default options. The possible key values include 

.. csv-table::
   :header: "Parameter", "Default", "Description"

   "timestep",               "the timestep that is parsed from the data",      "the timestep size of the simulation in hours"
   "cvct",               "1",      "the coefficient that penalizes phenotype conversion to the stationary phase"
   "cvcf",               "1",                "the coefficient that penalizes phenotype conversion from the stationary phase"
   "bcv",              "1",  "the highest fraction of species biomass that can convert phenotypes in a timestep"
   "cvmin",               "0",                "the lowest fraction of biomass that converts phenotypes in a single timestep"
   "v",              "0.4",  "the growth kinetics constant, which represents 1st-order kinetics"
   "carbon_sources",               "['cpd00136', 'cpd00179']",                "the ModelSEED IDs of the carbon sources in the media"
   "diffpos & diffneg",              "1",  "objective coefficients that correspond with the positive and negative differences between experimental and predicted bimoass values"

An example parameter dictionary may include the following::

.. code-block:: json

 {
            "cvct": 0.5,
            "cvcf": 0.5, 
            "v": 0.3,
            "carbon_sources": ["cpd00136"]
 }

- *zip_name* ``str``: the name of the export zip file to which content will be exported.
- *export_parameters* & *export_lp* ``bool``: specify whether the parameters and LP file will be exported.
                       
                       
----------------------
compute()
----------------------

The Linear Problem is simulated, and the primal values are parsed, optionally exported, and visualized as figures.

.. code-block:: python

 mscommfit.compute(graphs=[], zip_name=None)

- *graphs* ``list``: the graph specifications that specify which primal values will be graphed. Each list element describes a figure, with descriptive keys of ``trial``, ``content``, ``species``, ``phenotype``, and ``experimental_data``. The plots are all defined with time on the x-axis. The following graphs for the "B4" trial sample the range of possible figures that can be constructed through the API::

.. code-block:: json

 graph = [
    {
        "trial":"B4",
        "content": "biomass",
        "species": "ecoli",
        "phenotype": "*"
    },
    {
        "trial":"B4",
        "content": "biomass",
        "species": "pf",
        "phenotype": "*"
    },
    {
        "trial":"B4",
        "content": "OD",
        "experimental_data": true
    },
    {
        "trial":"B4",
        "content": "EX_cpd00179_e0",
        "species": "ecoli",
        "phenotype": "malt"
    }
 ]
 
The first and second figures will be plots of biomass for all phenotypes of *E. coli* and *P. fluorescens*, respectively. The third figure will overlay predicted community biomass and experimental OD biomass, which derives from converting OD signal to biomass from the model. The final figure displays the concentration of maltose, as its ModelSEED ID, for the maltose (``malt``) phenotype of *E. coli*.

- *zip_name* ``str``: the name of the export zip file to which content will be exported.
                       
                       
----------------------
graph()
----------------------

Primal values are visualized as figures.

.. code-block:: python

 mscommfit.compute(graphs=[], primal_values_filename=None, primal_values_zip_path=None, zip_name=None, data_timestep_hr=0.163)

- *graph* ``list``: the graph specifications that specify which primal values will be graphed, which is elaborated above for the ``compute`` function. 
- *primal_values_filename* ``str``: the name of the primal value JSON file ('primal_values.json')
- *primal_values_zip_path* ``str``: the path of the zip file that contains the primal values file
- *zip_name* ``str``: the name of the export zip file to which content will be exported.
- *data_timestep_hr* ``float``: the timestep value in hours of the data that is being graphed. This permits graphing primal values without previously simulating a model. The value is automatically overwritten by previously defined data timesteps in the ``MSCommFitting`` class object.

                       
----------------------
load_model()
----------------------

A JSON model file is imported.

.. code-block:: python

 mscommfit.load_model(mscomfit_json_path, zip_name=None, class_object=False)

- *mscomfit_json_path* ``str``: the path of the JSON model file that will be loaded and simulated. 
- *zip_name* ``str``: the path of the zip file that contains the JSON model file.
- *class_object* ``bool``: specifies whether the loaded model will be defined in the class object.
                       
**returns** *model* ``Optland.Model``: The model that is loaded via the .
 
----------------------
change_parameters()
----------------------

Primal values are visualized figures.

.. code-block:: python

 mscommfit.load_model(cvt=None, cvf=None, diff=None, vmax=None, mscomfit_json_path='mscommfitting.json', zip_name=None, class_object=False)

- *cvt*, *cvf*, *diff*, & *vmax* ``float`` or ``dict``: the parameter values that will replace existing values in the LP file. The parameters may be defined as either floats, which will be applied globally to all applicable instances in the model, or as dictionaries that defined values at specific times and possibly at specific trials for a certain time. The latter follows a dictionary structure of ``param["time"]["trial"]``, where the "trial" level can be omitted to applied a parameter value at every trial of a time. A default value can also be specified in the dictionary ``param["default"]`` that applies to times+trials that are not captured by the defined conditions.
- *mscomfit_json_path* ``str``: the path of the JSON model file that will be loaded and simulated.
- *zip_name* ``str``: the zipfile to which the edited LP JSON will be exported .


----------------------
Accessible content
----------------------

Several objects within the ``MSCommFitting`` class may be useful for subsequent post-processing or troubleshooting:

- *problem* ``Optlang.Model``: the LP model of the experimental system that is simulated.
- *carbon_conc* ``dict``: the media concentrations per substrate as defined in ``carbon_conc_series``.
- *variables* & *constraints* ``dict``: the complete collection of all variables and constraints that comprise the LP model.
