dFBA Package
--------------

++++++++++
dFBAPkg()
++++++++++

This class defines and executes dynamic FBA simulations of COBRA models:

.. code-block:: python

 from modelseedpy.fbapkg import dFBAPkg
 dfba = dFBAPkg(model, modelseed_db_path, solver = 'glpk', warnings = True, verbose = False, printing = False, jupyter = False)

- *model* ``cobra.core.model.Model``: the CobraKBase model that will be simulated. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *modelseed_db_path* ``str``: specifies the path to a local version of the ModelSEED Database.
- *solver* ``str``: specifies which linear programmating algorithm will be used to simulate the FBA model. The `glpk` solver is selected by default since it is free and universally accessible.
- *warnings*, *verbose*, & *printing* ``bool``: specifies whether simulation warnings, details and calculations, or results will be printed, respectively. These options are valuable for troubleshooting.
- *jupyter* ``bool``: specifies whether simulation is being conducted in a Jupyter notebook, in which case the printed DataFrames will be expressed with the ``display()`` function. 

           
----------------------
simulate()
----------------------

A cobrakbase model is simulated with the parameterized kinetics data over the defined time and conditions:

.. code-block:: python

 dfba.simulate(kinetics_path = None, initial_concentrations_M: dict = {}, total_time = 200, timestep = 20, export_name = None, 
               export_directory = None, kinetics_data = {}, temperature = 25, p_h = 7, cellular_dry_mass_fg = 222, cellular_fL = 1, 
               figure_title = 'Metabolic perturbation', included_metabolites = [], labeled_plots = True, visualize = True, export = True)

- *kinetics_path* & *kinetics_data* ``str`` & ``dict``: either the path to a `JSON` file that can be imported or a dictionary argument that provide the kinetics data that will constrain the model in the simulation. The `JSON` structure in both means of providing the data, and possesses the following nesting: ``<reaction>`` -> ``<source>`` -> ``substituted_rate_law`` + ``initial_concentrations_M`` + optional keys (such as ``metadata``). The ``substituted_rate_law`` key must contain the mathematically valid rate law expression as a string that can be evaluated through the `eval() built-in function <https://pythongeeks.org/python-eval-function/>`_ of Python, where each metabolite in the rate law is represented by a **single letter variable** that is defined with a `ModelSEED Compound ID <https://modelseed.org/biochem/compounds>`_ in the ``met_id`` key. The concentrations of the ``initial_concentrations_M`` key must be provided in units of Molar and defined for the rate law variable letters.

.. code-block:: json

  {
    "R_3OAS140": {
        "source_1": {
            "substituted_rate_law": "(68.0*A*B)/(50.0*0.34*C+360.0*B+0.34*A+A*B*C)",
            "initial_concentrations_M": {
                "A": 0.0200,
                "C": 0.022,
                "B": 0.0014
            },
            "met_id": {
                "A": "cpd11468",
                "B": "cpd00067",
                "C": "cpd11492"
            }
        }
    },
    "rxn2": {
        "source_1": {
            "substituted_rate_law": "(A*B)/(50.0*0.34*C+3*B+0.34*A+C)",
            "initial_concentrations_M": {
                "A": 0.0200,
                "C": 0.022,
                "B": 0.0012
            },
            "met_id": {
                "A": "cpd11468",
                "B": "cpd00011",
                "C": "cpd11492"
            }
        }
    }
  }
       
The additional keys can provide provenance of the datum source:
            
.. code-block:: json

 {
    "2-Oxogluterate dehydrogenase": {
        "55199": {
            "RateLaw": "Vmax*S/(Km+S)",
            "initial_concentrations_M": {
                "S": 1.6e-08
            },
            "metadata": {
                "Buffer": "[50 mm Mops, 8 mm TCEP, 50 mm Mops, 8 mm TCEP]",
                "Enzyme Variant": "wildtype",
                "KineticMechanismType": "Michaelis-Menten",
                "Organism": "Pisum sativum",
                "Pathway": null,
                "Product": "NADH;H+;Oxidized N-alpha-(benzyloxycarbonyl)-N-omega-(D,L-1,2-dithiolane-3-pentanoyl)-L-lysine",
                "Publication": "Neuburger M, Polidori AM, Pi\u00e8tre E, Faure M, Jourdain A, Bourguignon J, Pucci B, Douce R: Interaction between the lipoamide-containing H-protein and the lipoamide dehydrogenase (L-protein) of the glycine decarboxylase multienzyme system. 1. Biochemical studies., Eur J Biochem 2000 (267) , 2882-9",
                "Temperature": "30.0",
                "annotations": {
                    "ECNumber": "1.8.1.4",
                    "KeggReactionID": null,
                    "PubMedID": 10806385.0,
                    "SabioReactionID": 13969
                },
                "pH": "7.5",
                "reaction_string": " <-> Nicotinamide adenine dinucleotide-reduced + H+"
            },
            "substituted_parameters": {
                "Km": {
                    "comment": "-",
                    "deviat.": "10",
                    "end val.": "-",
                    "species": "N-alpha-(benzyloxycarbonyl)-N-omega-(D,L-1,2-dithiolane-3-pentanoyl)-L-lysine",
                    "start val.": "170.0",
                    "type": "Km",
                    "unit": "\u00b5M"
                },
                "Vmax": {
                    "comment": "-",
                    "deviat.": "7",
                    "end val.": "-",
                    "species": "-",
                    "start val.": "90.0",
                    "type": "Vmax",
                    "unit": "nmol/min"
                }
            },
            "substituted_rate_law": "1.5000000000000002e-09*S/(0.00016999999999999999+S)",
            "variables_molar": {
                "Km": "0.00016999999999999999",
                "Vmax": "1.5000000000000002e-09"
            },
            "variables_name": {
                "Km": "N-alpha-(benzyloxycarbonyl)-N-omega-(D,L-1,2-dithiolane-3-pentanoyl)-L-lysine",
                "S": "N-alpha-(benzyloxycarbonyl)-N-omega-(D,L-1,2-dithiolane-3-pentanoyl)-L-lysine",
                "Vmax": "-"
            }
        }
    }
 }

- *initial_concentrations_M* ``dict``: specifies initial concentrations of the simulated metabolites that supplant values from the kinetics data. Every metabolite in this dictionary must be defined in the model, and the concentrations in units of molar must be assigned according to the metabolite's ModelSEED Compound ID:
           
.. code-block:: json

 {
    "cpd00002":0.0200, 
    "cpd00008":0.0014
 }
 
- *total_time* & *timestep* ``float``: specify the total time and the timstep of the simulation in minutes.
- *export_name* & *export_directory* ``str``: specify the folder name and directory to which the simulation content will be exported, where `None` defaults to a unique folder name in the current working directory.
- *temperature* & *p_h* ``float``: optionally specify the temperature and pH at which the simulation will occur, respective, which allows the most closely matched datum to be parameterized, where multiple datum exist for the same reaction.
- *cellular_dry_mass_fg* & *cellular_fL* ``float``: The `dry mass <https://doi.org/10.1101/2021.12.30.474524>`_ and `volume <https://doi.org/10.1128/AEM.00117-14>`_ of the simulated cell, in base units of femto- grams and liters, respectively.  These values can be sourced from literature, and the standard values may approximate prokaryotic cells.
- *figure_title*, *included_metabolites*, & *labeled_plots* ``str``, ``list``, & ``bool``: specify the title of the simulation Figure, the metabolites that will be plotted in the simulation Figure, and where each plot will be labeled with text to clarify its identity.
- *visualize* & *export* ``bool``: specifies whether the simulation results will be visually depicted or exported to a specified folder, respectively.



----------------------
Accessible content
----------------------

A multitude of values are stored within the ``dFBA`` object, and can be subsequently used in a workflow. The complete list of content within the ``dFBA`` object can be printed through the built-in ``dir()`` function in the following example sequence:

.. code-block:: python

 # conduct a dFBA simulation
 from dfbapy import dFBA
 dfba = dFBA(model)
 dfba.simulate(reaction_kinetics, None, total_time, timestep)
 
 # evaluate the dFBA simulation contents
 print(dir(dfba))

The following list highlights stored content in the ``dFBA`` object after a simulation:

- *model* ``cobra.core.model.Model``: The cobrakbase model that is simulated.
- *concentrations* & *fluxes* ``pandas.core.frame.DataFrame``: `Pandas DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/frame.html>`_ that contain the ``mM`` concentrations (or changes thereof) for each metabolite and ``mmol/g_(dw)/hr`` fluxes for each reaction, respectively.
- *kinetics_data* ``dict``: A dictionary of the kinetics data constrains the model.
- *timestep_value* ``float``: The simulation timestep in minutes.
- *compound_ids* ``dict``: A dictionary of all ModelSEED IDs with their names as values, which is loaded from the ModelSEED Database via the parameterized path. 
- *cell_dry_mass* & *cell_liters* ``float``: The mass and volume of the simulated cell.
- *changed* & *unchanged* ``set``: The exclusive sets of metabolites whose concentrations either changed or did not change over the simulation, respectively.
- *constrained* ``OrderedDict``: A dictionary with reaction names as the keys and their respective kinetic constraints as the values.
- *solutions* ``list``: A list of the Cobra solutions from each timestep that constitute the columns of the `fluxes` DataFrame.