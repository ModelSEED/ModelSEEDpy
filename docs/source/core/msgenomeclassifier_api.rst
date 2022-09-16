msgenomeclassifier
--------------------------------------

-------------------------------------------
load_classifier_from_folder()
-------------------------------------------

A function that loads the model and model features from a local folder:

.. code-block:: python

 media = load_classifier_from_folder(directory, filename)

- *directory* ``str``: The directory in which the model and features files are provided.
- *filename* ``str``: The basename of the file, where the model file is named ``<filename>``.pickle and the features file is named ``<filename>``\_features.json.

**returns** *msgenclass* ``modelseedpy.core.msgenomeclassifier.MSGenomeClassifier``: The MSGenomeClassifier object of the respective model and model features.

+++++++++++++++++++++
MSGenomeClassifier()
+++++++++++++++++++++

This class classifies a model and its features:

.. code-block:: python

 from modelseedpy.core import MSGenomeClassifier
 genclass = MSGenomeClassifier(model, model_features)

- *model* ``cobra.core.model.Model``: The CobraKBase model whose genome will be classified.
- *model_features* ``dict``: A descriptive dictionary of the investigated model.

-------------------------------------------
extract_features_from_genome()
-------------------------------------------

A function that assembles a unique list of features for the specified genome:

.. code-block:: python

 genome_features = genclass.extract_features_from_genome(genome, ontology_term)

- *genome* ``ModelSEED Genome``: The ModelSEED Genome that will be classified.
- *ontology_term* ``str``: The ontological criteria that will assess the genome.

**returns** *genome_features* ``dict``: A list of genome features ``value`` with the key of ``"genome"``.

------------------
classify()
------------------

A function that predicts FBA solutions based upon a set of genome features and indicators:

.. code-block:: python

 media = genclass.classify(genome, ontology_term='RAST')

- *genome* ``ModelSEED Genome``: The ModelSEED Genome that will be classified.
- *ontology_term* ``str``: The ontological criteria that will assess the genome.

**returns** *prediction* ``str``: The numerical prediction of the model based upon the set of genome features and indicators.
