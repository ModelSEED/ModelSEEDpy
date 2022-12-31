predict_phenotype functions
--------------------------------------

-------------------------------------------
get_functional_roles()
-------------------------------------------

A function that determines the set of genome feature that match the parameterized ontological term:

.. code-block:: python

 roles = get_functional_roles(genome, ontology_term)

- *genome* ``ModelSEED Genome``: The ModelSEED Genome that will be classified.
- *ontology_term* ``str``: The ontological criteria that will assess the genome.

**returns** *roles* ``set``: The set of genome features that match the parameterized ontological term.

-------------------------------------------
get_list_functional_roles_from_kbase()
-------------------------------------------

A function that determines the set of genome feature that match the parameterized ontological term:

.. code-block:: python

 roles = get_list_functional_roles_from_kbase(genome_ref, ws_client)

- *genome_ref* ``str``: The KBase genome reference that is used to extract content from a respective KBase workspace.
- *ws_client* ``KBase object``: The KBase client that is used to acquire the genome from the KBase workspace.

**returns** *list_functional_roles* ``list``: The list of functional roles in the loaded genome.

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
