msgenome
------------

+++++++++++++++++++++
MSFeature
+++++++++++++++++++++

A class that defines attributes of a genome feature:

.. code-block:: python

 msfeat = MSFeature(feature_id, sequence, description=None)

- *feature_id*, *sequence*, & *description* ``str``: The ID, sequence, and description of the genome feature.

------------------------------------
add_ontology_term()
------------------------------------

Returns the collection of genes for all roles of all complexes from the template reaction:

.. code-block:: python

 msfeat.add_ontology_term(ontology_term, value)

- *ontology_term* ``str``: The ontological term that will be added to the dictionary of ontological terms.
- *value* ``str``: The sequence that corresponds to the ontological term.

+++++++++++++++++++++
MSGenome
+++++++++++++++++++++

A class that edits and parses a genome:

.. code-block:: python

 msgen = MSGenome()

----------------------
from_fasta()
----------------------

``staticMethod`` Defines genome features from a FASTA file:

.. code-block:: python

 genome = msgen.from_fasta(filename, split='|', h_func=None)

- *filename* ``str``: The name of the FASTA file that will be parsed and populated in the genome features.
- *split* ``str``: The delimiter that separates the sequence ID from the sequence description.
- *h_func* ``function``: A custom function that parses the file line of a FASTA file into its sequence ID and descriptoin.

**Returns** *genome* ``modelseedpy.core.msgenome.MSGenome``: The genome whose features have been parsed from a fasta file.

------------------------------------
from_protein_sequences_hash()
------------------------------------

``staticMethod`` Defines genome features from a dictionary of protein sequences:

.. code-block:: python

 genome_class = msgen.from_protein_sequences_hash(sequences)
 
- *sequences* ``dict``: Protein sequences (``values``) for various sequence IDs (``keys``) that will be added to the genome features.

**Returns** *genome* ``modelseedpy.core.msgenome.MSGenome``: The genome whose features have been parsed from a fasta file.

-------------------
alias_hash()
-------------------

Returns the gene for each alias in each gene of the features:

.. code-block:: python

 alias_hash = msgen.alias_hash()
 
**Returns** *alias_hash* ``dict``: The collection of all alias-gene (``key``:``value``) pairs for each gene in the genome features.

-------------------
search_for_gene()
-------------------

Returns the sought gene based upon a query term of features or aliases:

.. code-block:: python

 gene = msgen.search_for_gene(query)
 
- *query* ``str``: The search term of a feature ID or gene alias.

**Returns** *gene* ``modelseedpy.core.msgenome.MSGenome``: The gene that matches the search term, where ``None`` signifies that no match was discerned.