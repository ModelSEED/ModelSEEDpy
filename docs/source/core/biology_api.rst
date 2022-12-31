biology
-------------------

+++++++++++++++++++++
BiologPlate()
+++++++++++++++++++++

This class assembles media from lists of compounds and HTML representations of the well media:

.. code-block:: python

 from modelseedpy.core import BiologPlate
 bioplate = BiologPlate(plate_id, rows, cols)

- *plate_id* ``str``: The ID of the plate that will be parsed to develop a simulation media.
- *rows* & *cols* ``list``: The rows and columns of the experimental system, respectively.

----------------------
add_base()
----------------------

A list of compounds are ascribed a value to construct a media:

.. code-block:: python

 bioplate.add_base(compounds, value)

- *compounds* ``list``: The compounds that will constitute a media.
- *value* ``float``: The value that will be assigned to each compound.

----------------------
get_media()
----------------------

The generate a media from a defined well ID in the set of all wells:

.. code-block:: python

 bioplate.get_media(well_id)

- *well_id* ``dict``: The collection of compounds and values that supplant default entries from the ``add_base`` function, where compounds and values are provided as the ``values`` of the "compounds" and value" ``keys``, respectively.

**returns** *media* ``dict``: The media that is created from the parameterized list of compounds and the respective value.

----------------------
_repr_html_()
----------------------

An HTML table is constructed for the plate ID, for all of the experimental wells:

.. code-block:: python

 bioplate._repr_html_()

----------------------
Accessible content
----------------------

The ``BilevelPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *id* ``str``: The plate ID that will be assessed.
- *rows* & *cols* ``list``: The rows and columns of the experimental system, respectively.
- *wells* & *base* ``dict``: Dictionaries of the experimental wells and added solution that will be constructed into media, respectively.


+++++++++++++++++++++
Biolog()
+++++++++++++++++++++

This class applies plates to the model medium:

.. code-block:: python

 from modelseedpy.fbapkg import Biolog
 biol = Biolog()

----------------------
add_plate()
----------------------

The parameterized plate is added to the collection of plates:

.. code-block:: python

 biol.add_plate(plate)

- *plate* ``Plate object``: A plate object of the experimental content.

----------------------
run_plates()
----------------------

Each well of each plate is defined as the model media, and the simulation results from that media are stored will the respective plate and well:

.. code-block:: python

 biol.run_plates(model, biomass=None, cmp='e')

- *model* ``cobra.core.model.Model``: The model whose medium will be updated with the composed media.
- *cmp* ``str``: The compartment of the exchange metabolite.

----------------------
Accessible content
----------------------

The ``BilevelPkg`` class contains a couple of accessible content that may be useful for subsequent post-processing or troubleshooting:

- *plates* ``dict``: The collection of plates (``values``) for all plate IDs (``keys``), will is updated with simulation results from each well media.
