Media Packages
--------------------------------------

+++++++++++++++++++++
MediaCompound()
+++++++++++++++++++++

This class instantiates a media compound for potential manipulations:

.. code-block:: python

 from modelseedpy.core import MediaCompound
 media_comp = MediaCompound(compound_id, lower_bound, upper_bound, concentration=None)

- *compound_id* ``str`` & *concentration* ``float``: The ID and concentration of the media compound. 
- *lower_bound* & *upper_bound* ``float``: The lower and upper bounds of the exchange reaction for the media compound, respectively. 

+++++++++++++++++++++
MSMedia()
+++++++++++++++++++++

This class instantiates a media for investigation:

.. code-block:: python

 from modelseedpy.core import MSMedia
 msmedia = MSMedia(media_id)

- *media_id* ``str``: The ID of the investigated media. 

-------------------------------------------
from_dict()
-------------------------------------------

A function that converts a media dictionary into a media object:

.. code-block:: python

 media = msmedia.from_dict(media_dictionary)

- *media_dictionary* ``dict``: A dictionary representation of the media that contains either a list of exchange bounds or the uptake magnitude (``value``) for each compound in the media.

**returns** *media* ``modelseedpy.core.msmedia.MSMedia``: The media that is constructed from the dictionary format.

-------------------------------------------
get_media_constraints()
-------------------------------------------

A function that assigns a compartment to each compound in the media:

.. code-block:: python

 media = msmedia.get_media_constraints(cmp='e0')

- *cmp* ``str``: The compartment suffix that will be appended to all compounds in the media, while ``cmp`` is not ``None``.

**returns** *media* ``dict``: A dictionary of lower and upper bounds (``values``) for each metabolite ID (``key``).

-----------------
merge()
-----------------

A function that expands the media compounds according to amendments from the aforementioned ``MediaCompound`` object:

.. code-block:: python

 msmedia.merge(media, overwrite_overlap=False)

- *media* ``modelseedpy.core.msmedia.MSMedia``: The media whose compound compartments have been defined to the parameterized compartment.
- *overwrite_overlap* ``bool``: specifies whether existing metabolite IDs will be redefined by the ``modelseedpy.core.msmedia.MediaCompound`` compound object.