MSModel Package
--------------------------------------

-------------------------------------------
get_reaction_constraints_from_direction()
-------------------------------------------

A function that converts direction symbols (">" or "<") to lower and upper bound, where any other value is returned as reversible bounds:

.. code-block:: python

 lower_bound, upper_bound = get_reaction_constraints_from_direction(direction)

- *direction* ``str``: The name of the metabolite that is will be returned by the function.

**returns** *lower_bound* & *upper_bound* ``float``: The lower and upper bounds of a reaction are deduced from the parameterized direction string.

-------------------------------------------
get_direction_from_constraints()
-------------------------------------------

A function that deduces a direction symbol (">", "<", or "=") from the lower and upper reaction bounds:

.. code-block:: python

 rxn_direction = get_direction_from_constraints(lower_bound, upper_bound)

- *lower_bound* & *upper_bound* ``float``: The flux boundaries that are used to deduce the representative reaction direction symbol.

**returns** *rxn_direction* ``float``: The representative direction symbol for the reaction that is described by the lower and upper flux bounds.

-------------------------------------------
get_gpr_string()
-------------------------------------------

A function that constructs a GRP string from the parameterized iterable collection of GRPs:

.. code-block:: python

 gpr_string = get_gpr_string(gpr)

- *gpr* ``list``: The collection of GPRs that will be assembled into a string via the function.

**returns** *gpr_string* ``str``: An assembled string of all GPRs that are provided in the parameterized list.

-------------------------------------------
split_compartment_from_index()
-------------------------------------------

A function that splits an index from its associated compartment:

.. code-block:: python

 compartment, index = split_compartment_from_index(cmp_str)

- *cmp_str* ``str``: The compartment string that will be parsed.

**returns** *compartment* & *index* ``str``: The compartment and index of the parameterized compartment string, respectively.

-------------------------------------------
get_cmp_token()
-------------------------------------------

A function that filters a list of compartments to determine the compartment of greatest interest:

.. code-block:: python

 compartment_id = get_cmp_token(compartments)

- *compartments* ``list``: The collection of compartments that will be parsed.

**returns** *compartment_id* ``str``: The compartment of greatest interest among the parameterized collection of compartments.

-------------------------------------------
get_set()
-------------------------------------------

A function that filters a list of compartments to determine the compartment of greatest interest:

.. code-block:: python

 dnf_set = get_set(expression_string)

- *expression_string* ``str``: The expression string that will be parsed.

**returns** *dnf_set* ``set``: The set collection of all dnf objects from the parameterized expression string.

+++++++++++++++++++++
MSModel()
+++++++++++++++++++++

This class is a representation of ModelSEED models:

.. code-block:: python

 from modelseedpy.core import MSModel
 msmodel = MSModel(id_or_model=None, genome=None, template=None)

- *id_or_model* ``str || cobra.core.model.Model``: Either the ID of a model or the actual CobraKBase model object that will be represented through this class. The conversion from `standard COBRA models  <https://cobrapy.readthedocs.io/en/latest/autoapi/cobra/core/model/index.html>`_ to CobraKBase models is facilitated by the `cobrakbase` package. 
- *genome* ``genome object``: The representative object of the genome.
- *template* ``modelseedpy.core.mstemplate.MSTemplate``: The template of the represented model.

----------------------
template() & genome()
----------------------

**returns** *template* ``modelseedpy.core.mstemplate.MSTemplate`` & ``genome object`` : The template and genome that are associated with the model that is passed to the class instance.