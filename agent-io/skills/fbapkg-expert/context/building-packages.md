# Building Custom FBA Packages

## Package Architecture

FBA packages add variables and constraints to the COBRA model's solver (optlang). When you call `model.optimize()`, these constraints are active.

```
┌─────────────────────────────────────────────────────────────┐
│                    Your FBA Package                          │
│                                                              │
│  ┌─────────────────┐    ┌─────────────────┐                 │
│  │   Variables     │    │   Constraints   │                 │
│  │                 │    │                 │                 │
│  │ build_variable()│    │ build_constraint│                 │
│  └────────┬────────┘    └────────┬────────┘                 │
│           │                      │                           │
│           └──────────┬───────────┘                           │
│                      │                                       │
│                      ▼                                       │
│           ┌─────────────────────┐                           │
│           │   model.solver      │                           │
│           │   (optlang)         │                           │
│           └─────────────────────┘                           │
└─────────────────────────────────────────────────────────────┘
```

## Step-by-Step: Creating a Package

### Step 1: Define the Package Class

```python
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg

class MyCustomPkg(BaseFBAPkg):
    """
    My custom FBA package for [purpose].
    """

    def __init__(self, model):
        # Call parent constructor
        BaseFBAPkg.__init__(
            self,
            model,
            "my_custom",  # Package name (used for registration)
            # Variable types: {type_name: naming_scheme}
            {
                "myvar": "reaction",     # Named by reaction.id
                "auxvar": "none"         # Named by count
            },
            # Constraint types: {type_name: naming_scheme}
            {
                "myconst": "metabolite", # Named by metabolite.id
                "bound": "string"        # Named by provided string
            }
        )
        # Initialize package-specific state
        self.my_data = {}
```

### Step 2: Implement build_package()

```python
def build_package(self, parameters):
    # Validate parameters (required list, defaults dict)
    self.validate_parameters(
        parameters,
        ["required_param"],  # Must be provided
        {
            "optional_param": 10,      # Default values
            "another_param": "default"
        }
    )

    # Access validated parameters
    threshold = self.parameters["optional_param"]

    # Build variables
    for rxn in self.model.reactions:
        if some_condition(rxn):
            self.build_variable(
                "myvar",      # Type name
                0,            # Lower bound
                1000,         # Upper bound
                "continuous", # Variable type
                rxn           # COBRA object for naming
            )

    # Build constraints
    for met in self.model.metabolites:
        # Define coefficients: {variable: coefficient}
        coef = {}
        for var_name, var in self.variables["myvar"].items():
            coef[var] = 1.0

        self.build_constraint(
            "myconst",   # Type name
            0,           # Lower bound
            threshold,   # Upper bound
            coef,        # Coefficients
            met          # COBRA object for naming
        )
```

### Step 3: Naming Schemes

The naming scheme determines how variables/constraints are named:

| Scheme | cobra_obj Parameter | Resulting Name |
|--------|-------------------|----------------|
| `"none"` | Ignored | `"1_myvar"`, `"2_myvar"`, ... |
| `"string"` | String value | `"mystring_myvar"` |
| `"reaction"` | Reaction object | `"rxn00001_c0_myvar"` |
| `"metabolite"` | Metabolite object | `"cpd00001_c0_myvar"` |

### Step 4: Variable Types

```python
# Continuous variable (default)
self.build_variable("myvar", 0, 1000, "continuous", rxn)

# Binary variable (0 or 1)
self.build_variable("binvar", 0, 1, "binary", rxn)

# Integer variable
self.build_variable("intvar", 0, 10, "integer", rxn)
```

### Step 5: Constraint Coefficients

```python
# Constraint: sum(coef[i] * var[i]) between lb and ub
coef = {
    var1: 1.0,
    var2: -2.0,
    rxn.forward_variable: 1.0,
    rxn.reverse_variable: -1.0
}
self.build_constraint("myconst", 0, 100, coef, met)
```

## Complete Example: Reaction Count Package

This package limits the number of active reactions:

```python
from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from optlang.symbolics import Zero

class ReactionCountPkg(BaseFBAPkg):
    """
    Limits the total number of active reactions.
    """

    def __init__(self, model):
        BaseFBAPkg.__init__(
            self,
            model,
            "reaction_count",
            {"active": "reaction"},    # Binary per reaction
            {"total": "none"}          # Single constraint
        )

    def build_package(self, parameters):
        self.validate_parameters(
            parameters,
            [],
            {"max_reactions": 100}
        )

        max_rxns = self.parameters["max_reactions"]

        # Add binary variable for each reaction
        for rxn in self.model.reactions:
            if rxn.id.startswith("EX_"):
                continue  # Skip exchanges

            # Binary: 1 if reaction carries flux
            var = self.build_variable("active", 0, 1, "binary", rxn)

            # Link to flux: flux <= M * active
            M = 1000  # Big M
            self.build_constraint(
                "active_upper",
                -M,  # No lower bound
                0,   # Upper bound
                {
                    rxn.forward_variable: 1,
                    rxn.reverse_variable: 1,
                    var: -M
                },
                rxn
            )

        # Total active reactions <= max
        all_active = {v: 1 for v in self.variables["active"].values()}
        self.build_constraint("total", 0, max_rxns, all_active, "total")

# Usage:
pkg = pkgmgr.getpkg("ReactionCountPkg")
pkg.build_package({"max_reactions": 50})
solution = model.optimize()
```

## Advanced: Accessing Solver Directly

For complex operations, access optlang directly:

```python
def build_package(self, parameters):
    # Get solver interface
    solver = self.model.solver

    # Create variable manually
    from optlang import Variable
    my_var = Variable("custom_name", lb=0, ub=100, type="continuous")
    solver.add(my_var)

    # Create constraint manually
    from optlang import Constraint
    my_const = Constraint(
        my_var + rxn.flux_expression,
        lb=0,
        ub=100,
        name="custom_constraint"
    )
    solver.add(my_const)

    # Update solver
    solver.update()
```

## Package Registration

Packages self-register when instantiated. The registration happens in `BaseFBAPkg.__init__`:

```python
self.pkgmgr = MSPackageManager.get_pkg_mgr(model)
self.pkgmgr.addpkgobj(self)  # Registers package
```

For custom packages not in modelseedpy:

```python
# Add to available packages
pkgmgr.available_packages["MyCustomPkg"] = MyCustomPkg

# Now getpkg works
pkg = pkgmgr.getpkg("MyCustomPkg")
```

## Best Practices

1. **Clear before rebuild**: Call `self.clear()` if build_package may be called twice

2. **Use validate_parameters**: Provides defaults and required checking

3. **Track your objects**: Variables/constraints stored in `self.variables` and `self.constraints`

4. **Name consistently**: Use COBRA object IDs when possible

5. **Document parameters**: In docstring or class comments

6. **Handle empty models**: Check list lengths before iterating

7. **Update solver**: Call `self.model.solver.update()` after complex operations
