# FBA Packages Reference

## Core Infrastructure

### MSPackageManager
**File:** `fbapkg/mspackagemanager.py`

Central registry for FBA packages. Singleton per model.

```python
from modelseedpy.fbapkg import MSPackageManager

# Get or create manager
pkgmgr = MSPackageManager.get_pkg_mgr(model)

# List packages
pkgmgr.list_available_packages()  # All known packages
pkgmgr.list_active_packages()     # Currently loaded

# Get a package (creates if missing)
pkg = pkgmgr.getpkg("KBaseMediaPkg")

# Get without creating
pkg = pkgmgr.getpkg("KBaseMediaPkg", create_if_missing=False)
```

### BaseFBAPkg
**File:** `fbapkg/basefbapkg.py`

Base class all packages inherit from.

**Constructor Parameters:**
- `model` - cobra.Model or MSModelUtil
- `name` - Package name string
- `variable_types` - Dict mapping type names to naming schemes
- `constraint_types` - Dict mapping type names to naming schemes

**Key Methods:**
- `build_package(params)` - Override to add constraints
- `build_variable(type, lb, ub, vartype, cobra_obj)` - Create variable
- `build_constraint(type, lb, ub, coef, cobra_obj)` - Create constraint
- `clear()` - Remove all variables/constraints
- `validate_parameters(params, required, defaults)` - Check params

---

## Media and Exchange Packages

### KBaseMediaPkg
**File:** `fbapkg/kbasemediapkg.py`

Sets exchange reaction bounds based on media definition.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `media` | MSMedia | None | Media object |
| `default_uptake` | float | 0 | Default uptake bound |
| `default_excretion` | float | 100 | Default excretion bound |

**Example:**
```python
pkg = pkgmgr.getpkg("KBaseMediaPkg")
pkg.build_package({
    "media": media,
    "default_uptake": 0,
    "default_excretion": 100
})
# Or shorthand:
pkg.build_package(media)
```

### ElementUptakePkg
**File:** `fbapkg/elementuptakepkg.py`

Constrains total uptake of a specific element (e.g., carbon).

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `element` | str | "C" | Element to constrain |
| `max_uptake` | float | 10 | Maximum uptake rate |

---

## Gapfilling Packages

### GapfillingPkg
**File:** `fbapkg/gapfillingpkg.py` (~1200 lines)

MILP formulation for gapfilling. Adds reactions from templates and penalizes additions.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `default_gapfill_templates` | list | [] | Templates to add reactions from |
| `minimum_obj` | float | 0.01 | Minimum objective value |
| `reaction_scores` | dict | {} | Penalty scores per reaction |
| `blacklist` | list | [] | Reactions to exclude |
| `model_penalty` | float | 1 | Penalty for model reactions |
| `auto_sink` | list | [...] | Compounds to add sinks for |

**Variables Added:**
- `rmaxf` (reaction) - Max reverse flux
- `fmaxf` (reaction) - Max forward flux

**Constraints Added:**
- `rmaxfc` (reaction) - Reverse flux coupling
- `fmaxfc` (reaction) - Forward flux coupling

**Example:**
```python
pkg = pkgmgr.getpkg("GapfillingPkg")
pkg.build_package({
    "default_gapfill_templates": [template],
    "minimum_obj": 0.1,
    "reaction_scores": {"rxn00001": 0.5}
})
```

---

## Biomass Packages

### FlexibleBiomassPkg
**File:** `fbapkg/flexiblebiomasspkg.py`

Allows biomass composition to vary within bounds.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bio_rxn_id` | str | "bio1" | Biomass reaction ID |
| `flex_coefficient` | float | 0.1 | Flexibility (fraction) |
| `use_rna_class` | bool | True | Group RNA components |
| `use_protein_class` | bool | True | Group protein components |
| `use_dna_class` | bool | True | Group DNA components |

**Example:**
```python
pkg = pkgmgr.getpkg("FlexibleBiomassPkg")
pkg.build_package({
    "bio_rxn_id": "bio1",
    "flex_coefficient": 0.2  # 20% flexibility
})
```

---

## Thermodynamic Packages

### SimpleThermoPkg
**File:** `fbapkg/simplethermopkg.py`

Simple thermodynamic constraints (loopless FBA variant).

**Example:**
```python
pkg = pkgmgr.getpkg("SimpleThermoPkg")
pkg.build_package()
```

### FullThermoPkg
**File:** `fbapkg/fullthermopkg.py`

Full thermodynamic constraints with concentration variables.

**Variables Added:**
- `logconc` (metabolite) - Log concentration variables
- `dGrxn` (reaction) - Reaction Gibbs energy

---

## Reaction Control Packages

### ReactionUsePkg
**File:** `fbapkg/reactionusepkg.py`

Binary variables indicating whether reactions carry flux.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reaction_list` | list | [] | Reactions to add binaries for |

**Variables Added:**
- `use` (reaction) - Binary: 1 if reaction active

**Example:**
```python
pkg = pkgmgr.getpkg("ReactionUsePkg")
pkg.build_package({
    "reaction_list": model.reactions
})

# Access variables
for rxn_id, var in pkg.variables["use"].items():
    print(f"{rxn_id} active: {var.primal}")
```

### RevBinPkg
**File:** `fbapkg/revbinpkg.py`

Binary variables for reaction direction.

**Variables Added:**
- `revbin` (reaction) - Binary: 1 if forward, 0 if reverse

### ReactionActivationPkg
**File:** `fbapkg/reactionactivationpkg.py`

Activate/deactivate reactions based on expression data.

---

## Objective Packages

### ObjectivePkg
**File:** `fbapkg/objectivepkg.py`

Manage model objective function.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `objective` | str/Reaction | model.objective | Target reaction |
| `maximize` | bool | True | Maximize or minimize |

### ObjConstPkg
**File:** `fbapkg/objconstpkg.py`

Convert objective to constraint (for multi-objective).

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `objective_value` | float | - | Fix objective at value |

### TotalFluxPkg
**File:** `fbapkg/totalfluxpkg.py`

Minimize total flux (parsimonious FBA).

**Example:**
```python
pkg = pkgmgr.getpkg("TotalFluxPkg")
pkg.build_package()
# Now optimize minimizes total flux
```

### ChangeOptPkg
**File:** `fbapkg/changeoptpkg.py`

Change the solver/optimizer.

---

## Data Fitting Packages

### FluxFittingPkg
**File:** `fbapkg/fluxfittingpkg.py`

Fit model to measured flux data.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flux_data` | dict | {} | {rxn_id: measured_flux} |

### ProteomeFittingPkg
**File:** `fbapkg/proteomefittingpkg.py`

Fit model to proteome data.

### MetaboFBAPkg
**File:** `fbapkg/metabofbapkg.py`

Integrate metabolomics data.

---

## Utility Packages

### DrainFluxPkg
**File:** `fbapkg/drainfluxpkg.py`

Add drain reactions for specific metabolites.

**Parameters:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `metabolites` | list | [] | Metabolites to drain |

### ProblemReplicationPkg
**File:** `fbapkg/problemreplicationpkg.py`

Create multiple copies of the FBA problem.

### BilevelPkg
**File:** `fbapkg/bilevelpkg.py`

Bilevel optimization formulation.

---

## Package Interactions

### Packages That Work Together
- `KBaseMediaPkg` + any other (media is usually first)
- `GapfillingPkg` + `KBaseMediaPkg` (gapfilling needs media)
- `ReactionUsePkg` + `TotalFluxPkg` (minimize active reactions)
- `SimpleThermoPkg` or `FullThermoPkg` (not both)

### Order of Building
1. `KBaseMediaPkg` (set exchange bounds first)
2. Constraint packages (thermo, element uptake)
3. Objective packages
4. Analysis packages (gapfilling, fitting)

### Clearing Packages
```python
# Clear specific package
pkg.clear()

# Packages track their own variables/constraints
# clear() only removes what that package added
```
