# MSModelUtil Developer Guide

## Overview

`MSModelUtil` is the central utility wrapper class in ModelSEEDpy that encapsulates a COBRApy `Model` object and provides extensive model manipulation, analysis, and FBA-related functionality. It serves as the primary bridge between COBRApy models and ModelSEED-specific functionality.

**Location:** `modelseedpy/core/msmodelutl.py` (~2,000 lines)

## Architecture

### Design Pattern

MSModelUtil uses a **singleton-like caching pattern** where instances are cached by model object:

```python
class MSModelUtil:
    mdlutls = {}  # Static cache of MSModelUtil instances

    @staticmethod
    def get(model, create_if_missing=True):
        """Get or create MSModelUtil for a model"""
        if isinstance(model, MSModelUtil):
            return model
        if model in MSModelUtil.mdlutls:
            return MSModelUtil.mdlutls[model]
        elif create_if_missing:
            MSModelUtil.mdlutls[model] = MSModelUtil(model)
            return MSModelUtil.mdlutls[model]
        return None
```

This means you can safely call `MSModelUtil.get(model)` multiple times and always get the same instance.

### Core Dependencies

```
MSModelUtil
    ├── cobra.Model (wrapped object)
    ├── MSPackageManager (FBA constraint packages)
    ├── ModelSEEDBiochem (reaction/compound database)
    ├── FBAHelper (FBA utility functions)
    └── MSATPCorrection (lazy-loaded for ATP analysis)
```

### Key Instance Attributes

```python
self.model              # The wrapped cobra.Model
self.pkgmgr             # MSPackageManager for this model
self.wsid               # KBase workspace ID (if applicable)
self.atputl             # MSATPCorrection instance (lazy-loaded)
self.gfutl              # MSGapfill reference (set by gapfiller)
self.metabolite_hash    # Metabolite lookup cache
self.search_metabolite_hash  # Fuzzy search cache
self.test_objective     # Current test objective value
self.reaction_scores    # Gapfilling reaction scores
self.integrated_gapfillings  # List of integrated gapfilling solutions
self.attributes         # Model metadata dictionary
self.atp_tests          # Cached ATP test conditions
self.reliability_scores # Reaction reliability scores
```

---

## API Reference

### Initialization & Factory Methods

#### `MSModelUtil(model)`
Create a new MSModelUtil wrapping a cobra.Model.

```python
from modelseedpy.core.msmodelutl import MSModelUtil
import cobra

model = cobra.io.load_json_model("my_model.json")
mdlutl = MSModelUtil(model)
```

#### `MSModelUtil.get(model, create_if_missing=True)` [static]
Get or create an MSModelUtil instance. Preferred method for obtaining instances.

```python
# These are equivalent and return the same instance:
mdlutl1 = MSModelUtil.get(model)
mdlutl2 = MSModelUtil.get(model)
assert mdlutl1 is mdlutl2  # True

# Also accepts MSModelUtil directly (returns it unchanged)
mdlutl3 = MSModelUtil.get(mdlutl1)
assert mdlutl3 is mdlutl1  # True
```

#### `MSModelUtil.from_cobrapy(filename)` [static]
Load a model from a file and wrap it in MSModelUtil.

```python
# Supports .json and .xml/.sbml files
mdlutl = MSModelUtil.from_cobrapy("model.json")
mdlutl = MSModelUtil.from_cobrapy("model.xml")
```

#### `MSModelUtil.build_from_kbase_json_file(filename, kbaseapi)` [static]
Load a model from KBase JSON format.

```python
from modelseedpy.core.kbaseapi import KBaseAPI
kbaseapi = KBaseAPI()
mdlutl = MSModelUtil.build_from_kbase_json_file("kbase_model.json", kbaseapi)
```

---

### I/O Methods

#### `save_model(filename, format="json")`
Save the model to a file.

```python
mdlutl.save_model("output.json", format="json")
mdlutl.save_model("output.xml", format="xml")
```

#### `printlp(model=None, path="", filename="debug", print=False)`
Write the LP formulation to a file for debugging.

```python
mdlutl.printlp(print=True)  # Writes debug.lp
mdlutl.printlp(path="/tmp", filename="mymodel", print=True)
```

#### `print_solutions(solution_hash, filename="reaction_solutions.csv")`
Export multiple FBA solutions to CSV.

```python
solutions = {
    "glucose": model.optimize(),
    "acetate": model.optimize()  # after changing media
}
mdlutl.print_solutions(solutions, "flux_comparison.csv")
```

---

### Metabolite Search & Lookup

#### `find_met(name, compartment=None)`
Find metabolites by name, ID, or annotation. Returns a list of matching metabolites.

```python
# Find by ModelSEED ID
mets = mdlutl.find_met("cpd00001")  # Water

# Find by name
mets = mdlutl.find_met("glucose")

# Find in specific compartment
mets = mdlutl.find_met("cpd00001", "c0")  # Cytosolic water
mets = mdlutl.find_met("cpd00001", "e0")  # Extracellular water
```

#### `metabolite_msid(metabolite)` [static]
Extract the ModelSEED compound ID from a metabolite.

```python
msid = MSModelUtil.metabolite_msid(met)  # Returns "cpd00001" or None
```

#### `reaction_msid(reaction)` [static]
Extract the ModelSEED reaction ID from a reaction.

```python
msid = MSModelUtil.reaction_msid(rxn)  # Returns "rxn00001" or None
```

#### `msid_hash()`
Create a dictionary mapping ModelSEED IDs to metabolite lists.

```python
id_hash = mdlutl.msid_hash()
# id_hash["cpd00001"] = [<Metabolite cpd00001_c0>, <Metabolite cpd00001_e0>]
```

#### `build_metabolite_hash()`
Build internal metabolite lookup caches. Called automatically by `find_met()`.

```python
mdlutl.build_metabolite_hash()
# Now mdlutl.metabolite_hash and mdlutl.search_metabolite_hash are populated
```

#### `search_name(name)` [static]
Normalize a name for fuzzy searching (lowercase, strip compartment suffix, remove non-alphanumeric).

```python
MSModelUtil.search_name("D-Glucose_c0")  # Returns "dglucose"
```

---

### Reaction Search & Analysis

#### `rxn_hash()`
Create a dictionary mapping reaction stoichiometry strings to reactions.

```python
rxn_hash = mdlutl.rxn_hash()
# rxn_hash["cpd00001_c0+cpd00002_c0=cpd00003_c0"] = [<Reaction>, 1]
```

#### `find_reaction(stoichiometry)`
Find a reaction by its stoichiometry.

```python
stoich = {met1: -1, met2: -1, met3: 1}
result = mdlutl.find_reaction(stoich)
# Returns [reaction, direction] or None
```

#### `stoichiometry_to_string(stoichiometry)` [static]
Convert stoichiometry dict to canonical string representation.

```python
strings = MSModelUtil.stoichiometry_to_string(rxn.metabolites)
# Returns ["reactants=products", "products=reactants"]
```

#### `exchange_list()`
Get all exchange reactions (EX_ or EXF prefixed).

```python
exchanges = mdlutl.exchange_list()
```

#### `exchange_hash()`
Create a dictionary mapping metabolites to their exchange reactions.

```python
ex_hash = mdlutl.exchange_hash()
# ex_hash[<Metabolite cpd00001_e0>] = <Reaction EX_cpd00001_e0>
```

#### `nonexchange_reaction_count()`
Count non-exchange reactions that have non-zero bounds.

```python
count = mdlutl.nonexchange_reaction_count()
```

#### `is_core(rxn)`
Check if a reaction is a core metabolic reaction.

```python
if mdlutl.is_core("rxn00001_c0"):
    print("This is a core reaction")
```

---

### Exchange & Transport Management

#### `add_exchanges_for_metabolites(cpds, uptake=0, excretion=0, prefix="EX_", prefix_name="Exchange for ")`
Add exchange reactions for metabolites.

```python
# Add uptake-only exchanges
mdlutl.add_exchanges_for_metabolites(mets, uptake=1000, excretion=0)

# Add bidirectional exchanges
mdlutl.add_exchanges_for_metabolites(mets, uptake=1000, excretion=1000)
```

#### `add_transport_and_exchange_for_metabolite(met, direction="=", prefix="trans", override=False)`
Add a charge-balanced transport reaction and corresponding exchange.

```python
# Add bidirectional transport for a cytosolic metabolite
transport = mdlutl.add_transport_and_exchange_for_metabolite(
    met_c0, direction="=", prefix="trans"
)
```

#### `add_missing_exchanges(media)`
Add exchange reactions for media compounds that don't have them.

```python
missing = mdlutl.add_missing_exchanges(my_media)
# Returns list of compound IDs that needed exchanges added
```

---

### Media & FBA Configuration

#### `set_media(media)`
Set the model's growth media.

```python
from modelseedpy.core.msmedia import MSMedia

# From MSMedia object
mdlutl.set_media(my_media)

# From dictionary
mdlutl.set_media({"cpd00001": 1000, "cpd00007": 1000})
```

#### `set_objective_from_phenotype(phenotype, missing_transporters=[], create_missing_compounds=False)`
Configure the model objective based on a phenotype type.

```python
# For growth phenotypes, sets biomass objective
# For uptake/excretion phenotypes, sets appropriate exchange objectives
obj_str = mdlutl.set_objective_from_phenotype(phenotype)
```

---

### FBA Testing & Condition Management

#### `apply_test_condition(condition, model=None)`
Apply a test condition (media, objective, direction) to the model.

```python
condition = {
    "media": my_media,
    "objective": "bio1",
    "is_max_threshold": True,
    "threshold": 0.1
}
mdlutl.apply_test_condition(condition)
```

#### `test_single_condition(condition, apply_condition=True, model=None, report_atp_loop_reactions=False, analyze_failures=False, rxn_list=[])`
Test if a model meets a condition's threshold.

```python
passed = mdlutl.test_single_condition(condition)
# Returns True if threshold is NOT exceeded (for is_max_threshold=True)
```

#### `test_condition_list(condition_list, model=None, positive_growth=[], rxn_list=[])`
Test multiple conditions. Returns True only if ALL pass.

```python
all_passed = mdlutl.test_condition_list(conditions)
```

---

### Gapfilling Support

#### `test_solution(solution, targets, medias, thresholds=[0.1], remove_unneeded_reactions=False, do_not_remove_list=[])`
Test if gapfilling solution reactions are needed.

```python
# Solution format: {"new": {rxn_id: direction}, "reversed": {rxn_id: direction}}
# Or: list of [rxn_id, direction, label]
unneeded = mdlutl.test_solution(
    solution,
    targets=["bio1"],
    medias=[glucose_media],
    thresholds=[0.1]
)
```

#### `add_gapfilling(solution)`
Record an integrated gapfilling solution.

```python
mdlutl.add_gapfilling({
    "new": {"rxn00001_c0": ">"},
    "reversed": {"rxn00002_c0": "<"},
    "media": media,
    "target": "bio1",
    "minobjective": 0.1,
    "binary_check": True
})
```

#### `convert_solution_to_list(solution)`
Convert dictionary solution format to list format.

```python
solution_list = mdlutl.convert_solution_to_list(solution)
# Returns [[rxn_id, direction, "new"|"reversed"], ...]
```

---

### Reaction Expansion Testing

These methods are used for binary/linear search to find minimal reaction sets.

#### `reaction_expansion_test(reaction_list, condition_list, binary_search=True, attribute_label="gf_filter", positive_growth=[], resort_by_score=True, active_reaction_sets=[])`
Test which reactions in a list can be removed while still meeting conditions.

```python
filtered = mdlutl.reaction_expansion_test(
    reaction_list=[[rxn, ">"], [rxn2, "<"], ...],
    condition_list=conditions,
    binary_search=True
)
# Returns reactions that were filtered out
```

#### `binary_expansion_test(reaction_list, condition, currmodel, depth=0, positive_growth=[])`
Binary search variant of expansion testing.

#### `linear_expansion_test(reaction_list, condition, currmodel, positive_growth=[])`
Linear (one-by-one) variant of expansion testing.

---

### ATP Correction

#### `get_atputl(atp_media_filename=None, core_template=None, gapfilling_delta=0, max_gapfilling=0, forced_media=[], remake_atputil=False)`
Get or create the MSATPCorrection utility.

```python
atputl = mdlutl.get_atputl(core_template=template)
```

#### `get_atp_tests(core_template=None, atp_media_filename=None, recompute=False, remake_atputil=False)`
Get ATP test conditions.

```python
tests = mdlutl.get_atp_tests(core_template=template)
# Returns list of condition dicts with media, threshold, objective
```

---

### Reliability Scoring

#### `assign_reliability_scores_to_reactions(active_reaction_sets=[])`
Calculate reliability scores for all reactions based on biochemistry data.

```python
scores = mdlutl.assign_reliability_scores_to_reactions()
# scores[rxn_id][">"] = forward score
# scores[rxn_id]["<"] = reverse score
```

Scoring considers:
- Mass/charge balance status
- Delta G values
- Compound completeness
- ATP production direction
- Transported charge

---

### Biomass Analysis

#### `evaluate_biomass_reaction_mass(biomass_rxn_id, normalize=False)`
Calculate the mass balance of a biomass reaction.

```python
result = mdlutl.evaluate_biomass_reaction_mass("bio1")
# Returns {"ATP": atp_coefficient, "Total": total_mass}
```

#### `find_unproducible_biomass_compounds(target_rxn="bio1", ko_list=None)`
Find biomass compounds that cannot be produced.

```python
# Without knockouts
unproducible = mdlutl.find_unproducible_biomass_compounds()

# With knockouts to test sensitivity
ko_results = mdlutl.find_unproducible_biomass_compounds(
    ko_list=[["rxn00001_c0", ">"], ["rxn00002_c0", "<"]]
)
```

---

### Minimal Reaction Set Analysis

#### `analyze_minimal_reaction_set(solution, label, print_output=True)`
Analyze a minimal reaction set for alternatives and coupled reactions.

```python
output = mdlutl.analyze_minimal_reaction_set(fba_solution, "my_analysis")
# Writes CSV to nboutput/rxn_analysis/<label>-min_rxn_set_analysis.csv
```

---

### Model Attributes

#### `get_attributes(key=None, default=None)`
Get model attributes (metadata dictionary).

```python
# Get all attributes
attrs = mdlutl.get_attributes()

# Get specific attribute with default
pathways = mdlutl.get_attributes("pathways", {})
```

#### `save_attributes(value=None, key=None)`
Save attributes back to model (for KBase models).

```python
mdlutl.save_attributes({"my_key": "my_value"}, "my_key")
mdlutl.save_attributes()  # Sync to model.computed_attributes
```

---

### ModelSEED Reaction Addition

#### `add_ms_reaction(rxn_dict, compartment_trans=["c0", "e0"])`
Add reactions from ModelSEED database.

```python
# Add reactions with compartment mapping
reactions = mdlutl.add_ms_reaction({
    "rxn00001": "c0",
    "rxn00002": "c0"
})
```

#### `add_atp_hydrolysis(compartment)`
Add an ATP hydrolysis reaction.

```python
result = mdlutl.add_atp_hydrolysis("c0")
# Returns {"reaction": rxn, "direction": ">", "new": True/False}
```

---

### KBase Integration

#### `convert_cobra_compound_to_kbcompound(cpd, kbmodel, add_to_model=1)`
Convert a COBRA metabolite to KBase format.

#### `convert_cobra_reaction_to_kbreaction(rxn, kbmodel, cpd_hash, direction="=", add_to_model=1, reaction_genes=None)`
Convert a COBRA reaction to KBase format.

#### `create_kb_gapfilling_data(kbmodel, atpmedia_ws="94026")`
Create KBase gapfilling metadata from integrated gapfillings.

---

### Utility Methods

#### `parse_id(object)` [static]
Parse a ModelSEED-style ID into components.

```python
result = MSModelUtil.parse_id(metabolite)
# Returns (base_id, compartment, index) or None
# e.g., ("cpd00001", "c", "0") for cpd00001_c0
```

#### `compute_flux_values_from_variables()`
Get flux values from model solver variables.

```python
fluxes = mdlutl.compute_flux_values_from_variables()
# fluxes[rxn_id] = {"forward": float, "reverse": float}
```

#### `compare_reactions(reaction_list, filename)`
Compare reactions by their metabolite participation.

```python
mdlutl.compare_reactions([rxn1, rxn2, rxn3], "comparison.csv")
```

---

## Common Usage Patterns

### Pattern 1: Load and Analyze a Model

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msmedia import MSMedia

# Load model
mdlutl = MSModelUtil.from_cobrapy("my_model.json")

# Set media
media = MSMedia.from_dict({"EX_cpd00027_e0": 10})  # Glucose
mdlutl.set_media(media)

# Run FBA
solution = mdlutl.model.optimize()
print(f"Growth: {solution.objective_value}")
```

### Pattern 2: Gapfill and Test

```python
from modelseedpy.core.msgapfill import MSGapfill

# Create gapfiller
gapfill = MSGapfill(mdlutl, default_target="bio1")

# Run gapfilling
solution = gapfill.run_gapfilling(media, target="bio1")

# Test if solution reactions are all needed
unneeded = mdlutl.test_solution(
    solution,
    targets=["bio1"],
    medias=[media],
    thresholds=[0.1]
)
```

### Pattern 3: ATP Correction

```python
from modelseedpy.core.mstemplate import MSTemplateBuilder

# Get core template
template = MSTemplateBuilder.build_core_template()

# Get ATP tests
tests = mdlutl.get_atp_tests(core_template=template)

# Run tests
for test in tests:
    passed = mdlutl.test_single_condition(test)
    print(f"{test['media'].id}: {'PASS' if passed else 'FAIL'}")
```

### Pattern 4: Find and Add Reactions

```python
# Find a metabolite
glucose_list = mdlutl.find_met("glucose", "c0")
if glucose_list:
    glucose = glucose_list[0]

    # Add exchange if missing
    if glucose not in mdlutl.exchange_hash():
        mdlutl.add_exchanges_for_metabolites([glucose], uptake=10, excretion=0)

# Add a transport reaction
mdlutl.add_transport_and_exchange_for_metabolite(glucose, direction=">")
```

---

## Module Integration Map

```
┌─────────────────────────────────────────────────────────────────┐
│                         MSModelUtil                              │
│                    (You are here)                                │
└───────────────────────────┬─────────────────────────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
        ▼                   ▼                   ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────────────┐
│    MSFBA      │   │  MSGapfill    │   │   MSPackageManager    │
│  (FBA runner) │   │  (Gapfilling) │   │  (Constraint pkgs)    │
└───────────────┘   └───────────────┘   └───────────────────────┘
        │                   │                   │
        └───────────────────┼───────────────────┘
                            │
        ┌───────────────────┼───────────────────┐
        │                   │                   │
        ▼                   ▼                   ▼
┌───────────────┐   ┌───────────────┐   ┌───────────────────────┐
│  MSMedia      │   │MSATPCorrection│   │   ModelSEEDBiochem    │
│  (Media def)  │   │ (ATP tests)   │   │  (Reaction database)  │
└───────────────┘   └───────────────┘   └───────────────────────┘
```

---

## Debugging Tips

1. **Enable debug logging:**
   ```python
   import logging
   logging.getLogger("modelseedpy.core.msmodelutl").setLevel(logging.DEBUG)
   ```

2. **Print LP file for solver issues:**
   ```python
   mdlutl.printlp(print=True, filename="debug_problem")
   ```

3. **Check metabolite hash is built:**
   ```python
   if mdlutl.metabolite_hash is None:
       mdlutl.build_metabolite_hash()
   ```

4. **Verify MSModelUtil caching:**
   ```python
   print(f"Cached models: {len(MSModelUtil.mdlutls)}")
   ```

---

## Version History

- Part of ModelSEEDpy since initial release
- Continuously enhanced with new analysis methods
- Core API remains stable for backward compatibility
