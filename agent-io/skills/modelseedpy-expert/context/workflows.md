# ModelSEEDpy Common Workflows

## Workflow 1: Load and Analyze an Existing Model

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msmedia import MSMedia

# Load model from file
mdlutl = MSModelUtil.from_cobrapy("model.json")
# Or wrap an existing COBRApy model
# mdlutl = MSModelUtil.get(cobra_model)

# Inspect model
print(f"Reactions: {len(mdlutl.model.reactions)}")
print(f"Metabolites: {len(mdlutl.model.metabolites)}")

# Find specific metabolites
glucose_list = mdlutl.find_met("glucose", "c0")
if glucose_list:
    glucose = glucose_list[0]
    print(f"Found glucose: {glucose.id}")

# Set up media
media = MSMedia.from_dict({
    "EX_cpd00027_e0": 10,   # Glucose
    "EX_cpd00001_e0": 1000, # Water
    "EX_cpd00009_e0": 1000, # Phosphate
    # ... other nutrients
})

# Ensure exchanges exist
mdlutl.add_missing_exchanges(media)
mdlutl.set_media(media)

# Run FBA
solution = mdlutl.model.optimize()
print(f"Growth rate: {solution.objective_value}")
```

## Workflow 2: Gapfill a Non-Growing Model

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.msmedia import MSMedia

# Load model
mdlutl = MSModelUtil.from_cobrapy("draft_model.json")

# Define media
media = MSMedia.from_dict({"EX_cpd00027_e0": 10})
mdlutl.add_missing_exchanges(media)

# Check if model grows (probably not if draft)
mdlutl.set_media(media)
sol = mdlutl.model.optimize()
print(f"Pre-gapfill growth: {sol.objective_value}")

# Create gapfiller
gapfill = MSGapfill(
    mdlutl,
    default_target="bio1",
    minimum_obj=0.1  # Minimum required growth
)

# Run gapfilling
solution = gapfill.run_gapfilling(
    media=media,
    target="bio1"
)

print(f"Gapfilling solution: {solution}")

# Test which reactions are truly needed
unneeded = mdlutl.test_solution(
    solution,
    targets=["bio1"],
    medias=[media],
    thresholds=[0.1],
    remove_unneeded_reactions=True
)

# Record the gapfilling
mdlutl.add_gapfilling(solution)

# Verify growth
sol = mdlutl.model.optimize()
print(f"Post-gapfill growth: {sol.objective_value}")

# Save model
mdlutl.save_model("gapfilled_model.json")
```

## Workflow 3: ATP-Aware Gapfilling

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgapfill import MSGapfill
from modelseedpy.core.mstemplate import MSTemplateBuilder

# Load model
mdlutl = MSModelUtil.from_cobrapy("model.json")

# Get core template for ATP tests
template = MSTemplateBuilder.build_core_template()

# Get ATP test conditions (prevents ATP loops)
atp_tests = mdlutl.get_atp_tests(core_template=template)

# Create gapfiller with ATP constraints
gapfill = MSGapfill(mdlutl, default_target="bio1")

# Run ATP-constrained gapfilling
solution = gapfill.run_gapfilling(
    media=media,
    target="bio1",
    atp_tests=atp_tests  # Prevents solutions that produce free ATP
)
```

## Workflow 4: Test Growth on Multiple Media

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msmedia import MSMedia

mdlutl = MSModelUtil.from_cobrapy("model.json")

# Define test conditions
conditions = [
    {
        "media": MSMedia.from_dict({"EX_cpd00027_e0": 10}),  # Glucose
        "objective": "bio1",
        "is_max_threshold": False,  # Must grow ABOVE threshold
        "threshold": 0.1
    },
    {
        "media": MSMedia.from_dict({"EX_cpd00029_e0": 10}),  # Acetate
        "objective": "bio1",
        "is_max_threshold": False,
        "threshold": 0.05
    },
    {
        "media": MSMedia.from_dict({"EX_cpd00036_e0": 10}),  # Succinate
        "objective": "bio1",
        "is_max_threshold": False,
        "threshold": 0.05
    }
]

# Add missing exchanges for all media
for cond in conditions:
    mdlutl.add_missing_exchanges(cond["media"])

# Test all conditions
results = {}
for i, cond in enumerate(conditions):
    passed = mdlutl.test_single_condition(cond)
    media_name = f"condition_{i}"
    results[media_name] = passed
    print(f"{media_name}: {'PASS' if passed else 'FAIL'}")

# Or use batch testing
all_passed = mdlutl.test_condition_list(conditions)
print(f"All conditions passed: {all_passed}")
```

## Workflow 5: Build Model from Genome

```python
from modelseedpy.core.msbuilder import MSBuilder
from modelseedpy.core.msgenome import MSGenome
from modelseedpy.core.mstemplate import MSTemplateBuilder

# Load genome
genome = MSGenome.from_fasta("genome.fasta")
# Or from annotation
# genome = MSGenome.from_rast(annotation_data)

# Get template
template = MSTemplateBuilder.build_template("GramNegative")

# Build model
builder = MSBuilder(genome, template)
model = builder.build()

# Wrap in MSModelUtil for further operations
mdlutl = MSModelUtil.get(model)
print(f"Built model with {len(model.reactions)} reactions")
```

## Workflow 6: Community Modeling

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.community.mscommunity import MSCommunity

# Load individual models
model1 = MSModelUtil.from_cobrapy("species1.json").model
model2 = MSModelUtil.from_cobrapy("species2.json").model
model3 = MSModelUtil.from_cobrapy("species3.json").model

# Create community
community = MSCommunity(
    member_models=[model1, model2, model3],
    ids=["sp1", "sp2", "sp3"]
)

# Run community FBA
result = community.run_fba()

# Get individual contributions
for member in community.members:
    print(f"{member.id}: {member.growth_rate}")
```

## Workflow 7: Add Custom FBA Constraints

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg import MSPackageManager

mdlutl = MSModelUtil.from_cobrapy("model.json")
pkgmgr = mdlutl.pkgmgr

# Get reaction use package (binary variables for reaction on/off)
rxn_use_pkg = pkgmgr.getpkg("ReactionUsePkg")
rxn_use_pkg.build_package({
    "reaction_list": mdlutl.model.reactions
})

# Get total flux package (minimize total flux)
total_flux_pkg = pkgmgr.getpkg("TotalFluxPkg")
total_flux_pkg.build_package()

# Get thermodynamic package
thermo_pkg = pkgmgr.getpkg("SimpleThermoPkg")
thermo_pkg.build_package()

# Run FBA with all constraints active
solution = mdlutl.model.optimize()
```

## Workflow 8: Flexible Biomass Analysis

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg import MSPackageManager

mdlutl = MSModelUtil.from_cobrapy("model.json")

# Get flexible biomass package
flex_bio_pkg = mdlutl.pkgmgr.getpkg("FlexibleBiomassPkg")

# Build with flexibility parameters
flex_bio_pkg.build_package({
    "bio_rxn_id": "bio1",
    "flex_coefficient": 0.1,  # Allow 10% flexibility
    "use_rna_class": True,
    "use_protein_class": True
})

# Now biomass composition can vary within bounds
solution = mdlutl.model.optimize()
```

## Workflow 9: Compare Multiple Solutions

```python
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msmedia import MSMedia

mdlutl = MSModelUtil.from_cobrapy("model.json")

# Run FBA on different media and collect solutions
solutions = {}

media_glucose = MSMedia.from_dict({"EX_cpd00027_e0": 10})
mdlutl.add_missing_exchanges(media_glucose)
mdlutl.set_media(media_glucose)
solutions["glucose"] = mdlutl.model.optimize()

media_acetate = MSMedia.from_dict({"EX_cpd00029_e0": 10})
mdlutl.add_missing_exchanges(media_acetate)
mdlutl.set_media(media_acetate)
solutions["acetate"] = mdlutl.model.optimize()

# Export comparison to CSV
mdlutl.print_solutions(solutions, "flux_comparison.csv")
```

## Workflow 10: Debugging - Find Unproducible Biomass Components

```python
from modelseedpy.core.msmodelutl import MSModelUtil

mdlutl = MSModelUtil.from_cobrapy("model.json")

# Set up media
mdlutl.set_media(media)

# Find biomass components that can't be produced
unproducible = mdlutl.find_unproducible_biomass_compounds(
    target_rxn="bio1"
)

for met in unproducible:
    print(f"Cannot produce: {met.id} - {met.name}")

# Check sensitivity to specific reaction knockouts
ko_results = mdlutl.find_unproducible_biomass_compounds(
    target_rxn="bio1",
    ko_list=[["rxn00001_c0", ">"], ["rxn00002_c0", "<"]]
)
```
