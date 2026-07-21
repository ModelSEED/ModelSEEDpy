# ModelSEEDpy Architecture

## Overview

ModelSEEDpy is a Python package for metabolic model reconstruction, gapfilling, and analysis. It builds on COBRApy and integrates with the ModelSEED biochemistry database and KBase platform.

## Module Hierarchy

```
┌─────────────────────────────────────────────────────────────────┐
│                        User Code                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│  MSGapfill    │     │  MSCommunity  │     │  MSBuilder    │
│  (Gapfilling) │     │  (Community)  │     │  (Model build)│
└───────────────┘     └───────────────┘     └───────────────┘
        │                     │                     │
        └─────────────────────┼─────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      MSModelUtil                                 │
│              (Central Model Wrapper - core/msmodelutl.py)        │
│                                                                  │
│  • Wraps cobra.Model                                            │
│  • Provides metabolite/reaction search                          │
│  • Manages media, exchanges, tests                              │
│  • Coordinates with other components                            │
└─────────────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐     ┌───────────────┐     ┌───────────────┐
│MSPackageManager│    │ MSATPCorrection│    │ModelSEEDBiochem│
│ (FBA Packages) │    │ (ATP Analysis) │    │ (Biochem DB)  │
└───────────────┘     └───────────────┘     └───────────────┘
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│                     FBA Packages (fbapkg/)                       │
│                                                                  │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────┐             │
│  │GapfillingPkg │ │KBaseMediaPkg │ │FlexBiomassPkg│  ...        │
│  └──────────────┘ └──────────────┘ └──────────────┘             │
│                                                                  │
│  All inherit from BaseFBAPkg                                    │
│  Add variables/constraints to model.solver                      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      COBRApy Model                               │
│                   (cobra.Model object)                           │
│                                                                  │
│  • Reactions, Metabolites, Genes                                │
│  • model.solver (optlang)                                       │
│  • model.optimize()                                             │
└─────────────────────────────────────────────────────────────────┘
```

## Core Modules (modelseedpy/core/)

### MSModelUtil (msmodelutl.py) ~2000 lines
**The central hub for model operations.**

Key responsibilities:
- Wrap and extend cobra.Model
- Metabolite/reaction search and lookup
- Exchange and transport management
- Media configuration
- FBA testing and condition management
- Gapfilling support methods
- Integration with all other components

```python
from modelseedpy.core.msmodelutl import MSModelUtil
mdlutl = MSModelUtil.get(model)  # Singleton access
```

### MSGapfill (msgapfill.py) ~1200 lines
**Automated model gapfilling.**

Key features:
- Multi-media gapfilling
- ATP-aware gapfilling
- Binary/linear reaction filtering
- Solution testing and validation

```python
from modelseedpy.core.msgapfill import MSGapfill
gapfill = MSGapfill(mdlutl, default_target="bio1")
solution = gapfill.run_gapfilling(media, target="bio1")
```

### MSATPCorrection (msatpcorrection.py)
**ATP production analysis and correction.**

Prevents models from producing ATP without valid biochemistry.

```python
atputl = mdlutl.get_atputl(core_template=template)
atp_tests = mdlutl.get_atp_tests()
```

### MSFBA (msfba.py)
**Higher-level FBA execution with reporting.**

```python
from modelseedpy.core.msfba import MSFBA
fba = MSFBA(mdlutl)
result = fba.run_fba()
```

### MSMedia (msmedia.py)
**Growth media definitions.**

```python
from modelseedpy.core.msmedia import MSMedia
media = MSMedia.from_dict({"EX_cpd00027_e0": 10})
media = MSMedia.from_file("media.tsv")
```

### MSBuilder (msbuilder.py)
**Model construction from genome annotations.**

```python
from modelseedpy.core.msbuilder import MSBuilder
builder = MSBuilder(genome, template)
model = builder.build()
```

### MSTemplate (mstemplate.py)
**Model templates for reconstruction.**

Templates define which reactions can be added during reconstruction and their properties.

### MSGrowthPhenotypes (msgrowthphenotypes.py)
**Phenotype testing and comparison.**

Test model predictions against experimental growth data.

## FBA Packages (modelseedpy/fbapkg/)

### MSPackageManager
**Central registry for FBA packages.**

```python
from modelseedpy.fbapkg import MSPackageManager
pkgmgr = MSPackageManager.get_pkg_mgr(model)  # Singleton

# List available packages
pkgmgr.list_available_packages()

# Get or create a package
pkg = pkgmgr.getpkg("GapfillingPkg")
```

### BaseFBAPkg
**Base class for all FBA packages.**

All packages inherit from this and implement:
- `build_package(params)` - Add constraints/variables
- `clear()` - Remove constraints/variables

### Key Packages

| Package | Purpose |
|---------|---------|
| `GapfillingPkg` | Gapfilling MILP formulation |
| `KBaseMediaPkg` | Media exchange constraints |
| `FlexibleBiomassPkg` | Flexible biomass composition |
| `SimpleThermoPkg` | Simple thermodynamic constraints |
| `FullThermoPkg` | Full thermodynamic constraints |
| `ReactionUsePkg` | Binary reaction usage variables |
| `RevBinPkg` | Reversibility binary variables |
| `ObjectivePkg` | Objective function management |
| `TotalFluxPkg` | Total flux minimization |
| `BilevelPkg` | Bilevel optimization |

## Community Module (modelseedpy/community/)

### MSCommunity (mscommunity.py)
**Multi-species community modeling.**

```python
from modelseedpy.community.mscommunity import MSCommunity
community = MSCommunity(member_models=[m1, m2, m3])
```

### MSSteadyCom (mssteadycom.py)
**SteadyCom algorithm for community FBA.**

Computes steady-state community compositions.

## Biochemistry Module (modelseedpy/biochem/)

### ModelSEEDBiochem (modelseed_biochem.py)
**Access to ModelSEED reaction/compound database.**

```python
from modelseedpy.biochem import ModelSEEDBiochem
biochem = ModelSEEDBiochem.get()
reaction = biochem.get_reaction("rxn00001")
compound = biochem.get_compound("cpd00001")
```

## Key Design Patterns

### 1. Singleton/Cache Pattern
Used by MSModelUtil, MSPackageManager, ModelSEEDBiochem:

```python
# Same instance returned for same model
mdlutl1 = MSModelUtil.get(model)
mdlutl2 = MSModelUtil.get(model)
assert mdlutl1 is mdlutl2
```

### 2. Model/Utility Acceptance
All high-level classes accept either raw model or utility:

```python
def __init__(self, model_or_mdlutl):
    self.mdlutl = MSModelUtil.get(model_or_mdlutl)
    self.model = self.mdlutl.model
```

### 3. Package Registration
FBA packages self-register with MSPackageManager:

```python
class MyPkg(BaseFBAPkg):
    def __init__(self, model):
        super().__init__(model, "MyPkg", ...)
        # BaseFBAPkg.__init__ calls pkgmgr.addpkgobj(self)
```

### 4. Lazy Loading
Heavy components loaded on demand:

```python
# MSATPCorrection created only when needed
atputl = mdlutl.get_atputl()  # Creates if missing
```

## Data Flow Example: Gapfilling

```
User Request: "Gapfill model on glucose media"
                    │
                    ▼
            ┌───────────────┐
            │   MSGapfill   │
            │               │
            │ 1. Get media  │
            │ 2. Setup FBA  │
            │ 3. Run MILP   │
            │ 4. Filter     │
            └───────────────┘
                    │
        ┌───────────┼───────────┐
        │           │           │
        ▼           ▼           ▼
┌─────────────┐ ┌─────────┐ ┌─────────────┐
│MSModelUtil  │ │GapfillPkg│ │KBaseMediaPkg│
│             │ │          │ │             │
│set_media()  │ │build_pkg │ │build_pkg    │
│test_soln()  │ │(MILP)    │ │(bounds)     │
└─────────────┘ └─────────┘ └─────────────┘
        │           │           │
        └───────────┼───────────┘
                    │
                    ▼
            ┌───────────────┐
            │  cobra.Model  │
            │               │
            │ .solver       │
            │ .optimize()   │
            └───────────────┘
```
