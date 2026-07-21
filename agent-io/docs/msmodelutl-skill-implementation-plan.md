# MSModelUtil Expert Skill Implementation Plan

## Overview

This document describes how to implement a Claude Code skill/command that provides expert-level assistance for the `MSModelUtil` class from ModelSEEDpy. The skill should be implemented in your `claude_commands` repository and be invokable from any project context.

## Goal

Create a skill that allows you to:
1. Get expert help using the MSModelUtil API
2. Get suggestions for improving MSModelUtil code
3. Understand how MSModelUtil integrates with other ModelSEEDpy components
4. Debug issues in code that uses MSModelUtil

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    Your claude_commands Repo                     │
│                                                                  │
│  ┌────────────────────────────────────────────────────────────┐ │
│  │  /msmodelutl-expert                                        │ │
│  │  (Skill Definition)                                         │ │
│  │                                                             │ │
│  │  ├── prompt.md          # Skill prompt template             │ │
│  │  └── context/           # Reference documentation           │ │
│  │      ├── api-reference.md                                   │ │
│  │      ├── integration-map.md                                 │ │
│  │      └── patterns.md                                        │ │
│  └────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ invokes
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                     Any Project Context                          │
│                                                                  │
│  User: /msmodelutl-expert How do I add an exchange reaction?    │
│                                                                  │
│  Claude: [Reads skill context, provides expert answer]          │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Skill Implementation

### File Structure

In your `claude_commands` repository, create:

```
commands/
└── msmodelutl-expert/
    ├── skill.md              # Main skill definition
    └── context/
        ├── api-summary.md    # Condensed API reference
        ├── patterns.md       # Common usage patterns
        └── integration.md    # How it integrates with other modules
```

### Skill Definition (skill.md)

```markdown
---
name: msmodelutl-expert
description: Expert assistance for ModelSEEDpy's MSModelUtil class
version: 1.0.0
---

# MSModelUtil Expert

You are an expert on the MSModelUtil class from ModelSEEDpy. You have deep knowledge of:

1. **The MSModelUtil API** - All 55+ methods, their parameters, return values, and usage
2. **Integration patterns** - How MSModelUtil connects with MSGapfill, MSFBA, MSPackageManager, etc.
3. **Best practices** - Efficient ways to use the API, common pitfalls to avoid
4. **Debugging** - How to diagnose issues in code using MSModelUtil

## Your Knowledge Base

<context>
{{include:context/api-summary.md}}
</context>

<context>
{{include:context/patterns.md}}
</context>

<context>
{{include:context/integration.md}}
</context>

## Guidelines

When helping users:

1. **Be specific** - Reference exact method names, parameters, and return types
2. **Show examples** - Provide working code snippets
3. **Explain integration** - Show how methods connect to other ModelSEEDpy components
4. **Warn about pitfalls** - Mention common mistakes and how to avoid them

## Response Format

For API questions:
```
### Method: `method_name(params)`

**Purpose:** Brief description

**Parameters:**
- `param1` (type): Description
- `param2` (type, optional): Description

**Returns:** Description of return value

**Example:**
\`\`\`python
# Working example
\`\`\`

**Related methods:** List of related methods
```

For "how do I" questions:
```
### Approach

Brief explanation of the approach.

**Step 1:** Description
\`\`\`python
code
\`\`\`

**Step 2:** Description
\`\`\`python
code
\`\`\`

**Notes:** Any important considerations
```
```

### Context Files

#### context/api-summary.md

Create a condensed version of the API (key methods only, not the full 500-line doc):

```markdown
# MSModelUtil API Summary

## Core Concepts

- **Singleton pattern**: Use `MSModelUtil.get(model)` to get/create instances
- **Wraps cobra.Model**: Access via `mdlutl.model`
- **Integrates with MSPackageManager**: Access via `mdlutl.pkgmgr`

## Essential Methods

### Factory/Initialization
- `MSModelUtil.get(model)` - Get or create instance (PREFERRED)
- `MSModelUtil.from_cobrapy(filename)` - Load from file
- `MSModelUtil(model)` - Direct construction

### Metabolite Search
- `find_met(name, compartment=None)` - Find metabolites by name/ID
- `msid_hash()` - Get ModelSEED ID to metabolite mapping
- `metabolite_msid(met)` [static] - Extract ModelSEED ID from metabolite

### Reaction Operations
- `rxn_hash()` - Get stoichiometry to reaction mapping
- `find_reaction(stoichiometry)` - Find reaction by stoichiometry
- `exchange_list()` / `exchange_hash()` - Get exchange reactions
- `is_core(rxn)` - Check if reaction is core metabolism

### Exchange/Transport
- `add_exchanges_for_metabolites(cpds, uptake, excretion)` - Add exchanges
- `add_transport_and_exchange_for_metabolite(met, direction)` - Add transport
- `add_missing_exchanges(media)` - Fill media gaps

### Media/FBA
- `set_media(media)` - Configure growth media
- `apply_test_condition(condition)` - Apply test constraints
- `test_single_condition(condition)` - Run single test
- `test_condition_list(conditions)` - Run multiple tests

### Gapfilling Support
- `test_solution(solution, targets, medias, thresholds)` - Validate solutions
- `add_gapfilling(solution)` - Record integrated gapfilling
- `reaction_expansion_test(rxn_list, conditions)` - Find minimal sets

### ATP Correction
- `get_atputl()` - Get ATP correction utility
- `get_atp_tests()` - Get ATP test conditions

### Model Editing
- `add_ms_reaction(rxn_dict)` - Add ModelSEED reactions
- `add_atp_hydrolysis(compartment)` - Add ATP hydrolysis
- `get_attributes()` / `save_attributes()` - Model metadata

### Analysis
- `assign_reliability_scores_to_reactions()` - Score reactions
- `find_unproducible_biomass_compounds()` - Biomass sensitivity
- `analyze_minimal_reaction_set(solution, label)` - Alternative analysis
```

#### context/patterns.md

```markdown
# Common MSModelUtil Patterns

## Pattern 1: Safe Instance Access
```python
# Always use get() for consistent instance access
mdlutl = MSModelUtil.get(model)  # Works with model or mdlutl

# Functions should accept either
def my_function(model_or_mdlutl):
    mdlutl = MSModelUtil.get(model_or_mdlutl)
    model = mdlutl.model
```

## Pattern 2: Find and Operate on Metabolites
```python
# Always handle empty results
mets = mdlutl.find_met("glucose", "c0")
if mets:
    glucose = mets[0]
    # Do something with glucose
else:
    # Handle not found
```

## Pattern 3: Add Exchanges for Media
```python
# Before setting media, ensure exchanges exist
missing = mdlutl.add_missing_exchanges(media)
if missing:
    print(f"Added exchanges for: {missing}")
mdlutl.set_media(media)
```

## Pattern 4: Test Growth Conditions
```python
condition = {
    "media": media,
    "objective": "bio1",
    "is_max_threshold": True,  # True = must be BELOW threshold
    "threshold": 0.1
}
mdlutl.apply_test_condition(condition)
passed = mdlutl.test_single_condition(condition, apply_condition=False)
```

## Pattern 5: Gapfill and Validate
```python
# After gapfilling
solution = gapfiller.run_gapfilling(media, target="bio1")

# Test which reactions are actually needed
unneeded = mdlutl.test_solution(
    solution,
    targets=["bio1"],
    medias=[media],
    thresholds=[0.1],
    remove_unneeded_reactions=True  # Actually remove them
)
```

## Common Mistakes

1. **Not using get()**: Creating multiple MSModelUtil instances for same model
2. **Ignoring empty find_met results**: Always check if list is empty
3. **Forgetting build_metabolite_hash()**: Called automatically by find_met, but cached
4. **Wrong threshold interpretation**: is_max_threshold=True means FAIL if >= threshold
```

#### context/integration.md

```markdown
# MSModelUtil Integration Map

## Key Relationships

### MSModelUtil ↔ MSGapfill
- MSGapfill takes MSModelUtil in constructor
- Sets `mdlutl.gfutl = self` for bidirectional access
- Uses `mdlutl.test_solution()` for solution validation
- Uses `mdlutl.reaction_expansion_test()` for minimal solutions

### MSModelUtil ↔ MSPackageManager
- Created automatically: `self.pkgmgr = MSPackageManager.get_pkg_mgr(model)`
- Used for media: `self.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)`
- All FBA packages access model through MSPackageManager

### MSModelUtil ↔ MSATPCorrection
- Lazy-loaded via `get_atputl()`
- Sets `self.atputl` for caching
- Uses ATP tests for gapfilling constraints

### MSModelUtil ↔ ModelSEEDBiochem
- Used in `add_ms_reaction()` for reaction data
- Used in `assign_reliability_scores_to_reactions()` for scoring

### MSModelUtil ↔ MSFBA
- MSFBA wraps model_or_mdlutl input
- Uses MSModelUtil for consistent access

## Dependency Chain

```
User Code
    │
    ▼
MSGapfill / MSFBA / MSCommunity
    │
    ▼
MSModelUtil  ◄──────────────────┐
    │                           │
    ├── MSPackageManager ───────┤
    │       │                   │
    │       ▼                   │
    │   FBA Packages            │
    │                           │
    ├── MSATPCorrection ────────┤
    │                           │
    └── ModelSEEDBiochem        │
            │                   │
            └───────────────────┘
```
```

## Installation

### Option 1: User-level Skill

Place in your Claude Code user commands directory:
```
~/.claude/commands/msmodelutl-expert/
    ├── skill.md
    └── context/
        ├── api-summary.md
        ├── patterns.md
        └── integration.md
```

### Option 2: Project-level Skill (in claude_commands repo)

If your `claude_commands` repo is set up as a source for Claude commands:
```
claude_commands/
└── commands/
    └── msmodelutl-expert/
        └── [files as above]
```

## Usage

Once installed, invoke from any project:

```
/msmodelutl-expert How do I add an exchange reaction for glucose?

/msmodelutl-expert What's the difference between find_met and msid_hash?

/msmodelutl-expert Debug this code that's failing to set media correctly

/msmodelutl-expert How does MSModelUtil integrate with gapfilling?
```

## Keeping Documentation Current

When you update `msmodelutl.py`:

1. Update the full developer guide in ModelSEEDpy:
   - `agent-io/docs/msmodelutl-developer-guide.md`

2. Update the condensed skill context in claude_commands:
   - `commands/msmodelutl-expert/context/api-summary.md`

Consider creating a script or hook that syncs key changes.

## Alternative: Dynamic Documentation Loading

Instead of static context files, you could have the skill:

1. Read the current `msmodelutl.py` file dynamically
2. Read the developer guide from ModelSEEDpy repo
3. Parse and understand the current state

This requires the skill to know the path to ModelSEEDpy, which could be:
- Configured in the skill
- Passed as a parameter
- Discovered via environment variable

Example skill prompt with dynamic loading:
```markdown
Before answering, read the current MSModelUtil implementation and documentation:

1. Read: /Users/chenry/Dropbox/Projects/ModelSEEDpy/modelseedpy/core/msmodelutl.py
2. Read: /Users/chenry/Dropbox/Projects/ModelSEEDpy/agent-io/docs/msmodelutl-developer-guide.md

Use this current information to answer the user's question.
```

## Summary

This plan provides:
1. **Static context approach** - Pre-written documentation embedded in skill
2. **Dynamic loading approach** - Read current files at invocation time
3. **File structure** - How to organize the skill in your commands repo
4. **Context content** - What documentation to include

The static approach is faster but requires maintenance. The dynamic approach is always current but adds latency. A hybrid (static for patterns, dynamic for API details) may be optimal.
