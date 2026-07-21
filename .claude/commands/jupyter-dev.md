---
name: Jupyter Development
description: Create and manage Jupyter notebooks for data analysis and development
scope: universal
---

# Command: jupyter-dev

## Purpose

Develop Jupyter notebooks following a standardized workflow that emphasizes:
- Organized directory structure with data, models, and output segregation
- Independent, self-contained cells that can run in any order
- Centralized utilities and imports via util.py
- Intermediate data caching for debugging and efficiency
- Clear markdown documentation preceding each code cell

## Command Type

`jupyter-dev`

## Input

You will receive a request file containing:
- Notebook development task description
- Project name (for util.py configuration)
- Specific analysis or computation requirements
- Input data files (optional)
- User preferences (optional)

## Project Structure

All notebooks must follow this directory structure:

```
notebooks/
├── util.py                  # Centralized utilities and imports
├── <notebook-name>.ipynb    # Notebook files
├── data/                    # Input data (experimental, omics, expression data)
├── datacache/               # JSON output from util.save() function
├── genomes/                 # Genome files
├── models/                  # COBRA/COBRApy models
└── nboutput/                # Non-JSON output (TSV, Excel, tables, etc.)
```

### Directory Purposes

- **notebooks/**: Root directory containing all notebooks and util.py
- **data/**: All input data files (experimental data, omics data, expression data)
- **datacache/**: Intermediate JSON data saved via util.save() for cell independence
- **genomes/**: Genome files only
- **models/**: COBRA/COBRApy model files only
- **nboutput/**: Non-JSON output files (TSV, Excel, tables, plots, etc.)

## util.py Structure

The util.py file must follow this template:

```python
import sys
import os
import json
from os import path

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = os.path.dirname(os.path.dirname(script_dir))
folder_name = os.path.basename(script_dir)

print(base_dir+"/KBUtilLib/src")
sys.path = [base_dir+"/KBUtilLib/src",base_dir+"/cobrakbase",base_dir+"/ModelSEEDpy/"] + sys.path

# Import utilities with error handling
from kbutillib import NotebookUtils

import hashlib
import pandas as pd
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression

class NotebookUtil(NotebookUtils):
    def __init__(self,**kwargs):
        super().__init__(
            notebook_folder=script_dir,
            name="<PROJECT_NAME>",
            user="chenry",
            retries=5,
            proxy_port=None,
            **kwargs
        )

    # PLACE ALL UTILITY FUNCTIONS NEEDED FOR NOTEBOOKS HERE

# Initialize the NotebookUtil instance
util = NotebookUtil()
```

### Key Points for util.py

1. **Replace `<PROJECT_NAME>`** with the actual project name from user input
2. **Add all imports** needed by notebooks to this file
3. **Add all utility functions** as methods of the NotebookUtil class
4. **Keep it centralized**: All shared code goes here, not in notebook cells

## Notebook Cell Design Pattern

Every notebook must follow this strict cell pattern:

### 1. Markdown Cell (Always First)
```markdown
## [Step Name/Purpose]

[Explanation of what this code cell does and why]
- Key objective
- Input data used
- Output data produced
- Any important notes
```

### 2. Code Cell (Always Second)
```python
%run util.py

# Load required data from previous steps
data1 = util.load('data1_name')
data2 = util.load('data2_name')

# Perform analysis/computation
result = some_analysis(data1, data2)

# Save intermediate results for cell independence
util.save('result_name', result)
```

### Critical Cell Design Rules

1. **Every code cell starts with**: `%run util.py`
   - This instantiates the util class
   - This loads all imports
   - This ensures cell independence

2. **Load data at cell start**: Use `util.load('data_name')` for any data from previous cells
   - Only load what this cell needs
   - Data comes from datacache/ directory

3. **Save data at cell end**: Use `util.save('data_name', data)` for outputs
   - Save all intermediate results that other cells might need
   - Only JSON-serializable data structures
   - Saved to datacache/ directory

4. **Cell independence**: Each cell should run independently
   - Don't rely on variables from previous cells without loading them
   - Don't assume cells run in order
   - Enable debugging by re-running individual cells

5. **Markdown precedes code**: Every code cell has a markdown cell explaining it
   - What the cell does
   - Why it's needed
   - What data it uses and produces

## Process

### Phase 1: Setup Project Structure

1. **Check for notebooks/ Directory**
   - If `notebooks/` doesn't exist, create it
   - If it exists, verify subdirectories

2. **Create Required Subdirectories**
   - Create `notebooks/data/` if missing
   - Create `notebooks/datacache/` if missing
   - Create `notebooks/genomes/` if missing
   - Create `notebooks/models/` if missing
   - Create `notebooks/nboutput/` if missing

3. **Create or Validate util.py**
   - If `notebooks/util.py` doesn't exist, create it from template
   - Replace `<PROJECT_NAME>` with actual project name
   - If util.py exists, verify it has the NotebookUtil class
   - Document whether created or validated

### Phase 2: Understand Requirements

4. **Analyze Task Description**
   - Identify the scientific/analytical goal
   - Determine required input data
   - Identify computation steps needed
   - Plan logical cell breakdown
   - Determine what utility functions might be needed

5. **Plan Notebook Structure**
   - Break task into logical steps (cells)
   - Identify data flow between cells
   - Determine what gets saved/loaded at each step
   - Plan utility functions for util.py
   - Document the planned structure

### Phase 3: Develop Utility Functions

6. **Add Utility Functions to util.py**
   - Add any custom functions needed by notebooks
   - Add imports required for these functions
   - Add functions as methods to NotebookUtil class
   - Document each function with docstrings
   - Keep functions general and reusable

### Phase 4: Create/Modify Notebook

7. **Create Notebook Cells**
   - For each logical step:
     - Create markdown cell explaining the step
     - Create code cell with proper pattern:
       - Start with `%run util.py`
       - Load required data with util.load()
       - Perform computation
       - Save results with util.save()
   - Follow cell independence principles
   - Add clear variable names and comments

8. **Organize Data Files**
   - Move/reference input data to `notebooks/data/`
   - Reference genome files from `notebooks/genomes/`
   - Reference model files from `notebooks/models/`
   - Save non-JSON output to `notebooks/nboutput/`
   - Let util.save() handle datacache/ automatically

### Phase 5: Validate and Document

9. **Verify Notebook Standards**
   - Every code cell starts with `%run util.py`
   - Every code cell has preceding markdown explanation
   - Data dependencies use util.load()
   - Results saved with util.save()
   - Cells can run independently
   - All files in correct directories

10. **Create Summary Documentation**
    - Document notebook purpose and workflow
    - List required input data and locations
    - Describe each major step
    - Note any manual setup required
    - Include example usage

### Phase 6: Save Structured Output

11. **Save JSON Tracking File**
    - Document all files created/modified
    - List all utility functions added
    - Describe notebook cell structure
    - Note any issues or edge cases
    - Include completion status

## JSON Output Schema

The command execution tracking file must follow this structure:

```json
{
  "command_type": "jupyter-dev",
  "status": "complete | incomplete | user_query | error",
  "session_id": "string",
  "parent_session_id": "string | null",
  "session_summary": "Brief summary of notebook development work",

  "project": {
    "name": "string - project name used in util.py",
    "notebook_name": "string - name of notebook file",
    "purpose": "string - what this notebook does"
  },

  "structure": {
    "directories_created": ["data", "datacache", "genomes", "models", "nboutput"],
    "util_py_status": "created | existed | modified",
    "notebook_path": "notebooks/<notebook-name>.ipynb"
  },

  "notebook_cells": [
    {
      "cell_number": 1,
      "type": "markdown | code",
      "purpose": "Description of what this cell does",
      "data_loaded": ["data1", "data2"],
      "data_saved": ["result1"]
    }
  ],

  "utility_functions": [
    {
      "name": "function_name",
      "purpose": "What this utility function does",
      "added_to_util_py": true
    }
  ],

  "files": {
    "created": [
      {
        "path": "notebooks/util.py",
        "purpose": "Centralized utilities and imports",
        "type": "code"
      }
    ],
    "modified": [
      {
        "path": "notebooks/analysis.ipynb",
        "changes": "Added 5 cells for data loading and analysis"
      }
    ],
    "data_files": [
      {
        "path": "notebooks/data/experimental_data.csv",
        "purpose": "Input experimental data",
        "type": "input"
      }
    ]
  },

  "artifacts": {
    "notebook_filename": "notebooks/<notebook-name>.ipynb",
    "util_py_path": "notebooks/util.py",
    "cell_count": 10,
    "utility_function_count": 3
  },

  "validation": {
    "all_cells_have_markdown": true,
    "all_cells_start_with_run_util": true,
    "data_loading_uses_util_load": true,
    "data_saving_uses_util_save": true,
    "cells_independent": true,
    "files_in_correct_directories": true
  },

  "comments": [
    "Created notebook structure with 5 analysis steps",
    "Added 3 utility functions for data processing",
    "All cells follow independence pattern with util.load/save",
    "Input data placed in notebooks/data/",
    "Output tables saved to notebooks/nboutput/"
  ],

  "queries_for_user": [],

  "errors": []
}
```

## Command JSON Output Requirements

Your command execution JSON output must include:

**Required Fields:**
- `command_type`: "jupyter-dev"
- `status`: "complete", "user_query", or "error"
- `session_id`: Session ID for this execution
- `session_summary`: Brief summary of notebook development
- `project`: Project name and notebook details
- `structure`: Directory and util.py status
- `files`: All files created, modified, or referenced
- `artifacts`: Paths to notebook and util.py
- `validation`: Checklist confirming standards followed
- `comments`: Notes about development process

**For user_query status:**
- `queries_for_user`: Questions needing clarification
- `context`: Save partial work and notebook state

**Example Comments:**
- "Created notebooks directory structure with all required subdirectories"
- "Generated util.py with project name 'MetabolicAnalysis'"
- "Created notebook with 8 cells following independence pattern"
- "Added 4 utility functions for COBRA model manipulation"
- "All intermediate results saved to datacache/ for cell independence"
- "Placed genome files in genomes/, model files in models/"

## Design Principles

### Cell Independence Philosophy

The notebook design prioritizes **cell independence** for several critical reasons:

1. **Debugging Efficiency**: Re-run individual cells without executing entire notebook
2. **Time Savings**: Skip expensive computations by loading cached results
3. **Error Recovery**: Recover from failures without losing all progress
4. **Experimentation**: Test variations by modifying single cells
5. **Collaboration**: Others can understand and modify individual steps

### Implementation Strategy

- **util.load()** and **util.save()** create checkpoints
- **datacache/** stores intermediate results as JSON
- **%run util.py** ensures consistent environment
- **Markdown cells** provide context for each step

### When to Save Data

Save data when:
- Results took significant time to compute
- Data will be used by multiple subsequent cells
- Intermediate results are worth preserving
- Enabling cell re-runs would save time

Don't save data when:
- Quick computations (< 1 second)
- Data only used in next cell
- Data is not JSON-serializable (save to nboutput/ instead)

## Utility Function Guidelines

Add functions to util.py when:
- Code is used by multiple cells
- Complex operations that need documentation
- Interactions with external systems (APIs, databases)
- Data transformations used repeatedly
- Model-specific operations

Keep in notebooks when:
- Code is cell-specific analysis
- One-time exploratory code
- Visualization/plotting specific to that cell
- Simple operations that don't need abstraction

## Quality Checklist

Before marking complete, verify:
- ✅ notebooks/ directory exists with all 5 subdirectories
- ✅ util.py exists and has correct project name
- ✅ util.py contains NotebookUtil class with needed functions
- ✅ Every code cell starts with `%run util.py`
- ✅ Every code cell has preceding markdown explanation
- ✅ Data dependencies use util.load()
- ✅ Results saved with util.save() where appropriate
- ✅ Cells can run independently (tested)
- ✅ Input data in data/ directory
- ✅ Models in models/ directory
- ✅ Genomes in genomes/ directory
- ✅ Non-JSON output in nboutput/ directory
- ✅ JSON output handled by util.save() to datacache/
- ✅ Markdown cells explain reasoning and purpose
- ✅ All imports in util.py, not scattered in cells
- ✅ Utility functions documented with docstrings

## Error Handling

Handle these scenarios gracefully:

1. **Missing Dependencies**: If KBUtilLib or ModelSEEDpy not available, note in errors
2. **Existing Files**: Don't overwrite util.py if it already exists; validate instead
3. **Non-JSON Data**: Guide user to save to nboutput/ and load manually
4. **Complex Analysis**: Break into multiple cells for independence
5. **Long-Running Cells**: Emphasize saving intermediate results

## Privacy and Security Considerations

- Don't include API keys or credentials in util.py or notebooks
- Use environment variables or config files for sensitive data
- Document if manual credential setup is needed
- Don't log sensitive data in datacache/ files
- Note if data files contain sensitive information

## Example Workflow

For a typical metabolic modeling notebook:

1. **Cell 1**: Load genome data from genomes/
   - Markdown: Explain which genome and why
   - Code: Load, parse, save processed genome data

2. **Cell 2**: Load COBRA model from models/
   - Markdown: Explain model selection and purpose
   - Code: Load model, save to datacache

3. **Cell 3**: Load experimental data from data/
   - Markdown: Describe experimental conditions
   - Code: Load CSV, process, save data structure

4. **Cell 4**: Run flux balance analysis
   - Markdown: Explain FBA parameters and objectives
   - Code: Load model, run FBA, save results

5. **Cell 5**: Generate result tables
   - Markdown: Describe what tables show
   - Code: Load FBA results, create tables, save to nboutput/

Each cell independent, each with clear purpose, each properly cached.
