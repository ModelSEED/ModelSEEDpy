---
name: Create New Project
description: Bootstrap a new project with directory structure, git, and Claude configuration
scope: universal
---

# Command: create-new-project

## Purpose

Create a new project with complete setup including Cursor workspace configuration, virtual environment management, Claude commands installation, and optional git repository initialization. This command orchestrates multiple setup steps to create a fully configured development environment.

## Command Type

`create-new-project`

## Input

You will receive a request file containing:
- Project name (required)
- Project directory path (required - can be relative or absolute)
- Project type (optional: python, javascript, typescript, jupyter, multi-language)
- Initialize git repository (optional: boolean, default true)
- Python version for venv (optional: e.g., "3.11", default to system python3)
- Additional workspace folders (optional)
- Workspace settings preferences (optional)

## Process

### Phase 1: Create Project Directory

1. **Setup Project Directory**
   - Create project directory at specified path if it doesn't exist
   - Convert to absolute path for consistency
   - Verify write permissions
   - Document the project path

2. **Validate Project Name**
   - Use provided project name or derive from directory name
   - Sanitize for use in filenames (remove special characters)
   - Check for conflicts with existing projects
   - Document the final project name

### Phase 2: Initialize Git Repository (Optional)

3. **Git Initialization**
   - If git initialization requested (default: true):
     - Run: git init
     - Create .gitignore file with common patterns
     - Create initial commit with project structure
     - Document git initialization status
   - If git initialization skipped:
     - Note in comments why it was skipped
     - Continue with setup

4. **Create .gitignore**
   - Add common patterns based on project type:
     - Python: venv/, __pycache__/, *.pyc, .pytest_cache/, *.egg-info/
     - JavaScript/Node: node_modules/, dist/, .cache/
     - Jupyter: .ipynb_checkpoints/, notebooks/datacache/
     - General: .DS_Store, .vscode/, *.swp, *.swo
     - Claude: .claude/commands/ (managed by claude-commands)
   - Keep: .claude/CLAUDE.md, .claude/settings.local.json
   - Document .gitignore creation

### Phase 3: Setup Virtual Environment with venvman

5. **Register with venvman**
   - Run: venvman add PROJECT_NAME PROJECT_PATH
   - If Python version specified, create venv with that version
   - If not specified, use system default python3
   - Document venvman registration
   - Note the virtual environment path

6. **Activate and Setup Python Environment**
   - If project type is Python or Jupyter:
     - Install basic dependencies (pip, setuptools, wheel)
     - Create requirements.txt if it doesn't exist
     - Document Python setup
   - If not Python project:
     - Note that venv was created but is optional

### Phase 4: Install Claude Commands

7. **Register with claude-commands**
   - Run: claude-commands addproject PROJECT_PATH
   - This installs SYSTEM-PROMPT.md to .claude/CLAUDE.md
   - Installs all command files to .claude/commands/
   - Document claude-commands registration
   - Count and list installed commands

### Phase 5: Create Cursor Workspace

8. **Generate Workspace File**
   - Create workspace file: EXCLAMATION + project-name.code-workspace
   - Include current directory as primary folder
   - Configure settings based on project type
   - Add file exclusions (venv/, node_modules/, __pycache__, etc.)
   - Add search exclusions for performance
   - Document workspace creation

9. **Configure Project-Specific Settings**
   - Python projects: Black formatter, pytest, type checking
   - JavaScript/TypeScript: Prettier, ESLint
   - Jupyter: Notebook settings, output limits
   - Add extension recommendations
   - Document all settings configured

### Phase 6: Create Project Structure

10. **Create Standard Directories**
    - ALWAYS create: agent-io/ directory for Claude command tracking files
    - Based on project type, create:
      - Python: src/, tests/, docs/
      - Jupyter: notebooks/, notebooks/data/, notebooks/datacache/, notebooks/genomes/, notebooks/models/, notebooks/nboutput/, notebooks/util.py
      - JavaScript: src/, tests/, dist/
      - General: docs/, README.md
    - Document directory structure created

11. **Create Initial Files**
    - README.md with project name and description
    - requirements.txt (for Python projects)
    - package.json (for JavaScript projects)
    - For Jupyter: notebooks/util.py with NotebookUtil template
    - Document files created

### Phase 7: Finalize Setup

12. **Create Initial Git Commit (if git enabled)**
    - Stage all created files
    - Create commit: "Initial project setup: PROJECT_NAME"
    - Include setup details in commit message
    - Document commit creation

13. **Generate Setup Summary**
    - List all tools registered (venvman, claude-commands)
    - List all files and directories created
    - Provide next steps for user
    - Document complete setup status

### Phase 8: Save Structured Output

14. **Save JSON Tracking File**
    - IMPORTANT: Save all agent-io output to the NEW project directory, NOT the current working directory
    - Create agent-io/ directory in the new project if it doesn't exist
    - Save tracking JSON to: NEW_PROJECT_PATH/agent-io/create-new-project-session-SESSIONID.json
    - Document all setup steps completed
    - List all files and directories created
    - Record all command executions
    - Note any errors or warnings
    - Include completion status

## JSON Output Schema

```json
{
  "command_type": "create-new-project",
  "status": "complete | incomplete | user_query | error",
  "session_id": "string",
  "parent_session_id": "string | null",
  "session_summary": "Brief summary of project creation",

  "project": {
    "name": "string - project name",
    "path": "string - absolute path to project",
    "type": "python | javascript | typescript | jupyter | multi-language | other"
  },

  "git": {
    "initialized": true,
    "initial_commit": true,
    "commit_hash": "string - git commit hash",
    "gitignore_created": true
  },

  "venvman": {
    "registered": true,
    "command_run": "venvman add PROJECT_NAME PROJECT_PATH",
    "venv_path": "string - path to virtual environment",
    "python_version": "3.11"
  },

  "claude_commands": {
    "registered": true,
    "command_run": "claude-commands addproject .",
    "commands_installed": 5,
    "system_prompt_installed": true,
    "commands_list": ["create-prd", "doc-code-for-dev", "doc-code-usage", "jupyter-dev", "cursor-setup"]
  },

  "workspace": {
    "filename": "string - workspace file with ! prefix",
    "path": "string - absolute path to workspace file",
    "folders_count": 1,
    "settings_configured": true,
    "extensions_recommended": ["ms-python.python", "ms-toolsai.jupyter"]
  },

  "directories_created": [
    "agent-io/",
    "src/",
    "tests/",
    "docs/",
    "notebooks/",
    "notebooks/data/",
    "notebooks/datacache/",
    "notebooks/genomes/",
    "notebooks/models/",
    "notebooks/nboutput/"
  ],

  "files": {
    "created": [
      {
        "path": "!ProjectName.code-workspace",
        "purpose": "Cursor workspace configuration",
        "type": "config"
      },
      {
        "path": ".gitignore",
        "purpose": "Git ignore patterns",
        "type": "config"
      },
      {
        "path": "README.md",
        "purpose": "Project documentation",
        "type": "documentation"
      },
      {
        "path": "requirements.txt",
        "purpose": "Python dependencies",
        "type": "config"
      },
      {
        "path": ".claude/CLAUDE.md",
        "purpose": "Claude system prompt",
        "type": "documentation"
      },
      {
        "path": "agent-io/create-new-project-session-SESSIONID.json",
        "purpose": "Claude command execution tracking for this session",
        "type": "tracking"
      }
    ],
    "modified": []
  },

  "artifacts": {
    "project_path": "absolute path to project",
    "workspace_file": "path to workspace file",
    "readme_file": "path to README.md",
    "tracking_file": "agent-io/create-new-project-session-SESSIONID.json"
  },

  "next_steps": [
    "Open workspace: code !ProjectName.code-workspace",
    "Activate venv: venvman activate ProjectName",
    "Install dependencies: pip install -r requirements.txt",
    "Start developing!"
  ],

  "comments": [
    "Created project directory at /path/to/project",
    "Created agent-io/ directory for Claude command tracking",
    "Initialized git repository with initial commit",
    "Registered with venvman using Python 3.11",
    "Installed 5 Claude commands to .claude/commands/",
    "Created Cursor workspace with Python settings",
    "Created standard Python project structure (src/, tests/, docs/)",
    "Generated README.md and requirements.txt",
    "Saved tracking JSON to NEW_PROJECT_PATH/agent-io/"
  ],

  "queries_for_user": [],

  "errors": []
}
```

## Command JSON Output Requirements

**Required Fields:**
- `command_type`: "create-new-project"
- `status`: "complete", "user_query", or "error"
- `session_id`: Session ID for this execution
- `session_summary`: Brief summary of project creation
- `project`: Project details (name, path, type)
- `git`: Git initialization status
- `venvman`: Virtual environment registration
- `claude_commands`: Claude commands registration
- `workspace`: Cursor workspace details
- `directories_created`: List of directories created
- `files`: All files created
- `artifacts`: Key file paths
- `next_steps`: User guidance for next actions
- `comments`: Detailed notes about setup process

**For user_query status:**
- `queries_for_user`: Questions needing clarification
- `context`: Save partial setup state

**Example Comments:**
- "Created new project 'MetabolicModeling' at ~/Projects/MetabolicModeling"
- "Initialized git repository with initial commit (abc123f)"
- "Registered with venvman using Python 3.11 at ~/Projects/MetabolicModeling/venv"
- "Installed 5 Claude commands to .claude/commands/"
- "Created Cursor workspace: !MetabolicModeling.code-workspace"
- "Created Jupyter notebook structure with util.py template"
- "Generated .gitignore with Python and Jupyter patterns"

## .gitignore Template

### Python Projects
```
# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
venv/
env/
ENV/
.venv
pip-log.txt
pip-delete-this-directory.txt
.pytest_cache/
.coverage
htmlcov/
*.egg-info/
dist/
build/

# Jupyter
.ipynb_checkpoints/
notebooks/datacache/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db

# Claude (commands are managed by claude-commands)
.claude/commands/

# Agent-IO (Claude command tracking - keep in git for project history)
# agent-io/ is intentionally tracked

# Keep these
!.claude/CLAUDE.md
!.claude/settings.local.json
```

### JavaScript/Node Projects
```
# Node
node_modules/
npm-debug.log*
yarn-debug.log*
yarn-error.log*
.pnpm-debug.log*
dist/
build/
.cache/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db

# Claude
.claude/commands/

# Agent-IO (keep in git)
# agent-io/ is intentionally tracked

# Keep these
!.claude/CLAUDE.md
!.claude/settings.local.json
```

### Jupyter Projects
```
# Jupyter
.ipynb_checkpoints/
notebooks/datacache/

# Python
__pycache__/
*.py[cod]
venv/
*.egg-info/

# Data (keep structure, ignore large files)
notebooks/data/*.csv
notebooks/data/*.tsv
notebooks/data/*.xlsx
notebooks/genomes/*.fasta
notebooks/genomes/*.gbk
notebooks/models/*.xml
notebooks/models/*.json
notebooks/nboutput/*

# Keep these data directory files
!notebooks/data/.gitkeep
!notebooks/genomes/.gitkeep
!notebooks/models/.gitkeep

# IDE
.vscode/
.DS_Store

# Claude
.claude/commands/

# Agent-IO (keep in git)
# agent-io/ is intentionally tracked

# Keep these
!.claude/CLAUDE.md
```

## README.md Template

```markdown
# PROJECT_NAME

[Brief project description]

## Setup

This project was created with the `create-new-project` Claude command.

### Prerequisites

- Python 3.11+ (or appropriate version)
- venvman for virtual environment management
- claude-commands for Claude Code integration

### Installation

1. Activate the virtual environment:
   ```bash
   venvman activate PROJECT_NAME
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

### Development

Open the Cursor workspace:
```bash
code !PROJECT_NAME.code-workspace
```

### Project Structure

- `agent-io/` - Claude command execution tracking and session history
- `src/` - Source code
- `tests/` - Test files
- `docs/` - Documentation
- `notebooks/` - Jupyter notebooks (if applicable)
- `.claude/` - Claude Code configuration (commands managed by claude-commands)

### Claude Code Integration

This project includes Claude Code integration:
- Command tracking stored in `agent-io/` for project history
- Commands automatically installed to `.claude/commands/` (managed by claude-commands)
- Update commands: `claude-commands update`

## License

[Add license information]
```

## Jupyter util.py Template

For Jupyter projects, create notebooks/util.py:

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
            name="PROJECT_NAME",
            user="chenry",
            retries=5,
            proxy_port=None,
            **kwargs
        )

    # PLACE ALL UTILITY FUNCTIONS NEEDED FOR NOTEBOOKS HERE

# Initialize the NotebookUtil instance
util = NotebookUtil()
```

## Quality Checklist

Before marking complete, verify:
- ✅ Project directory created at specified path
- ✅ agent-io/ directory created in NEW project directory
- ✅ Git repository initialized (if requested)
- ✅ .gitignore created with appropriate patterns (agent-io/ kept in git)
- ✅ Initial git commit created (if git enabled)
- ✅ Registered with venvman successfully
- ✅ Virtual environment created with correct Python version
- ✅ Registered with claude-commands successfully
- ✅ Claude commands and SYSTEM-PROMPT installed to .claude/
- ✅ Cursor workspace file created with exclamation prefix
- ✅ Workspace settings configured for project type
- ✅ Standard directory structure created
- ✅ README.md generated with project info
- ✅ requirements.txt or package.json created (if applicable)
- ✅ For Jupyter: notebooks/util.py created with project name
- ✅ All setup steps documented in comments
- ✅ Tracking JSON saved to NEW_PROJECT_PATH/agent-io/ directory
- ✅ Next steps provided for user

## Error Handling

Handle these scenarios gracefully:

1. **Directory Already Exists**: Ask user whether to use existing or create new name
2. **Git Not Installed**: Skip git initialization, note in comments
3. **venvman Not Found**: Note error, continue with other setup steps
4. **claude-commands Not Found**: Note error, continue with other setup steps
5. **Permission Issues**: Document error and suggest manual fix
6. **Invalid Project Name**: Sanitize name and notify user of changes
7. **Python Version Not Available**: Fall back to system default, note in comments

## Command Execution Order

Critical: Execute commands in this exact order to avoid conflicts:

1. Create project directory
2. Change to project directory
3. Create agent-io/ directory
4. Initialize git (optional)
5. Create .gitignore
6. Register with venvman
7. Register with claude-commands
8. Create workspace file
9. Create directory structure (including agent-io/)
10. Create initial files
11. Create initial git commit (if enabled)
12. Save tracking file to NEW_PROJECT_PATH/agent-io/

## Integration Notes

### venvman Integration
- venvman stores virtual environments centrally
- Command: `venvman add PROJECT_NAME PROJECT_PATH`
- Activate with: `venvman activate PROJECT_NAME`
- List all: `venvman list`

### claude-commands Integration
- Installs commands to .claude/commands/
- Updates can be pulled with: `claude-commands update`
- List tracked projects: `claude-commands list`

### Cursor Workspace
- Workspace file appears at top of directory (! prefix)
- Open with: `code !PROJECT_NAME.code-workspace`
- Settings are project-specific and version-controlled

## Privacy and Security Considerations

- Don't include API keys or credentials in generated files
- .gitignore should exclude sensitive data directories
- README template should not expose internal paths
- Virtual environment paths are local, not in git
- .claude/commands/ excluded from git (managed by claude-commands)
- Keep .claude/CLAUDE.md in git for project-specific settings

## Next Steps After Project Creation

Provide users with clear next steps:

1. **Open Workspace**
   ```bash
   code !PROJECT_NAME.code-workspace
   ```

2. **Activate Virtual Environment**
   ```bash
   venvman activate PROJECT_NAME
   ```

3. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
   # or
   npm install
   ```

4. **Start Development**
   - Begin coding in src/
   - Write tests in tests/
   - Document in docs/
   - For Jupyter: Create notebooks in notebooks/

5. **Commit Changes**
   ```bash
   git add .
   git commit -m "Add initial implementation"
   ```
