---
name: Cursor Setup
description: Configure Cursor IDE with Claude Code integration and project settings
scope: universal
---

# Command: cursor-setup

## Purpose

Create a Cursor workspace file for the current project directory, enabling multi-root workspace features, custom settings, and organized project management in Cursor IDE.

## Command Type

`cursor-setup`

## Input

You will receive a request file containing:
- Project name (required)
- Additional workspace folders to include (optional)
- Workspace-specific settings (optional)
- Extensions to recommend (optional)

## Process

### Phase 1: Gather Project Information

1. **Determine Project Name**
   - Use project name from input request
   - If not provided, derive from current directory name
   - Sanitize name for filename use (remove special characters)
   - Document the project name

2. **Identify Project Structure**
   - Examine current directory structure
   - Identify key folders (src, tests, docs, etc.)
   - Note any existing configuration files (.vscode, .cursor, etc.)
   - Document project type (Python, Node.js, multi-language, etc.)

### Phase 2: Create Workspace File

3. **Generate Workspace Configuration**
   - Create workspace file with naming pattern: EXCLAMATION-project-name.code-workspace
   - The exclamation mark prefix ensures the file appears at top of directory listings
   - Include current directory as primary folder
   - Add any additional folders specified in request
   - Configure workspace settings appropriate for project type

4. **Configure Workspace Settings**
   - Add workspace-level settings for:
     - File associations
     - Editor preferences
     - Language-specific settings
     - Search exclusions
     - Extension recommendations
   - Preserve any existing settings from .vscode/settings.json
   - Document all settings added

### Phase 3: Register with ClaudeCommands

5. **Add Project to ClaudeCommands Database**
   - Run the command: claude-commands addproject .
   - This registers the project directory in the ClaudeCommands tracking system
   - Installs the latest Claude commands and SYSTEM-PROMPT.md to the project
   - Document the registration in comments
   - If the command fails, note the error but continue with workspace setup

### Phase 4: Validate and Document

6. **Validate Workspace File**
   - Verify JSON structure is valid
   - Ensure all paths are relative to workspace file location
   - Check that workspace file can be opened in Cursor
   - Document workspace structure

7. **Create Documentation**
   - Document workspace file location
   - Explain workspace structure
   - List any workspace-specific settings
   - Provide usage instructions

### Phase 5: Save Structured Output

8. **Save JSON Tracking File**
   - Document workspace file creation
   - List all settings configured
   - Note any issues or recommendations
   - Include completion status

## Workspace File Template

The workspace file should follow this structure:

```json
{
  "folders": [
    {
      "path": ".",
      "name": "<Project Name>"
    }
  ],
  "settings": {
    "files.exclude": {
      "**/__pycache__": true,
      "**/*.pyc": true,
      "**/.pytest_cache": true,
      "**/.DS_Store": true,
      "**/node_modules": true,
      "**/.git": false
    },
    "search.exclude": {
      "**/__pycache__": true,
      "**/*.pyc": true,
      "**/node_modules": true,
      "**/.git": true
    },
    "files.watcherExclude": {
      "**/__pycache__/**": true,
      "**/node_modules/**": true
    }
  },
  "extensions": {
    "recommendations": []
  }
}
```

### Workspace Settings by Project Type

**Python Projects:**
```json
{
  "python.analysis.typeCheckingMode": "basic",
  "python.analysis.autoImportCompletions": true,
  "[python]": {
    "editor.defaultFormatter": "ms-python.black-formatter",
    "editor.formatOnSave": true,
    "editor.codeActionsOnSave": {
      "source.organizeImports": true
    }
  },
  "files.exclude": {
    "**/__pycache__": true,
    "**/*.pyc": true,
    "**/.pytest_cache": true
  }
}
```

**Node.js/JavaScript Projects:**
```json
{
  "[javascript]": {
    "editor.defaultFormatter": "esbenp.prettier-vscode",
    "editor.formatOnSave": true
  },
  "[typescript]": {
    "editor.defaultFormatter": "esbenp.prettier-vscode",
    "editor.formatOnSave": true
  },
  "files.exclude": {
    "**/node_modules": true,
    "**/dist": true,
    "**/.cache": true
  }
}
```

**Jupyter Notebook Projects:**
```json
{
  "jupyter.notebookFileRoot": "${workspaceFolder}/notebooks",
  "notebook.output.textLineLimit": 500,
  "[python]": {
    "editor.defaultFormatter": "ms-python.black-formatter"
  },
  "files.exclude": {
    "**/.ipynb_checkpoints": true,
    "**/__pycache__": true
  }
}
```

## JSON Output Schema

```json
{
  "command_type": "cursor-setup",
  "status": "complete | incomplete | user_query | error",
  "session_id": "string",
  "parent_session_id": "string | null",
  "session_summary": "Brief summary of workspace setup",

  "project": {
    "name": "string - project name",
    "type": "python | javascript | typescript | jupyter | multi-language | other",
    "workspace_filename": "string - filename with ! prefix"
  },

  "workspace": {
    "folders": [
      {
        "path": "string - relative path",
        "name": "string - folder display name"
      }
    ],
    "settings_count": 10,
    "extensions_recommended": 3
  },

  "claude_commands": {
    "registered": true,
    "command_run": "claude-commands addproject .",
    "commands_installed": 5,
    "system_prompt_installed": true
  },

  "files": {
    "created": [
      {
        "path": "string - workspace file with ! prefix",
        "purpose": "Cursor workspace configuration",
        "type": "config"
      }
    ],
    "modified": []
  },

  "artifacts": {
    "workspace_filename": "string - workspace file with ! prefix",
    "workspace_path": "absolute path to workspace file"
  },

  "comments": [
    "Created workspace file with name prefix '!' for top sorting",
    "Configured Python-specific settings for project",
    "Added file exclusions for __pycache__ and .pyc files",
    "Workspace can be opened in Cursor via File > Open Workspace",
    "Registered project with ClaudeCommands database",
    "Installed 5 Claude commands to .claude/commands/"
  ],

  "queries_for_user": [],

  "errors": []
}
```

## Command JSON Output Requirements

**Required Fields:**
- `command_type`: "cursor-setup"
- `status`: "complete", "user_query", or "error"
- `session_id`: Session ID for this execution
- `session_summary`: Brief summary of workspace creation
- `project`: Project name and workspace details
- `workspace`: Configuration details
- `claude_commands`: Registration status with ClaudeCommands database
- `files`: Workspace file created
- `artifacts`: Path to workspace file
- `comments`: Notes about workspace configuration

**For user_query status:**
- `queries_for_user`: Questions about project structure or preferences
- `context`: Save partial workspace configuration

**Example Comments:**
- "Created workspace file with exclamation prefix for top sorting"
- "Configured Python development settings with Black formatter"
- "Added exclusions for common Python cache directories"
- "Included notebooks/ folder as additional workspace folder"
- "Recommended extensions: Python, Jupyter, Black Formatter"

## Workspace File Naming Convention

The workspace file must be named with an exclamation mark prefix followed by the project name and .code-workspace extension.

Format: EXCLAMATION + project-name + .code-workspace

**Why the exclamation mark prefix?**
- Ensures workspace file appears at top of alphabetical directory listings
- Makes workspace file easy to find and identify
- Common convention for important configuration files
- Visual indicator of workspace root file

**Examples:**
- Exclamation mark + MetabolicModeling.code-workspace
- Exclamation mark + ClaudeCommands.code-workspace
- Exclamation mark + WebsiteRedesign.code-workspace

## Quality Checklist

Before marking complete, verify:
- ✅ Workspace file created with exclamation mark prefix in filename
- ✅ JSON structure is valid and properly formatted
- ✅ Current directory included as primary folder
- ✅ Workspace settings appropriate for project type
- ✅ File exclusions configured to hide build artifacts
- ✅ Search exclusions configured for better performance
- ✅ Extension recommendations included (if applicable)
- ✅ All paths are relative to workspace file location
- ✅ Workspace file can be opened in Cursor
- ✅ Project registered with ClaudeCommands (claude-commands addproject .)
- ✅ Claude commands and SYSTEM-PROMPT installed to .claude/ directory
- ✅ Documentation includes usage instructions

## Error Handling

Handle these scenarios gracefully:

1. **No Project Name**: Use current directory name as fallback
2. **Existing Workspace File**: Ask user whether to overwrite or merge
3. **Invalid Characters in Name**: Sanitize project name for filename
4. **Unknown Project Type**: Use generic workspace template
5. **Permission Issues**: Document if unable to write file
6. **ClaudeCommands Not Found**: Note error in comments, continue with workspace setup

## Usage Instructions

After creating workspace file, users can:

1. **Open Workspace in Cursor**
   - File > Open Workspace from File
   - Select the workspace file (begins with exclamation mark)
   - Or double-click the workspace file

2. **Benefits of Workspace**
   - Consistent settings across team members
   - Multi-root folder support
   - Workspace-specific extensions
   - Organized project structure
   - Easy project switching

3. **Customization**
   - Edit workspace file to add more folders
   - Add custom tasks and launch configurations
   - Configure language-specific settings
   - Add extension recommendations

## Advanced Workspace Features

Optionally include these advanced features:

**Tasks Configuration:**
```json
{
  "tasks": {
    "version": "2.0.0",
    "tasks": [
      {
        "label": "Run Tests",
        "type": "shell",
        "command": "pytest",
        "group": "test"
      }
    ]
  }
}
```

**Launch Configurations:**
```json
{
  "launch": {
    "version": "0.2.0",
    "configurations": [
      {
        "name": "Python: Current File",
        "type": "python",
        "request": "launch",
        "program": "${file}"
      }
    ]
  }
}
```

## Privacy and Security Considerations

- Don't include absolute paths that expose user directory structure
- Use relative paths for all folder references
- Don't include API keys or credentials in workspace settings
- Don't commit sensitive workspace settings to version control
- Use workspace file for team-shared settings only
