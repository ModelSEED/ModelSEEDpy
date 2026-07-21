# venvman CLI Reference

## Environment Commands

### list
List all virtual environments in `$VIRTUAL_ENVIRONMENT_DIRECTORY`.

```bash
venvman list
```

**Output:** One environment name per line (e.g., `myapp-py3.12`)

---

### create
Create a new virtual environment and generate `activate.sh` in the project directory.

```bash
venvman create --project <name> --dir <path> [--python <version>] [--force] [--install-deps]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--project` | Yes | Project name (becomes part of env folder name) |
| `--dir` | Yes | Path to project directory |
| `--python` | No | Python version (e.g., `3.12`). Uses pyenv or system python |
| `--force` | No | Replace existing activate.sh |
| `--install-deps` | No | Install from requirements.txt/pyproject.toml |

**Examples:**
```bash
# Basic
venvman create --project myapp --dir ~/projects/myapp

# With Python version
venvman create --project myapp --dir ~/projects/myapp --python 3.12

# With dependencies
venvman create --project myapp --dir . --python 3.12 --install-deps
```

---

### delete
Delete a virtual environment from centralized storage.

```bash
venvman delete --project <name>
venvman delete --env <exact-name>
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--project` | Either | Project name (prompts if multiple versions) |
| `--env` | Either | Exact environment name (e.g., `myapp-py3.12`) |

**Note:** Prompts for confirmation. Does NOT remove symlinks from project directories.

---

### info
Display information about a virtual environment.

```bash
venvman info --project <name>
venvman info --env <exact-name>
```

**Output:**
```
Environment: myapp-py3.12
Path:        /Users/user/VirtualEnvironments/myapp-py3.12
Python:      3.12
Interpreter: /usr/bin/python3.12
Size:        45.2 MB
Created:     2025-01-13 10:30:00
```

---

## Configuration Commands

### setenv
Set `VIRTUAL_ENVIRONMENT_DIRECTORY` in shell configuration.

```bash
venvman setenv <directory>
```

**What it does:**
1. Creates directory if it doesn't exist
2. Updates `~/.bash_profile` or `~/.bashrc`
3. Sets both `VIRTUAL_ENVIRONMENT_DIRECTORY` and `VENVMAN_DIRECTORY`

**Example:**
```bash
venvman setenv ~/VirtualEnvironments
source ~/.bash_profile
```

---

### set_home
Set environment directory with optional migration of existing environments.

```bash
venvman set_home <directory>
```

**What it does:**
1. Prompts to migrate existing environments
2. Copies environments to new location
3. Optionally deletes old environments
4. Updates shell configuration

---

## Project Tracking Commands

### listprojects
List all tracked projects with status.

```bash
venvman listprojects
```

**Output:**
```
Tracked projects (5):

  [ok] myapp
      Path: /Users/user/projects/myapp
      Venv: myapp-py3.12
      Deps: 2025-01-13T10:30:00

  [!] oldproject
      Path: /Users/user/projects/old
      Venv: oldproject-py3.11
      Warning: Project path not found
```

Status `[ok]` = path and venv exist, `[!]` = issue detected

---

### addproject
Add a project directory to tracking.

```bash
venvman addproject <directory> [--project <name>] [--venv <venv-name>]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `directory` | Yes | Path to project directory |
| `--project` | No | Project name (defaults to directory name) |
| `--venv` | No | Virtual environment subdirectory name |

**Examples:**
```bash
# Auto-detect project name from directory
venvman addproject ~/projects/myapp

# Specify project name
venvman addproject ~/projects/myapp --project my-app

# Link to existing environment
venvman addproject ~/projects/myapp --venv myapp-py3.12
```

---

### removeproject
Remove a project from tracking (does not delete files).

```bash
venvman removeproject <project-name>
```

---

### bootstrap
Import existing environments into tracking.

```bash
venvman bootstrap
```

**What it does:**
1. Scans `$VIRTUAL_ENVIRONMENT_DIRECTORY` for environment directories
2. Parses names matching `<project>-py<version>` pattern
3. Adds to projects.json with `path: null`

**Note:** After bootstrap, use `addproject` to set project paths.

---

### update
Regenerate `activate.sh` scripts in all tracked projects.

```bash
venvman update
```

**What it does:**
1. Iterates through tracked projects
2. Removes `.venv` symlinks (no longer needed)
3. Regenerates `activate.sh` with current template

---

### installdeps
Install dependencies from requirements.txt.

```bash
venvman installdeps --project <name>
```

**Requirements:**
- Project must be tracked
- `VIRTUAL_ENVIRONMENT_DIRECTORY` must be set
- `requirements.txt` must exist in project directory

---

### help
Display full README documentation.

```bash
venvman help
```
