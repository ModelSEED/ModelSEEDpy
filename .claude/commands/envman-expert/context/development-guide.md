# venvman Development Guide

## Architecture Overview

```
venvman.py (single file CLI)
    │
    ├── Environment Storage
    │   └── $VIRTUAL_ENVIRONMENT_DIRECTORY/<project>-py<version>/
    │
    ├── Project Tracking
    │   └── data/projects.json
    │
    └── Per-Project Files
        └── activate.sh (generated)
```

## Key Design Decisions

### Single-File CLI
The entire tool is in `venvman.py` for simplicity. No package structure, no dependencies beyond Python stdlib.

### Environment Variables over Symlinks
The current design uses `VIRTUAL_ENVIRONMENT_DIRECTORY` environment variable instead of `.venv` symlinks. This makes `activate.sh` portable across machines.

### Project Tracking
Projects are tracked in `data/projects.json` to enable bulk operations like `update` and to remember the venv-to-project mapping.

## Code Organization

### Entry Point
```python
def main():
    parser = argparse.ArgumentParser(...)
    sub = parser.add_subparsers(dest="cmd", required=True)
    # Register commands...
    args = parser.parse_args()
    args.func(args)
```

### Command Handler Pattern
Each command is a function that receives the parsed `args`:
```python
def command_name(args):
    """Docstring describing the command."""
    # 1. Validate inputs
    # 2. Load state if needed
    # 3. Perform operation
    # 4. Save state if needed
    # 5. Print output
```

### Core Helper Functions

```python
# Paths
def script_dir() -> Path
def data_dir() -> Path
def projects_file() -> Path
def venv_home() -> Path

# Project tracking
def load_projects() -> Dict[str, dict]
def save_projects(projects: Dict[str, dict]) -> None

# Python resolution
def find_python(pyver: str | None) -> Path | None
def python_version_str(python_bin: Path) -> str

# File operations
def run(cmd: list[str]) -> subprocess.CompletedProcess
def ensure_symlink(link: Path, target: Path, force: bool)
def write_activate_sh(repo_dir: Path, venv_subdir: str)
def install_dependencies(env_dir: Path, repo_dir: Path)

# Shell config
def update_shell_rc(var_name: str, new_directory: str) -> bool
```

## Adding a New Command

### Step 1: Create Handler Function
```python
def my_new_command(args):
    """
    Description of what the command does.

    Args:
        args: Parsed command-line arguments
    """
    # Your implementation
    print("Done!")
```

### Step 2: Register in main()
```python
# In main(), after other command registrations:
p_mycommand = sub.add_parser("mycommand", help="Short description")
p_mycommand.add_argument("--option", required=True, help="Option description")
p_mycommand.add_argument("--flag", action="store_true", help="Boolean flag")
p_mycommand.set_defaults(func=my_new_command)
```

### Step 3: Update README
Add documentation for the new command in README.md.

## Common Modifications

### Adding a Field to Project Tracking

1. **Update save location** (in create_env or relevant command):
```python
projects[args.project] = {
    "path": str(repo_dir),
    "venv_subdir": venv_subdir,
    "last_deps_install": last_deps_install,
    "new_field": new_value  # Add here
}
```

2. **Update display** (in listprojects):
```python
new_field = info.get("new_field")
if new_field:
    print(f"      NewField: {new_field}")
```

### Modifying activate.sh Template

Edit the string in `write_activate_sh()`:
```python
def write_activate_sh(repo_dir: Path, venv_subdir: str):
    script = repo_dir / "activate.sh"
    script_content = f'''#!/usr/bin/env bash
# Your modified template here
VENV_SUBDIR="{venv_subdir}"
# ...
'''
    script.write_text(script_content)
    script.chmod(0o755)
```

After modifying, run `venvman update` to regenerate all scripts.

### Adding Validation

```python
def my_command(args):
    # Input validation
    if not args.required_option:
        print("Error: --required-option is required", file=sys.stderr)
        sys.exit(1)

    # Path validation
    path = Path(args.dir).expanduser().resolve()
    if not path.exists():
        print(f"Error: Directory not found: {path}", file=sys.stderr)
        sys.exit(1)

    # Continue with operation...
```

## Error Handling Conventions

- Print errors to `sys.stderr`
- Use `sys.exit(1)` for fatal errors
- Include helpful context in error messages
- Suggest next steps when possible

```python
print(f"Error: Project '{args.project}' not found.", file=sys.stderr)
print("\nUse 'venvman listprojects' to see tracked projects.", file=sys.stderr)
sys.exit(1)
```

## Testing

### Manual Testing
```bash
# Run from the EnvironmentManager directory
cd ~/Dropbox/Projects/EnvironmentManager

# Test commands
python venvman.py list
python venvman.py listprojects
python venvman.py info --project SomeProject

# Test with temporary project
mkdir /tmp/test-project
python venvman.py create --project test --dir /tmp/test-project --python 3.12
python venvman.py delete --project test
```

### Testing activate.sh
```bash
# Test generation
python venvman.py update

# Test activation
cd /path/to/tracked/project
source activate.sh
which python  # Should show venv path
```

## File Format References

### projects.json
```json
{
  "ProjectName": {
    "path": "/absolute/path/to/project",
    "venv_subdir": "ProjectName-py3.12",
    "last_deps_install": "2025-01-13T10:30:00.000000"
  }
}
```

### dependencies.yaml (optional, per-project)
```yaml
dependencies:
  - name: LocalLibrary
    path: ../LocalLibrary
  - name: AnotherLib
    path: /absolute/path/to/lib
```

## Best Practices

1. **Keep it simple** - This is meant to be a straightforward tool
2. **No external dependencies** - Stdlib only
3. **Fail fast** - Validate early, exit on errors
4. **Be explicit** - Print what you're doing
5. **Preserve data** - Don't delete without confirmation
