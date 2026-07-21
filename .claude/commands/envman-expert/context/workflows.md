# venvman Common Workflows

## Initial Setup

### New Machine Setup
```bash
# 1. Clone EnvironmentManager (or have it synced via Dropbox)
cd ~/Dropbox/Projects/EnvironmentManager

# 2. Set the environment variable for where venvs will be stored
python venvman.py setenv ~/VirtualEnvironments

# 3. Restart shell or source the profile
source ~/.bash_profile

# 4. Verify it's set
echo $VIRTUAL_ENVIRONMENT_DIRECTORY
# Output: /Users/you/VirtualEnvironments

# 5. If you have existing environments (from sync or backup), bootstrap
python venvman.py bootstrap

# 6. Add wrapper to PATH (optional)
# Add to ~/.bashrc or ~/.zshrc:
alias venvman='python ~/Dropbox/Projects/EnvironmentManager/venvman.py'
```

### Setting Up venvman Alias
```bash
# Option 1: Alias in shell config
echo "alias venvman='python ~/Dropbox/Projects/EnvironmentManager/venvman.py'" >> ~/.bashrc

# Option 2: Wrapper script in ~/bin
cat > ~/bin/venvman << 'EOF'
#!/bin/bash
python ~/Dropbox/Projects/EnvironmentManager/venvman.py "$@"
EOF
chmod +x ~/bin/venvman
```

## Project Workflows

### New Project
```bash
# 1. Create project directory
mkdir ~/projects/myapp
cd ~/projects/myapp

# 2. Initialize git, create requirements.txt, etc.
git init
echo "requests" > requirements.txt

# 3. Create environment with dependencies
venvman create --project myapp --dir . --python 3.12 --install-deps

# 4. Activate
source activate.sh

# 5. Start working
python -c "import requests; print('Success!')"
```

### Existing Project (not yet tracked)
```bash
# If you already have a project with a venv elsewhere:
cd ~/projects/existing-app

# 1. Add to tracking (if venv already exists in VirtualEnvironments)
venvman addproject . --venv existing-app-py3.11

# 2. Or create new environment
venvman create --project existing-app --dir . --python 3.12

# 3. Install dependencies
venvman installdeps --project existing-app
```

### Clone Existing Tracked Project
```bash
# After cloning a repo that was using venvman:
cd ~/projects/cloned-repo

# 1. Check if environment exists
venvman info --project cloned-repo

# 2. If environment exists, just add the project
venvman addproject .

# 3. If not, create it
venvman create --project cloned-repo --dir . --python 3.12 --install-deps
```

## Multi-Python-Version Workflow

### Multiple Versions for Same Project
```bash
# Create Python 3.11 environment
venvman create --project myapp --dir ~/projects/myapp --python 3.11

# Create Python 3.12 environment (won't overwrite the 3.11 one)
venvman create --project myapp --dir ~/projects/myapp --python 3.12

# List both
venvman list | grep myapp
# myapp-py3.11
# myapp-py3.12

# Switch by editing activate.sh or recreating with desired version
venvman create --project myapp --dir ~/projects/myapp --python 3.11
```

## Maintenance Workflows

### Update All Projects
```bash
# Regenerate activate.sh in all tracked projects
venvman update

# This is useful after:
# - Updating venvman itself
# - Changing VIRTUAL_ENVIRONMENT_DIRECTORY
# - Migrating to a new machine
```

### Clean Up Unused Environments
```bash
# 1. List all environments
venvman list

# 2. List tracked projects
venvman listprojects

# 3. Check info on suspicious ones
venvman info --env old-project-py3.10

# 4. Delete unused
venvman delete --env old-project-py3.10
```

### Migrate to New Storage Location
```bash
# 1. Set new home (with migration)
venvman set_home /new/path/to/VirtualEnvironments

# 2. Follow prompts to:
#    - Copy environments to new location
#    - Delete old environments
#    - Update shell config

# 3. Restart shell
source ~/.bash_profile

# 4. Update all project scripts
venvman update
```

## Syncing Across Machines

### Via Dropbox/Cloud Sync

The `data/projects.json` syncs automatically. Environments do NOT sync (too large).

**On each machine:**
```bash
# 1. Set environment variable (same on all machines)
venvman setenv ~/VirtualEnvironments

# 2. Create environments locally
# For each project you need:
venvman create --project ProjectName --dir /path/to/project --python 3.12 --install-deps
```

### Recreating Environments from Tracking
```bash
# View tracked projects
venvman listprojects

# For projects with [!] status (missing venv):
venvman create --project ProjectName --dir /path/to/project --python 3.12 --install-deps
```

## Troubleshooting Workflows

### Fix Broken Project
```bash
# 1. Check status
venvman listprojects
# Look for [!] markers

# 2. If path is wrong, update it
venvman addproject /correct/path --project myapp

# 3. If venv is missing, recreate
venvman create --project myapp --dir /path/to/project --python 3.12

# 4. Reinstall dependencies
venvman installdeps --project myapp
```

### Reset activate.sh
```bash
# If activate.sh is corrupted or outdated:
venvman create --project myapp --dir /path/to/project --python 3.12

# Or for all projects:
venvman update
```

### Environment Variable Not Set
```bash
# Check current value
echo $VIRTUAL_ENVIRONMENT_DIRECTORY

# If empty, set it
venvman setenv ~/VirtualEnvironments

# Then reload shell
exec $SHELL
# or
source ~/.bash_profile
```

## CI/CD Integration

### GitHub Actions Example
```yaml
# In .github/workflows/test.yml
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: '3.12'

      # Don't use venvman in CI - just use standard venv
      - name: Create virtualenv
        run: python -m venv .venv

      - name: Install dependencies
        run: |
          source .venv/bin/activate
          pip install -r requirements.txt
```

**Note:** venvman is designed for local development, not CI. In CI, use standard `python -m venv` or actions/setup-python.
