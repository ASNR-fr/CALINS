# Git Hooks

This folder contains the project's Git hook scripts.

## Installation

After cloning the repository, run the appropriate setup script **once**:

### On Linux/Mac:
```bash
bash git-hooks/setup-hooks.sh
```

### On Windows (PowerShell):
```powershell
.\git-hooks\setup-hooks.ps1
```

This command configures Git to use hooks from this folder.

## Available Hooks

### `pre-commit`
Runs automatically **before each commit**.

**Features:**
- Automatically generates the version number
- Format: `{baseVersion}.dev0+{hash}.{date}`
- Updates `calins/src/version.py`
- Includes the previous commit hash
