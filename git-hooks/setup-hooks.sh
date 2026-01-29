#!/bin/bash
# Setup script for git hooks
# Run this once after cloning the repo to configure git hooks

echo "Setting up git hooks..."
git config core.hooksPath git-hooks
echo "✅ Git hooks configured successfully!"
echo "The pre-push hook will now run automatically before each push."
