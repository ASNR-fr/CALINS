# Setup script for git hooks (Windows PowerShell)
# Run this once after cloning the repo to configure git hooks

Write-Host "Setting up git hooks..." -ForegroundColor Green
git config core.hooksPath git-hooks
Write-Host "✅ Git hooks configured successfully!" -ForegroundColor Green
Write-Host "The pre-push hook will now run automatically before each push."
