# Git Hooks

Ce dossier contient les scripts git hooks du projet.

## Installation

Après avoir cloné le repo, exécutez **une seule fois** le script de configuration approprié :

### Sur Linux/Mac :
```bash
bash git-hooks/setup-hooks.sh
```

### Sur Windows (PowerShell) :
```powershell
.\git-hooks\setup-hooks.ps1
```

Cette commande configure Git pour utiliser les hooks depuis ce dossier.

## Hooks disponibles

### `pre-commit`
S'exécute automatiquement **avant chaque commit**.

**Fonctionnalité :**
- Génère automatiquement le numéro de version
- Format : `{baseVersion}.dev0+{hash}.{date}`
- Met à jour `calins/src/version.py`
- Inclut la version dans le commit (pas de commit supplémentaire)

## Notes

- Les hooks doivent être exécutables (sur Linux/Mac, c'est automatique)
- Cette configuration est locale à votre clone du repo
- La configuration persiste même après des `git pull`
