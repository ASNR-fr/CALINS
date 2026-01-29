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

### `pre-push`
S'exécute automatiquement **avant chaque push**.

**Fonctionnalité :**
- Génère automatiquement le numéro de version avec le hash correct
- Format : `{baseVersion}.dev0+{hash}.{date}`
- Met à jour `calins/src/version.py`
- **Amend** le dernier commit avec les changements de version (pas de commit supplémentaire)
- Tout se fait automatiquement, vous ne voyez rien

## Notes

- Les hooks doivent être exécutables (sur Linux/Mac, c'est automatique)
- Cette configuration est locale à votre clone du repo
- La configuration persiste même après des `git pull`
