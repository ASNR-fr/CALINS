# Installation

## Prerequisites

- **Python**: Version 3.9 or higher
- **Operating System**: Windows, macOS, or Linux
- **Dependencies**: Listed in `requirements.txt` (automatically installed)
  - numpy ≥ 1.24.2
  - pandas ≥ 2.0.0
  - plotly ≥ 5.14.0
  - matplotlib ≥ 3.8.3
  - tabulate ≥ 0.9.0
  - scipy ≥ 1.13.0

## Installing CALINS

### Method 1: From GitHub (Latest Version)

1. **Clone or download the repository**:
   ```bash
   git clone https://github.com/ASNR-fr/CALINS.git
   cd CALINS
   ```

2. **Install the package**:
   ```bash
   pip install .
   ```

### Method 2: From ZIP Archive

If you have downloaded a compressed archive:

```bash
pip install calins-archive.zip
```

### Method 3: Development Installation

For development or contributing, install in editable mode:

```bash
pip install -e .
```

## Verification

Verify that CALINS is installed correctly:

```python
import calins as cl
print(cl.__name__)  # Should print: calins
```

You can also check for auto-completion in your Python IDE/environment to confirm the installation.