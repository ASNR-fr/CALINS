# Loading Covariance Data

The `NDCovariances` class provides a unified interface for loading and managing nuclear data covariance matrices from various file formats. This is the **recommended method** for handling covariance data in CALINS.

## Basic Usage

### Auto-detect Format

The easiest way to load covariance data is to use automatic format detection:

```python
import calins as cl

# Load covariance data with auto-detection
cov_data = cl.NDCovariances(input_path='path/to/covariance', format='auto')
```

### Specify Format Explicitly

You can also specify the format explicitly for better control:

```python
import calins as cl

# Load SCALE COVERX binary format
cov_data = cl.NDCovariances(input_path='scale_44g', format='coverx')

# Load SCALE COVERX text format
cov_data = cl.NDCovariances(input_path='scale_44g.txt', format='coverx_text')

# Load COMAC format
cov_data = cl.NDCovariances(input_path='path/to/comac_folder', format='comac')

# Load GENDF format
cov_data = cl.NDCovariances(input_path='gendf_file', format='gendf')

# Load from Excel (previously exported)
cov_data = cl.NDCovariances(input_path='covariance.xlsx', format='xlsx')
```

Supported formats: `'coverx'`, `'coverx_text'`, `'comac'`, `'gendf'`, `'xlsx'`, `'auto'`

## Accessing Data and Metadata

Once loaded, you can access the covariance data and metadata:

```python
import calins as cl

cov_data = cl.NDCovariances(input_path='scale_44g', format='coverx')

# Access the DataFrame containing covariance data
cov_df = cov_data.cov_dataf

# Access metadata
print(f"Format: {cov_data.format}")
print(f"Energy groups: {cov_data.group_nb}")
print(f"Energy bins: {cov_data.e_bins}")
print(f"Number of iso-reac pairs: {len(cov_data.iso_reac_list)}")

# Display first few iso-reac pairs
print("Iso-reac pairs present in covariance data:")
for iso, reac in cov_data.iso_reac_list[:5]:
    print(f"  Isotope: {iso}, Reaction: {reac}")
```

## Exporting Covariance Data

### Export to Excel

```python
import calins as cl

cov_data = cl.NDCovariances(input_path='scale_44g', format='coverx')

# Export to Excel
cov_data.write_xlsx('covariance_export.xlsx')
```

### Export to COVERX Text Format

```python
import calins as cl

cov_data = cl.NDCovariances(input_path='scale_44g', format='coverx')

# Export to COVERX text format
cov_data.write_txt('covariance_export.txt')
```

## Round-trip Export/Import

The Excel format supports full round-trip export/import cycles, preserving all data:

```python
import calins as cl

# Load original data
cov_data = cl.NDCovariances(input_path='scale_44g', format='coverx')

# Export to Excel
cov_data.write_xlsx('covariance_backup.xlsx')

# Later, reload from Excel
cov_reloaded = cl.NDCovariances(input_path='covariance_backup.xlsx', format='xlsx')

# Verify data integrity
assert cov_reloaded.group_nb == cov_data.group_nb
assert len(cov_reloaded.iso_reac_list) == len(cov_data.iso_reac_list)
```

## Using with CALINS Functions

The `NDCovariances` object can be used directly with all CALINS functions:

```python
import calins as cl

# Load covariance data
cov_data = cl.NDCovariances(input_path='scale_44g', format='auto')

# Calculate uncertainty
uncertainty = cl.calcul_uncertainty(
    study_case='path/to/case.sdf',
    cov_data=cov_data  # Use NDCovariances object directly
)

# Calculate Ck similarity index
case1 = cl.Case('path/to/case1.sdf')
case2 = cl.Case('path/to/case2.sdf')
Ck = cl.calcul_Ck(case1, case2, cov_data)

# Perform data assimilation
assimilation = cl.Assimilation(
    study_case='path/to/study.sdf',
    benchmarks_list=['bench1.sdf', 'bench2.sdf'],
    cov_data=cov_data  # NDCovariances object (recommended)
)
```