# CALINS

### **CAL**culations and **I**nvestigations on **N**uclear data uncertainties and **S**ensitivities

[![Python Version](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD--3--Clause-blue.svg)](LICENSE)

> A Python package for computing nuclear data uncertainty propagation, sensitivity analyses, and data assimilation using the Generalized Linear Least Squares Method (GLLSM).

**Developed by**: French Authority for Nuclear Safety and Radiation Protection (ASNR)

---

## About

CALINS is a specialized Python package designed for nuclear safety and criticality calculations. It enables:
- **Uncertainty propagation** from nuclear data to integral responses
- **Sensitivity analysis** of nuclear systems to cross-section variations
- **Data assimilation** using the Generalized Linear Least Squares Method (GLLSM)
- **Bias estimation** for calculated responses based on experimental benchmarks

The package is particularly useful for criticality safety analyses where understanding uncertainties and biases in calculated k<sub>eff</sub> values is essential.

## Features

- ✅ **Uncertainty Calculation**: Propagate nuclear data uncertainties using the sandwich formula
- ✅ **Sensitivity Analysis**: Process and visualize sensitivity profiles from SDF files
- ✅ **Similarity Indices**: Calculate E, C<sub>k</sub>, G, and SSR indices between cases
- ✅ **GLLSM Assimilation**: Assimilate experimental benchmark data to reduce uncertainties
- ✅ **Multiple Covariance Formats**: Support for SCALE (COVERX binary/text), COMAC, GENDF, and Excel (xlsx) formats, with auto-detection
- ✅ **Interactive Visualizations**: Generate HTML reports with Plotly graphs
- ✅ **Chi-squared Filtering**: Automatic filtering of inconsistent benchmark cases

<div align="center"><img src="./images/Assim.PNG" alt="Assim" style="width:70%;border:1px solid #1d1d53ff"></div>
<div align="left"><img src="./images/Bench_list.PNG" alt="Benchmark List" style="width:70%;border:1px solid #1d1d53ff;"></div>
<div align="right"><img src="./images/Case_sensi.PNG" alt="Case Sensi" style="width:70%;border:1px solid #1d1d53ff;"></div>
<div align="left"><img src="./images/Sensi_heatmap.PNG" alt="Sensi Heatmap" style="width:70%;border:1px solid #1d1d53ff;"></div>
<div align="right"><img src="./images/Iso_reac_list.PNG" alt="Iso Reac List" style="width:70%;border:1px solid #1d1d53ff;"></div>
<div align="left"><img src="./images/Contributions.PNG" alt="Contributions" style="width:70%;border:1px solid #1d1d53ff;"></div>