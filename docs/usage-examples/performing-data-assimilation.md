# Performing Data Assimilation

Assimilate experimental benchmark data using GLLSM:

```python
import calins as cl

# Define application case
appl_case = cl.Case('path/to/appl_case.sdf')

# Define benchmark cases
benchmarks = [
    cl.Case('path/to/benchmark1.sdf'),
    cl.Case('path/to/benchmark2.sdf'),
    cl.Case('path/to/benchmark3.sdf')
]

# Load covariance data (recommended: use NDCovariances object)
cov_data = cl.NDCovariances(input_path='path/to/covariance', format='auto')  # auto-detect format

# Perform assimilation
assimilation = cl.Assimilation(
    appl_case=appl_case,
    benchmarks_list=benchmarks,
    cov_data=cov_data,  # NDCovariances object
    reac_list=[2, 4, 16, 18, 101, 452, 1018],
    targetted_chi2=1.2,  # Optional: chi-squared filtering threshold
    Ck_threshold=0.8,    # Optional: C_k similarity threshold
    output_html_path='assimilation_results.html'
)

# Access results
print(f"Prior uncertainty: {assimilation.prior_uncertainty.value} pcm")
print(f"Posterior uncertainty: {assimilation.post_uncertainty.value} pcm")
print(f"Posterior bias: {assimilation.bias.value} pcm")
print(f"Chi-squared prior: {assimilation.prior_chi2}")
print(f"Chi-squared post: {assimilation.post_chi2}")

# Access validation methods results (each is a dictionary)
# USL = 1 - CM - MOS, where MOS is the Margin of Subcriticality (default 0.05)
print(f"MOS (Margin of Subcriticality): {assimilation.MOS}")

# GLLSM USL
print(f"GLLSM USL: {assimilation.USL_gllsm['USL']:.5f}")
print(f"  Calculational Margin: {assimilation.USL_gllsm['calculational_margin']*1e5:.0f} pcm")
print(f"  Coverage factor K: {assimilation.USL_gllsm['K']:.4f}")

# Parametric USL (assumes normal distribution of C/E)
print(f"Parametric USL: {assimilation.USL_parametric['USL']:.5f}")
print(f"  Bias beta: {assimilation.USL_parametric['beta']:.6f}")
print(f"  Normality test passed: {assimilation.USL_parametric['normality_passed']}")

# Nonparametric USL (distribution-free, based on worst-case C/E)
# Note: USL and calculational_margin can be None if CNP <= 0.4 (not enough benchmarks)
print(f"Nonparametric USL: {assimilation.USL_nonparametric['USL']}")
print(f"  CNP: {assimilation.USL_nonparametric['CNP']:.4f}")
```

## Accounting for experimental correlations

When benchmarks share common experimental biases (e.g. same fission chamber, same facility), their experimental uncertainties are correlated. Use `expe_correlations` to pass a symmetric correlation matrix:

```python
import numpy as np

# Correlation matrix (N x N), values in [0, 1], must be symmetric.
# IMPORTANT: rows/columns must follow the same order as benchmarks_list.
expe_corr = np.array([
    [1.0, 0.3, 0.1],   # bench1 vs bench1, bench2, bench3
    [0.3, 1.0, 0.2],   # bench2 vs ...
    [0.1, 0.2, 1.0]    # bench3 vs ...
])

assimilation = cl.Assimilation(
    appl_case=appl_case,
    benchmarks_list=benchmarks,
    cov_data=cov_data,
    expe_correlations=expe_corr,  # populates off-diagonal of C_bench
    output_html_path='assimilation_correlated.html'
)
```