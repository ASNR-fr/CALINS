# Performing Data Assimilation

Assimilate experimental benchmark data using GLLSM:

```python
import calins as cl

# Define study case
study_case = cl.Case('path/to/study_case.sdf')

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
    study_case=study_case,
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
```