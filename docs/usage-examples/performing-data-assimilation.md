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

# Load covariance data
cov_df = cl.format_scale_binary_to_dataframe('path/to/covariance')

# Perform assimilation
assimilation = cl.Assimilation(
    study_case=study_case,
    benchmark_list=benchmarks,
    cov_data=cov_df,
    chi2_target=1.5,  # Optional: chi-squared filtering threshold
    Ck_target=0.7,    # Optional: C_k similarity threshold
    output_html_path='assimilation_results.html'
)

# Access results
print(f"Prior uncertainty: {assimilation.prior_uncertainty.value} pcm")
print(f"Posterior uncertainty: {assimilation.post_uncertainty.value} pcm")
print(f"Posterior bias: {assimilation.bias.value} pcm")
print(f"Chi-squared: {assimilation.chi2_initial} → {assimilation.chi2_final}")
```