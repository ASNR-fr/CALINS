# Computing Similarity Indices

Evaluate similarity between different cases:

```python
import calins as cl

# Load two cases to compare
case1 = cl.Case('path/to/case1.sdf')
case2 = cl.Case('path/to/case2.sdf')

# Load covariance data
cov_df = cl.format_scale_binary_to_dataframe('path/to/covariance')

# Calculate E similarity index (0 to 1)
E_index = cl.calcul_E(case1, case2)
print(f"E similarity index: {E_index}")

# Calculate C_k index (weighted by covariances)
Ck_index = cl.calcul_Ck(case1, case2, cov_df)
print(f"C_k similarity index: {Ck_index}")

# Calculate SS overlap index (Shared Sensitivity)
SS_index = cl.calcul_SS(study_case=case1, bench_case=case2, reference=case1)
print(f"SS overlap index: {SS_index}")
```