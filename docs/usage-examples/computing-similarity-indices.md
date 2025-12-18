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
E_index = cl.calcul_E(case1, case2, reac_list=[2, 4, 16, 18, 102, 103, 452, 1018], exclude_iso="Pu239")
print(f"E similarity index: {E_index}")

# Calculate C_k index (weighted by covariances)
Ck_index = cl.calcul_Ck(case1, case2, cov_df, reac_list=[2, 4, 16, 18, 101, 452, 1018], exclude_iso="Pu239")
print(f"C_k similarity index: {Ck_index}")

# Calculate SS overlap index (Shared Sensitivity)
SSR_index = cl.calcul_SSR(study_case=case1, bench_case=case2, reference=case1, reac_list=[2, 4, 16, 18, 101, 452, 1018])
print(f"SS overlap index: {SSR_index}")
```