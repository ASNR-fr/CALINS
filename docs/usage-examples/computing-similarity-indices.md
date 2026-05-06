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
SSR_index = cl.calcul_SSR(appl_case=case1, bench_case=case2, reference=case1, reac_list=[2, 4, 16, 18, 101, 452, 1018])
print(f"SS overlap index: {SSR_index}")
```

## Partial decomposition and sub-list indices

All indices support `return_partial=True` (for E, Ck, G) or `return_decomposition=True` (for SSR), returning a DataFrame with each isotope-reaction contribution.

```python
# Partial decomposition per iso-reac
E_result = cl.calcul_E(case1, case2, reac_list=reac_list, return_partial=True)
# E_result = {"total_value": float, "partial": DataFrame[ISO, REAC, ISO_NAME, REAC_NAME, SIMILARITY]}

Ck_result = cl.calcul_Ck(case1, case2, cov_df, reac_list=reac_list, return_partial=True)
G_result = cl.calcul_G(case1, case2, reac_list=reac_list, return_partial=True)

# SSR decomposition (SIMILARITY_NUMERATOR, SIMILARITY_DENOMINATOR per iso-reac)
SSR_result = cl.calcul_SSR(appl_case=case1, bench_case=case2, reac_list=reac_list, return_decomposition=True)
SSR_decomp = SSR_result['decomposition']

# Compute partial SSR for a user-defined sub-list
sublist_actinides = [(94239, 18), (94239, 102),  # Pu239 fission & capture
                    (94240, 18), (94240, 102),   # Pu240 fission & capture
                    (92238, 18), (92238, 102)]   # U238 fission & capture

sub_df = SSR_decomp[SSR_decomp.apply(lambda row: (row["ISO"], row["REAC"]) in sublist_actinides, axis=1)]
partial_numerator = sub_df["SIMILARITY_NUMERATOR"].sum()
partial_denominator = sub_df["SIMILARITY_DENOMINATOR"].sum()
SSR_partial = partial_numerator / partial_denominator if partial_denominator != 0 else 0.0
print(f"SSR partial (actinides): {SSR_partial:.4f}")
```