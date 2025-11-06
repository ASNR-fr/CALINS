# Calculating Uncertainty

Calculate nuclear data uncertainty using the sandwich formula:

```python
import calins as cl

# Load covariance matrix (example: SCALE format)
scale_file_path = 'path/to/scale_44g'
cov_df = cl.format_scale_binary_to_dataframe(input_path=scale_file_path)

# Load sensitivity data
sensi_file_path = 'path/to/sensitivity.sdf'

# Calculate a priori uncertainty
uncertainty = cl.calcul_uncertainty(
    study_case=sensi_file_path,  # Can also be a Case object
    cov_data=cov_df,
    output_html_path='uncertainty_report.html'
)

print(f"A priori uncertainty: {uncertainty.value} pcm")
print(f"Calculated k_eff: {uncertainty.resp_calc}")

# Export decomposition to Excel
uncertainty.decomposition.to_excel('uncertainty_decomposition.xlsx')
```