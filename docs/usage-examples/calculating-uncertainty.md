# Calculating Uncertainty

Calculate nuclear data uncertainty using the sandwich formula:

```python
import calins as cl

# Load covariance matrix (recommended: use NDCovariances object)
scale_file_path = 'path/to/scale_44g'
cov_data = cl.NDCovariances(input_path=scale_file_path, format='coverx')  # or format='auto'

# Load sensitivity data
sensi_file_path = 'path/to/sensitivity.sdf'

# Calculate a priori uncertainty
uncertainty = cl.calcul_uncertainty(
    study_case=sensi_file_path,  # Can also be a Case object
    cov_data=cov_data,  # NDCovariances object
    output_html_path='uncertainty_report.html'
)

print(f"A priori uncertainty: {uncertainty.value} pcm")
print(f"Calculated k_eff: {uncertainty.resp_calc}")

# Export decomposition to Excel
uncertainty.decomposition.to_excel('uncertainty_decomposition.xlsx')

# Optional: Export/import covariance data in Excel format
cov_data.write_xlsx('covariance_data.xlsx')
# Later, reload from Excel:
cov_data_reloaded = cl.NDCovariances(input_path='covariance_data.xlsx', format='xlsx')
```