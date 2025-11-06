# Creating a Case Object

A `Case` object represents a study case or benchmark with sensitivity data:

```python
import calins as cl

# Path to your sensitivity file (.sdf format)
sensi_file_path = 'path/to/sensitivity.sdf'

# Create Case object
# occurrences_rule: how to handle multiple occurrences of same iso-reac pair
#   Options: "first", "last", "sum" (default)
my_case = cl.Case(sdf_path=sensi_file_path, occurrences_rule="sum")

# Generate HTML visualization of sensitivities
my_case.plot_case_sensi(
    output_html_path='case_sensitivities.html',
    plotting_unit='pcm'  # Options: 'pcm' or 'relative'
)

# Access case attributes
print(f"Case name: {my_case.casename}")
print(f"Energy groups: {my_case.group_nb}")
print(f"Calculated k_eff: {my_case.resp_calc} ± {my_case.sigma_resp_calc}")
if my_case.resp_expe:
    print(f"Experimental k_eff: {my_case.resp_expe} ± {my_case.sigma_resp_expe}")
```