# Sensitivity Analysis and Case Screening

CALINS provides advanced capabilities for sensitivity analysis and case screening, allowing users to identify the most sensitive cases and analyze sensitivity patterns across multiple benchmarks.

These functions are particularly useful for:
- **Benchmark selection**: Identifying the most sensitive cases, targeting specific nuclear data
- **Pattern identification**: Spotting systematic sensitivity trends across benchmark suites and understanding which energy ranges drive uncertainties
- **Isotope comparison**: Comparing sensitivity patterns between different isotopes
- **Benchmark characterization**: Understanding the physical characteristics of different experiments
- **Cross-section priority**: Determining which nuclear data require the most attention

## Finding Sensitive Cases

Use `find_sensitive_cases()` to identify cases with significant sensitivities to specific isotope-reaction pairs:

```python
import calins as cl
import pickle

# Path to a pickle containing a dictionary of Case objects
pickle_db = "/PATH/TO/DATABASE/cases.pickle" 

with open(pickle_db, "rb") as pickle_file:
    dict_pickle = pickle.load(pickle_file)

# Selecting cases to consider for analysis - here all MIX-COMP-THERM except series 006
filtered_cases = cl.filter_cases(dict_pickle, to_keep=["MIX-COMP-THERM"], to_exclude=["MIX-COMP-THERM-006"]) 

# Searching for cases sensitive to the fission cross-sections of U-235 and Pu-239 in the thermal domain,
# based on integral sensitivities, not normalized in lethargy.
# Only bins with a sensitivity ≥ 2e-3, representing more than 10% of the integral sensitivity, are considered sensitive.
sensitive_cases = cl.find_sensitive_cases(
    dict_pickle,
    filtered_cases,
    isotopes=[92235, 94239],
    reactions=[18, 18],
    energy_regions=["thermal", "thermal"],
    sensitivity_threshold=2e-3,
    fraction_threshold = 0.10,
    integrate_data=True,
    normalize_lethargy = False,
)

# Displaying in the terminal the 50 most sensitive cases according to the previous search criteria
cl.display_most_sensitive_cases(dict_sensitive = sensitive_cases, max_number_of_cases_displayed=50)
```

## Creating Sensitivity Heatmaps

Use `plot_sensitivity_heatmap()` to create visual heatmaps showing how many cases show the same sensitivity patterns across multiple energy regions:

![](../images/Sensi_heatmap.PNG)

```python
import calins as cl
import pickle

# Path to a pickle containing a dictionary of Case objects
pickle_path = "/PATH/TO/DATABASE/database_dict.pickle"

with open(pickle_path, "rb") as pickle_file:
    dict_pickle = pickle.load(pickle_file)

# Selecting cases to consider for analysis - here all LEU-COMP-THERM except series 001
list_cases = cl.filter_cases(dict_pickle, to_keep=[], to_exclude=[]) 

# Creating a heatmap of cases sensitive to cross-sections of SCATTER (0), TOTAL (1), ELASTIC (2), N,GAMMA (102)
# Only bins with a sensitivity ≥ 2e-3 and representing more than 10% of the integral sensitivity will be considered sensitive.
cl.plot_sensitivity_heatmap(
    dict_data=dict_pickle,
    cases=list_cases,
    isotopes=[17035, 17037],
    reactions=[0, 1, 2, 102],
    sensitivity_threshold=2e-3,
    fraction_threshold = 0.10,
    user_energy_bins=None,
    title="Cases pickle db sensitivity heatmap for chlorine isotopes",
    namefig=None,
    save_logfile=True,
)
```