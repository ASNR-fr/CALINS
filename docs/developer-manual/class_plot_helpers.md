# Developer Manual - class_plot_helpers.py

`class_plot_helpers.py` contains the shared plotting workflows used by the CALINS domain classes. It sits between the object-specific plotting methods in `classes.py` and the generic Plotly/HTML primitives in `plots.py`.

This module is an internal implementation module. Its functions are intended to keep the public methods of `Case`, `Uncertainty`, and `Assimilation` concise and consistent; users should use those class methods rather than calling these helpers directly.

## Module boundaries

The plotting responsibilities are split as follows:

- `plots.py` provides generic Plotly figure factories and HTML components. It has no knowledge of CALINS objects or their data model.
- `class_plot_helpers.py` prepares CALINS-specific vectors, isotope-reaction selections, units, and reusable sets of figures.
- `classes.py` exposes the public plotting/export methods and assembles the complete HTML reports, including their tables, tabs, and report-specific content.

This separation makes it possible to reuse a figure workflow across multiple classes without moving high-level report composition outside the owning class.

## Output validation and display units

- `validate_image_output_path(output_path, format)` validates a static-image format (`png`, `jpg`, `jpeg`, or `pdf`), appends the corresponding extension when needed, and checks that the destination can be created.
- `validate_html_output_path(output_html_path)` appends the `.html` extension when necessary and checks that the destination can be created. `None` is preserved to support display-only interactive figures.
- `resolve_plotting_unit(plotting_unit, resp_calc, squared=False)` returns the numerical conversion factor and display label for either relative values or pcm. When `squared=True`, the factor is squared for variance-like quantities.

## Domain-data preparation

These helpers centralize operations that are common to plots based on isotope-reaction pairs.

- `normalize_filter_iso_reac_list(...)` normalizes isotope/reaction identifiers, applies optional isotope and reaction filters, and returns both the selected pairs and the normalized requested pairs.
- `align_trace_styles(...)` keeps user-provided colors and dash styles aligned with the pairs remaining after filtering.
- `select_top_iso_reac_pairs(...)` ranks candidate pairs by the absolute value of a decomposition column and returns the requested number of leading contributors.
- `select_top_iso_reac_pairs_by_vector(...)` ranks candidate pairs by the absolute value of their vector integrals, for quantities such as the assimilation nuclear-data adjustment that have no decomposition table.
- `build_vector_from_decomp(...)` extracts and concatenates the energy-group slices associated with selected pairs from a decomposition vector.
- `build_lethargy_vector(...)` builds a vector for profiles expressed per unit lethargy. Missing data can either propagate as an error or be replaced with zero-valued energy-group slices.

## Interactive HTML report figure builders

The following functions produce Plotly figures already configured for their inclusion in an interactive CALINS HTML report. They return figures; the calling class decides how to organize and write the final report.

- `build_interactive_integral_profile_figures(...)` creates the paired group-wise integral and energy-profile figures for one vector.
- `build_interactive_sensitivity_figures(...)` creates sensitivity integral, energy-profile, and per-unit-lethargy profile figures.
- `build_interactive_covariance_figures(...)` creates covariance-integral, diagonal variance-profile, and row-sum covariance-profile figures. Sparse covariance matrices are supported.

## Static publication-image export

`export_static_plotly_figure(...)` applies the shared publication layout (white background, Times New Roman font, grids, margins, and legend style) to a Plotly figure, exports it to an image, and optionally displays it. The public class methods remain responsible for selecting the domain data and figure type before calling this helper.
