"""Shared plotting workflows for CALINS domain objects.

``plots`` provides generic Plotly and HTML primitives. This module prepares
class-specific data and combines those primitives into reusable interactive
report blocks or static-export workflows. It does not compose complete HTML
reports; that remains the responsibility of the calling class.
"""

import numpy as np

from . import errors, methods, plots


IMAGE_FORMATS = ["png", "jpg", "jpeg", "pdf"]


# ---------------------------------------------------------------------------
# Output validation and display units
# ---------------------------------------------------------------------------

def validate_image_output_path(output_path, format):
	if format not in IMAGE_FORMATS:
		raise errors.UserInputError("Format must be either 'png', 'jpg', 'jpeg' or 'pdf'")
	if output_path is None:
		raise errors.UserInputError("output_path must be provided")
	if not output_path.endswith(f".{format}"):
		output_path = output_path + f".{format}"
	with open(output_path, "w") as f:
		None
	return output_path


def validate_html_output_path(output_html_path):
	if output_html_path is None:
		return None
	if not output_html_path.endswith(".html"):
		output_html_path = output_html_path + ".html"
	with open(output_html_path, "w") as f:
		None
	return output_html_path


def resolve_plotting_unit(plotting_unit, resp_calc, squared=False):
	if plotting_unit not in ["relative", "pcm"]:
		raise errors.UserInputError("Unit must be either 'relative' or 'pcm'")

	if plotting_unit == "relative":
		factor = 1.0
		unit_str = "%[resp]"
	else:
		factor = resp_calc * 1e5
		unit_str = "pcm"

	if squared:
		factor = factor**2

	return factor, unit_str


# ---------------------------------------------------------------------------
# Domain-data preparation
# ---------------------------------------------------------------------------

def normalize_filter_iso_reac_list(
	reference_iso_reac_list,
	iso_reac_list=None,
	iso_list=None,
	reac_list=None,
	require_target=False,
	default_to_all=False,
):
	normalized_input_iso_reac_list = None

	if iso_reac_list in [None, []]:
		if default_to_all or (require_target and iso_list not in [None, []]):
			final_iso_reac_list = list(reference_iso_reac_list)
		elif require_target:
			raise errors.UserInputError("You should target for specific isotopes with 'iso_reac_list' or 'iso_list'")
		else:
			final_iso_reac_list = []
	else:
		normalized_input_iso_reac_list = [methods.normalize_iso_reac((iso, reac)) for iso, reac in iso_reac_list]
		final_iso_reac_list = [(iso, reac) for (iso, reac) in normalized_input_iso_reac_list if (iso, reac) in reference_iso_reac_list]

	if iso_list not in [None, []]:
		iso_list = [methods.convert_iso_string_to_id(iso) if isinstance(iso, str) else int(iso) for iso in iso_list]
		final_iso_reac_list = [(iso, reac) for (iso, reac) in final_iso_reac_list if iso in iso_list]

	if reac_list not in [None, []]:
		reac_list = [int(reac) if isinstance(reac, int) else int(methods.reac_trad_inv[str(reac).upper()]) for reac in reac_list]
		final_iso_reac_list = [(iso, reac) for (iso, reac) in final_iso_reac_list if reac in reac_list]

	return final_iso_reac_list, normalized_input_iso_reac_list


def align_trace_styles(final_iso_reac_list, requested_iso_reac_list=None, traces_colors=None, dashed_traces=None):
	if traces_colors not in [None, []]:
		ref_list = requested_iso_reac_list if requested_iso_reac_list not in [None, []] else final_iso_reac_list
		if len(traces_colors) != len(ref_list):
			raise errors.UserInputError("The length of traces_colors must be equal to the length of iso_reac_list")
		if requested_iso_reac_list not in [None, []]:
			traces_colors = [c for i, c in enumerate(traces_colors) if ref_list[i] in final_iso_reac_list]

	if dashed_traces not in [None, []]:
		ref_list = requested_iso_reac_list if requested_iso_reac_list not in [None, []] else final_iso_reac_list
		if len(dashed_traces) != len(ref_list):
			raise errors.UserInputError("The length of dashed_traces must be equal to the length of iso_reac_list")
		if requested_iso_reac_list not in [None, []]:
			dashed_traces = [d for i, d in enumerate(dashed_traces) if ref_list[i] in final_iso_reac_list]

	return traces_colors, dashed_traces


def select_top_iso_reac_pairs(decomposition, candidate_iso_reac_list, contrib_col, nb_top_iso_reac):
	if nb_top_iso_reac is None:
		return list(candidate_iso_reac_list)

	selected_pairs = {(int(iso), int(reac)) for iso, reac in candidate_iso_reac_list}
	sel_df = decomposition[
		decomposition[["ISO", "REAC"]].apply(lambda row: (int(row["ISO"]), int(row["REAC"])) in selected_pairs, axis=1)
	].copy()
	sel_df["contrib_abs"] = sel_df[contrib_col].apply(lambda x: abs(float(x)))
	sel_df = sel_df.sort_values("contrib_abs", ascending=False)
	top_df = sel_df.head(nb_top_iso_reac)

	return [(int(iso), int(reac)) for iso, reac in zip(top_df["ISO"], top_df["REAC"])]


def select_top_iso_reac_pairs_by_vector(vector, candidate_iso_reac_list, full_iso_reac_list, group_nb, nb_top_iso_reac):
	if nb_top_iso_reac is None:
		return list(candidate_iso_reac_list)

	vector = np.asarray(vector)
	contributions = []
	for iso_reac in candidate_iso_reac_list:
		idx = full_iso_reac_list.index(iso_reac)
		integral = np.sum(vector[idx * group_nb : (idx + 1) * group_nb])
		contributions.append((abs(integral), iso_reac))

	contributions.sort(key=lambda item: item[0], reverse=True)
	return [iso_reac for _, iso_reac in contributions[:nb_top_iso_reac]]


def build_vector_from_decomp(decomp_vec, full_iso_reac_list, selected_iso_reac_list, group_nb):
	vector = []
	for iso, reac in selected_iso_reac_list:
		idx = full_iso_reac_list.index((iso, reac))
		vector.extend(decomp_vec[idx * group_nb : (idx + 1) * group_nb])
	return vector


def build_lethargy_vector(reference_iso_reac_list, group_nb, getter, fallback_to_zeros=False):
	vector = []
	for iso, reac in reference_iso_reac_list:
		try:
			vector += getter(iso, reac)
		except errors.MissingDataError:
			if not fallback_to_zeros:
				raise
			vector += [0 for i in range(group_nb)]
	return np.array(vector)


# ---------------------------------------------------------------------------
# Interactive HTML report figure builders
# ---------------------------------------------------------------------------

def build_interactive_integral_profile_figures(
	vector,
	iso_reac_list,
	group_nb,
	e_bins,
	factor,
	integral_title,
	profile_title,
	integral_yaxis_title,
	profile_yaxis_title,
):
	integral_fig = plots.plot_integrals_per_iso_reac(
		vector=vector,
		iso_reac_list=iso_reac_list,
		group_nb=group_nb,
		factor=factor,
		title=integral_title,
		yaxis_title=integral_yaxis_title,
	)
	plots.apply_interactive_report_layout(integral_fig)

	profile_fig = plots.plot_profiles_per_iso_reac(
		vector=vector,
		iso_reac_list=iso_reac_list,
		e_bins=e_bins,
		factor=factor,
		title=profile_title,
		yaxis_title=profile_yaxis_title,
	)
	plots.apply_interactive_report_layout(profile_fig)

	return integral_fig, profile_fig


def build_interactive_sensitivity_figures(vector, lethargy_vector, iso_reac_list, group_nb, e_bins, factor, unit_str, title_prefix=None):
	if title_prefix is None:
		integral_title = "Sensitivities (group-wise integrals)"
		profile_title = "Sensitivities"
	else:
		integral_title = f"{title_prefix} sensitivities (integrals)"
		profile_title = f"{title_prefix} sensitivities"

	integral_fig, profile_fig = build_interactive_integral_profile_figures(
		vector=vector,
		iso_reac_list=iso_reac_list,
		group_nb=group_nb,
		e_bins=e_bins,
		factor=factor,
		integral_title=integral_title,
		profile_title=profile_title,
		integral_yaxis_title=f"Sensitivity ({unit_str}/%[ND] - group-wise integral)",
		profile_yaxis_title=f"Sensitivity ({unit_str}/%[ND])",
	)

	lethargy_fig = plots.plot_profiles_per_iso_reac(
		vector=lethargy_vector,
		iso_reac_list=iso_reac_list,
		e_bins=e_bins,
		factor=factor,
		title="Sensitivities per unit lethargy",
		yaxis_title=f"Sensitivity ({unit_str}/%[ND] per unit lethargy)",
	)
	plots.apply_interactive_report_layout(lethargy_fig)

	return integral_fig, profile_fig, lethargy_fig


def build_interactive_covariance_figures(cov_mat, iso_reac_list, group_nb, e_bins, integral_title, profile_title, covariance_profile_title, integral_yaxis_title, profile_yaxis_title, covariance_profile_yaxis_title):
	integral_fig = plots.plot_matrix_integrals_per_iso_reac(
		cov_mat=cov_mat,
		iso_reac_list=iso_reac_list,
		group_nb=group_nb,
		title=integral_title,
		yaxis_title=integral_yaxis_title,
	)
	plots.apply_interactive_report_layout(integral_fig)

	diagonal_vector = np.diag(cov_mat.toarray()) if hasattr(cov_mat, "toarray") else np.diag(cov_mat)
	profile_fig = plots.plot_profiles_per_iso_reac(
		vector=diagonal_vector,
		iso_reac_list=iso_reac_list,
		e_bins=e_bins,
		title=profile_title,
		yaxis_title=profile_yaxis_title,
	)
	plots.apply_interactive_report_layout(profile_fig)

	covariance_vector = np.asarray(cov_mat.sum(axis=1)).ravel()
	covariance_profile_fig = plots.plot_profiles_per_iso_reac(
		vector=covariance_vector,
		iso_reac_list=iso_reac_list,
		e_bins=e_bins,
		title=covariance_profile_title,
		yaxis_title=covariance_profile_yaxis_title,
	)
	plots.apply_interactive_report_layout(covariance_profile_fig)

	return integral_fig, profile_fig, covariance_profile_fig


# ---------------------------------------------------------------------------
# Static publication-image export
# ---------------------------------------------------------------------------

def export_static_plotly_figure(fig, output_path, format, width, height, show):
	fig.update_layout(
		template="plotly_white",
		paper_bgcolor="rgb(255, 255, 255)",
		font=dict(family="Times New Roman", size=14),
		margin=dict(l=80, r=30, t=70, b=70),
		width=width,
		height=height,
		xaxis=dict(showgrid=True, gridcolor="lightgrey", minor=dict(showgrid=True, gridcolor="#f0f0f0")),
		yaxis=dict(showgrid=True, gridcolor="lightgrey"),
		legend=dict(bordercolor="black", borderwidth=1, font=dict(size=11)),
	)
	fig.write_image(output_path, scale=2, format=format)
	if show:
		fig.show()
	return fig
