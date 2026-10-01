import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os, re
import matplotlib.pyplot as plt


from . import errors, methods

global HTML_intro
global HTML_tab
global HTML_end

color_map_reac = {
    "ELASTIC": "darkcyan",
    "INELASTIC": "indianred",
    "NXN": "forestgreen",
    "FISSION": "darksalmon",
    "CAPTURE": "darkgoldenrod",
    "N,GAMMA": "chartreuse",
    "NUBAR": "saddlebrown",
}


def plot_integrals_per_iso_reac(vector, iso_reac_list, group_nb, factor=1.0, output_html_path: str = None, show=False, title="", yaxis_title=None, xaxis_title=None):

    if len(vector) != group_nb * len(iso_reac_list):
        raise errors.DimError(
            f"The dimmensions of the sensi vector, the isotop/reaction list and the number of groups is not consistent to plot\n {len(vector)} != {group_nb} x {len(iso_reac_list)}"
        )

    dikt = {"isotope": [], "reaction": [], "value": []}

    for i, (iso, reac) in enumerate(iso_reac_list):

        current_iso = methods.convert_iso_id_to_string(iso)
        current_reac = methods.reac_trad[str(reac)]
        sensi = sum(vector[i * group_nb : (i + 1) * group_nb]) * factor

        dikt["isotope"].append(current_iso)
        dikt["reaction"].append(current_reac)
        dikt["value"].append(sensi)

    df = pd.DataFrame(dikt)
    df = df.sort_values(by="value", ascending=False, key=lambda col: abs(col))

    fig = px.bar(df, x="isotope", y="value", hover_data=["reaction"], color="reaction", color_discrete_map=color_map_reac, title=title)
    fig.update_xaxes(categoryorder="sum descending", title=xaxis_title if xaxis_title else "Isotope")
    fig.update_yaxes(title=yaxis_title if yaxis_title else title)

    if show:
        fig.show()

    if output_html_path != None:
        if not output_html_path.endswith(".html"):
            output_html_path = output_html_path + ".html"

        fig.write_html(output_html_path)

    return fig


def plot_profiles_per_iso_reac(
    vector,
    iso_reac_list,
    e_bins: list,
    factor=1.0,
    output_html_path: str = None,
    show=False,
    title="",
    yaxis_title=None,
    xaxis_title=None,
    traces_colors=None,
    dashed_traces=None,
):

    group_nb = len(e_bins) - 1

    if len(vector) != group_nb * len(iso_reac_list):
        raise errors.DimError(
            f"The dimmensions of the sensi vector, the isotope/reaction list and the number of groups is not consistent to plot\n {len(vector)} != {group_nb} x {len(iso_reac_list)}"
        )

    if traces_colors not in [None, []] and len(traces_colors) != len(iso_reac_list):
        raise errors.UserInputError("The length of traces_colors must be equal to the length of iso_reac_list")

    if dashed_traces not in [None, []] and len(dashed_traces) != len(iso_reac_list):
        raise errors.UserInputError("The length of dashed_traces must be equal to the length of iso_reac_list")

    e_bins = sorted(e_bins, reverse=True)

    fig = go.Figure()

    vector = [val * factor for val in vector]
    vector_iso_reac = [vector[i * group_nb : (i + 1) * group_nb] for i in range(len(iso_reac_list))]
    dikt_sensi = {i_r: vec for i_r, vec in zip(iso_reac_list, vector_iso_reac)}
    dikt_sensi = dict(
        sorted(
            dikt_sensi.items(),
            key=lambda item: methods.convert_iso_id_to_string(item[0][0]) + methods.reac_trad.get(str(item[0][1]), f"REAC_{item[0][1]}"),
        )
    )

    style_by_iso_reac = {}
    for idx, iso_reac in enumerate(iso_reac_list):
        color = traces_colors[idx] if traces_colors not in [None, []] else None
        dashed = bool(dashed_traces[idx]) if dashed_traces not in [None, []] else False
        style_by_iso_reac[iso_reac] = {"color": color, "dash": "3px,2px" if dashed else "solid"}

    for (iso, reac), vec in dikt_sensi.items():

        sensis = list(vec)
        sensis.append(0)

        iso_str = methods.convert_iso_id_to_string(iso)
        reac_str = methods.reac_trad.get(str(reac), f"REAC_{reac}")
        line_style = {"dash": style_by_iso_reac[(iso, reac)]["dash"]}
        if style_by_iso_reac[(iso, reac)]["color"] is not None:
            line_style["color"] = style_by_iso_reac[(iso, reac)]["color"]

        trace = go.Scatter(x=e_bins, y=sensis, name=f"{iso_str} - {reac_str}", line_shape="hv", line=line_style, showlegend=True)
        fig.add_trace(trace)

    fig.update_layout(title=title)
    fig.update_xaxes(title=xaxis_title if xaxis_title else "Energy (eV)", type="log", tickformat=".0e")
    fig.update_yaxes(title=yaxis_title if yaxis_title else title)


    if show:
        fig.show()

    if output_html_path != None:
        if not output_html_path.endswith(".html"):
            output_html_path = output_html_path + ".html"

        fig.write_html(output_html_path)

    return fig


def plot_matrix_integrals_per_iso_reac(cov_mat, iso_reac_list, group_nb, title, output_html_path: str = None, show=False, yaxis_title=None):

    if np.shape(cov_mat) != (group_nb * len(iso_reac_list), group_nb * len(iso_reac_list)):
        raise errors.DimError("The dimmensions of the covar matrix, the isotope/reaction list and the number of groups is not consistent to plot")

    # Sum along axis 0, reshape by group_nb, then sum along axis 1
    if hasattr(cov_mat, "sum"):
        var_ints = np.asarray(cov_mat.sum(axis=0)).flatten()
    else:
        var_ints = np.sum(cov_mat, axis=0)
    var_ints = np.sum(var_ints.reshape(len(iso_reac_list), group_nb), axis=1)

    iso_strings = [methods.convert_iso_id_to_string(iso) for iso, _ in iso_reac_list]
    reac_strings = [methods.reac_trad.get(str(reac), f"REAC_{reac}") for _, reac in iso_reac_list]

    dikt = {"isotope": iso_strings, "reaction": reac_strings, "var_and_covar_integral": var_ints.tolist()}

    df = pd.DataFrame(dikt)

    title = title

    fig = px.bar(
        df, x="isotope", y="var_and_covar_integral", hover_data=["reaction"], title=title, color="reaction", color_discrete_map=color_map_reac
    )
    fig.update_xaxes(categoryorder="sum descending", title="Isotope")
    if yaxis_title == None:
        fig.update_yaxes(title=title)
    else:
        fig.update_yaxes(title=yaxis_title)

    if show:
        fig.show()

    if output_html_path != None:
        if not output_html_path.endswith(".html"):
            output_html_path = output_html_path + ".html"

        fig.write_html(output_html_path)

    return fig


def plot_submatrix(
    cov_mat,
    iso_reac_list,
    group_nb,
    iso_reac_pair_to_plot: tuple,
    title,
    color_scale: list = None,
    scale_min=None,
    output_html_path: str = None,
    show=False,
):

    ((iso1, reac1), (iso2, reac2)) = iso_reac_pair_to_plot

    try:
        iso1 = methods.convert_iso_string_to_id(iso1)
        iso2 = methods.convert_iso_string_to_id(iso2)
    except:
        iso1 = int(iso1)
        iso2 = int(iso2)

    try:
        reac1 = int(methods.reac_trad_inv[str(reac1).upper()])
        reac2 = int(methods.reac_trad_inv[str(reac2).upper()])
    except:
        reac1 = int(reac1)
        reac2 = int(reac2)

    iso_reac_idx_h = iso_reac_list.index((iso1, reac1))
    iso_reac_idx_v = iso_reac_list.index((iso2, reac2))

    if iso_reac_idx_h > len(iso_reac_list) + 1 or iso_reac_idx_v > len(iso_reac_list) + 1:
        raise errors.DimError("The iso-reac index you're targetting is too high for the iso-reac list")

    sub_mat = cov_mat[iso_reac_idx_v * group_nb : (iso_reac_idx_v + 1) * group_nb, iso_reac_idx_h * group_nb : (iso_reac_idx_h + 1) * group_nb]
    sub_mat = sub_mat.toarray()

    iso_reac_str_H = (
        methods.convert_iso_id_to_string(iso_reac_list[iso_reac_idx_h][0]) + " " + methods.reac_trad[str(iso_reac_list[iso_reac_idx_h][1])]
    )
    iso_reac_str_V = (
        methods.convert_iso_id_to_string(iso_reac_list[iso_reac_idx_v][0]) + " " + methods.reac_trad[str(iso_reac_list[iso_reac_idx_v][1])]
    )

    if color_scale != None:
        zmin, zmax = np.amin(color_scale), np.amax(color_scale)
    else:
        zmin, zmax = np.amin(sub_mat), np.amax(sub_mat)

    fig = px.imshow(
        sub_mat,
        title=f"{title} <span style='color:grey'>{iso_reac_str_H} / {iso_reac_str_V}     </span> \n <span style='color:black'>{group_nb} Energy groups</span>",
        zmax=zmax,
        zmin=zmin,
    )
    fig.update_xaxes(
        {
            "title": {"text": iso_reac_str_H, "standoff": 2},
            "side": "top",
            "tickmode": "array",
            "tickvals": [0, group_nb - 1],
            "ticktext": ["20 MeV", "0 MeV"],
        }
    )
    fig.update_yaxes(
        {
            "title": {"text": iso_reac_str_V, "standoff": 2},
            "side": "top",
            "tickmode": "array",
            "tickvals": [0, group_nb - 1],
            "ticktext": ["20 MeV", "0 MeV"],
        }
    )

    if show:
        fig.show()

    if output_html_path != None:
        if not output_html_path.endswith(".html"):
            output_html_path = output_html_path + ".html"

        fig.write_html(output_html_path)

    if np.sum(sub_mat) == 0:
        return None, None
    else:
        return fig, (zmin, zmax)


def plot_specific_covariance_matrix(
    covariance_block,
    e_bins,
    title="",
    covariance_horizontal=None,
    covariance_vertical=None,
    horizontal_label=None,
    vertical_label=None,
):

    group_nb = covariance_block.shape[0]
    is_cross_covariance = covariance_vertical is not None

    std_horizontal = np.sqrt(np.diag(covariance_horizontal if covariance_horizontal is not None else covariance_block)) * 100
    std_vertical = np.sqrt(np.diag(covariance_vertical if covariance_vertical is not None else covariance_block)) * 100

    correlation_matrix = np.zeros_like(covariance_block)
    for row_idx in range(group_nb):
        for col_idx in range(group_nb):
            if std_vertical[row_idx] > 0 and std_horizontal[col_idx] > 0:
                correlation_matrix[row_idx, col_idx] = covariance_block[row_idx, col_idx] / (
                    std_vertical[row_idx] / 100 * std_horizontal[col_idx] / 100
                )
    
    if len(e_bins) == group_nb : e_bins.append(1E-6)
    if len(e_bins) != group_nb + 1:
        raise errors.DimError(f"The length of the energy bins list must be equal to the number of groups + 1. {len(e_bins)} != {group_nb} + 1")
    e_bins_sorted = np.sort(e_bins)
    log_edges = np.log10(e_bins_sorted)


    std_horizontal_ascending = std_horizontal[::-1]
    std_vertical_ascending = std_vertical[::-1]
    correlation_matrix_ascending = correlation_matrix[::-1, ::-1]

    if is_cross_covariance:
        fig_width, fig_height = 10.0, 9.0
        matrix_x0, matrix_y0, matrix_h = 0.09, 0.08, 0.62
        matrix_w = matrix_h * fig_height / fig_width
        top_y0, top_h = matrix_y0 + matrix_h + 0.05, 0.13
    else:
        fig_width, fig_height = 8.0, 9.0
        matrix_x0, matrix_y0, matrix_h = 0.10, 0.08, 0.65
        matrix_w = matrix_h * fig_height / fig_width
        top_y0, top_h = matrix_y0 + matrix_h + 0.04, 0.13

    fig = plt.figure(figsize=(fig_width, fig_height))

    ax_top = fig.add_axes([matrix_x0, top_y0, matrix_w, top_h])
    ax_top.step(log_edges[:-1], std_horizontal_ascending, where="post", color="blue", linewidth=1.5)
    ax_top.hlines(std_horizontal_ascending[-1], log_edges[-2], log_edges[-1], color="blue", linewidth=1.5)
    ax_top.set_xlim(log_edges[0], log_edges[-1])
    ax_top.set_ylabel("Δσ/σ (%)")
    if is_cross_covariance and horizontal_label is not None:
        ax_top.set_title(horizontal_label, fontsize=12)
    else:
        ax_top.set_title(title, fontsize=14)
    ax_top.tick_params(axis="x", labelbottom=False)
    ax_top.grid(True, alpha=0.3)

    ax_matrix = fig.add_axes([matrix_x0, matrix_y0, matrix_w, matrix_h], sharex=ax_top)
    x_mesh, y_mesh = np.meshgrid(log_edges, log_edges)
    pcm = ax_matrix.pcolormesh(x_mesh, y_mesh, correlation_matrix_ascending, cmap="RdYlGn", vmin=-1, vmax=1, shading="flat")
    ax_matrix.set_xlim(log_edges[0], log_edges[-1])
    ax_matrix.set_ylim(log_edges[0], log_edges[-1])
    ax_matrix.set_aspect("equal")

    if is_cross_covariance and horizontal_label is not None:
        ax_matrix.set_xlabel(f"Energy (eV) — {horizontal_label}")
    else:
        ax_matrix.set_xlabel("Energy (eV)")
    if is_cross_covariance and vertical_label is not None:
        ax_matrix.set_ylabel(f"Energy (eV) — {vertical_label}")
    else:
        ax_matrix.set_ylabel("Energy (eV)")

    tick_values = np.arange(np.ceil(log_edges[0]), np.floor(log_edges[-1]) + 1)
    tick_labels = [f"$10^{{{int(val)}}}$" for val in tick_values]
    ax_top.set_xticks(tick_values)
    ax_matrix.set_xticks(tick_values)
    ax_matrix.set_xticklabels(tick_labels)
    ax_matrix.set_yticks(tick_values)
    ax_matrix.set_yticklabels(tick_labels)

    if is_cross_covariance:
        right_x0 = matrix_x0 + matrix_w + 0.015
        right_w = 0.09
        ax_right = fig.add_axes([right_x0, matrix_y0, right_w, matrix_h], sharey=ax_matrix)
        x_points, y_points = [], []
        for group_idx in range(group_nb):
            y_points += [log_edges[group_idx], log_edges[group_idx + 1]]
            x_points += [std_vertical_ascending[group_idx], std_vertical_ascending[group_idx]]
        ax_right.plot(x_points, y_points, color="blue", linewidth=1.5)
        ax_right.set_ylim(log_edges[0], log_edges[-1])
        ax_right.set_xlabel("Δσ/σ (%)")
        ax_right.tick_params(axis="y", labelleft=False)
        ax_right.grid(True, alpha=0.3)
        if vertical_label is not None:
            ax_right.set_title(vertical_label, fontsize=10)

        cax = fig.add_axes([right_x0 + right_w + 0.05, matrix_y0, 0.025, matrix_h])
    else:
        cax = fig.add_axes([matrix_x0 + matrix_w + 0.02, matrix_y0, 0.025, matrix_h])

    plt.colorbar(pcm, cax=cax, label="Correlation")
    return fig


def html_setup():

    with open(os.path.join(os.path.dirname(__file__), "html_outputfile", "html.model"), "r") as f:
        lines = f.readlines()
    global HTML_intro
    global HTML_tab
    global HTML_end
    HTML_intro, HTML_tab, HTML_end = [], [], []
    idx = 0
    for line in lines:
        if re.search("@", line):
            idx += 1
            continue
        if idx == 1:
            HTML_intro.append(line)
        elif idx == 2:
            HTML_tab.append(line)
        elif idx == 3:
            HTML_end.append(line)


def create_html_tabs(names: list = []):

    txt = HTML_tab[0]
    for name in names:
        txt += HTML_tab[1].replace("$$Name", name)
    txt += HTML_tab[2]

    return txt


def create_html_tip(txt):
    """Return a small '?' icon with a hover tooltip."""
    safe = txt.replace("'", "&#39;").replace('"', "&quot;")
    return f'<span class="tip">?<span class="tiptext">{safe}</span></span>'


def apply_interactive_report_layout(fig, **overrides):
    """Apply the standard layout for a Plotly figure embedded in a CALINS HTML report."""
    layout = {
        "height": 500,
        "width": 1200,
        "font_size": 14,
        "font_family":"Times New Roman",
        "template": "plotly_white",
        "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
        "xaxis_showgrid":True,
        "xaxis_gridcolor":"lightgrey",
        "xaxis_minor":dict(showgrid=True, gridcolor='#f0f0f0'),
        "yaxis_showgrid":True,
        "yaxis_gridcolor":"lightgrey",
        "legend":dict(
            bordercolor='black',
            borderwidth=1,
            font=dict(size=11))
    }
    layout.update(overrides)
    fig.update_layout(layout)
    return fig

def create_html_table(headers=None, lines=None, color_per_lines=None, table_attrs=None, header_style=None, centered=True):
    if color_per_lines is None:
        color_per_lines = ["white" for i in range(len(lines[0]))]

    table = ""
    if centered:
        table += "<h1> </h1>"
        table += '<div style="display: flex; align-items: center; justify-content: center;">\n'

    if table_attrs is None:
        table += '<table style="font-size:14px;" border="0" bordercolor="#363636" bgcolor="#e9d4c9">\n'
    else:
        table += f"<table {table_attrs}>\n"

    # Create the table's column headers
    if headers:
        if header_style:
            table += f"  <tr style='{header_style}'>\n"
        else:
            table += "  <tr>\n"
        for column in headers:
            table += f"    <th>{column}</th>\n"
        table += "  </tr>\n"

    # Create the table's row data
    for r in range(len(lines[0])):
        table += "  <tr>\n"
        for c in range(len(lines)):
            table += f'    <td bgcolor="{color_per_lines[r]}">{lines[c][r]}</td>\n'
        table += "  </tr>\n"

    table += "</table>"
    if centered:
        table += "</div>"
        table += "<h1> </h1>"

    return table


html_setup()
