import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os, re


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


def plot_integrals_per_iso_reac(vector, iso_reac_list, group_nb, factor=1.0, output_html_path: str = None, show=False, title="", yaxis_title=None):

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


def plot_profiles_per_iso_reac(vector, iso_reac_list, e_bins: list, factor=1.0, output_html_path: str = None, show=False, title="", yaxis_title=None):

    group_nb = len(e_bins) - 1

    if len(vector) != group_nb * len(iso_reac_list):
        raise errors.DimError(
            f"The dimmensions of the sensi vector, the isotope/reaction list and the number of groups is not consistent to plot\n {len(vector)} != {group_nb} x {len(iso_reac_list)}"
        )

    e_bins.sort(reverse=True)

    fig = go.Figure()

    iso_reac_list = [(methods.convert_iso_id_to_string(iso), methods.reac_trad[str(reac)]) for iso, reac in iso_reac_list]

    vector = [val * factor for val in vector]
    vector_iso_reac = [vector[i * group_nb : (i + 1) * group_nb] for i in range(len(iso_reac_list))]
    dikt_sensi = {i_r: vec for i_r, vec in zip(iso_reac_list, vector_iso_reac)}
    dikt_sensi = dict(sorted(dikt_sensi.items(), key=lambda item: item[0][0] + item[0][1]))

    for iso_reac, vec in dikt_sensi.items():

        sensis = list(vec)
        sensis.append(0)

        trace = go.Scatter(x=e_bins, y=sensis, name=f"{iso_reac[0]} - {iso_reac[1]}", line_shape="hv")
        fig.add_trace(trace)

    fig.update_layout(title=title)
    fig.update_xaxes(title="Energy (eV)", type="log", tickformat="e")
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


def plot_matrix_integrals_per_iso_reac(cov_mat, iso_reac_list, group_nb, title, output_html_path: str = None, show=False, yaxis_title=None):

    if np.shape(cov_mat) != (group_nb * len(iso_reac_list), group_nb * len(iso_reac_list)):
        raise errors.DimError("The dimmensions of the covar matrix, the isotope/reaction list and the number of groups is not consistent to plot")

    dikt = {"isotope": [], "reaction": [], "var_and_covar_integral": []}
    for i, (iso, reac) in enumerate(iso_reac_list):

        current_iso = methods.convert_iso_id_to_string(iso)
        current_reac = methods.reac_trad[str(reac)]
        var_int = np.sum(cov_mat[:, i * group_nb : (i + 1) * group_nb])

        dikt["isotope"].append(current_iso)
        dikt["reaction"].append(current_reac)
        dikt["var_and_covar_integral"].append(var_int)

    df = pd.DataFrame(dikt)

    # title = title + " - Var-covar profile (integral on every var-covar involved for each Isotope/Reaction)"
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

    # Penser à gerer un dataframe en input aussi

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


def create_html_table(headers=None, lines=None, color_per_lines=None):
    if color_per_lines == None:
        color_per_lines = ["white" for i in range(len(lines[0]))]

    table = "<h1> </h1>"

    table += '<div style="display: flex; align-items: center; justify-content: center;">\n'
    table += '<table style="font-size:14px;" border="0" bordercolor="#363636" bgcolor="#e9d4c9">\n'

    # Create the table's column headers
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
    table += "</div>"
    table += "<h1> </h1>"

    return table


html_setup()
