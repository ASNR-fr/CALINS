import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from tabulate import tabulate
import os, copy, re, warnings
import plotly.offline as po
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
from math import ceil, sqrt
from importlib.metadata import version

from . import methods, errors, plots
from logs import log_exec, warn, write_and_print

global HTML_intro
global HTML_tab
global HTML_end

HTML_intro, HTML_tab, HTML_end = plots.HTML_intro.copy(), plots.HTML_tab.copy(), plots.HTML_end.copy()

global plotlyjs_fig_include
global OFFLINE_EXPORT
global OFFLINE_CONTROL
OFFLINE_EXPORT = False
OFFLINE_CONTROL = False
plotlyjs_fig_include = "cdn"

global plot_height
global plot_width
global block_margin
global font_size
global plotly_template

plot_height, plot_width, block_margin, font_size, plotly_template = 500, 1200, 20, 10, "plotly_white"


__pkg_name__ = methods.__package__.split(".")[0]
try:
    __version__ = version(__pkg_name__)
except:
    __version__ = f"no-version-found"

global pkg_fullname
pkg_fullname = f"{__pkg_name__} v. {__version__}"


def offline_control():
    global OFFLINE_CONTROL
    global plotlyjs_fig_include
    global HTML_intro
    if OFFLINE_EXPORT and not OFFLINE_CONTROL:
        plotlyjs_fig_include = False
        plotly_js = "<script type='text/javascript'>\n"
        plotly_js += po.get_plotlyjs() + "\n"
        plotly_js += "</script>\n"
        for i, line in enumerate(HTML_intro):
            if re.search("<head>", line):
                HTML_intro[i] = line + plotly_js
                OFFLINE_CONTROL = True
    else:
        HTML_intro = plots.HTML_intro
        plotlyjs_fig_include = "cdn"

    return None


class Case:
    """
    A class to represent a case for sensitivity analysis.

    Attributes
    -----------
    sdf_path : str
        Path to the SDF file.
    casename : str
        Name of the case.
    resp_calc : float
        Calculated response. Read from SDF file if not provided (if possible).
    sigma_resp_calc : float
        Uncertainty in the calculated response. Read from SDF file if not provided (if possible).
    resp_expe : float
        Experimental response. Read from SDF file if not provided (if possible).
    sigma_resp_expe : float
        Uncertainty in the experimental response. Read from SDF file if not provided (if possible).
    sensitivities : DataFrame
        Sensitivity data (containing standard deviations if the argument std is activated). Output of the function format_sensi_to_dataframe().
    group_nb : int
        Number of energy groups.
    e_bins : list
        Energy bins.
    iso_reac_list : list
        List of isotope-reaction pairs.
    iso_list : list
        List of isotopes present in the case.
    reac_list : list
        List of reactions present in the case.

    Methods
    --------
    export_to_html(output_html_path: str, plotting_unit="pcm", show=False):
        Exports the sensitivity data for the case to an HTML file with interactive plots.
    condense_sensi(output_ebins: list):
        Condenses the sensitivity coefficients and std to a specified energy binning (if compatible).
    create_sdf(output_sdf_path="", header=""):
        Creates an SDF file from the parsed, and/or condensed sensitivity coefficients.
    get_sensitivity_values(isotope, reaction):
        Retrieves sensitivity values for a given isotope and reaction.
    get_integral_sensitivity_value(isotope, reaction, absolute_value=False):
        Retrieves the integral sensitivity value for a given isotope and reaction.
    get_normalized_lethargy_sensitivity(isotope, reaction):
        Returns the normalized lethargy sensitivity for a given isotope and reaction.
    filter_energy_bins(energy_region):
        Filters energy bins based on the specified energy region.
    get_filtered_sensitivity_by_energy_region(isotope, reaction, energy_region, normalize_lethargy):
        Filters sensitivity values for a given isotope and reaction based on the specified energy region.
    """

    @log_exec()
    def __init__(
        self,
        sdf_path,
        resp_calc=None,
        sigma_resp_calc=None,
        resp_expe=None,
        sigma_resp_expe=None,
        filter_min_ratio_sensi_abs=0.0,
        occurrences_rule="sum",
        mcnp=False,
        std=False,
    ) -> None:
        """
        Initializes the Case object.

        Parameters
        -----------
        sdf_path : str
            [Required] Path to the SDF file.
        resp_calc : float, optional
            Calculated response. Read from SDF file if not provided (if possible).
        sigma_resp_calc : float, optional
            Uncertainty in the calculated response. Read from SDF file if not provided (if possible).
        resp_expe : float, optional
            Experimental response. Read from SDF file if not provided (if possible).
        sigma_resp_expe : float, optional
            Uncertainty in the experimental response. Read from SDF file if not provided (if possible).
        filter_min_ratio_sensi_abs : float, optional
            Minimum ratio of the absolute sensitivity integral (for one isotope-reaction pair) to the sum of all absolute sensitivities. The default is 0.0.
        occurences_rule : str, optional
            Rule for handling occurrences of same isotope-reaction in one SDf file. Either "first", "last" or "sum", to chosse which one to consider. Defaults to "sum".
        mcnp : bool, optional
            Flag to be activated when parsing MCNP format of SDF file.
        std : bool, optional
            Flag to be activated for parsing standard deviations of sensitivity coefficients.
        """
        self.sdf_path = sdf_path
        self.casename = os.path.basename(sdf_path)
        self.resp_calc = resp_calc
        self.sigma_resp_calc = sigma_resp_calc
        self.resp_expe = resp_expe
        self.sigma_resp_expe = sigma_resp_expe

        self.sensitivities = methods.format_sensi_to_dataframe(
            input_sdf_path=sdf_path, filter_min_ratio_sensi_abs=filter_min_ratio_sensi_abs, occurrences_rule=occurrences_rule, mcnp=mcnp, std=std
        )

        self.read_group_number()
        if self.group_nb != len(self.sensitivities["SENSI"][0]):
            raise errors.DimError(
                f"The energy group number in the header of the sdf files is not equal to the first iso-reac sensitivity vector length - Check your sdf file : {self.casename}"
            )

        self.read_energy_bins()

        self.iso_reac_list = [(iso, reac) for (iso, reac) in zip(self.sensitivities["ISO"].to_list(), self.sensitivities["REAC"].to_list())]
        self.iso_list = list(set([iso for (iso, reac) in self.iso_reac_list]))
        self.reac_list = list(set([reac for (iso, reac) in self.iso_reac_list]))

        if self.resp_calc == None or self.sigma_resp_calc == None:

            self.read_resp_calc()

        if self.resp_expe == None or self.sigma_resp_expe == None:

            self.read_resp_expe()

    def read_group_number(self):
        with open(self.sdf_path, "r") as f:
            lines_sensi = f.readlines()

        for i, line in enumerate(lines_sensi):
            if re.search("number of neutron groups", line):
                self.group_nb = int(line.split()[0])
                break

            if i > 50:
                raise errors.EmptyParsingError(
                    f'The energy groups number is not written before the first 50 lines of your sdf file - Check your sdf file : {self.casename} - Has to be written before the sentence "number of neutron groups"'
                )

    def read_energy_bins(self):
        with open(self.sdf_path, "r") as f:
            lines_sensi = f.readlines()

        for i, line in enumerate(lines_sensi):
            if re.search("energy boundaries", line):
                self.e_bins = [float(x) for x in "".join(lines_sensi[i + 1 : i + 1 + ceil((self.group_nb + 1) / 5)]).split()]
                break

            if i > 500:
                raise errors.EmptyParsingError(
                    f"The energy bins are not written before the first 500 lines of your sdf file - Check your sdf file : {self.casename}"
                )

    def read_resp_calc(self):
        with open(self.sdf_path, "r") as f:
            lines_sensi = f.readlines()

        try:
            if lines_sensi[3].split()[1] == "+/-":
                self.resp_calc = float(lines_sensi[3].split()[0])
                self.sigma_resp_calc = float(lines_sensi[3].split()[2])
            elif re.search("k-eff from the forward case", lines_sensi[3]):
                self.resp_calc = float(lines_sensi[3].split()[0])
                self.sigma_resp_calc = 0

            if self.resp_calc > 2 or self.resp_calc < 0 or self.sigma_resp_calc > 2 or self.sigma_resp_calc < 0:
                raise errors.EmptyParsingError(
                    f"The sdf case file {self.casename} doesn't have the right format for simulated resp values \n\
                    The format should be on line nb 4 as : [resp calc] +/- [sigma resp calc] ..."
                )
        except:
            warn(f"No calculated response was extracted from the file {self.sdf_path}.")
            None

    def read_resp_expe(self):
        with open(self.sdf_path, "r") as f:
            lines_sensi = f.readlines()

        try:
            if re.search("keff expe :", lines_sensi[0]):
                self.resp_expe = float(lines_sensi[0].split()[-3])
                self.sigma_resp_expe = float(lines_sensi[0].split()[-1])
            elif len(lines_sensi[3].split()) in [5, 10]:
                self.resp_expe = float(lines_sensi[3].split()[3])
                self.sigma_resp_expe = float(lines_sensi[3].split()[4])

            if self.resp_expe > 2 or self.resp_expe < 0 or self.sigma_resp_expe > 2 or self.sigma_resp_expe < 0:
                raise errors.EmptyParsingError(
                    f"The sdf case file {self.casename} doesn't have the right format for simulated resp values \n\
                    The format should be on line nb 1 as : [resp calc] +/- [sigma resp calc]"
                )
        except:
            warn(f"No experimental response was extracted from the file {self.sdf_path}.")
            None

    @log_exec()
    def export_to_html(self, output_html_path: str, plotting_unit="pcm", show=False):
        """
        Exports the sensitivity data for the case to an HTML file with interactive plots.

        Parameters
        -----------
        output_html_path : str
            [Required] Path to save the output HTML file.
        plotting_unit : str, optional
            Unit for the sensitivity coefficients when plotted. Can be "relative" or "pcm". Defaults to "relative".
        show : bool, optional
            Flag to display the plot (if possible).
        """
        if not output_html_path.endswith(".html"):
            output_html_path = output_html_path + ".html"
        with open(output_html_path, "w") as f:
            None
        if plotting_unit not in ["relative", "pcm"]:
            raise errors.UserInputError("Unit must be either 'relative' or 'pcm'")

        # --------------------------------
        text_intro = ""
        text_intro += f"<h1><center>{pkg_fullname}</center></h1>"
        text_intro += f"<h2>Study case : {self.casename}</h2>"

        # --------------------------------
        headers = [
            "Calculated response of study case",
            "Energy groups number",
        ]
        headers = ["<b>" + h + ":</b>" for h in headers]

        results = [
            f"{self.resp_calc} +/- {self.sigma_resp_calc}",
            self.group_nb,
        ]

        table_res = plots.create_html_table(headers=[], lines=[headers, results])

        if plotting_unit == "relative":
            val_factor = 1.0
            unit_str = "%[resp]"
        elif plotting_unit == "pcm":
            val_factor = self.resp_calc * 1e5
            unit_str = "pcm"

        [case_vec], iso_reac_list = methods.make_sensi_vectors(cases_list=[self])
        trace_sensi_integrals = plots.plot_integrals_per_iso_reac(
            vector=case_vec,
            iso_reac_list=iso_reac_list,
            group_nb=self.group_nb,
            factor=val_factor,
            show=show,
            title="Sensitivities (group-wise integrals)",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND] - group-wise integral)",
        )
        trace_sensi_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        case_vec_per_unit_lethargy = []
        for iso, reac in iso_reac_list:
            case_vec_per_unit_lethargy += self.get_normalized_lethargy_sensitivity(iso, reac)

        case_vec_per_unit_lethargy = np.array(case_vec_per_unit_lethargy)

        trace_sensi_profiles = plots.plot_profiles_per_iso_reac(
            vector=case_vec,
            iso_reac_list=iso_reac_list,
            e_bins=self.e_bins,
            factor=val_factor,
            show=show,
            title="Sensitivities",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND])",
        )
        trace_sensi_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        trace_sensi_profiles_lethargy = plots.plot_profiles_per_iso_reac(
            vector=case_vec_per_unit_lethargy,
            iso_reac_list=iso_reac_list,
            e_bins=self.e_bins,
            factor=val_factor,
            show=show,
            title="Sensitivities per unit lethargy",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND] per unit lethargy)",
        )
        trace_sensi_profiles_lethargy.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        with open(output_html_path, "a", encoding="utf-8") as f:

            offline_control()
            f.writelines(HTML_intro)
            f.write(text_intro)
            f.write(table_res)
            f.write(plots.create_html_tabs(["Study case sensitivities"]))

            f.write('<div id="Study case sensitivities" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#a8b0af, white)">\n' + "<br>\n")
            f.write(trace_sensi_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_sensi_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "\n<br>\n")
            f.write(trace_sensi_profiles_lethargy.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + " \n")
            f.write("</section>\n")
            f.write("</div>\n")

            f.writelines(HTML_end)

    @log_exec()
    def condense_sensi(self, output_ebins: list):
        """
        Condenses the sensitivity coefficients and std to a specified energy binning (if compatible).

        Parameters
        -----------
        output_ebins : list
            List of energy bins to condense the sensitivities to.

        """

        energy_bins = self.e_bins
        for o_ebin in output_ebins:
            try:
                bound = energy_bins.index(o_ebin)
                found = True
            except:
                found = False
                for i, ebin in enumerate(energy_bins):
                    if abs(o_ebin - ebin) / o_ebin < 0.0001:
                        found = True
                        warn(
                            f"The ouput energy binning for condensation have a binning that is (<)0.01% different from input binning : input {ebin}, ouput {o_ebin}"
                        )
                        break
            if not found:
                raise errors.UserInputError(
                    f"The two energy binings are not compatibles for compression\n\
                    Output ebin {o_ebin} not in input ebins {energy_bins}"
                )

        pass_line = 0

        for i, row in self.sensitivities.iterrows():
            sensis = row["SENSI"]
            self.sensitivities.at[i, "SENSI"] = methods.condense_binning(sensis, i_ebins=energy_bins, o_ebins=output_ebins)

        self.e_bins = output_ebins
        self.group_nb = len(output_ebins) - 1

    def create_sdf(self, output_sdf_path="", header=""):
        """
        Creates an SDF file from the parsed, and/or condensed sensitivity coefficients.

        Parameters
        -----------
        output_sdf_path : str, optional
            Path to save the output SDF file.
        header : str, optional
            Header for the SDF file.
        """
        title = os.path.basename(self.sdf_path) + " - " + header

        sensi_nuclide = {}

        keff, keff_std = self.resp_calc, self.sigma_resp_calc

        input_energy_bins = self.e_bins

        input_energy_bins = sorted(input_energy_bins, reverse=True)
        input_energy_bins = ["{:.6E}".format(j) for j in input_energy_bins]

        for i, row in self.sensitivities.iterrows():

            iso = row["ISO"]

            reac = row["REAC"]

            means = row["SENSI"]
            stds = [0.0 for i in range(len(means))]

            integral = np.sum(means)
            integral_std = 0

            integral_abs = np.sum(np.abs(means))
            if integral > 0:
                integral_oppose = np.sum([m for m in means if m < 0])
                sigma_oppose = np.sqrt(np.sum([stds[i] ** 2 for i in range(len(means)) if means[i] < 0]))
            elif integral < 0:
                integral_oppose = np.sum([m for m in means if m > 0])
                sigma_oppose = np.sqrt(np.sum([stds[i] ** 2 for i in range(len(means)) if means[i] > 0]))
            else:
                integral_oppose = 0
                sigma_oppose = 0

            means = ["{:.6E}".format(j) for j in means]
            stds = ["{:.6E}".format(j) for j in stds]
            integral = "{:.6E}".format(integral)
            integral_std = "{:.6E}".format(integral_std)
            integral_abs = "{:.6E}".format(integral_abs)
            integral_oppose = "{:.6E}".format(integral_oppose)
            sigma_oppose = "{:.6E}".format(sigma_oppose)

            sensi_nuclide[f"{iso} {reac}"] = {
                "means": means,
                "stds": stds,
                "integral": integral,
                "integral_std": integral_std,
                "integral_abs": integral_abs,
                "integral_oppose": integral_oppose,
                "sigma_oppose": sigma_oppose,
            }

        # CONSTRUCTION OF SDF TEXT
        # INTRO
        text = f"{title} - {self.group_nb}gr{', keff expe : ' + str(self.resp_expe) + ' +/- ' + str(self.sigma_resp_expe) if self.resp_expe is not None else ''}\n"
        text += f"       {self.group_nb:3} number of neutron groups\n"
        text += f"       {len(sensi_nuclide):3}   number of sensitivity profiles         {len(sensi_nuclide):3} are region integrated\n"
        text += f"  {keff:.6f} +/-   {keff_std:.6f} k-eff from the forward case\n"
        text += "energy boundaries:\n"

        # ENERGY BINS
        for idx, energy in enumerate(input_energy_bins):
            if idx > 0 and idx % 5 == 0:
                text += "\n"
            text += f"{energy: >14}"
        text += "\n"

        # SENSIBILITIES and STDS
        for iso_reac, sensi_nuclide in sensi_nuclide.items():
            iso_idx = iso_reac.split()[0]
            reac_idx = iso_reac.split()[1]
            iso_trad_str = methods.convert_iso_id_to_string(iso_idx).lower()
            reac_trad_str = methods.reac_trad.get(str(int(reac_idx)), f"REAC_{int(reac_idx)}").lower()
            integral_abs = sensi_nuclide["means"]

            text += f"{iso_trad_str: <13}{reac_trad_str: <17}{iso_idx: >5}{reac_idx: >7}\n"
            text += "      0      0\n"
            text += "  0.000000E+00  0.000000E+00      0      0\n"
            text += f"{sensi_nuclide['integral']: >14}{sensi_nuclide['integral_std']: >14}{sensi_nuclide['integral_abs']: >14}{sensi_nuclide['integral_oppose']: >14}{sensi_nuclide['sigma_oppose']: >14}\n"

            for idx, sensi in enumerate(sensi_nuclide["means"]):
                if idx > 0 and idx % 5 == 0:
                    text += "\n"
                text += f"{sensi: >14}"
            text += "\n"

            for idx, std in enumerate(sensi_nuclide["stds"]):
                if idx > 0 and idx % 5 == 0:
                    text += "\n"
                text += f"{std: >14}"
            text += "\n"

        if output_sdf_path == "":
            output_sdf_path = self.sdf_path + f".condensed-{self.group_nb}g.sdf"
        with open(output_sdf_path, "w+") as f:
            f.write(text)

    def get_sensitivity_values(self, isotope, reaction):
        """
        Retrieve sensitivity values for a given isotope and reaction.

        Parameters
        ----------
        isotope : str or int
            The isotope of interest, either as a string (e.g., 'U235') or as an integer ID.
        reaction : int
            The reaction type ID.

        Returns
        -------
        list
            A list of sensitivity values for the specified isotope and reaction.
        """

        if isinstance(isotope, str):
            isotope = methods.convert_iso_string_to_id(isotope)
        elif not isinstance(isotope, int):
            raise errors.UserInputError(f"isotope must be a str or an int and was provided as {type(isotope).__name__}")

        if not isinstance(reaction, int):
            raise errors.UserInputError(f"reaction must be an int, but was provided as {type(reaction).__name__}")

        try:
            df = self.sensitivities
            filtered_df = df[(df["ISO"] == isotope) & (df["REAC"] == reaction)]
            sensi_values = filtered_df["SENSI"].values[0]
        except:
            raise errors.MissingDataError(f"Couldn't find any data for isotope {isotope} and reaction {reaction}")

        return sensi_values

    def get_integral_sensitivity_value(self, isotope, reaction, absolute_value=False):
        """
        Retrieve the integral sensitivity value for a given isotope and reaction.

        Parameters
        ----------
        isotope : str or int
            The isotope of interest, either as a string (e.g., 'U235') or as an integer ID.
        reaction : int
            The reaction type ID.
        absolute_value : bool
            boolean to recover absolute integral sensitivity value

        Returns
        -------
        list
            An integral sensitivity value for the specified isotope and reaction.
        """

        if isinstance(isotope, str):
            isotope = methods.convert_iso_string_to_id(isotope)
        elif not isinstance(isotope, int):
            raise errors.UserInputError(f"isotope must be a str or an int and was provided as {type(isotope).__name__}")

        if not isinstance(reaction, int):
            raise errors.UserInputError(f"reaction must be an int, but was provided as {type(reaction).__name__}")

        try:
            df = self.sensitivities
            filtered_df = df[(df["ISO"] == isotope) & (df["REAC"] == reaction)]
            if not absolute_value:
                sensi_value = filtered_df["SENSI_INTEGRAL"].values[0]
            else:
                sensi_value = filtered_df["SENSI_INTEGRAL_ABS"].values[0]
        except:
            raise errors.MissingDataError(f"Couldn't find any data for isotope {isotope} and reaction {reaction}")
        return sensi_value

    def get_normalized_lethargy_sensitivity(self, isotope, reaction):
        """
        Return the normalized lethargy sensitivity for a given isotope and reaction.

        Parameters
        -----------
        isotope : str or int
            The isotope identifier, which can be either a string (e.g., "U235") or an integer (e.g., 92235).
            If a string is provided, it will be converted to the corresponding integer ID.

        reaction : int
            The reaction identifier as an integer.

        Returns
        --------
        list
            A list of normalized lethargy sensitivity values.
        """

        if isinstance(isotope, str):
            isotope = methods.convert_iso_string_to_id(isotope)
        elif not isinstance(isotope, int):
            raise errors.UserInputError(f"isotope must be a str or an int and was provided as {type(isotope).__name__}")

        if not isinstance(reaction, int):
            raise errors.UserInputError(f"reaction must be an int, but was provided as {type(reaction).__name__}")

        sensi_values = self.get_sensitivity_values(isotope, reaction)
        lethargy_widths = np.diff(np.log(self.e_bins))

        if len(sensi_values) != len(lethargy_widths):
            raise errors.DimError("The length of sensi_values and lethargy_widths do not match.")

        # Using the absolute value of lethargy widths to avoid division by negative values (in case of decreasing binning)
        normalized_lethargy_sensi_values = sensi_values / np.abs(lethargy_widths)

        return list(normalized_lethargy_sensi_values)

    def filter_energy_bins(self, energy_region):
        """
        Filter energy bins based on the specified energy region.

        Parameters
        -----------
        energy_region : str or list
            If str, should be one of 'thermal', 'epithermal', 'fast' or 'all', which correspond
            to predefined energy bounds (in eV):
                - 'thermal': (0.625, float('-inf'))
                - 'epithermal': (1e5, 0.625)
                - 'fast': (float('-inf'), 1e5)
                - 'all' : (float('inf'), float('-inf'))
            If list, should contain two numerical values representing the energy bounds [Emax, Emin] in eV.

        Returns
        --------
        list of int
            Indices of the energy bins that fall within the specified energy region.

        list of float
            Energy bin values that fall within the specified energy region.
        """

        standard_bounds = {
            "thermal": (0.625, float("-inf")),
            "epithermal": (1e5, 0.625),
            "fast": (float("inf"), 1e5),
            "all": (float("inf"), float("-inf")),
        }

        if isinstance(energy_region, str):
            try:
                boundaries = standard_bounds[energy_region]
            except KeyError:
                raise errors.UserInputError("Invalid energy_region. Must be 'thermal', 'epithermal', 'fast' or 'all'.")
        elif isinstance(energy_region, list) and len(energy_region) == 2:
            boundaries = energy_region
        else:
            raise errors.UserInputError("energy_region must be a str ('thermal', 'epithermal', 'fast' or 'all') or a list of two values [Emax, Emin]")

        max_energy, min_energy = boundaries

        if min_energy >= max_energy:
            raise errors.UserInputError("Emax should be greater than Emin")

        filtered_indices = [
            i for i, e in enumerate(self.e_bins) if (min_energy is None or e >= min_energy) and (max_energy is None or e <= max_energy)
        ]
        filtered_ebins = [self.e_bins[i] for i in filtered_indices]

        return filtered_indices, filtered_ebins

    def get_filtered_sensitivity_by_energy_region(self, isotope, reaction, energy_region, normalize_lethargy):
        """
        Filters sensitivity values for a given isotope and reaction based on the specified energy region.

        Parameters
        ----------
        isotope : str or int
            The isotope identifier, which can be a string or an integer.
        reaction : int
            The reaction identifier, provided as an integer.
        energy_region : str
            The energy region to filter the sensitivity values.
        normalize_lethargy : bool
            Flag to indicate whether to use lethargy-normalized sensitivity values.

        Returns
        -------
        list
            A list of filtered sensitivity values corresponding to the specified energy region.
        """
        indices, ebins = self.filter_energy_bins(energy_region)
        sensi_values = (
            self.get_normalized_lethargy_sensitivity(isotope, reaction) if normalize_lethargy else self.get_sensitivity_values(isotope, reaction)
        )
        filtered_sensi_values = [sensi_values[index] for index in indices if index != indices[-1]]

        return filtered_sensi_values


class NDCovariances:
    """
    A class to handle nuclear data covariance matrices from various file formats.

    Attributes
    -----------
    input_path : str
        Path to the covariance file.
    format : str
        Format of the covariance file (coverx, coverx_text, comac, gendf, xlsx).
    cov_dataf : pd.DataFrame
        Covariance data in DataFrame format.
    e_bins : list
        Energy bins boundaries.
    group_nb : int
        Number of energy groups.
    header : str
        Header information from the covariance file.
    iso_reac_list : list of tuples
        List of (isotope_ID, reaction_ID) pairs present in the covariance data.
    iso_list : list
        List of isotopes present in the covariance data.
    reac_list : list
        List of reactions present in the covariance data.

    Methods
    --------
    write_xlsx(output_path):
        Writes the covariance data to an Excel file.
    write_txt(output_path, header):
        Writes the covariance data to a text file in COVERX format.
    """

    @log_exec()
    def __init__(self, input_path: str, format: str = "auto"):
        """
        Initializes the NDCovariances object.

        Parameters
        -----------
        input_path : str
            [Required] Path to the covariance file.
        format : str
            Format of the covariance file. Can be "auto", "coverx", "coverx_text", "comac", "gendf", or "xlsx".
            Defaults to "auto".
            The "xlsx" format allows re-importing Excel files previously exported with write_xlsx().
            String representations of lists in the Excel file are automatically converted back to lists of floats.
        """
        self.input_path = input_path
        self.format = format.lower()

        if self.format == "auto":
            self.format = self._detect_format()
        if self.format == "coverx":
            self.cov_dataf, self.e_bins, self.group_nb, self.header = methods.format_scale_binary_to_dataframe(input_path)
        elif self.format in ["coverx_text", "coverx_txt"]:
            self.cov_dataf, self.e_bins, self.group_nb, self.header = methods.format_scale_txt_to_dataframe(input_path)
        elif self.format == "comac":
            self.cov_dataf, self.e_bins, self.group_nb, self.header = methods.format_comac_to_dataframe(input_path)
        elif self.format == "gendf":
            self.cov_dataf, self.e_bins, self.group_nb, self.header = methods.format_gendf_to_dataframe(input_path)
        elif self.format == "xlsx":
            self.cov_dataf, self.e_bins, self.group_nb, self.header = methods.format_xlsx_to_dataframe(input_path)
        else:
            raise errors.UserInputError(f"Format must be either 'coverx', 'coverx_txt', 'comac', 'gendf', or 'xlsx', but was provided as {format}")

        cov_dataf = self.cov_dataf

        # Extract list of (isotope, reaction) pairs present in covariance data
        iso_reac_cov_H = cov_dataf.apply(lambda x: (x["ISO_H"], x["REAC_H"]), axis=1).to_list()
        iso_reac_cov_V = cov_dataf.apply(lambda x: (x["ISO_V"], x["REAC_V"]), axis=1).to_list()
        all_iso_reac = list(set(iso_reac_cov_H) | set(iso_reac_cov_V))
        # Filter out energy bins entries (0, 0)
        self.iso_reac_list = [(iso, reac) for iso, reac in all_iso_reac if not (iso == 0)]
        self.iso_list = list(set([iso for (iso, reac) in self.iso_reac_list]))
        self.reac_list = list(set([reac for (iso, reac) in self.iso_reac_list]))

    def _detect_format(self) -> str:
        """
        Detects the format of the covariance file based on its content.

        Parameters
        -----------
        input_path : str
            Path to the covariance file.

        Returns
        --------
        str
            Detected format of the covariance file.
        """
        if self.input_path.endswith(".xlsx"):
            return "xlsx"
        elif "scale" in self.input_path or "SCALE" in self.input_path:
            if "txt" in self.input_path:
                return "coverx_text"
            return "coverx"
        elif "comac" in self.input_path or "COMAC" in self.input_path:
            return "comac"
        elif "gendf" in self.input_path or "GENDF" in self.input_path or "endf" in self.input_path or "ENDF" in self.input_path:
            return "gendf"
        else:
            raise errors.UserInputError(
                "Could not detect the format of the covariance file. Please specify the format explicitly among 'coverx', 'coverx_text', 'comac', or 'gendf'."
            )

    def write_xlsx(self, output_path: str):
        """
        Writes the covariance data to an Excel file.

        The file can be re-imported later with NDCovariances(path, format='xlsx').
        String representations of lists will be automatically converted back to lists.

        Parameters
        -----------
        output_path : str
            [Required] Path to save the output Excel file.
        """
        if output_path.endswith(".xlsx"):
            self.cov_dataf.to_excel(output_path)
        elif os.path.isdir(output_path):
            self.cov_dataf.to_excel(os.path.join(output_path, os.path.basename(self.input_path) + ".xlsx"))
        else:
            self.cov_dataf.to_excel(output_path + ".xlsx")

    def write_txt(self, output_path: str, header: str = None):
        """
        Writes the covariance data to a text file in COVERX format.

        Parameters
        -----------
        output_path : str
            [Required] Path to save the output text file.
        header : str
            Header to include in the output file.
        """
        if header == None:
            header = self.header

        with open(output_path, "w") as output_file:
            output_file.write(header)

            if self.e_bins != []:
                for idx, e_bin in enumerate(self.e_bins):
                    if idx > 0 and idx % 5 == 0:
                        e_bin = "{:.7E}".format(e_bin)
                        output_file.write("\n")
                    output_file.write(f"{e_bin: >15}")
                output_file.write("\n")
            else:
                output_file.write(f"The number of groups is {self.group_nb}. Covariance data file doesn't have its energy binning embeded.\n")

            for idx, row in self.cov_dataf.iterrows():
                if row["ISO_H"] == 0 or row["ISO_V"] == 0:
                    continue

                output_file.write(f"{row['ISO_H']: >12}{row['REAC_H']: >12}{row['ISO_V']: >12}{row['REAC_V']: >12}{'1': >12}\n")

                for std_x in row["STD"]:
                    for i, std in enumerate(std_x):
                        if i > 0 and i % 5 == 0:
                            output_file.write("\n")
                        std = "{:.7E}".format(std)
                        output_file.write(f"{std: >15}")
                    output_file.write("\n")


class Assimilation:
    """
    A class to compute the GLLSM assimilation process.

    Attributes
    -----------
    bench_cases : list of Case
        List of benchmark cases.
    bench_list : DataFrame
        DataFrame containing calculated data of the benchmark cases, after assimilation.
    study_case : Case
        Study case.
    cov_data : NDCovariances
        Covariance data object.
    bias : Bias
        Bias object of the study case, after assimilation process.
    prior_uncertainty : Uncertainty
        A priori uncertainty of the study case (as an Unceratinty object), before assimilation.
    post_uncertainty : Uncertainty
        A posteriori uncertainty of the study case (as an Unceratinty object), after assimilation.
    output_html_path : str
        Path to save the output HTML file.
    plotting_unit : str
        Unit for the sensitivity coefficients and uncertainty/bias contributions when plotted. Can be "relative" or "pcm". Defaults to "relative".
    Ck_threshold : float
        Threshold for the Ck value to filter out benchmark cases.
    targetted_chi2 : float
        Targetted chi2 value to to filter out benchmark cases.
    filtering_chi2_method : "delta" or "diagonal"
        Method for filtering chi2 values.
    iso_reac_list : list of tuples if ints
        List of isotope-reaction pairs (iso_ID, reac_ID) to consider for assimilation.
    reac_list : list of int or str, optional
        List of reactions considered for assimilation.
    iso_list : list of int or str, optional
        List of isotopes considered for assimilation.
    exclude_iso : list of int or str, optional
        List of isotopes to exclude from assimilation.
    exclude_reac : list of int or str, optional
        List of reactions to exclude from assimilation.
    isotopes_to_detail : list of int or str, optional
        List of isotopes to provide detailed variances-covariances (as heatmaps) in the HTML outputfile.
    [+ explicit internal attributes...]

    Methods
    --------
    export_to_html(output_html_path, plotting_unit="pcm", isotopes_to_detail=[]):
        Exports the results of the assimilation process to an HTML file with interactive plots.
    plot_study_case_sensi(output_html_path, show=False):
        Plots the sensitivity data for the study case in HTML.
    [+ explicit internal methods...]
    """

    @log_exec()
    def __init__(
        self,
        benchmarks_list,
        cov_data,
        study_case=None,
        iso_reac_list: list = None,
        reac_list: list = None,
        iso_list: list = None,
        exclude_iso: list = None,
        exclude_reac: list = None,
        Ck_threshold=None,
        targetted_chi2=None,
        filtering_chi2_method="delta",
        output_html_path=None,
        plotting_unit="pcm",
        isotopes_to_detail=[],
    ) -> None:
        """
        Initializes the Assimilation object.

        Parameters
        -----------
        benchmarks_list : list
            [Required] List of benchmark cases or sdf paths.
        study_case : Case or sdf path
            [Required] Study case.
        cov_data : NDCovariances
            [Required] Covariance data object.
        output_html_path : str
            Path to save the output HTML file.
        plotting_unit : str
            Unit for the sensitivity coefficients and uncertainty/bias contributions when plotted. Can be "relative" or "pcm". Defaults to "relative".
        Ck_threshold : float, optional
            Threshold for the Ck value to filter out benchmark cases.
        targetted_chi2 : float, optional
            Targetted chi2 value to filter out benchmark cases.
        filtering_chi2_method : "delta" or "diagonal"
            Method for filtering based on chi2 values. Defaults to "delta".
        iso_reac_list : list of tuples (int ISO, int REAC), optional
            List of isotope-reaction pairs (iso_ID, reac_ID) to consider for assimilation.
        reac_list : list of int or str, optional
            List of reactions to consider for assimilation.
        iso_list : list of int or str, optional
            List of isotopes to consider for assimilation.
        exclude_iso : list of int or str, optional
            List of isotopes to exclude from assimilation.
        exclude_reac : list of int or str, optional
            List of reactions to exclude from assimilation.
        isotopes_to_detail : list of int or str, optional
            List of isotopes to provide detailed variances-covariances (as heatmaps) in the HTML outputfile.
        """

        self.bench_cases = benchmarks_list[:]

        # Store the original covariance data object
        self.cov_data = cov_data

        for i, bench_data in enumerate(self.bench_cases):
            if isinstance(bench_data, (Path, str)):

                self.bench_cases[i] = Case(sdf_path=bench_data)

            elif isinstance(bench_data, Case):
                None

            else:
                raise errors.SensInputError(
                    f"Wrong sensitivity type for {i+1}th case of the benchmarks list : {bench_data} - Choose case object, Path, or string"
                )

            if self.bench_cases[i].resp_expe == None or self.bench_cases[i].sigma_resp_expe == None:
                raise errors.SensInputError(
                    f"The case {self.bench_cases[i].casename} is used as a experimental benchmark case but doesn't have an experimental response value/sigma - Check you sdf file or your case construction (arguments)"
                )

        if isinstance(study_case, (Path, str)):
            study_case = Case(sdf_path=study_case)

        elif isinstance(study_case, Case):
            study_case = study_case
        elif study_case == None:
            None
        else:
            raise errors.SensInputError("Wrong sensitivity type for study case data - Choose case, Path, or string")

        if iso_reac_list != None:
            for i, pair in enumerate(iso_reac_list):
                if not isinstance(pair, tuple):
                    raise errors.UserInputError(f"The isotope-reaction pair '{pair}' should be a tuple.")

        self.output_html_path = output_html_path

        self.Ck_threshold = Ck_threshold

        self.targetted_chi2 = targetted_chi2

        if filtering_chi2_method not in ["delta", "diagonal"]:
            raise Exception('The filtering chi2 method available currently are "delta" or "diagonal"')

        self.filtering_chi2_method = filtering_chi2_method

        self.bench_list = pd.DataFrame({"PATH": [bench_case.sdf_path for bench_case in self.bench_cases]})
        bench_resp_expe, bench_sigma_expe, bench_resp_calc, bench_sigma_calc, bench_EC_prior = [], [], [], [], []
        for bench in self.bench_cases:
            bench_resp_expe.append(bench.resp_expe)
            bench_sigma_expe.append(bench.sigma_resp_expe)
            bench_resp_calc.append(bench.resp_calc)
            bench_sigma_calc.append(bench.sigma_resp_calc)
            bench_EC_prior.append(round((bench.resp_expe - bench.resp_calc) * 1e5))

        (
            self.bench_list["RESP EXPE"],
            self.bench_list["SIGMA RESP EXPE"],
            self.bench_list["RESP CALC"],
            self.bench_list["SIGMA RESP CALC"],
            self.bench_list["E - C_PRIOR (pcm)"],
        ) = (bench_resp_expe, bench_sigma_expe, bench_resp_calc, bench_sigma_calc, bench_EC_prior)

        self.study_case = study_case

        self.make_sensimat_covmat_casevec_expemat_and_deltaCE(
            iso_reac_list=iso_reac_list,
            reac_list=reac_list,
            iso_list=iso_list,
            exclude_iso=exclude_iso,
            exclude_reac=exclude_reac,
        )

        self.check_dimmensions()

        methods.check_correspondences(sensi_vec=self.case_vec, cov_mat=self.cov_mat, iso_reac_list=self.iso_reac_list, group_nb=self.group_nb)

        if self.study_case != None:
            self.e_bins = self.study_case.e_bins
        else:
            self.e_bins = self.bench_cases[0].e_bins

        self.calcul_prior_uncertainty()
        self.calcul_matrix_assimilation()
        self.calcul_bias()
        self.calcul_post_chi2()
        self.calcul_post_uncertainty()
        self.export_to_html(output_html_path=self.output_html_path, plotting_unit=plotting_unit, isotopes_to_detail=isotopes_to_detail)

    def make_sensimat_covmat_casevec_expemat_and_deltaCE(
        self, iso_reac_list: list = None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None, exclude_reac: list = None
    ):
        """Create sensitivity matrix, covariance matrix, case vector, experimental matrix, and delta C-E.

        Parameters
        ----------
        iso_reac_list : list of tuples (int ISO, int REAC), optional
            List of (isotope, reaction) tuples to consider
        reac_list : list, optional
            List of reactions to consider
        iso_list : list, optional
            List of isotopes to consider
        exclude_iso : list, optional
            List of isotopes to exclude
        exclude_reac : list, optional
            List of reactions to exclude
        """
        # Construct the sensitivity vectors and covariances matrix from the iso_reac_list given or the union of every case's iso_reac_list
        if self.study_case != None:
            sensi_vecs, cov_mat, iso_reac_list = methods.make_sensi_vectors_and_cov_matrix(
                cases_list=[*self.bench_cases, self.study_case],
                cov_data=self.cov_data,
                operation="union",
                iso_reac_list=iso_reac_list,
                reac_list=reac_list,
                iso_list=iso_list,
                exclude_iso=exclude_iso,
                exclude_reac=exclude_reac,
            )
            benchmark_vecs = sensi_vecs[: len(self.bench_cases)]
            case_vec = sensi_vecs[-1]
        else:
            sensi_vecs, cov_mat, iso_reac_list = methods.make_sensi_vectors_and_cov_matrix(
                cases_list=[*self.bench_cases],
                cov_data=self.cov_data,
                operation="union",
                iso_reac_list=iso_reac_list,
                reac_list=reac_list,
                iso_list=iso_list,
                exclude_iso=exclude_iso,
                exclude_reac=exclude_reac,
            )
            benchmark_vecs = sensi_vecs
            case_vec = []

        # Construct the sensitivity benchmarks matrix, the experimental response matrix, and the C-E vector
        sensi_matrix = []
        for benchmark_vec in benchmark_vecs:
            sensi_matrix.append(benchmark_vec)

        self.bench_sensi_mat = np.array(sensi_matrix)
        self.case_vec = np.array(case_vec)
        self.cov_mat = cov_mat
        self.iso_reac_list = iso_reac_list
        self.iso_list = list(set([iso for (iso, reac) in self.iso_reac_list]))
        self.reac_list = list(set([reac for (iso, reac) in self.iso_reac_list]))

        expe_mat = np.zeros((len(self.bench_cases), len(self.bench_cases)), dtype=float)
        delta_CE_vec = np.zeros(len(self.bench_cases), dtype=float)

        for idx, bench_case in enumerate(self.bench_cases):
            # expe_mat[idx, idx] = bench_case.sigma_resp_expe**2 / bench_case.resp_expe**2
            expe_mat[idx, idx] = bench_case.sigma_resp_expe**2 / bench_case.resp_calc**2
            # delta_CE_vec[idx] = (bench_case.resp_calc - bench_case.resp_expe) / bench_case.resp_expe
            delta_CE_vec[idx] = (bench_case.resp_calc - bench_case.resp_expe) / bench_case.resp_calc

        self.expe_mat = np.array(expe_mat)
        self.delta_CE_vec = np.array(delta_CE_vec)

    def check_dimmensions(self):
        if (
            not len(self.bench_cases)
            == np.shape(self.bench_sensi_mat)[0]
            == np.shape(self.expe_mat)[0]
            == np.shape(self.expe_mat)[1]
            == len(self.delta_CE_vec)
        ):

            raise errors.DimError("The vertical dimension of the benchmarks matrix is not equal to the length of the benchmark list")

        if self.study_case != None:
            if not np.shape(self.cov_mat)[0] == np.shape(self.cov_mat)[1] == len(self.case_vec) == np.shape(self.bench_sensi_mat)[1]:

                raise errors.DimError(
                    f"The sensitivities and uncertainties data doesn't have the same dimensions - Check your energy groups number - (or internal iso-reac variable error)\n\
                    DIMENSIONS : Cov Mat {np.shape(self.cov_mat)[0]} (should be square matrix) | Study Case Vec {len(self.case_vec)} | Bench Mat (vertical axe) {np.shape(self.bench_sensi_mat)[1]}"
                )

            if self.study_case.group_nb != int(np.shape(self.cov_mat)[0] / len(self.iso_reac_list)):
                raise errors.DimError(
                    f"The energy group number of your study case is not equal to the covariance matrix group number - Check your sdf file {self.study_case.casename}, and you covariance file"
                )
        else:
            if not np.shape(self.cov_mat)[0] == np.shape(self.cov_mat)[1] == np.shape(self.bench_sensi_mat)[1]:

                raise errors.DimError(
                    f"The sensitivities and uncertainties data doesn't have the same dimensions - Check your energy groups number - (or internal iso-reac variable error)\n\
                    DIMENSIONS : Cov Mat {np.shape(self.cov_mat)[0]} (should be square matrix) | Study Case Vec {len(self.case_vec)} | Bench Mat (vertical axe) {np.shape(self.bench_sensi_mat)[1]}"
                )

        if self.study_case != None:
            self.group_nb = self.study_case.group_nb
        else:
            self.group_nb = self.bench_cases[0].group_nb

    def calcul_matrix_assimilation(self):

        # Function to filter out experiments with high contributions to a priori global Chi2
        def filter_chi2(C_SexpW0SexpT_inv):

            individual_chi2 = []

            if self.filtering_chi2_method == "delta":
                for i in range(len(self.bench_cases)):

                    C_i = np.delete(self.expe_mat, (i), axis=0)
                    C_i = np.delete(C_i, (i), axis=1)
                    Sexp_i = np.delete(self.bench_sensi_mat, (i), axis=0)
                    delta_CE_i = np.delete(self.delta_CE_vec, (i))

                    W0SexpT_i = self.cov_mat @ Sexp_i.T
                    SexpW0SexpT_i = Sexp_i @ W0SexpT_i
                    C_SexpW0SexpT_i = np.add(C_i, SexpW0SexpT_i)

                    # INVERSION
                    try:
                        cholesky = np.linalg.inv(np.linalg.cholesky(C_SexpW0SexpT_i))
                        C_SexpW0SexpT_inv_i = cholesky.T @ cholesky

                    except:
                        write_and_print(
                            "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
                        )
                        raise Exception(
                            "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
                        )

                    individual_chi2.append(self.prior_chi2 * len(self.bench_cases) - delta_CE_i @ C_SexpW0SexpT_inv_i @ delta_CE_i.T)

            elif self.filtering_chi2_method == "diagonal":
                for i in range(len(self.bench_cases)):

                    individual_chi2.append((self.bench_cases[i].resp_calc - self.bench_cases[i].resp_expe) ** 2 * C_SexpW0SexpT_inv[i][i])

            self.bench_list["INDIVIDUAL CHI2"] = individual_chi2

            # --- Those arrays are changing over time, everytime a bench case is removed
            C_SexpW0SexpT_inv_i = None
            chi2_list, C_i, Sexp_i, delta_CE_i = (
                copy.deepcopy(individual_chi2),
                copy.deepcopy(self.expe_mat),
                copy.deepcopy(self.bench_sensi_mat),
                copy.deepcopy(self.delta_CE_vec),
            )
            for i in range(len(self.bench_cases)):

                if len(chi2_list) == 1:
                    raise errors.UserInputError(
                        "Your Chi2 filtering threshold is leading to the removal of every benchmark case of your list. Please, choose other benchmark cases, or increasing your threshold."
                    )
                indiv_chi2_idx = chi2_list.index(max(chi2_list))

                bench_path_to_remove = list(self.bench_list[self.bench_list["INDIVIDUAL CHI2"] == max(chi2_list)]["PATH"])[0]
                self.bench_list.at[self.bench_list.index[self.bench_list["INDIVIDUAL CHI2"] == max(chi2_list)][0], "REMOVED"] = True

                write_and_print(
                    f"{2*'    '}Assimilation - Removing benchmark {os.path.basename(bench_path_to_remove)} from assimilation ...\n - Indivudual Chi2 to remove : {max(chi2_list)}"
                )

                C_i = np.delete(C_i, (indiv_chi2_idx), axis=0)
                C_i = np.delete(C_i, (indiv_chi2_idx), axis=1)
                Sexp_i = np.delete(Sexp_i, (indiv_chi2_idx), axis=0)
                delta_CE_i = np.delete(delta_CE_i, (indiv_chi2_idx))
                chi2_list.pop(indiv_chi2_idx)

                W0SexpT_i = self.cov_mat @ Sexp_i.T
                SexpW0SexpT_i = Sexp_i @ W0SexpT_i
                C_SexpW0SexpT_i = np.add(C_i, SexpW0SexpT_i)

                try:
                    cholesky = np.linalg.inv(np.linalg.cholesky(C_SexpW0SexpT_i))
                    C_SexpW0SexpT_inv_i = cholesky.T @ cholesky
                except:
                    write_and_print(
                        "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
                    )
                    raise Exception(
                        "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
                    )

                chi2 = delta_CE_i @ C_SexpW0SexpT_inv_i @ delta_CE_i.T
                chi2 = chi2 / len(chi2_list)

                if chi2 < self.targetted_chi2:
                    self.prior_chi2 = chi2
                    break

                write_and_print(f"{2*'    '}Assimilation - Chi2 removal {i+1} result : {chi2}")

            self.bench_sensi_mat = Sexp_i
            self.expe_mat = C_i
            self.delta_CE_vec = delta_CE_i

            self.delta_mu = -1 * W0SexpT_i @ (C_SexpW0SexpT_inv_i @ self.delta_CE_vec)

            self.cov_mat_delta = csr_matrix((W0SexpT_i @ C_SexpW0SexpT_inv_i) @ (self.bench_sensi_mat @ self.cov_mat))

        # Function to filter out experiments with low Ck similarity coefficients
        def filter_Ck():

            C_SexpW0SexpT_inv_i = None
            C_i, Sexp_i, delta_CE_i = (
                copy.deepcopy(self.expe_mat),
                copy.deepcopy(self.bench_sensi_mat),
                copy.deepcopy(self.delta_CE_vec),
            )

            for i in range(len(self.bench_list)):

                if self.bench_list["Ck"][i] < self.Ck_threshold and self.bench_list.at[i, "REMOVED"] == False:
                    self.bench_list.at[i, "REMOVED"] = True

                    write_and_print(
                        f"{2*'    '}Assimilation - Removing benchmark {os.path.basename(self.bench_list['PATH'][i])} from assimilation ...\n - Ck to remove : {self.bench_list['Ck'][i]}"
                    )

                    indiv_ck_idx = i - sum(self.bench_list["REMOVED"][:i])

                    C_i = np.delete(C_i, (indiv_ck_idx), axis=0)
                    C_i = np.delete(C_i, (indiv_ck_idx), axis=1)
                    Sexp_i = np.delete(Sexp_i, (indiv_ck_idx), axis=0)
                    delta_CE_i = np.delete(delta_CE_i, (indiv_ck_idx))

            self.bench_sensi_mat = Sexp_i
            self.expe_mat = C_i
            self.delta_CE_vec = delta_CE_i

            W0SexpT = self.cov_mat @ self.bench_sensi_mat.T

            SexpW0SexpT = self.bench_sensi_mat @ W0SexpT

            C_SexpW0SexpT = np.add(self.expe_mat, SexpW0SexpT)

            cholesky = np.linalg.inv(np.linalg.cholesky(C_SexpW0SexpT))
            C_SexpW0SexpT_inv = cholesky.T @ cholesky

            self.delta_mu = -1 * W0SexpT @ (C_SexpW0SexpT_inv @ self.delta_CE_vec)

            self.cov_mat_delta = csr_matrix((W0SexpT @ C_SexpW0SexpT_inv) @ (self.bench_sensi_mat @ self.cov_mat))

        # Beginning of the assimilation process
        W0SexpT = self.cov_mat @ self.bench_sensi_mat.T

        SexpW0SexpT = self.bench_sensi_mat @ W0SexpT
        self.bench_list["UNC_PRIOR (pcm)"] = [(bench.resp_calc * sqrt(unc)) * 1e5 for bench, unc in zip(self.bench_cases, np.diag(SexpW0SexpT))]

        C_SexpW0SexpT = np.add(self.expe_mat, SexpW0SexpT)

        # Inversion
        try:
            cholesky = np.linalg.inv(np.linalg.cholesky(C_SexpW0SexpT))
            C_SexpW0SexpT_inv = cholesky.T @ cholesky

        except:
            # C_SexpW0SexpT_inv = np.linalg.pinv(C_SexpW0SexpT)
            write_and_print(
                "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
            )
            raise Exception(
                "ERROR : The covarariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
            )

        # Compute a priori global Chi2
        chi2 = self.delta_CE_vec @ C_SexpW0SexpT_inv @ self.delta_CE_vec.T
        chi2 = chi2 / len(self.bench_cases)
        self.prior_chi2 = chi2

        write_and_print(f"{2*'    '}Assimilation - initial Chi2 = {self.prior_chi2}")

        self.bench_list["INDIVIDUAL CHI2"] = [None for i in range(len(self.bench_cases))]
        self.bench_list["REMOVED"] = [False for i in range(len(self.bench_cases))]

        # Compute Ck similarity coefficients between the study case and the benchmark cases
        Ck_bench = []
        norm_sensi_case = self.case_vec @ self.cov_mat @ self.case_vec
        for i in range(len(self.bench_cases)):
            norm_sensi_bench = SexpW0SexpT[i, i]
            if norm_sensi_case == 0 or norm_sensi_bench == 0:
                Ck_bench.append(0)
            else:
                Ck_bench.append(sqrt(abs((self.case_vec @ self.cov_mat @ self.bench_sensi_mat[i, :]) ** 2 / (norm_sensi_case * norm_sensi_bench))))

        self.bench_list["Ck"] = Ck_bench

        # Filtering on chi2 and Ck if needed
        if self.targetted_chi2 != None and self.prior_chi2 > self.targetted_chi2:
            filter_chi2(C_SexpW0SexpT_inv=C_SexpW0SexpT_inv)

        if self.Ck_threshold != None:
            filter_Ck()

        if self.Ck_threshold == None and self.targetted_chi2 == None:

            self.delta_mu = -1 * W0SexpT @ (C_SexpW0SexpT_inv @ self.delta_CE_vec)

            self.cov_mat_delta = csr_matrix((W0SexpT @ C_SexpW0SexpT_inv) @ (self.bench_sensi_mat @ self.cov_mat))

        self.cov_mat_post = np.subtract(self.cov_mat, self.cov_mat_delta)

    @log_exec()
    def calcul_post_chi2(self):
        # Compute a posteriori global chi2
        len_bench_cases = len([i for i in self.bench_list["REMOVED"] if i == False])
        expe_mat_post = np.zeros((len_bench_cases, len_bench_cases), dtype=float)
        delta_CE_vec_post = np.zeros(len_bench_cases, dtype=float)
        b_present = 0
        for b in range(len(self.bench_list)):
            if self.bench_list["REMOVED"][b] == False:
                sigma_expe = self.bench_list["SIGMA RESP EXPE"].to_list()[b]
                resp_expe, EC_post = self.bench_list["RESP EXPE"].to_list()[b], self.bench_list["E - C_POST (pcm)"].to_list()[b]
                resp_calc_post = resp_expe - EC_post * 1e-5
                expe_mat_post[b_present, b_present] = sigma_expe**2 / resp_calc_post**2
                delta_CE_vec_post[b_present] = (resp_calc_post - self.bench_list["RESP EXPE"].to_list()[b]) / resp_calc_post
                b_present += 1

        W0SexpT = self.cov_mat_post @ self.bench_sensi_mat.T

        SexpW0SexpT = self.bench_sensi_mat @ W0SexpT

        C_SexpW0SexpT = np.add(expe_mat_post, SexpW0SexpT)

        # Inversion
        try:
            cholesky = np.linalg.inv(np.linalg.cholesky(C_SexpW0SexpT))
            C_SexpW0SexpT_inv = cholesky.T @ cholesky

        except:
            # C_SexpW0SexpT_inv = np.linalg.pinv(C_SexpW0SexpT)
            write_and_print(
                "ERROR : The covariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
            )
            raise Exception(
                "ERROR : The covariances data you choose leads to an non-invertible matrix (singular matrix). The assimimlation process is not possible."
            )

        # Calcul initial Chi2
        chi2 = delta_CE_vec_post @ C_SexpW0SexpT_inv @ delta_CE_vec_post.T
        chi2 = chi2 / len_bench_cases
        self.post_chi2 = chi2

    @log_exec()
    def calcul_bias(self):

        bench_biases = self.bench_sensi_mat @ self.delta_mu

        bench_EC_post = []
        idx_present = 0
        for b in range(len(self.bench_list)):
            bench = self.bench_cases[b]
            if self.bench_list["REMOVED"][b] == False:
                bench_EC_post.append(round((bench.resp_expe - (bench.resp_calc + bench_biases[idx_present] * bench.resp_calc)) * 1e5))
                idx_present += 1
            else:
                bench_EC_post.append("Not calculated")

        self.bench_list["E - C_POST (pcm)"] = bench_EC_post

        if self.study_case != None:
            self.bias = Bias(
                study_case=self.study_case,
                sensi_vec=self.case_vec,
                resp_calc=self.study_case.resp_calc,
                delta_mu=self.delta_mu,
                iso_reac_list=self.iso_reac_list,
            )

        return self.bias

    @log_exec()
    def calcul_post_uncertainty(self):

        W1SexpT = self.cov_mat_post @ self.bench_sensi_mat.T
        SexpW1SexpT = self.bench_sensi_mat @ W1SexpT

        bench_unc_post = []
        idx_present = 0
        for b in range(len(self.bench_list)):
            bench = self.bench_cases[b]
            if self.bench_list["REMOVED"][b] == False:
                bench_unc_post.append((bench.resp_calc * sqrt(SexpW1SexpT[idx_present, idx_present])) * 1e5)
                idx_present += 1
            else:
                bench_unc_post.append("Not calculated")
        self.bench_list["UNC_POST (pcm)"] = bench_unc_post

        if self.study_case != None:
            self.post_uncertainty = Uncertainty(
                study_case=self.study_case,
                cov_data=self.cov_data,
                sensi_vec=self.case_vec,
                resp_calc=self.study_case.resp_calc,
                cov_mat=self.cov_mat_post,
                iso_reac_list=self.iso_reac_list,
            )
            return self.post_uncertainty

    @log_exec()
    def calcul_prior_uncertainty(self):
        if self.study_case != None:
            self.prior_uncertainty = Uncertainty(
                study_case=self.study_case,
                cov_data=self.cov_data,
                sensi_vec=self.case_vec,
                resp_calc=self.study_case.resp_calc,
                cov_mat=self.cov_mat,
                iso_reac_list=self.iso_reac_list,
            )
            return self.prior_uncertainty

    def plot_study_case_sensi(self, output_html_path: str = None, show=False):
        """
        Plot the sensitivity of the study case.

        Parameters
        -----------
        output_html_path : str, optional
            Path to save the output HTML file.
        show : bool, optional
            Flag to display the plot. Defaults to False.
        """
        self.study_case.export_to_html(output_html_path=output_html_path, show=show)

    @log_exec()
    def export_to_html(self, output_html_path: str, plotting_unit="pcm", isotopes_to_detail=[]):
        """
        Export the results of the assimilation process to an HTML file.

        Parameters
        -----------
        output_html_path : str
            [Required] Path to save the output HTML file.
        plotting_unit : str, optional
            Unit for the sensitivity coefficients and uncertainty/bias contributions when plotted. Can be "relative" or "pcm". Defaults to "relative".
        isotopes_to_detail : list, optional
            List of isotopes (ID) to provide detailed variances-covariances (as heatmaps) in the HTML outputfile.
        """
        if output_html_path is not None:
            if not output_html_path.endswith(".html"):
                output_html_path = output_html_path + ".html"
            with open(output_html_path, "w") as f:
                None
        else:
            return None

        if not isinstance(isotopes_to_detail, list):
            isotopes_to_detail = [isotopes_to_detail]

        # --------------------------------
        text_intro = ""
        text_intro += f"<h1><center>{pkg_fullname}</center></h1>"
        text_intro += f"<h1>GLLSM assimilation results</h1>"
        if self.study_case != None:
            text_intro += f"<h2>Study case : {self.study_case.casename}</h2>"
        else:
            text_intro += f"<h2>No study case</h2>"

        # --------------------------------
        headers = [
            "Calculated response of study case",
            "Uncertainty a priori",
            "Uncertainty a posteriori",
            "Evaluated bias from GLLSM (C<sub>post</sub> - C<sub>prior</sub>)",
            "Chi2 a priori ",
            "Chi2 a posteriori",
            "Nb of benchmark cases removed",
            "Energy groups number",
        ]
        headers_b = ["<b>" + h + ":</b>" for h in headers]

        if self.study_case != None:
            results = [
                f"{self.study_case.resp_calc} +/- {self.study_case.sigma_resp_calc}",
                f"{int(round(self.prior_uncertainty.value))} pcm",
                f"{int(round(self.post_uncertainty.value))} pcm",
                f"{int(round(self.bias.value))} pcm ",
                round(self.prior_chi2, 5),
                round(self.post_chi2, 5),
                list(self.bench_list["REMOVED"]).count(True),
                self.group_nb,
            ]
        else:
            results = [
                "None",
                "None",
                "None",
                "None",
                round(self.prior_chi2, 5),
                round(self.post_chi2, 5),
                list(self.bench_list["REMOVED"]).count(True),
                self.group_nb,
            ]

        table_res = plots.create_html_table(headers=[], lines=[headers_b, results])
        write_and_print("\n" + tabulate(np.array([["Casename", *headers], [self.study_case.casename, *results]]).T))

        # --------------------------------
        bench_list_included = self.bench_list[self.bench_list["REMOVED"] == False]
        bench_list_custom = pd.concat(
            [
                pd.DataFrame(
                    {
                        "PATH": bench_list_included["PATH"],
                        "E - C (pcm)": bench_list_included["E - C_PRIOR (pcm)"],
                        "PRIOR/POST": "PRIOR",
                        "3SIGMA EXPE+ND (pcm)": 3
                        * np.sqrt(
                            (bench_list_included["SIGMA RESP EXPE"].astype(float) * 1e5) ** 2
                            + bench_list_included["UNC_PRIOR (pcm)"].astype(float) ** 2
                        ),
                    }
                ),
                pd.DataFrame(
                    {
                        "PATH": bench_list_included["PATH"],
                        "E - C (pcm)": bench_list_included["E - C_POST (pcm)"],
                        "PRIOR/POST": "POST",
                        "3SIGMA EXPE+ND (pcm)": 3
                        * np.sqrt(
                            (bench_list_included["SIGMA RESP EXPE"].astype(float) * 1e5) ** 2
                            + bench_list_included["UNC_POST (pcm)"].astype(float) ** 2
                        ),
                    }
                ),
            ],
            ignore_index=True,
        )

        trace_bias = px.scatter(
            bench_list_custom,
            x="PATH",
            y="E - C (pcm)",
            error_y="3SIGMA EXPE+ND (pcm)",
            color="PRIOR/POST",
            symbol="PRIOR/POST",
            hover_name="PATH",
            title="Benchmark cases bias progression",
        )
        trace_bias.update_yaxes(title_text="E - C (pcm) (3 sigma=sqrt(expe²+ND²))")
        trace_bias.update_traces(marker_size=8, error_y_thickness=0.5)
        trace_bias.update_xaxes(tickvals=bench_list_custom["PATH"], ticktext=[os.path.basename(path) for path in bench_list_custom["PATH"]])
        trace_bias.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        trace_sim = px.scatter(
            bench_list_included,
            x="Ck",
            y="E - C_PRIOR (pcm)",
            error_y=3
            * np.sqrt((bench_list_included["SIGMA RESP EXPE"].astype(float) * 1e5) ** 2 + bench_list_included["UNC_PRIOR (pcm)"].astype(float) ** 2),
            hover_name="PATH",
            title="Benchmark cases similarity with study case (Ck) vs bias (E-C)",
        )

        x_data = bench_list_included["Ck"].values
        y_data = bench_list_included["E - C_PRIOR (pcm)"].values

        coeffs = np.polyfit(x_data, y_data, 1)
        y_extrapolated = coeffs[0] * 1.0 + coeffs[1]

        x_trendline = np.linspace(min(x_data), 1.0, 100)
        y_trendline = coeffs[0] * x_trendline + coeffs[1]
        trace_sim.add_scatter(
            x=x_trendline, y=y_trendline, mode="lines", line=dict(dash="dash", color="gray"), name="Linear extrapolation", showlegend=True
        )

        trace_sim.add_scatter(
            x=[1.0],
            y=[y_extrapolated],
            mode="markers+text",
            marker=dict(size=10, color="gray", symbol="star"),
            text=[f"{y_extrapolated:.0f} pcm"],
            textposition="top center",
            name="Extrapolated bias to Ck=1",
            showlegend=True,
        )

        x_data_3s = x_data
        y_data_3s = (
            y_data
            + 3
            * np.sqrt(
                (bench_list_included["SIGMA RESP EXPE"].astype(float) * 1e5) ** 2 + bench_list_included["UNC_PRIOR (pcm)"].astype(float) ** 2
            ).values
        )

        coeffs_3s = np.polyfit(x_data_3s, y_data_3s, 1)
        y_extrapolated_3s = coeffs_3s[0] * 1.0 + coeffs_3s[1]

        x_trendline_3s = np.linspace(min(x_data_3s), 1.0, 100)
        y_trendline_3s = coeffs_3s[0] * x_trendline_3s + coeffs_3s[1]

        trace_sim.add_scatter(
            x=x_trendline_3s,
            y=y_trendline_3s,
            mode="lines",
            line=dict(dash="dash", color="green"),
            name="Conservative linear extrapolation of bias + 3*sigma",
            showlegend=True,
        )

        trace_sim.add_scatter(
            x=[1.0],
            y=[y_extrapolated_3s],
            mode="markers+text",
            marker=dict(size=10, color="green", symbol="star"),
            text=[f"{y_extrapolated_3s:.1f} pcm"],
            textposition="bottom center",
            name="Conservative extrapolated bias+3*sigma to Ck=1",
            showlegend=True,
        )

        trace_sim.update_yaxes(title_text="E - C (pcm) (3 sigma=sqrt(expe²+ND²))")
        trace_sim.update_traces(marker_size=8, error_y_thickness=0.5)
        trace_sim.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        headers = [
            "Benchmark cases",
            "Resp expe",
            "Resp calc (prior)",
            "Resp calc (post)",
            "ND uncertainty (prior) (pcm)",
            "ND uncertainty (post) (pcm)",
            "E - C<sub>prior</sub> (pcm)",
            "E - C<sub>post</sub> (pcm)",
            "Chi2 individual",
            "Ck",
        ]
        lines = [
            self.bench_list["PATH"],
            [f"{round(x.resp_expe, 5)} +/- {round(x.sigma_resp_expe, 5)}" for x in self.bench_cases],
            [f"{round(x.resp_calc, 5)} +/- {round(x.sigma_resp_calc, 5)}" for x in self.bench_cases],
            [
                f"{round(x_case.resp_expe - EC_post/1E5, 5)} +/- {round(x_case.sigma_resp_calc, 5)}" if type(EC_post) is not str else None
                for x_case, EC_post in zip(self.bench_cases, self.bench_list["E - C_POST (pcm)"])
            ],
            [int(round(x)) for x in self.bench_list["UNC_PRIOR (pcm)"]],
            [int(round(x)) if type(x) is not str else None for x in self.bench_list["UNC_POST (pcm)"]],
            [int(round(x)) for x in self.bench_list["E - C_PRIOR (pcm)"]],
            [int(round(x)) if type(x) is not str else None for x in self.bench_list["E - C_POST (pcm)"]],
            [round(chi2, 5) if chi2 != None else None for chi2 in self.bench_list["INDIVIDUAL CHI2"]],
            [round(ck, 4) if ck != None else None for ck in self.bench_list["Ck"]],
        ]
        table_bench_fig = plots.create_html_table(
            headers=headers,
            lines=lines,
            color_per_lines=["red" if rem else "white" for rem in self.bench_list["REMOVED"]],
        )

        lines[0] = [os.path.basename(p) for p in self.bench_list["PATH"]]
        write_and_print("\n" + tabulate(np.array(lines).T, headers=headers))

        # --------------------------------
        if plotting_unit == "relative":
            val_factor = 1.0
            unit_str = "%[resp]"
        elif plotting_unit == "pcm":
            val_factor = self.study_case.resp_calc * 1e5
            unit_str = "pcm"

        if self.study_case != None:

            trace_case_sensi_integrals = plots.plot_integrals_per_iso_reac(
                vector=self.case_vec,
                iso_reac_list=self.iso_reac_list,
                group_nb=self.group_nb,
                factor=val_factor,
                title="Study case sensitivities (integrals)",
                yaxis_title=f"Sensitivity ({unit_str}/%[ND] - group-wise integral)",
            )
            trace_case_sensi_integrals.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            trace_case_sensi_profiles = plots.plot_profiles_per_iso_reac(
                vector=self.case_vec,
                iso_reac_list=self.iso_reac_list,
                e_bins=self.e_bins,
                factor=val_factor,
                title="Study case sensitivities",
                yaxis_title=f"Sensitivity ({unit_str}/%[ND])",
            )
            trace_case_sensi_profiles.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            case_vec_per_unit_lethargy = []
            for iso, reac in self.iso_reac_list:
                if (iso, reac) in self.study_case.iso_reac_list:
                    case_vec_per_unit_lethargy += self.study_case.get_normalized_lethargy_sensitivity(iso, reac)
                else:
                    case_vec_per_unit_lethargy += [0 for i in range(self.group_nb)]
            case_vec_per_unit_lethargy = np.array(case_vec_per_unit_lethargy)

            trace_case_sensi_profiles_lethargy = plots.plot_profiles_per_iso_reac(
                vector=case_vec_per_unit_lethargy,
                iso_reac_list=self.iso_reac_list,
                e_bins=self.e_bins,
                factor=val_factor,
                title="Sensitivities per unit lethargy",
                yaxis_title=f"Sensitivity ({unit_str}/%[ND] per unit lethargy)",
            )
            trace_case_sensi_profiles_lethargy.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            # --------------------------------
            dec_prior_trace_integ = plots.plot_integrals_per_iso_reac(
                vector=self.prior_uncertainty.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                group_nb=self.group_nb,
                factor=val_factor**2,
                title="Integral contribution to relative a priori uncertainy² (covariances included)",
                yaxis_title=f"Unc² ({unit_str}²) (covariances included)",
            )

            dec_prior_trace_integ.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            dec_prior_trace = plots.plot_profiles_per_iso_reac(
                vector=self.prior_uncertainty.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                e_bins=self.e_bins,
                factor=val_factor**2,
                title="Contribution to relative a priori uncertainy² (covariances included)",
                yaxis_title=f"Unc² ({unit_str}²) (covariances included)",
            )
            dec_prior_trace.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            # --------------------------------
            dec_post_trace_integ = plots.plot_integrals_per_iso_reac(
                vector=self.post_uncertainty.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                group_nb=self.group_nb,
                factor=val_factor**2,
                title="Integral contribution to relative a posteriori uncertainy² (covariances included)",
                yaxis_title=f"Unc² ({unit_str}²)",
            )

            dec_post_trace_integ.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            dec_post_trace = plots.plot_profiles_per_iso_reac(
                vector=self.post_uncertainty.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                e_bins=self.e_bins,
                factor=val_factor**2,
                title="Contribution to relative a posteriori uncertainy² (covariances included)",
                yaxis_title=f"Unc² ({unit_str}²)",
            )
            dec_post_trace.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            # --------------------------------
            dec_bias_trace_integ = plots.plot_integrals_per_iso_reac(
                vector=self.bias.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                group_nb=self.group_nb,
                factor=val_factor,
                title="Integral contribution to relative bias",
                yaxis_title=f"Bias ({unit_str})",
            )
            dec_bias_trace_integ.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

            dec_bias_trace = plots.plot_profiles_per_iso_reac(
                vector=self.bias.decomp_vec,
                iso_reac_list=self.iso_reac_list,
                e_bins=self.e_bins,
                factor=val_factor,
                title="Contribution to relative bias",
                yaxis_title=f"Bias ({unit_str})",
            )
            dec_bias_trace.update_layout(
                {
                    "height": plot_height,
                    "width": plot_width,
                    "font_size": font_size,
                    "template": plotly_template,
                    "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                }
            )

        # --------------------------------
        trace_cov_prior_integrals = plots.plot_matrix_integrals_per_iso_reac(
            cov_mat=self.cov_mat,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            title="Variances-covariances matrix before assimilation (group-wise integrals - covariances included)",
            yaxis_title="Unc(%)² (group-wise integral - covariances included)",
        )
        trace_cov_prior_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_diag = np.diag(self.cov_mat.toarray())
        trace_matrix_profiles = plots.plot_profiles_per_iso_reac(
            vector=matrix_diag,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix before assimilation (diagonal elements)",
            yaxis_title="Unc(%)²",
        )
        trace_matrix_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_prof_covar = [np.sum(x) for x in self.cov_mat]
        trace_matrix_profiles_cov = plots.plot_profiles_per_iso_reac(
            vector=matrix_prof_covar,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix before assimilation (covariances included)",
            yaxis_title="Unc(%)² (covariances included)",
        )
        trace_matrix_profiles_cov.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        trace_cov_delta_integrals = plots.plot_matrix_integrals_per_iso_reac(
            cov_mat=-1 * self.cov_mat_delta,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            title="Variances-covariances matrix deltas (Cov_post - Cov_prior) after assimilation (group-wise integrals - covariances included)",
            yaxis_title="Variance & Covariance (group-wise integral - covariances included)",
        )
        trace_cov_delta_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_delta_diag = np.diag(-1 * self.cov_mat_delta.toarray())
        trace_matrix_delta_profiles = plots.plot_profiles_per_iso_reac(
            vector=matrix_delta_diag,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix deltas (Cov_post - Cov_prior) after assimilation (diagonal elements)",
            yaxis_title="Unc(%)²",
        )
        trace_matrix_delta_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_delta_prof_covar = [np.sum(x) for x in -1 * self.cov_mat_delta]
        trace_matrix_delta_profiles_cov = plots.plot_profiles_per_iso_reac(
            vector=matrix_delta_prof_covar,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix deltas (Cov_post - Cov_prior) after assimilation (covariances included)",
            yaxis_title="Unc(%)² (covariances included)",
        )
        trace_matrix_delta_profiles_cov.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        delta_mu_integrals = plots.plot_integrals_per_iso_reac(
            vector=self.delta_mu,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            title="Delta mu (group-wise integrals)",
            yaxis_title=f"Delta mu (group-wise integral) - Delta[ND]/ND",
        )
        delta_mu_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        delta_mu_profiles = plots.plot_profiles_per_iso_reac(
            vector=self.delta_mu, iso_reac_list=self.iso_reac_list, e_bins=self.e_bins, title="Delta mu", yaxis_title=f"Delta mu - Delta[ND]/ND"
        )
        delta_mu_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )
        # --------------------------------
        case_iso_reac = []
        if self.study_case != None:
            case_iso_reac = [(iso, reac) for iso, reac in self.study_case.iso_reac_list if reac != 1]
        benchs_iso_reac = [
            (iso, reac)
            for b in range(len(self.bench_cases))
            for iso, reac in self.bench_cases[b].iso_reac_list
            if reac != 1 and self.bench_list["REMOVED"][b] == False
        ]
        benchs_iso_reac = list(set(benchs_iso_reac))
        common_iso_reac = list(set(case_iso_reac).union(benchs_iso_reac))

        # Get iso_reac_list from cov_data
        cov_iso_reac = self.cov_data.iso_reac_list

        iso_str, reac_str, case_present, benchs_present, cov_present, calc_present = [], [], [], [], [], []
        for iso, reac in common_iso_reac:
            iso_str.append(methods.convert_iso_id_to_string(iso))
            reac_str.append(methods.reac_trad.get(str(reac), f"REAC_{reac}"))
            case_present.append(True if (iso, reac) in case_iso_reac else False)
            benchs_present.append(True if (iso, reac) in benchs_iso_reac else False)
            cov_present.append(
                "added*"
                if (iso, reac) in self.iso_reac_list and (iso, reac) not in cov_iso_reac
                else True if (iso, reac) in self.iso_reac_list else False
            )
            calc_present.append(True if (iso, reac) in self.iso_reac_list else False)

        iso_reac_df = {
            "Isotope": iso_str,
            "Reaction": reac_str,
            "Iso_ID": [iso for iso, reac in common_iso_reac],
            "Reac_ID": [reac for iso, reac in common_iso_reac],
            "Study_case": case_present,
            "Benchmark_cases": benchs_present,
            "Covariance_data": cov_present,
            "Used_in_calculation": calc_present,
        }
        iso_reac_df = pd.DataFrame(iso_reac_df)
        iso_reac_df.sort_values(by=["Isotope", "Reaction"], inplace=True)
        iso_reac_df.reset_index(drop=True, inplace=True)

        fill_colors = [
            "white",
            "white",
            "white",
            "white",
            ["green" if p else "red" for p in iso_reac_df["Study_case"]],
            ["green" if p else "red" for p in iso_reac_df["Benchmark_cases"]],
            ["orange" if p == "added*" else "green" if p == True else "red" for p in iso_reac_df["Covariance_data"]],
            ["green" if p else "red" for p in iso_reac_df["Used_in_calculation"]],
        ]

        table_iso_reac = go.Figure(
            data=[
                go.Table(
                    header=dict(values=iso_reac_df.columns, fill_color="#809684", align="center"),
                    cells=dict(values=iso_reac_df.values.T, fill_color=fill_colors, align="center"),
                )
            ]
        )
        table_iso_reac.update_layout(
            {"width": plot_width, "font_size": font_size, "paper_bgcolor": "rgba(255, 255, 255, 0.3)", "margin": dict(l=20, r=20, t=20, b=20)}
        )

        with open(output_html_path, "a", encoding="utf-8") as f:
            offline_control()
            f.writelines(HTML_intro)
            f.write(text_intro)
            f.write(table_res)

            if self.study_case != None:
                f.write(
                    plots.create_html_tabs(
                        [
                            "Benchmark list",
                            "Study case sensitivities",
                            "Uncertainty/Bias decomposition",
                            "Covariances matrix",
                            "Isotopes and reactions",
                        ]
                    )
                )
            else:
                f.write(plots.create_html_tabs(["Benchmark list", "Covariances matrix", "Isotopes and reactions"]))

            f.write('<div id="Benchmark list" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#d2978e, white) ; padding: 14px 100px">\n' + "<br>\n")
            f.write(trace_bias.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_sim.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(table_bench_fig)
            f.write("</section>\n")
            f.write("</div>\n")

            if self.study_case != None:
                f.write('<div id="Study case sensitivities" class="tabcontent">\n')
                f.write(f'<section style="background:linear-gradient(#a8b0af, white)">\n' + "<br>\n")
                f.write(trace_case_sensi_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(trace_case_sensi_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(trace_case_sensi_profiles_lethargy.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write("</section>\n")
                f.write("</div>\n")

                f.write('<div id="Uncertainty/Bias decomposition" class="tabcontent">\n')
                f.write(f'<section style="background:linear-gradient(#d2ebbb, white)">\n' + "<br>\n")
                f.write(dec_prior_trace_integ.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(dec_prior_trace.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(dec_post_trace_integ.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(dec_post_trace.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(dec_bias_trace_integ.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write(dec_bias_trace.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
                f.write("</section>\n")
                f.write("</div>\n")

            f.write('<div id="Covariances matrix" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#bcc1d4, white)">\n' + "<br>\n")
            f.write(trace_cov_prior_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_profiles_cov.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_cov_delta_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_delta_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_delta_profiles_cov.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(delta_mu_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(delta_mu_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")

            iso_reac_to_plot = []
            for isotope in isotopes_to_detail:

                if isinstance(isotope, str):
                    try:
                        isotope_id = methods.convert_iso_string_to_id(isotope)
                    except:
                        isotope_id = int(isotope)
                else:
                    isotope_id = int(isotope)

                iso_reac_to_plot.extend([iso_reac for iso_reac in self.iso_reac_list if iso_reac[0] == isotope_id])

            for iso_reac_1 in iso_reac_to_plot:
                for iso_reac_2 in iso_reac_to_plot:

                    write_div = False
                    trace_cov_prior_submat, color_scale = plots.plot_submatrix(
                        cov_mat=self.cov_mat,
                        iso_reac_list=self.iso_reac_list,
                        iso_reac_pair_to_plot=(iso_reac_1, iso_reac_2),
                        group_nb=self.group_nb,
                        title="Prior var-covar submatrix",
                    )

                    trace_cov_delta_submat, color_scale = plots.plot_submatrix(
                        cov_mat=self.cov_mat_delta,
                        iso_reac_list=self.iso_reac_list,
                        iso_reac_pair_to_plot=(iso_reac_1, iso_reac_2),
                        group_nb=self.group_nb,
                        color_scale=color_scale,
                        title="Var-covar delta submatrix",
                    )

                    if trace_cov_prior_submat != None:
                        write_div = True
                        trace_cov_prior_submat.update_layout(
                            {
                                "height": plot_height,
                                "width": plot_width,
                                "font_size": font_size,
                                "template": plotly_template,
                                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                            }
                        )
                        html_plot = trace_cov_prior_submat.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include)
                        f.write(f'<section style="border: 4px solid #4B0082; border-radius: 25px ; margin:{block_margin}; padding:10">\n')
                        f.write(html_plot + "<br>\n")

                    if trace_cov_delta_submat != None:
                        trace_cov_delta_submat.update_layout(
                            {
                                "height": plot_height,
                                "width": plot_width,
                                "font_size": font_size,
                                "template": plotly_template,
                                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                            }
                        )

                        html_plot = trace_cov_delta_submat.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include)
                        if not write_div:
                            f.write(f'<section style="border: 4px solid #4B0082; border-radius: 25px ; margin:{block_margin}; padding:10">\n')
                        f.write(html_plot + "<br>\n")

                    if trace_cov_delta_submat != None or trace_cov_prior_submat != None:
                        f.write("</section>\n")
            f.write("</section>\n")
            f.write("</div>\n")

            f.write('<div id="Isotopes and reactions" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#809684, white)">\n' + "<br>\n")
            f.write(f'<h2 style="text-align: center;">Isotopes and reactions in all cases :</h2>')
            f.write(table_iso_reac.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            if "added*" in iso_reac_df["Covariance_data"].values:
                f.write(
                    f"<h4>*added : this isotope-reaction is not present in the covariances data but its natural (int(iso_ID / 1000) * 1000) or rebased state (int(iso / 1000) * 1000 + int(str(iso_ID)[-2:])) was used instead.</h4>"
                )
            f.write("</section>\n")
            f.write("</div>\n")

            f.writelines(HTML_end)


class Uncertainty:
    """
    Class to calculate the uncertainty of a response based on the covariance matrix and sensitivity vector.
    The covariances matrix and sensitivity vectors are previously constructed and aligned when calling the function calcul_uncertainty().
    This class is the output of the function calcul_uncertainty().
    It contains the uncertainty value and the decomposition of the uncertainty into contributions from isotopes and reactions.

    Attributes
    -----------
    value : float
        The calculated uncertainty value in pcm.
    decomposition : pd.DataFrame
        A DataFrame containing the decomposition of the uncertainty into contributions from isotopes and reactions.
    study_case : Case
        The study case object containing the response and sensitivity data.
    cov_data : NDCovariances or Assimilation
        Covariance data object (NDCovariances or Assimilation object).
    __sensi_vec : np.ndarray
        The sensitivity vector.
    resp_calc : float
        The calculated response of the study case.
    __cov_mat : np.ndarray
        The covariance matrix.
    decomp_vec : list
        List of contributions to the uncertainty from every isotopes, reactions and energy group (last step of matrix multiplication, before summing and square rooting).
    iso_reac_list : list
        List of isotopes and reactions taken into account.
    iso_list : list
        List of isotopes taken into account.
    reac_list : list
        List of reactions taken into account.
    e_bins : list
        Energy bins.
    group_nb : int
        Number of energy groups.
    output_html_path : str, optional
        Path to save the output HTML file.

    Methods
    --------
    export_to_html(output_html_path, plotting_unit="pcm", isotopes_to_detail=[])
        Exports the results of the uncertainty calculation to an HTML file with interactive plots.
    """

    def __init__(
        self, study_case, cov_data, sensi_vec, resp_calc, cov_mat, iso_reac_list, output_html_path=None, plotting_unit="pcm", isotopes_to_detail=[]
    ) -> None:
        """
        Initialize the Uncertainty class.

        Parameters
        -----------
        study_case : Case
            The study case object containing the response and sensitivity data.
        cov_data : NDCovariances or Assimilation
            Covariance data as NDCovariances object or Assimilation object.
        sensi_vec : np.ndarray
            The sensitivity vector.
        resp_calc : float
            The calculated response of the study case.
        cov_mat : np.ndarray
            The covariance matrix.
        iso_reac_list : list
            List of isotopes and reactions taken into account.
        output_html_path : str, optional
            Path to save the output HTML file.
        plotting_unit : str, optional
            Unit for the sensitivity coefficients and uncertainty/bias contributions when plotted. Can be "relative" or "pcm". Defaults to "relative".
        isotopes_to_detail : list, optional
            List of isotopes (ID) to provide detailed variances-covariances (as heatmaps) in the HTML outputfile.
        """
        self.study_case = study_case
        self.cov_data = cov_data
        self.iso_reac_list = iso_reac_list
        self.iso_list = list(set([iso for iso, reac in iso_reac_list]))
        self.reac_list = list(set([reac for iso, reac in iso_reac_list]))
        self.resp_calc = resp_calc
        self.e_bins = study_case.e_bins
        self.group_nb = self.study_case.group_nb
        self.output_html_path = output_html_path

        dikt = {
            "ISO": [],
            "REAC": [],
            "ISO_NAME": [],
            "REAC_NAME": [],
            "CASE SENSIB INTEGRAL": [],
            "CONTRIBUTION TO RELATIVE UNC_SQUARED (COVAR WITH OTHER ISO-REAC INCLUDED)": [],
            "CONTRIBUTION INTEGRAL TO RELATIVE UNC SQUARED (COVAR WITH OTHER ISO-REAC INCLUDED)": [],
            "CONTRIBUTION INTEGRAL TO RELATIVE UNC SQUARED (FROM AUTO-CORRELATION ONLY)": [],
        }
        sum = 0
        decomp_vec = []
        for i, (iso, reac) in enumerate(iso_reac_list):

            sub_cov_mat = cov_mat[:, i * self.group_nb : (i + 1) * self.group_nb]

            sub_sensi_vec = sensi_vec[i * self.group_nb : (i + 1) * self.group_nb]

            unc_partial_covar = sensi_vec @ sub_cov_mat
            unc_partial_covar_detail = [x * y for x, y in zip(unc_partial_covar, sub_sensi_vec)]
            decomp_vec += unc_partial_covar_detail
            unc_partial_covar = unc_partial_covar @ sub_sensi_vec
            unc_partial_var = sub_sensi_vec @ sub_cov_mat[i * self.group_nb : (i + 1) * self.group_nb, :] @ sub_sensi_vec

            sum += unc_partial_covar

            sensi_df = study_case.sensitivities

            row = sensi_df[(sensi_df["ISO"] == iso) & (sensi_df["REAC"] == reac)]
            if len(row) != 0:
                integral = list(row["SENSI_INTEGRAL"])[0]
            else:
                integral = None

            dikt["ISO"].append(iso)
            dikt["REAC"].append(reac)
            dikt["ISO_NAME"].append(methods.convert_iso_id_to_string(iso))
            dikt["REAC_NAME"].append(methods.reac_trad.get(str(reac), f"REAC_{reac}"))
            dikt["CASE SENSIB INTEGRAL"].append(integral)
            dikt["CONTRIBUTION TO RELATIVE UNC_SQUARED (COVAR WITH OTHER ISO-REAC INCLUDED)"].append(unc_partial_covar_detail)
            dikt["CONTRIBUTION INTEGRAL TO RELATIVE UNC SQUARED (COVAR WITH OTHER ISO-REAC INCLUDED)"].append(unc_partial_covar)
            dikt["CONTRIBUTION INTEGRAL TO RELATIVE UNC SQUARED (FROM AUTO-CORRELATION ONLY)"].append(unc_partial_var)

        self.decomp_vec = decomp_vec
        self.decomposition = pd.DataFrame(dikt)

        self.value = sqrt(abs(sum)) * 1e5 * resp_calc

        self.__sensi_vec = sensi_vec
        self.__cov_mat = cov_mat

        if output_html_path is not None:
            self.export_to_html(output_html_path=self.output_html_path, plotting_unit=plotting_unit, isotopes_to_detail=isotopes_to_detail)

    @log_exec()
    def export_to_html(self, output_html_path: str, plotting_unit="pcm", isotopes_to_detail=[]):
        """
        Export the results of the uncertainty calculation to an HTML file.
        It contains the uncertainty value and the decomposition of the uncertainty into contributions from isotopes and reactions.

        Parameters
        -----------
        output_html_path : str
            [Required] Path to save the output HTML file.
        plotting_unit : str, optional
            Unit for the sensitivity coefficients and uncertainty/bias contributions when plotted. Can be "relative" or "pcm". Defaults to "relative".
        isotopes_to_detail : list, optional
            List of isotopes (ID) to provide detailed variances-covariances (as heatmaps) in the HTML outputfile.
        """
        if output_html_path is not None:
            if not output_html_path.endswith(".html"):
                output_html_path = output_html_path + ".html"
            with open(output_html_path, "w") as f:
                None
        else:
            return None

        if not isinstance(isotopes_to_detail, list):
            isotopes_to_detail = [isotopes_to_detail]

        # --------------------------------
        text_intro = ""
        text_intro += f"<h1><center>{pkg_fullname}</center></h1>"
        text_intro += f"<h1>Sandwich formula results</h1>"
        text_intro += f"<h2>Study case : {self.study_case.casename}</h2>"

        # --------------------------------
        headers = [
            "Calculated response of study case",
            "Uncertainty a priori",
            "Energy groups number",
        ]
        headers = ["<b>" + h + ":</b>" for h in headers]

        results = [
            f"{self.resp_calc} +/- {self.study_case.sigma_resp_calc}",
            f"{int(round(self.value))} pcm",
            self.group_nb,
        ]

        table_res = plots.create_html_table(headers=[], lines=[headers, results])

        # --------------------------------
        if plotting_unit == "relative":
            val_factor = 1.0
            unit_str = "%[resp]"
        elif plotting_unit == "pcm":
            val_factor = self.resp_calc * 1e5
            unit_str = "pcm"

        trace_case_sensi_integrals = plots.plot_integrals_per_iso_reac(
            vector=self.__sensi_vec,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            factor=val_factor,
            title="Study case sensitivities (integrals)",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND] - group-wise integral)",
        )
        trace_case_sensi_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        trace_case_sensi_profiles = plots.plot_profiles_per_iso_reac(
            vector=self.__sensi_vec,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            factor=val_factor,
            title="Study case sensitivities",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND])",
        )
        trace_case_sensi_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        case_vec_per_unit_lethargy = []
        for iso, reac in self.iso_reac_list:
            try:
                case_vec_per_unit_lethargy += self.study_case.get_normalized_lethargy_sensitivity(iso, reac)
            except errors.MissingDataError:
                case_vec_per_unit_lethargy += [0 for i in range(self.group_nb)]
        case_vec_per_unit_lethargy = np.array(case_vec_per_unit_lethargy)

        trace_case_sensi_profiles_lethargy = plots.plot_profiles_per_iso_reac(
            vector=case_vec_per_unit_lethargy,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            factor=val_factor,
            title="Sensitivities per unit lethargy",
            yaxis_title=f"Sensitivity ({unit_str}/%[ND] per unit lethargy)",
        )
        trace_case_sensi_profiles_lethargy.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        # --------------------------------
        dec_trace_integ = plots.plot_integrals_per_iso_reac(
            vector=self.decomp_vec,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            factor=val_factor**2,
            title="Integral contribution to relative uncertainy² (covariances included)",
            yaxis_title=f"Unc² ({unit_str}²) (covariances included)",
        )

        dec_trace_integ.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        dec_trace = plots.plot_profiles_per_iso_reac(
            vector=self.decomp_vec,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            factor=val_factor**2,
            title="Contribution to relative uncertainy² (covariances included)",
            yaxis_title=f"Unc² ({unit_str}²) (covariances included)",
        )
        dec_trace.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )
        # --------------------------------
        trace_cov_integrals = plots.plot_matrix_integrals_per_iso_reac(
            cov_mat=self.__cov_mat,
            iso_reac_list=self.iso_reac_list,
            group_nb=self.group_nb,
            title="Variances-covariances matrix (group-wise integrals - covariances included)",
            yaxis_title="Unc(%)² (group-wise integral - covariances included)",
        )
        trace_cov_integrals.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_diag = np.diag(self.__cov_mat.toarray())
        trace_matrix_profiles = plots.plot_profiles_per_iso_reac(
            vector=matrix_diag,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix (diagonal elements)",
            yaxis_title="Unc(%)²",
        )
        trace_matrix_profiles.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )

        matrix_prof_covar = [np.sum(x) for x in self.__cov_mat]
        trace_matrix_profiles_cov = plots.plot_profiles_per_iso_reac(
            vector=matrix_prof_covar,
            iso_reac_list=self.iso_reac_list,
            e_bins=self.e_bins,
            title="Variances-covariances matrix (covariances included)",
            yaxis_title="Unc(%)² (covariances included)",
        )
        trace_matrix_profiles_cov.update_layout(
            {
                "height": plot_height,
                "width": plot_width,
                "font_size": font_size,
                "template": plotly_template,
                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
            }
        )
        # --------------------------------
        case_iso = [iso for iso, reac in self.study_case.iso_reac_list if reac != 1]
        case_reac = [reac for iso, reac in self.study_case.iso_reac_list if reac != 1]
        case_iso_str = [methods.convert_iso_id_to_string(iso) for iso in case_iso]
        case_reac_str = [methods.reac_trad.get(str(reac), f"REAC_{reac}") for reac in case_reac]

        # Get iso_reac_list from cov_data
        if isinstance(self.cov_data, NDCovariances):
            cov_iso_reac_list = self.cov_data.iso_reac_list
        elif isinstance(self.cov_data, Assimilation):
            cov_iso_reac_list = self.cov_data.cov_data.iso_reac_list

        cov_present = [True if (iso, reac) in cov_iso_reac_list else False for iso, reac in zip(case_iso, case_reac)]
        calc_present = [True if (iso, reac) in self.iso_reac_list else False for iso, reac in zip(case_iso, case_reac)]

        headers = [["Isotope"], ["Reaction"], ["Iso_ID"], ["Reac_ID"], ["Covariance_data"], ["Used_in_calculation"]]

        results = [
            case_iso_str,
            case_reac_str,
            case_iso,
            case_reac,
            cov_present,
            calc_present,
        ]

        fill_colors = [
            "white",
            "white",
            "white",
            "white",
            ["green" if p else "red" for p in cov_present],
            ["green" if p else "red" for p in calc_present],
        ]

        table_iso_reac = go.Figure(
            data=[
                go.Table(
                    header=dict(values=headers, fill_color="#809684", align="center"),
                    cells=dict(values=results, fill_color=fill_colors, align="center"),
                )
            ]
        )
        table_iso_reac.update_layout({"font_size": font_size, "paper_bgcolor": "rgba(255, 255, 255, 0.3)", "margin": dict(l=20, r=20, t=20, b=20)})

        with open(output_html_path, "a", encoding="utf-8") as f:

            offline_control()
            f.writelines(HTML_intro)
            f.write(text_intro)
            f.write(table_res)
            f.write(plots.create_html_tabs(["Study case sensitivities", "Uncertainty decomposition", "Covariances matrix", "Isotopes and reactions"]))

            f.write('<div id="Study case sensitivities" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#a8b0af, white)">\n' + "<br>\n")
            f.write(trace_case_sensi_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_case_sensi_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "\n<br>\n")
            f.write(trace_case_sensi_profiles_lethargy.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + " \n")
            f.write("</section>\n")
            f.write("</div>\n")

            f.write('<div id="Uncertainty decomposition" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#d2ebbb, white)">\n' + "<br>\n")
            f.write(dec_trace_integ.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(dec_trace.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write("</section>\n")
            f.write("</div>\n")

            f.write('<div id="Covariances matrix" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#bcc1d4, white)">\n' + "<br>\n")
            f.write(trace_cov_integrals.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_profiles.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write(trace_matrix_profiles_cov.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")

            iso_reac_to_plot = []
            for isotope in isotopes_to_detail:

                try:
                    isotope_id = methods.convert_iso_string_to_id(isotope)
                except:
                    isotope_id = int(isotope)

                iso_reac_to_plot.extend([iso_reac for iso_reac in self.iso_reac_list if iso_reac[0] == isotope_id])

            for iso_reac_1 in iso_reac_to_plot:
                for iso_reac_2 in iso_reac_to_plot:

                    write_div = False
                    trace_cov_submat, color_scale = plots.plot_submatrix(
                        cov_mat=self.__cov_mat,
                        iso_reac_list=self.iso_reac_list,
                        iso_reac_pair_to_plot=(iso_reac_1, iso_reac_2),
                        group_nb=self.group_nb,
                        title="Prior var-covar submatrix",
                    )

                    if trace_cov_submat != None:
                        trace_cov_submat.update_layout(
                            {
                                "height": plot_height,
                                "width": plot_width,
                                "font_size": font_size,
                                "template": plotly_template,
                                "paper_bgcolor": "rgba(255, 255, 255, 0.8)",
                            }
                        )
                        f.write(trace_cov_submat.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")

            f.write("</section>\n")
            f.write("</div>\n")

            f.write('<div id="Isotopes and reactions" class="tabcontent">\n')
            f.write(f'<section style="background:linear-gradient(#809684, white)">\n' + "<br>\n")
            f.write(f'<h2 style="text-align: center;">Isotopes and reactions in study case :</h2>')
            f.write(table_iso_reac.to_html(full_html=False, include_plotlyjs=plotlyjs_fig_include) + "<br>\n")
            f.write("</section>\n")
            f.write("</div>\n")

            f.writelines(HTML_end)


class Bias:
    """
    Class to calculate the bias of a response based on the sensitivity vector and delta m (global nuclear data variation due to assimilation).
    The delta mu and sensitivity vectors are previously constructed and aligned when calling the function calcul_bias() (inside the class Assimilation).
    This class is the output of the function calcul_bias().
    It contains the bias value and the decomposition of the bias into contributions from isotopes and reactions.

    Attributes
    -----------
    value : float
        The calculated bias value in pcm.
    decomposition : pd.DataFrame
        DataFrame containing the decomposition of the bias into contributions from isotopes and reactions.
    decomp_vec : list
        List of contributions to the bias from every isotopes, reactions and energy group (last step of vectors multiplication, before summing).
    study_case : Case
        The study case object containing the response and sensitivity data.
    resp_calc : float
        The calculated response of the study case.
    iso_reac_list : list
        List of isotopes and reactions taken into account.
    iso_list : list
        List of isotopes taken into account.
    reac_list : list
        List of reactions taken into account.
    e_bins : list
        Energy bins.
    group_nb : int
        Number of energy groups.
    output_html_path : str, optional
        Path to save the output HTML file.
    """

    def __init__(self, study_case, sensi_vec, resp_calc, delta_mu, iso_reac_list) -> None:

        self.study_case = study_case
        self.iso_reac_list = iso_reac_list
        self.iso_list = list(set([iso for iso, reac in iso_reac_list]))
        self.reac_list = list(set([reac for iso, reac in iso_reac_list]))
        self.resp_calc = resp_calc
        self.e_bins = study_case.e_bins
        self.group_nb = self.study_case.group_nb

        dikt = {
            "ISO": [],
            "REAC": [],
            "ISO_NAME": [],
            "REAC_NAME": [],
            "CASE SENSIB INTEGRAL": [],
            "CONTRIBUTION TO RELATIVE BIAS": [],
            "CONTRIBUTION INTEGRAL TO RELATIVE BIAS": [],
        }
        sum = 0
        decomp_vec = []
        for i, (iso, reac) in enumerate(iso_reac_list):

            sub_delta_mu = delta_mu[i * self.group_nb : (i + 1) * self.group_nb]

            sub_sensi_vec = sensi_vec[i * self.group_nb : (i + 1) * self.group_nb]

            bias_partial_detail = [x * y for x, y in zip(sub_delta_mu, sub_sensi_vec)]
            decomp_vec += bias_partial_detail

            bias_partial = sub_sensi_vec @ sub_delta_mu

            sum += bias_partial

            sensi_df = study_case.sensitivities

            row = sensi_df[(sensi_df["ISO"] == iso) & (sensi_df["REAC"] == reac)]
            if len(row) != 0:
                integral = list(row["SENSI_INTEGRAL"])[0]
            else:
                integral = None

            if bias_partial != 0.0:

                dikt["ISO"].append(iso)
                dikt["REAC"].append(reac)
                dikt["ISO_NAME"].append(methods.convert_iso_id_to_string(iso))
                dikt["REAC_NAME"].append(methods.reac_trad.get(str(reac), f"REAC_{reac}"))
                dikt["CASE SENSIB INTEGRAL"].append(np.sum(sub_sensi_vec))
                dikt["CONTRIBUTION TO RELATIVE BIAS"].append(bias_partial_detail)
                dikt["CONTRIBUTION INTEGRAL TO RELATIVE BIAS"].append(bias_partial)

        self.decomp_vec = decomp_vec
        self.decomposition = pd.DataFrame(dikt)

        self.value = sum * 1e5 * resp_calc
