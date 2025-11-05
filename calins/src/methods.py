import os, re, struct, json
import numpy as np
import pandas as pd
from scipy.sparse import lil_matrix
from pathlib import Path
from math import ceil, sqrt
import matplotlib.pyplot as plt

from . import classes, errors
from logs import log_exec, warn, write_and_print

# ---Loading of id's for reactions and isotopes---
global reac_trad
global reac_trad_inv
global iso_z
global iso_z_inv

reac_trad = json.load(open(os.path.join(os.path.dirname(__file__), "isotopes_reactions/reac_id.json"), "r"))
reac_trad_inv = {val: key for key, val in reac_trad.items()}

iso_z = json.load(open(os.path.join(os.path.dirname(__file__), "isotopes_reactions/iso_id.json"), "r"))
iso_z_inv = {val: key for key, val in iso_z.items()}
# ------------------------------------------------


@log_exec()
def format_comac_to_dataframe(input_path: str, output_path: str = None):
    """
    Function to parse the COMAC (CEA) covariance matrix folder and extract the covariance values into a DataFrame, that can be passes to other functions.

    Parameters
    ----------
    input_path : str
        [Required] Path to the input covariance matrix folder.
    output_path : str, optional
        Path to store the DataFrame as an Excel file. The default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the covariance values. This DataFrame can be passed to other functions.
    """

    reac_trad_comac = {"2": "ELASTIC", "4": "INELASTIC", "16": "NXN", "18": "FISSION", "101": "CAPTURE", "452": "NU", "455": "NU_D"}
    reac_trad_comac_inv = {val: key for key, val in reac_trad_comac.items()}

    dikt_cov = {"ISO_H": [], "REAC_H": [], "ISO_V": [], "REAC_V": [], "STD": []}

    cov_root = os.path.join(input_path)

    for file in os.listdir(cov_root):

        with open(os.path.join(cov_root, file), "r") as f:
            lines = f.readlines()

        group_nb = len(lines[1].split())  # ! Important assumption !

        for i_line, line in enumerate(lines):

            if line.split()[0] in reac_trad_comac_inv and line.split()[2] in reac_trad_comac_inv:
                dikt_cov["REAC_V"].append(int(reac_trad_comac_inv[line.split()[0]]))
                dikt_cov["REAC_H"].append(int(reac_trad_comac_inv[line.split()[2]]))

                try:
                    for idx, direction in enumerate(["V", "H"], start=1):
                        Z = re.findall("[A-Z]+", line.split()[idx * 2 - 1].upper())
                        A = re.findall("[0-9]+", line.split()[idx * 2 - 1])

                        if len(A) != 0:
                            iso_idx = iso_z[Z[0]] * 1000 + int(A[0])
                        else:
                            iso_idx = iso_z[Z[0]] * 1000

                        dikt_cov[f"ISO_{direction}"].append(iso_idx)

                except:
                    warn(f"Isotopes in line '{line}' could not be identified")
                    if len(dikt_cov["ISO_V"]) == len(dikt_cov["ISO_H"]) + 1:
                        dikt_cov["ISO_V"].pop(-1)
                    dikt_cov["REAC_V"].pop()
                    dikt_cov["REAC_H"].pop()
                    continue

                stds_v = [float(x) / 100 for x in lines[i_line + 1].split()]
                stds_h = [float(x) / 100 for x in lines[i_line + 2].split()]

                corr_vals = [[float(val) for val in x.split()] for x in lines[i_line + 3 : i_line + group_nb + 3]]

                std_vals = corr_vals[:]
                for line_idx, corr_row in enumerate(corr_vals):
                    for col_idx, val in enumerate(corr_row):

                        std_vals[line_idx][col_idx] = round(val * stds_h[col_idx] * stds_v[line_idx], 13)

                dikt_cov["STD"].append(std_vals)

    cov_df = pd.DataFrame(dikt_cov)

    if cov_df.empty:
        raise errors.EmptyParsingError(f"No data was extracted from your var-covar matrix file : {input_path}")

    if output_path != None:
        if output_path.endswith(".xlsx"):
            cov_df.to_excel(output_path)
        elif os.path.isdir(output_path):
            cov_df.to_excel(os.path.join(output_path, "COV_MAT_COMAC.xlsx"))
        else:
            cov_df.to_excel(output_path + ".xlsx")

    return cov_df


@log_exec()
def format_scale_corr_to_dataframe(input_path: str, output_path: str = None):
    """
    Function to parse the SCALE correlation matrix text file and extract the correlation values into a DataFrame, that can be passes to other functions.

    Parameters
    ----------
    input_path : str
        [Required] Path to the input correlation matrix text file.
    output_path : str, optional
        Path to store the DataFrame as an Excel file. The default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the correlation values. This DataFrame can be passed to other functions.
    """
    corr_root = os.path.join(input_path)

    with open(corr_root, "r") as f:
        lines_corr = f.readlines()

    dikt_corr = {"ISO_H": [], "REAC_H": [], "ISO_V": [], "REAC_V": [], "CORR": []}

    group_nb = int(lines_corr[1].split()[1])

    corr_vals_block_len = ceil(group_nb**2 / 6)
    iso_reac_block_dist = ceil((group_nb + 1) / 6) + 1

    for i, line in enumerate(lines_corr):
        if line.startswith(" 9d") and i != 1:
            iso_reac_line = lines_corr[i - iso_reac_block_dist]
            if iso_reac_line.split()[2] not in ["1", "101"] and iso_reac_line.split()[4] not in ["1", "101"]:
                dikt_corr["ISO_H"].append(int(iso_reac_line.split()[1]))
                dikt_corr["REAC_H"].append(int(iso_reac_line.split()[2]))
                dikt_corr["ISO_V"].append(int(iso_reac_line.split()[3]))
                dikt_corr["REAC_V"].append(int(iso_reac_line.split()[4]))

                corr_val = [float(x) for x in "".join(lines_corr[i : i + corr_vals_block_len]).split()[1:]]
                corr_val = np.reshape(corr_val, ((group_nb, group_nb)))
                dikt_corr["CORR"].append(corr_val.tolist())

    corr_df = pd.DataFrame(dikt_corr)

    if corr_df.empty:
        raise errors.EmptyParsingError(f"No data was extracted from your correlation matrix file : {input_path}")

    if output_path != None:
        if output_path.endswith(".xlsx"):
            corr_df.to_excel(output_path)
        elif os.path.isdir(output_path):
            corr_df.to_excel(os.path.join(output_path, "CORR_MAT_SCALE.xlsx"))
        else:
            corr_df.to_excel(output_path + ".xlsx")

    return corr_df


@log_exec()
def format_scale_binary_to_dataframe(input_path: str, group_nb: int, output_xlsx_path: str = None, output_txt_path: str = None):
    """
    Function to parse the SCALE covariance matrix binary file and extract the covariance values into a DataFrame, that can be passes to other functions.

    Parameters
    ----------
    input_path : str
        [Required] Path to the input covariance matrix text file.
    group_nb : int
        [Required] Number of energy groups in the covariance matrix.
    output_xlsx_path : str, optional
        Path to store the DataFrame as an Excel file. The default is None.
    output_txt_path : str, optional
        Path to store the covariance matrix as a text file. The default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the covariance values. This DataFrame can be passed to other functions.
    """

    if group_nb not in [44, 56]:
        raise errors.UserInputError(
            f"The group number {group_nb} of the SCALE matrix should be 44 or 56. It determines the byte order for int and float numbers. Use 44 for big endian and 56 for normal, if your matrix has a different group number."
        )

    big_endian = ""
    if group_nb == 44:
        big_endian = ">"

    single_word = 6  # Byte-size for words
    single_number = 4  # Byte-size for float and integer

    dikt_cov = {"ISO_H": [], "REAC_H": [], "ISO_V": [], "REAC_V": [], "STD": []}

    with open(input_path, "rb") as input_file:

        input_file.read(single_number)

        # Read the file identification
        hname = input_file.read(single_word).decode("ascii")
        huse = input_file.read(2 * single_word).decode("ascii")
        ivers = struct.unpack(f"{big_endian}i", input_file.read(single_number))[0]

        input_file.read(2 * single_number)

        # Read the file control
        (ngroup, nngrup, nggrup, ntype, nmmp, nmtrix, nholl) = struct.unpack(f"{big_endian}7i", input_file.read(7 * single_number))

        input_file.read(2 * single_number)

        # Read the file description
        desc = input_file.read(nholl * single_word).decode("ascii")

        input_file.read(2 * single_number)

        # Read the neutron group boundaries
        ebc = []
        for i in range((nngrup + 1)):

            ebc.append(struct.unpack(f"{big_endian}f", input_file.read(4))[0])

        ebc = ["{:.7E}".format(j) for j in ebc]

        # Read the material-reaction control
        iso_reac = {"matid": [], "mtid": [], "mwgt": [], "xs": [], "err": []}

        input_file.read(2 * single_number)

        for i in range(nmmp):
            matid, mtid, mwgt = struct.unpack(f"{big_endian}3i", input_file.read(3 * single_number))

            if matid == 1801:
                matid = 8001001
            if matid == 1802:
                matid = 8001002

            # if i != 0:
            iso_reac["matid"].append(matid)
            iso_reac["mtid"].append(mtid)
            iso_reac["mwgt"].append(mwgt)

        # Read the material-reaction cross sections and error files
        for i in range(nmmp):
            xs, err = [], []
            input_file.read(2 * single_number)

            for j in range(ngroup):
                err_micro = struct.unpack(f"{big_endian}f", input_file.read(single_number))[0]
                err.append(err_micro)
            for j in range(ngroup):
                xs_micro = struct.unpack(f"{big_endian}f", input_file.read(single_number))[0]
                xs.append(xs_micro)

            iso_reac["xs"].append(xs)
            iso_reac["err"].append(err)

        # iso_reac = pd.DataFrame(iso_reac)

        input_file.read(single_number)

        # Read the matrix control
        for i, _ in enumerate(range(nmtrix)):

            input_file.read(single_number)

            mat1, mt1, mat2, mt2, nblock = struct.unpack(f"{big_endian}5i", input_file.read(5 * single_number))

            if mat1 == 1801:
                mat1 = 8001001
            if mat1 == 1802:
                mat1 = 8001002
            if mat2 == 1801:
                mat2 = 8001001
            if mat2 == 1802:
                mat2 = 8001002

            jband, ijj, lgpr = [], [], []
            for j in range(ngroup):
                jband_micro, ijj_micro = struct.unpack(f"{big_endian}2i", input_file.read(2 * single_number))
                jband.append(jband_micro)
                ijj.append(ijj_micro)

            for j in range(nblock):
                lgpr_micro = struct.unpack(f"{big_endian}i", input_file.read(single_number))[0]
                lgpr.append(lgpr_micro)

            input_file.read(2 * single_number)

            cov = np.zeros((nngrup, nngrup))
            for j in range(ngroup):

                i = 1 + j - ijj[j]
                while i < jband[j] + 1 + j - ijj[j]:
                    cov[j][i] = struct.unpack(f"{big_endian}f", input_file.read(single_number))[0]
                    i += 1

            input_file.read(single_number)

            if np.sum(cov) == 0.0:
                continue

            dikt_cov["ISO_H"].append(mat1)
            dikt_cov["REAC_H"].append(mt1)
            dikt_cov["ISO_V"].append(mat2)
            dikt_cov["REAC_V"].append(mt2)
            dikt_cov["STD"].append(cov.tolist())

        if output_txt_path != None:

            with open(output_txt_path, "w") as output_file:
                output_file.write(f"{hname} {huse} {ivers} {desc}\n")
                output_file.write(f"{ngroup: >12}{ntype: >12}{nmtrix: >12}{nmmp: >12}\n")

                for idx, energy in enumerate(ebc):
                    if idx > 0 and idx % 5 == 0:
                        output_file.write("\n")
                    output_file.write(f"{energy: >15}")
                output_file.write("\n")

                for idx, iso in enumerate(dikt_cov["ISO_H"]):

                    output_file.write(
                        f"{dikt_cov['ISO_H'][idx]: >12}{dikt_cov['REAC_H'][idx]: >12}{dikt_cov['ISO_V'][idx]: >12}{dikt_cov['REAC_V'][idx]: >12}{'1': >12}\n"
                    )

                    for std_x in dikt_cov["STD"][idx]:
                        for i, std in enumerate(std_x):
                            if i > 0 and i % 5 == 0:
                                output_file.write("\n")
                            std = "{:.7E}".format(std)
                            output_file.write(f"{std: >15}")
                        output_file.write("\n")

    cov_df = pd.DataFrame(dikt_cov)

    if cov_df.empty:
        raise errors.EmptyParsingError(f"No data was extracted from your var-covar matrix file : {input_path}")

    if output_xlsx_path != None:
        if output_xlsx_path.endswith(".xlsx"):
            cov_df.to_excel(output_xlsx_path)
        elif os.path.isdir(output_xlsx_path):
            cov_df.to_excel(os.path.join(output_xlsx_path, os.path.basename(input_path) + ".xlsx"))
        else:
            cov_df.to_excel(output_xlsx_path + ".xlsx")

    return cov_df


@log_exec()
def format_scale_txt_to_dataframe(input_path: str, output_path: str = None):
    """
    Function to parse the SCALE covariance matrix text file and extract the covariance values into a DataFrame, that can be passes to other functions.

    Parameters
    ----------
    input_path : str
        [Required] Path to the input covariance matrix text file.
    output_path : str, optional
        Path to store the DataFrame as an Excel file. The default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the covariance values. This DataFrame can be passed to other functions.
    """
    cov_root = os.path.join(input_path)

    with open(cov_root, "r") as f:
        lines_cov = f.readlines()

    dikt_cov = {"ISO_H": [], "REAC_H": [], "ISO_V": [], "REAC_V": [], "STD": []}

    group_nb = int(lines_cov[1].split()[0])

    block_len = ceil(group_nb / 5)

    for i, line in enumerate(lines_cov):
        # if line.startswith('     ') and i != 1 :
        if i > 2 and len(line.split()) == 5:
            cov_val = []
            if line.split()[1] != "1" and line.split()[3] != "1" and line.split()[1].upper() in reac_trad and line.split()[3].upper() in reac_trad:
                dikt_cov["ISO_H"].append(int(line.split()[0]))
                dikt_cov["REAC_H"].append(int(line.split()[1]))
                dikt_cov["ISO_V"].append(int(line.split()[2]))
                dikt_cov["REAC_V"].append(int(line.split()[3]))

                std_nb = 0
                j = 0
                while std_nb < group_nb * group_nb:
                    for val in lines_cov[i + 1 + j].split():
                        cov_val.append(float(val))
                    std_nb += len(lines_cov[i + 1 + j].split())

                    j += 1

                dikt_cov["STD"].append(np.reshape(cov_val, (group_nb, group_nb)).tolist())

    cov_df = pd.DataFrame(dikt_cov)

    if cov_df.empty:
        raise errors.EmptyParsingError(f"No data was extracted from your var-covar matrix file : {input_path}")

    if output_path != None:
        if output_path.endswith(".xlsx"):
            cov_df.to_excel(output_path)
        elif os.path.isdir(output_path):
            cov_df.to_excel(os.path.join(output_path, "COV_MAT_SCALE.xlsx"))
        else:
            cov_df.to_excel(output_path + ".xlsx")

    return cov_df


def condense_binning(values, i_ebins, o_ebins, std=False):
    sorted_input, sorted_input_reverse = sorted(i_ebins[:], reverse=False), sorted(i_ebins[:], reverse=True)
    if sorted_input == i_ebins:
        reverse_sort = False
    elif sorted_input_reverse == i_ebins:
        reverse_sort = True
    else:
        raise ValueError("The energy binning of the sdf file is not sorted")

    o_ebins = sorted(o_ebins, reverse=reverse_sort)

    output_values = []

    boundaries_idx = []
    for o_ebin in o_ebins:
        try:
            boundaries_idx.append(i_ebins.index(o_ebin))
        except:
            found = False
            for i, ebin in enumerate(i_ebins):
                if abs(o_ebin - ebin) / o_ebin < 0.0001:
                    boundaries_idx.append(i)
                    found = True
                    break
            if not found:
                raise errors.UserInputError(
                    f"The two energy binings are not compatibles for compression\n\
                    Output ebin {o_ebin} not in input ebins {i_ebins}"
                )
    for k in range(0, len(boundaries_idx) - 1):
        if std == False:
            output_values.append(sum(values[boundaries_idx[k] : boundaries_idx[k + 1]]))
        else:
            output_values.append(np.sqrt(sum(np.array(values[boundaries_idx[k] : boundaries_idx[k + 1]]) ** 2)))

    return output_values


@log_exec()
def condense_sdf(input_sdf_path: str, output_ebins: list, output_sdf_path: str):
    """
    Function to condense the SDF file from its energy binning to a wider compatible energy binning, and export it as a new SDF file.

    Parameters
    ----------
    input_sdf_path : str
        Path to the input SDF file.
    output_ebins : list
        List of energy bins for the output SDF file. Not necessarily sorted.
    output_sdf_path : str
        Path to write the new SDF file.
    """
    with open(input_sdf_path, "r") as f:
        lines_sensi = f.readlines()
        if output_sdf_path != None:
            lines_output = lines_sensi[:]

    group_nb = int(lines_sensi[1].split()[0])

    if lines_sensi[3].split()[1] == "+/-":
        formatting = "TSUNAMI-B"
        headers_length = 4
        headers_line_nb = 4
        block_length = headers_line_nb - 1 + ceil(group_nb / 5) * 2
    else:
        formatting = "TSUNAMI-A"
        headers_length = 6
        headers_line_nb = 2
        block_length = headers_line_nb - 1 + ceil(group_nb / 5)

    for i, line in enumerate(lines_sensi):
        if line.startswith("energy boundaries"):
            energy_bins = [float(x) for x in "".join(lines_sensi[i + 1 : i + 1 + ceil((group_nb + 1) / 5)]).split()]
            energy_bins = sorted(energy_bins, reverse=True)

            ngroup = len(output_ebins) - 1

            lines_output[1] = f"{ngroup: >10} number of neutron groups\n"

            text = ""
            for idx, energy in enumerate(output_ebins):
                energy = f"{energy:.6E}"
                if idx > 0 and idx % 5 == 0:
                    text += "\n"
                text += f"{energy: >14}"
            text += "\n"

            for idx in range(ceil((group_nb + 1) / 5)):
                if idx == 0:
                    lines_output[i + 1 + idx] = text
                else:
                    lines_output[i + 1 + idx] = ""

            break

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

    def detect_header(i):

        if not (
            len(lines_sensi[i].split()) == headers_length
            and i > (4 + ceil((group_nb + 1) / 5))
            and i < len(lines_sensi) - block_length
            and lines_sensi[i].split()[3].replace('-', '').isdigit()
        ):
            return False

        if formatting == "TSUNAMI-B":
            if not len(lines_sensi[i + 1].split()) == 2 and len(lines_sensi[i + 2].split()) == 4:
                return False
        elif formatting == "TSUNAMI-A":
            if not len(lines_sensi[i + 1]) == 3:
                return False

        return True

    pass_line = 0

    for i, line in enumerate(lines_sensi):

        if pass_line != 0:
            pass_line -= 1
            continue

        if detect_header(i):

            pass_line = block_length

            sensis = [float(x) for x in "".join(lines_sensi[i + headers_line_nb : i + headers_line_nb + ceil(group_nb / 5)]).split()]
            stds = [
                float(x)
                for x in "".join(lines_sensi[i + headers_line_nb + ceil(group_nb / 5) : i + headers_line_nb + ceil(group_nb / 5) * 2]).split()
            ]

            sensis = condense_binning(sensis, i_ebins=energy_bins, o_ebins=output_ebins)
            stds = condense_binning(stds, i_ebins=energy_bins, o_ebins=output_ebins)

            sensis, stds = [f"{sensi:.6E}" for sensi in sensis], [f"{std:.6E}" for std in stds]

            text = ""
            for idx, sensi in enumerate(sensis):
                if idx > 0 and idx % 5 == 0:
                    text += "\n"
                text += f"{sensi: >14}"
            text += "\n"

            for idx, std in enumerate(stds):
                if idx > 0 and idx % 5 == 0:
                    text += "\n"
                text += f"{std: >14}"
            text += "\n"

            for idx in range(ceil(group_nb / 5) * 2):
                if idx == 0:
                    lines_output[i + headers_line_nb + idx] = text
                else:
                    lines_output[i + headers_line_nb + idx] = ""

    with open(output_sdf_path, "w") as f:
        f.writelines(lines_output)

    return None


@log_exec()
def format_sensi_to_dataframe(
    input_sdf_path: str, output_path: str = None, filter_min_ratio_sensi_abs=0.0, occurrences_rule="sum", mcnp=False, std=False
):
    """
    SDF PARSER
    Function to parse the SDF file from TSUNAMI A or B formats and extract the sensitivities (and std) of the isotopes and reactions.

    Parameters
    ----------
    input_sdf_path : str
        Path to the input SDF file.
    output_path : str, optional
        Path to store the DataFrame as an Excel file. The default is None.
    filter_min_ratio_sensi_abs : float, optional
        Minimum ratio of the absolute sensitivity integral (for one isotope-reaction pair) to the sum of all absolute sensitivities. The default is 0.0.
    occurences_rule : str, optional
        Rule for handling occurrences of same isotope-reaction in one SDf file. Either "first", "last" or "sum", to chosse which one to consider. Defaults to "sum".
    mcnp : bool, optional
        Flag to be activated when parsing MCNP format of SDF file.
    std : bool, optional
        Flag to be activated for parsing standard deviations of sensitivity coefficients.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the sensitivities (and standard deviations if asked). This DataFrame is included inside the class object Case.
    """

    dikt_sensi = {"ISO": [], "REAC": [], "SENSI": [], "SENSI_INTEGRAL": [], "STD_INTEGRAL": [], "SENSI_INTEGRAL_ABS": []}
    if std:
        dikt_sensi["SENSI_STD"] = []

    iso_reac_micro_region = []

    if occurrences_rule not in ["sum", "first", "last"]:
        errors.UserInputError(f"The occurrence_rule argument to dictate how to parse mutliple iso-reac must be 'sum'/'first'/'last'.")

    with open(input_sdf_path, "r") as f:
        lines_sensi = f.readlines()

    group_nb = int(lines_sensi[1].split()[0])

    if lines_sensi[3].split()[1] == "+/-":
        formatting = "TSUNAMI-B"
        headers_length = 4
        headers_line_nb = 4
        block_length = headers_line_nb - 1 + ceil(group_nb / 5) * 2
    else:
        formatting = "TSUNAMI-A"
        headers_length = 6
        headers_line_nb = 2
        block_length = headers_line_nb - 1 + ceil(group_nb / 5)

    def detect_header(i):

        if not (
            len(lines_sensi[i].split()) == headers_length
            and i > (4 + ceil((group_nb + 1) / 5))
            and i < len(lines_sensi) - block_length
            and lines_sensi[i].split()[3].replace('-', '').isdigit()
        ):
            return False

        if formatting == "TSUNAMI-B":
            if not len(lines_sensi[i + 1].split()) == 2 and len(lines_sensi[i + 2].split()) == 4:
                return False
        elif formatting == "TSUNAMI-A":
            if not len(lines_sensi[i + 1]) == 3:
                return False

        return True

    sensi_profiles_nb = int(lines_sensi[2].split()[0])
    count_sensi_profiles = 0
    pass_line = 0

    for i, line in enumerate(lines_sensi):

        if pass_line != 0:
            pass_line -= 1
            continue

        if detect_header(i):

            pass_line = block_length
            count_sensi_profiles += 1

            if (
                (formatting == "TSUNAMI-B" and lines_sensi[i + 1].split()[0] == "0") or (formatting == "TSUNAMI-A" and line.split()[-2] == "0")
            ) and line.split()[3] in reac_trad:

                # Let's read the block of sensitivities values in one go and then pass all the lines until the next block (with the help of pass_line)
                sensi = [float(x) for x in "".join(lines_sensi[i + headers_line_nb : i + headers_line_nb + ceil(group_nb / 5)]).split()]
                if std:
                    sensi_std = [
                        float(x)
                        for x in "".join(lines_sensi[i + headers_line_nb + ceil(group_nb / 5) : i + headers_line_nb + ceil(group_nb / 5) * 2]).split()
                    ]

                sensi_integ = float(lines_sensi[i + 3].split()[0])
                std_integ = float(lines_sensi[i + 3].split()[1])
                # sensi_integ_abs = float( lines_sensi[i+3].split()[2] )
                sensi_integ_abs = np.sum(np.abs(sensi))

                reac = int(line.split()[3])
                iso = int(line.split()[2])
                if mcnp:
                    try:
                        iso = int(line.split()[0].split(".")[0])
                    except:
                        if line.split()[0].split(".")[0] in ["h-ch2", "h-ch2"]:
                            iso = 1001
                        else:
                            raise errors.SensInputError(f"The isotope {line.split()[0]} wasn't able to be identified.")
                    if reac == -2:
                        reac = 101
                    elif reac < 0:
                        warn(
                            f"The sdf file from MCNP calculation has a the following reaction not taken into account : REAC ID {reac}, ISO ID {iso}."
                        )
                        break

                occurrence_idx = None
                for i in range(len(dikt_sensi["ISO"])):
                    if dikt_sensi["ISO"][i] == iso and dikt_sensi["REAC"][i] == reac:
                        occurrence_idx = i
                        break

                if occurrence_idx != None:
                    if occurrences_rule == "sum":
                        if sensi_integ_abs != 0.0:
                            dikt_sensi["SENSI"][occurrence_idx] = (np.array(dikt_sensi["SENSI"][occurrence_idx]) + np.array(sensi)).tolist()
                            dikt_sensi["SENSI_INTEGRAL"][occurrence_idx] = dikt_sensi["SENSI_INTEGRAL"][occurrence_idx] + sensi_integ
                            dikt_sensi["SENSI_INTEGRAL_ABS"][occurrence_idx] = np.sum(np.abs(dikt_sensi["SENSI"][occurrence_idx]))
                            if isinstance(dikt_sensi["STD_INTEGRAL"][occurrence_idx], float):
                                dikt_sensi["STD_INTEGRAL"][occurrence_idx] = [dikt_sensi["STD_INTEGRAL"][occurrence_idx], std_integ]
                                if std:
                                    dikt_sensi["SENSI_STD"][occurrence_idx] = [dikt_sensi["SENSI_STD"][occurrence_idx], sensi_std]
                            elif isinstance(dikt_sensi["STD_INTEGRAL"][occurrence_idx], list):
                                dikt_sensi["STD_INTEGRAL"][occurrence_idx].append(std_integ)
                                if std:
                                    dikt_sensi["SENSI_STD"][occurrence_idx].append(sensi_std)

                        warn(
                            f"The sensitivities for the isotope {line.split()[0]} / {iso} and the reaction {line.split()[1]} / {reac} is being sum up with its previous encounter.{' But the current values are zeros anyway' if sensi_integ_abs == 0.0 else ''}"
                        )

                    elif occurrences_rule == "last":
                        dikt_sensi["SENSI"][occurrence_idx] = sensi
                        dikt_sensi["SENSI_INTEGRAL"][occurrence_idx] = sensi_integ
                        if std:
                            dikt_sensi["SENSI_STD"][occurrence_idx] = sensi_std
                        dikt_sensi["STD_INTEGRAL"][occurrence_idx] = std_integ
                        dikt_sensi["SENSI_INTEGRAL_ABS"][occurrence_idx] = sensi_integ_abs
                        dikt_sensi["ISO"][occurrence_idx] = iso
                        dikt_sensi["REAC"][occurrence_idx] = reac

                        warn(
                            f"The previous encounter for the sensitivities of the isotope {line.split()[0]} / {iso} and the reaction {line.split()[1]} / {reac} is being ignored. New values have just been found (arg occurrences_rule:'last').{' But be carefull because its all zeros !' if sensi_integ_abs == 0.0 else ''}"
                        )

                    elif occurrences_rule == "first":
                        warn(
                            f"The sensitivities for the isotope {line.split()[0]} / {iso} and the reaction {line.split()[1]} / {reac} is being ignored. Its first encounter has already been stored (arg occurrences_rule:'first').{' But the current values are zeros anyway' if sensi_integ_abs == 0.0 else ''}"
                        )

                else:
                    dikt_sensi["SENSI"].append(sensi)
                    dikt_sensi["SENSI_INTEGRAL"].append(sensi_integ)
                    if std:
                        dikt_sensi["SENSI_STD"].append(sensi_std)
                    dikt_sensi["STD_INTEGRAL"].append(std_integ)
                    dikt_sensi["SENSI_INTEGRAL_ABS"].append(sensi_integ_abs)
                    dikt_sensi["ISO"].append(iso)
                    dikt_sensi["REAC"].append(reac)

                    if sensi_integ_abs == 0.0:
                        warn(
                            f"The sensitivities for the isotope {line.split()[0]} / {iso} and the reaction {line.split()[1]} / {reac} are all zeros."
                        )

            elif line.split()[3] in reac_trad:
                reac = int(line.split()[3])
                iso = int(line.split()[2])
                if mcnp:
                    try:
                        iso = int(line.split()[0].split(".")[0])
                    except:
                        if line.split()[0].split(".")[0] in ["h-ch2", "h-ch2"]:
                            iso = 1001
                        else:
                            raise errors.SensInputError(f"The isotope {line.split()[0]} wasn't able to be identified.")
                    if reac == -2:
                        reac = 101
                    elif reac < 0:
                        warn(
                            f"The sdf file from MCNP calculation has a the following reaction not taken into account : REAC ID {reac}, ISO ID {iso}."
                        )
                        break

                if (iso, reac) not in iso_reac_micro_region:
                    iso_reac_micro_region.append((iso, reac))

        if count_sensi_profiles == sensi_profiles_nb:
            break

    for iso_reac in iso_reac_micro_region:
        if iso_reac not in zip(dikt_sensi["ISO"], dikt_sensi["REAC"]):
            warn(
                f"The isotope-reaction {iso_reac[0]} - {iso_reac[1]} has sensitivity data only for sub-regions of the case and not for the sum of all regions."
            )

    for i in range(len(dikt_sensi["STD_INTEGRAL"])):
        if isinstance(dikt_sensi["STD_INTEGRAL"][i], list):
            sum_integ = 0
            sum_std = np.zeros(group_nb)
            for j, std_integ in enumerate(dikt_sensi["STD_INTEGRAL"][i]):
                sum_integ += std_integ**2
                if std:
                    sum_std += np.array(dikt_sensi["SENSI_STD"][i][j]) ** 2

            dikt_sensi["STD_INTEGRAL"][i] = sqrt(sum_integ)
            if std:
                dikt_sensi["SENSI_STD"][i] = np.sqrt(sum_std).tolist()

    sensi_df = pd.DataFrame(dikt_sensi)

    sensi_df = sensi_df[sensi_df["REAC"] != 0]

    sum_integ_abs = np.sum(sensi_df[sensi_df["REAC"] != 1]["SENSI_INTEGRAL_ABS"].to_list())

    if filter_min_ratio_sensi_abs != 0.0:
        sensi_df = sensi_df[sensi_df["SENSI_INTEGRAL_ABS"] > sum_integ_abs * filter_min_ratio_sensi_abs]

    sensi_df = sensi_df.sort_values(by=["ISO", "REAC"]).reset_index(drop=True)

    if sensi_df.empty:
        raise errors.EmptyParsingError(f"No data was extracted from your sdf file : {input_path}")

    if output_path != None:
        if output_path.endswith(".xlsx"):
            sensi_df.to_excel(output_path)
        elif os.path.isdir(output_path):
            xls_name = os.path.basename(input_path)[:-4] + ".sensi.xlsx"
            sensi_df.to_excel(os.path.join(output_path, xls_name))

    return sensi_df


# ----- Build the covariance matrix according to the 'iso_reac_list'
def make_cov_matrix(cov_dataf: pd.DataFrame, iso_reac_list: list):

    # ----- Loop over col_idx and line_idx coordinates (iso-reac_Horizontal and iso-reac_Vertical) and insert a submatrix (or its transpose) if the input covariance matrix contains values for this iso-reac pair
    group_nb = len(list(cov_dataf["STD"])[0])

    cov_inter_iso = []
    for idx in cov_dataf.index:
        if cov_dataf["ISO_H"][idx] != cov_dataf["ISO_V"][idx]:
            cov_inter_iso.append((cov_dataf["ISO_H"][idx], cov_dataf["ISO_V"][idx]))

    # Create new submatrices for versatile isotopes (e.g., H-poly / 1901) and isotopes not present but associated with their natural isotope
    iso_reac_cov_H, iso_reac_cov_V = (
        cov_dataf.apply(lambda x: (x["ISO_H"], x["REAC_H"]), axis=1).to_list(),
        cov_dataf.apply(lambda x: (x["ISO_V"], x["REAC_V"]), axis=1).to_list(),
    )
    iso_reac_cov = list(set(iso_reac_cov_H) | set(iso_reac_cov_V))
    iso_cov = [iso for iso, reac in iso_reac_cov]

    for iso in list(set([iso for iso, reac in iso_reac_list])):
        # Handle isotopes not present in the matrix, if their natural isotope is present
        natural_idx = int(iso / 1000) * 1000
        if iso not in iso_cov and natural_idx in iso_cov:
            for i, row in cov_dataf[(cov_dataf["ISO_H"] == natural_idx) | (cov_dataf["ISO_V"] == natural_idx)].iterrows():
                if row["ISO_H"] == natural_idx:
                    row["ISO_H"] = iso
                if row["ISO_V"] == natural_idx:
                    row["ISO_V"] = iso
                cov_dataf = pd.concat([cov_dataf, pd.DataFrame([row])], ignore_index=True)
                warn(
                    f"The isotope {iso} is not present in the var-covar matrix file. It will be replaced by its natural isotope {natural_idx} for the construction of the covariance matrix."
                )

        # Handle unconventional isotopes (e.g., H-poly / 1901) if their standard isotope (e.g., H-1 / 1001) is present
        rebased_idx = int(iso / 1000) * 1000 + int(str(iso)[-2:])
        if iso not in iso_cov and rebased_idx in iso_cov:
            for i, row in cov_dataf[(cov_dataf["ISO_H"] == rebased_idx) | (cov_dataf["ISO_V"] == rebased_idx)].iterrows():
                if row["ISO_H"] == rebased_idx:
                    row["ISO_H"] = iso
                if row["ISO_V"] == rebased_idx:
                    row["ISO_V"] = iso
                cov_dataf = pd.concat([cov_dataf, pd.DataFrame([row])], ignore_index=True)
                warn(
                    f"The unconventional isotope {iso} is not present in the var-covar matrix file. It will be replaced by its base isotope {rebased_idx} for the construction of the covariance matrix."
                )

    iso_reac_cov_H, iso_reac_cov_V = (
        cov_dataf.apply(lambda x: (x["ISO_H"], x["REAC_H"]), axis=1).to_list(),
        cov_dataf.apply(lambda x: (x["ISO_V"], x["REAC_V"]), axis=1).to_list(),
    )
    iso_reac_cov = list(set(iso_reac_cov_H) | set(iso_reac_cov_V))

    iso_reac_inter = list(set(iso_reac_cov) & set(iso_reac_list))
    if len(iso_reac_inter) == 0:
        raise errors.EmptyParsingError(
            f"No iso-reac found for var-covar matrix construction\n\
        Iso-reac searched for the construction : {iso_reac_list}"
        )

    iso_reac_notfound = list(set(iso_reac_list) - set(iso_reac_inter))
    iso_reac_notfound = [convert_iso_id_to_string(iso) + " " + str(reac) for iso, reac in iso_reac_notfound]
    if len(iso_reac_notfound) > 0:
        warn(
            f"The following iso-reac are not present in the var-covar matrix file - they are ignored for the construction of the covariance matrix : \n{list(iso_reac_notfound)}"
        )

    cov_dim = len(iso_reac_inter) * group_nb
    cov_mat = lil_matrix((cov_dim, cov_dim), dtype=float)

    for col_idx, (iso_H, reac_H) in enumerate(iso_reac_inter):
        for line_idx, (iso_V, reac_V) in enumerate(iso_reac_inter):

            # iso_H,reac_H,iso_V,reac_V = str(iso_H), str(reac_H), str(iso_V), str(reac_V)

            # ----- If the horizontal and vertical isotopes are not part of the present inter-isotope covariances, skip
            if iso_H != iso_V and (iso_H, iso_V) not in cov_inter_iso and (iso_V, iso_H) not in cov_inter_iso:
                continue

            row = cov_dataf[
                (cov_dataf["ISO_H"] == iso_H) & (cov_dataf["REAC_H"] == reac_H) & (cov_dataf["ISO_V"] == iso_V) & (cov_dataf["REAC_V"] == reac_V)
            ]
            if len(row) != 0:
                cov_mat[line_idx * group_nb : (line_idx + 1) * group_nb, col_idx * group_nb : (col_idx + 1) * group_nb] = np.array(
                    list(row["STD"])[0]
                )

            else:

                row = cov_dataf[
                    (cov_dataf["ISO_H"] == iso_V) & (cov_dataf["REAC_H"] == reac_V) & (cov_dataf["ISO_V"] == iso_H) & (cov_dataf["REAC_V"] == reac_H)
                ]
                if len(row) != 0:
                    cov_mat[line_idx * group_nb : (line_idx + 1) * group_nb, col_idx * group_nb : (col_idx + 1) * group_nb] = np.array(
                        list(row["STD"])[0]
                    ).T

    return cov_mat, iso_reac_inter


def get_common_iso_reac_list(
    cases_list: list, iso_reac_list: list = None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None, operation="union"
):
    default_reac_list = [int(x) for x in reac_trad.keys() if x != "1"]

    if reac_list == None:
        reac_list = default_reac_list

    if iso_list != None:
        for i, iso in enumerate(iso_list):
            if not isinstance(iso, int):
                iso_list[i] = convert_iso_string_to_id(iso)
    if exclude_iso != None:
        for i, iso in enumerate(exclude_iso):
            if not isinstance(iso, int):
                exclude_iso[i] = convert_iso_string_to_id(iso)

    if iso_reac_list == None:
        iso_reac_lists = []
        for sensi_case in cases_list:
            iso_reac_lists.append(sensi_case.iso_reac_list)

        if len(cases_list) == 1:
            common_list = sensi_case.iso_reac_list
        elif operation == "intersection":
            common_list = list(set(iso_reac_lists[0]).intersection(*iso_reac_lists))
        elif operation == "union":
            common_list = list(set(iso_reac_lists[0]).union(*iso_reac_lists))
        else:
            raise ValueError("The 'operation' argument can only be 'intersection' or 'union'")

        common_list = sorted(common_list)

        common_list = [(iso, reac) for (iso, reac) in common_list if reac in reac_list]

        if iso_list != None:
            common_list = [(iso, reac) for (iso, reac) in common_list if iso in iso_list]
        if exclude_iso != None:
            common_list = [(iso, reac) for (iso, reac) in common_list if iso not in exclude_iso]

    else:
        for i, (iso, reac) in enumerate(iso_reac_list):
            if not isinstance(iso, int):
                iso_reac_list[i] = (convert_iso_string_to_id(iso), reac)
            if not isinstance(reac, int):
                raise errors.UserInputError(f"The reaction '{reac}' should be an integer.")
        common_list = iso_reac_list

    if len(common_list) == 0:
        raise errors.EmptyParsingError(
            f"The iso-reac list in common between the sensitivities vectors is empty \n\
        Cases : {[case.casename for case in cases_list]} {reac_list}"
        )

    return common_list


def make_sensi_vectors(
    cases_list: list, iso_reac_list: list = None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None, operation="union"
):

    common_list = get_common_iso_reac_list(**locals())

    # Check dimensions :
    group_nb = None
    for case_i in cases_list:
        if group_nb == None:
            group_nb = case_i.group_nb

        if case_i.group_nb != group_nb:
            raise errors.DimError(
                f"The energy groups number of your sensitivities is not consistent between your case 1 {cases_list[0].casename} ({cases_list[0].group_nb} groups number), and your case 2 {case_i.casename} ({case_i.group_nb} groups number)"
            )

    # ----- Construction of sensitivity vectors by following the iso-reac list 'common_list'
    sensi_vecs_output = []
    found_element = False
    for sensi_case in cases_list:
        group_nb = sensi_case.group_nb
        sensi_df = sensi_case.sensitivities
        sensi_vec = []
        iso_reac_notfound = []
        for iso, reac in common_list:
            # iso, reac = str(iso), str(reac)
            row = sensi_df[(sensi_df["ISO"] == iso) & (sensi_df["REAC"] == reac)]
            if len(row) != 0:
                sensi_vec += list(row["SENSI"])[0]
                found_element = True
            else:
                sensi_vec += [0.0] * group_nb
                iso_reac_notfound.append((iso, reac))

        iso_reac_notfound = [convert_iso_id_to_string(iso) + " " + str(reac) for iso, reac in iso_reac_notfound]
        if len(iso_reac_notfound) > 0:
            warn(
                f"No sensitivity data found for the following isotopes-reactions in the case {sensi_case.casename} - zeros were inserted for the vector construction : \n{list(set(iso_reac_notfound))}"
            )

        sensi_vecs_output.append(np.array(sensi_vec))

    if not found_element and iso_reac_list != None:
        raise errors.EmptyParsingError(
            f"No iso-reac from the given list has been found inside the sensitivities vectors \n\
        Given list : {iso_reac_list}"
        )

    return (sensi_vecs_output, common_list)


def make_sensi_vectors_and_cov_matrix(
    cases_list: list,
    cov_dataf: pd.DataFrame,
    iso_reac_list: list = None,
    reac_list=None,
    iso_list: list = None,
    exclude_iso: list = None,
    operation="union",
):

    common_list = get_common_iso_reac_list(
        cases_list=cases_list, iso_reac_list=iso_reac_list, reac_list=reac_list, iso_list=iso_list, exclude_iso=exclude_iso, operation=operation
    )

    cov_mat, iso_reac_list_present = make_cov_matrix(cov_dataf=cov_dataf, iso_reac_list=common_list)

    (sensi_vecs_output, common_list) = make_sensi_vectors(cases_list=cases_list, iso_reac_list=iso_reac_list_present, operation=operation)

    return sensi_vecs_output, cov_mat, common_list


def convert_iso_id_to_string(iso_id):
    """
    Convert an isotope ID to its corresponding string representation.
    It uses the dictionary iso_z_inv to map the isotope atomic number to its corresponding string.

    Parameters
    ----------
    iso_id : int or str
        The isotope ID to convert.

    Returns
    -------
    str
        The corresponding isotope string representation.
    """
    if isinstance(iso_id, str):
        iso_id = int(iso_id)

    try:
        iso_str = iso_z_inv[int(iso_id / 1000)] + str(int(str(iso_id)[-3:]))
    except:
        warn(f"Isotope ID '{iso_id}' is not assigned to any know isotope name, or is not a type str or int")
        iso_str = str(iso_id)

    return iso_str


def convert_iso_string_to_id(iso_str):
    """
    Convert an isotope string to its corresponding ID.
    The string can be in the format 'A-1' or 'A1' (e.g., 'U-235', 'U235').
    It uses the dictionary iso_z to map the isotope to its corresponding atomic number.

    Parameters
    ----------
    iso_str : str
        The isotope string to convert.

    Returns
    -------
    int
        The corresponding isotope ID.
    """

    try:
        [(z, a)] = re.findall("([A-Za-z]+)-?([0-9]+)", iso_str)
        z = z.upper()
        z = iso_z[z]
        id = z * 1000 + int(a)

    except:
        raise errors.UserInputError(
            f"Isotope '{iso_str}' is not assigned to any know isotope name, or should be an int if its an isotope ID, or is not a type str if its an isotope name."
        )

    return id


# ====================================================================================================


@log_exec()
def calcul_E(case_1, case_2, return_iso_reac_list=False, iso_reac_list=None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None):
    """
    Calculate the E similarity coefficient between two cases.

    Parameters
    ----------
    case_1 : str, Path, or Case object
        The first case for which the E is calculated.
    case_2 : str, Path, or Case object
        The second case for which the E is calculated.
    return_iso_reac_list : bool, optional
        Flag to return the iso-reac list with the E value. Default is False.
    iso_reac_list : list, optional
        The list of iso-reac pairs to consider. If None, all iso-reac pairs are used.
    reac_list : list, optional
        The list of reactions to consider. If None, all reactions are used.
    iso_list : list, optional
        The list of isotopes to consider. If None, all isotopes are used.
    exclude_iso : list, optional
        The list of isotopes to exclude. If None, no isotopes are excluded.

    Returns
    -------
    float
        The E similarity coefficient between the two cases.
    or
    tuple (float, list)
        If return_iso_reac_list is True, returns a tuple with the E value and the iso-reac list.
    """
    cases = []
    for case_i in [case_1, case_2]:
        if isinstance(case_i, (Path, str)):

            cases.append(classes.Case(sdf_path=case_i))

        elif isinstance(case_i, classes.Case):
            cases.append(case_i)

        else:
            raise TypeError(f"Wrong sensitivity type for {case_i} - Choose case object, Path, or string")

    [sensi_vec1, sensi_vec2], iso_reac_list = make_sensi_vectors(
        cases_list=cases, operation="union", iso_reac_list=iso_reac_list, reac_list=reac_list, iso_list=iso_list, exclude_iso=exclude_iso
    )

    E = abs(sensi_vec1 @ sensi_vec2 / (np.linalg.norm(sensi_vec1) * np.linalg.norm(sensi_vec2)))

    if return_iso_reac_list:
        return E, iso_reac_list
    else:
        return E


@log_exec()
def calcul_SS(
    study_case, bench_case, return_iso_reac_list=False, iso_reac_list=None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None
):
    """
    Calculate the SS similarity coefficient between two cases.

    Parameters
    ----------
    case_1 : str, Path, or Case object
        The first case for which the SS is calculated.
    case_2 : str, Path, or Case object
        The second case for which the SS is calculated.
    return_iso_reac_list : bool, optional
        Flag to return the iso-reac list with the SS value. Default is False.
    iso_reac_list : list, optional
        The list of iso-reac pairs to consider. If None, all iso-reac pairs are used.
    reac_list : list, optional
        The list of reactions to consider. If None, all reactions are used.
    iso_list : list, optional
        The list of isotopes to consider. If None, all isotopes are used.
    exclude_iso : list, optional
        The list of isotopes to exclude. If None, no isotopes are excluded.

    Returns
    -------
    float
        The SS similarity coefficient between the two cases.
    or
    tuple (float, list)
        If return_iso_reac_list is True, returns a tuple with the SS value and the iso-reac list.
    """
    if isinstance(study_case, (Path, str)):

        study_case = classes.Case(sdf_path=study_case)

    elif not isinstance(study_case, classes.Case):
        raise TypeError(f"Wrong sensitivity type for {study_case}- Choose case object, Path, or string")

    if isinstance(bench_case, (Path, str)):

        bench_case = classes.Case(sdf_path=bench_case)

    elif not isinstance(bench_case, classes.Case):
        raise TypeError(f"Wrong sensitivity type {bench_case}- Choose case object, Path, or string")

    [study_vec, bench_vec], iso_reac_list = make_sensi_vectors(
        cases_list=[study_case, bench_case],
        operation="union",
        iso_reac_list=iso_reac_list,
        reac_list=reac_list,
        iso_list=iso_list,
        exclude_iso=exclude_iso,
    )

    integ_abs_study_case = 0
    SS_sum = 0
    for val_study, val_bench in zip(study_vec, bench_vec):
        integ_abs_study_case += abs(val_study)
        if val_study * val_bench > 0:
            SS_sum += min([abs(val_study), abs(val_bench)])

    if integ_abs_study_case == 0:
        SS = 0
    else:
        SS = SS_sum / integ_abs_study_case

    if return_iso_reac_list:
        return SS, iso_reac_list
    else:
        return SS


@log_exec()
def calcul_G(
    study_case, bench_case, return_iso_reac_list=False, iso_reac_list=None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None
):
    """
    Calculate the G similarity coefficient between two cases.

    Parameters
    ----------
    case_1 : str, Path, or Case object
        The first case for which the G is calculated.
    case_2 : str, Path, or Case object
        The second case for which the G is calculated.
    return_iso_reac_list : bool, optional
        Flag to return the iso-reac list with the G value. Default is False.
    iso_reac_list : list, optional
        The list of iso-reac pairs to consider. If None, all iso-reac pairs are used.
    reac_list : list, optional
        The list of reactions to consider. If None, all reactions are used.
    iso_list : list, optional
        The list of isotopes to consider. If None, all isotopes are used.
    exclude_iso : list, optional
        The list of isotopes to exclude. If None, no isotopes are excluded.

    Returns
    -------
    float
        The G similarity coefficient between the two cases.
    or
    tuple (float, list)
        If return_iso_reac_list is True, returns a tuple with the G value and the iso-reac list.
    """
    if isinstance(study_case, (Path, str)):

        study_case = classes.Case(sdf_path=study_case)

    elif not isinstance(study_case, classes.Case):
        raise TypeError(f"Wrong sensitivity type for {study_case}- Choose case object, Path, or string")

    if isinstance(bench_case, (Path, str)):

        bench_case = classes.Case(sdf_path=bench_case)

    elif not isinstance(bench_case, classes.Case):
        raise TypeError(f"Wrong sensitivity type for {bench_case} - Choose case object, Path, or string")

    [study_vec, bench_vec], iso_reac_list = make_sensi_vectors(
        cases_list=[study_case, bench_case],
        operation="union",
        iso_reac_list=iso_reac_list,
        reac_list=reac_list,
        iso_list=iso_list,
        exclude_iso=exclude_iso,
    )

    integ_study_case = 0
    G_sum = 0
    for val_study, val_bench in zip(study_vec, bench_vec):
        integ_study_case += val_study
        if val_study * val_bench > 0:
            if abs(val_study) >= abs(val_bench):
                G_sum += val_study - val_bench

        else:
            G_sum += val_study

    G_coef = 1 - (G_sum / integ_study_case)

    if return_iso_reac_list:
        return G_coef, iso_reac_list
    else:
        return G_coef


@log_exec()
def calcul_Ck(
    case_1, case_2, cov_data, return_iso_reac_list=False, iso_reac_list=None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None
):
    """
    Calculate the Ck similarity coefficient between two cases using a covariance DataFrame.

    Parameters
    ----------
    case_1 : str, Path, or Case object
        The first case for which the Ck is calculated.
    case_2 : str, Path, or Case object
        The second case for which the Ck is calculated.
    cov_data : pd.DataFrame or Assimilation object
        The covariance data as a DataFrame or inside an Assimilation object.
    return_iso_reac_list : bool, optional
        Flag to return the iso-reac list with the Ck value. Default is False.
    iso_reac_list : list, optional
        The list of iso-reac pairs to consider. If None, all iso-reac pairs are used.
    reac_list : list, optional
        The list of reactions to consider. If None, all reactions are used.
    iso_list : list, optional
        The list of isotopes to consider. If None, all isotopes are used.
    exclude_iso : list, optional
        The list of isotopes to exclude. If None, no isotopes are excluded.

    Returns
    -------
    float
        The Ck similarity coefficient between the two cases.
    or
    tuple (float, list)
        If return_iso_reac_list is True, returns a tuple with the Ck value and the iso-reac list.
    """

    cases = []
    for case_i in [case_1, case_2]:
        if isinstance(case_i, (Path, str)):

            cases.append(classes.Case(sdf_path=case_i))

        elif isinstance(case_i, classes.Case):
            cases.append(case_i)

        else:
            raise TypeError(f"Wrong sensitivity type for {case_i} - Choose case object, Path, or string")

    if isinstance(cov_data, pd.DataFrame):

        [sensi_vec1, sensi_vec2], cov_mat, iso_reac_list = make_sensi_vectors_and_cov_matrix(
            cases_list=cases,
            cov_dataf=cov_data,
            iso_reac_list=iso_reac_list,
            reac_list=reac_list,
            iso_list=iso_list,
            exclude_iso=exclude_iso,
            operation="union",
        )
        check_dimmensions(casename=cases[0].casename, sensi_vec=sensi_vec1, cov_mat=cov_mat, iso_reac_list=iso_reac_list)

    elif isinstance(cov_data, classes.Assimilation):

        if iso_reac_list != None or reac_list != None or iso_list != None or exclude_iso != None:
            warn(
                f"The list of iso-reac taking into account is fixed by the object Assimilation in argument. If you need to restrict the list of iso-reac, insert a cov-matrix (Dataframe format) as argument 'cov_data'."
            )

        iso_reac_list = cov_data.iso_reac_list

        [sensi_vec1, sensi_vec2], dump = make_sensi_vectors(cases_list=cases, iso_reac_list=iso_reac_list)

        cov_mat = np.substract(cov_data.cov_mat, cov_data.cov_mat_delta)

        check_dimmensions(casename=cases[0].casename, sensi_vec=sensi_vec1, cov_mat=cov_mat, iso_reac_list=iso_reac_list)

    if (sensi_vec1 @ cov_mat @ sensi_vec1) * (sensi_vec2 @ cov_mat @ sensi_vec2) == 0:
        Ck = 0
    else:
        Ck = ((sensi_vec1 @ cov_mat @ sensi_vec2) ** 2) / ((sensi_vec1 @ cov_mat @ sensi_vec1) * (sensi_vec2 @ cov_mat @ sensi_vec2))
        Ck = sqrt(abs(Ck))

    if return_iso_reac_list:
        return Ck, iso_reac_list
    else:
        return Ck


@log_exec()
def calcul_uncertainty(
    study_case, cov_data, iso_reac_list=None, reac_list: list = None, iso_list: list = None, exclude_iso: list = None, output_html_path=None
):
    """
    Calculate the uncertainty of a study case using a covariances Dataframe.

    Parameters
    ----------
    study_case : str, Path, or Case object
        The study case for which the uncertainty is calculated.
    cov_data : pd.DataFrame or Assimilation object
        The covariance data as a DataFrame or inside an Assimilation object.
    iso_reac_list : list, optional
        The list of iso-reac pairs to consider. If None, all iso-reac pairs are used.
    reac_list : list, optional
        The list of reactions to consider. If None, all reactions are used.
    iso_list : list, optional
        The list of isotopes to consider. If None, all isotopes are used.
    exclude_iso : list, optional
        The list of isotopes to exclude. If None, no isotopes are excluded.

    Returns
    -------
    Uncertainty object
        The uncertainty object containing the uncertainty value, decomposition, and a function to store the results inside a HTML file.
    """
    if isinstance(study_case, (Path, str)):

        study_case = classes.Case(sdf_path=study_case)

    elif isinstance(study_case, classes.Case):
        study_case = study_case

    else:
        raise TypeError(f"Wrong sensitivity type for {study_case}- Choose case, Path, or string")

    if isinstance(cov_data, pd.DataFrame):

        [sensi_vec], cov_mat, iso_reac_list = make_sensi_vectors_and_cov_matrix(
            cases_list=[study_case], cov_dataf=cov_data, iso_reac_list=iso_reac_list, reac_list=reac_list, iso_list=iso_list, exclude_iso=exclude_iso
        )

    elif isinstance(cov_data, classes.Assimilation):

        if iso_reac_list != None or reac_list != None or iso_list != None or exclude_iso != None:
            warn(
                f"The list of iso-reac taking into account is fixed by the object Assimilation in argument. If you need to restrict the list of iso-reac, insert a cov-matrix (Dataframe format) as argument 'cov_data'."
            )

        iso_reac_list = cov_data.iso_reac_list

        [sensi_vec], dump = make_sensi_vectors(cases_list=[study_case], iso_reac_list=iso_reac_list)

        cov_mat = np.subtract(cov_data.cov_mat, cov_data.cov_mat_delta)

    check_dimmensions(casename=study_case.casename, sensi_vec=sensi_vec, cov_mat=cov_mat, iso_reac_list=iso_reac_list)
    check_correspondences(sensi_vec=sensi_vec, cov_mat=cov_mat, iso_reac_list=iso_reac_list, group_nb=study_case.group_nb)

    unc = classes.Uncertainty(
        study_case=study_case,
        cov_data=cov_data,
        sensi_vec=sensi_vec,
        resp_calc=study_case.resp_calc,
        cov_mat=cov_mat,
        iso_reac_list=iso_reac_list,
        output_html_path=output_html_path,
    )

    return unc


@log_exec()
def calcul_bias(study_case, assimilation):
    """
    Calculate the bias of a study case using an Assimilation object.

    Parameters
    ----------
    study_case : str, Path, or Case object
        The study case for which the bias is calculated.
    assimilation : Assimilation object
        The assimilation object containing the covariance matrix and other parameters.

    Returns
    -------
    Bias object
        The bias object containing the bias value and decomposition.
    """

    if not isinstance(assimilation, classes.Assimilation):

        raise TypeError(f"Wrong Assimilation object type")

    if isinstance(study_case, (Path, str)):

        study_case = classes.Case(sdf_path=study_case)

    elif isinstance(study_case, classes.Case):
        study_case = study_case

    else:
        raise TypeError(f"Wrong sensitivity type for {study_case}- Choose case, Path, or string")

    [sensi_vec], dump = make_sensi_vectors(cases_list=[study_case], iso_reac_list=assimilation.iso_reac_list)

    check_correspondences(sensi_vec=sensi_vec, cov_mat=assimilation.cov_mat, iso_reac_list=assimilation.iso_reac_list, group_nb=assimilation.group_nb)

    bias_post = classes.Bias(
        study_case=study_case,
        sensi_vec=sensi_vec,
        resp_calc=study_case.resp_calc,
        delta_mu=assimilation.delta_mu,
        iso_reac_list=assimilation.iso_reac_list,
    )

    return bias_post


def check_dimmensions(casename, sensi_vec, cov_mat, iso_reac_list):

    if len(sensi_vec) != np.shape(cov_mat)[0]:
        raise errors.DimError(
            f"The energy group number of your study case is not equal to the covariance matrix group number - Check your sdf file {casename}, and you covariance file\n\
            DIMENSIONS : Cov Mat {np.shape(cov_mat)} (should be square matrix) | Study Case Vec {len(sensi_vec)}\n\
            GROUP NB : Cov Mat {int(np.shape(cov_mat)[0]/len(iso_reac_list))} (should be square matrix) | Study Case Vec {int(len(sensi_vec)/len(iso_reac_list))}"
        )


def check_correspondences(sensi_vec, cov_mat, iso_reac_list, group_nb):

    dikt = {"iso-reac": [], "sensi-integral-abs": [], "cov-integral": []}
    for i, iso_reac in enumerate(iso_reac_list):
        dikt["iso-reac"].append(iso_reac)
        dikt["sensi-integral-abs"].append(np.sum(np.abs(sensi_vec[i * group_nb : (i + 1) * group_nb])))
        dikt["cov-integral"].append(np.sum(cov_mat[i * group_nb : (i + 1) * group_nb, i * group_nb : (i + 1) * group_nb]))

    df = pd.DataFrame(dikt)
    df = df.sort_values(by="sensi-integral-abs", ascending=False, key=lambda col: abs(col))

    nb_iso_reac_check = 100
    i = 0
    while i < nb_iso_reac_check and i < len(df):
        if df.iloc[i, 2] == 0.0:
            warn(
                f"The isotope-reaction {convert_iso_id_to_string(df.iloc[i, 0][0])} - {reac_trad[str(df.iloc[i, 0][1])]} doesn't have data inside the covariances matrix you selected ; it is the {i+1}-th isotope-reaction with the highest absolute integral sensitivity."
            )

        i += 1


def compare_sensitivities(
    case1,
    case2,
    sigma_nb_target=4,
    tol_on_total_sensi_percent=0.1,
    tol_between_sensi_percent=0.0,
    results_format="terminal",
):
    """Compare integral sensitivities for each isotope-reaction and display values outside sigma_target sigma"""

    if isinstance(case1, (Path, str)):
        case1 = classes.Case(sdf_path=case1)
    elif isinstance(case1, classes.Case):
        pass
    else:
        raise errors.SensInputError("bad type for case1")

    if isinstance(case2, (Path, str)):
        case2 = classes.Case(sdf_path=case2)
    elif isinstance(case2, classes.Case):
        pass
    else:
        raise errors.SensInputError("bad type for case2")

    if not (type(sigma_nb_target) == int or type(sigma_nb_target) == float):
        raise errors.SensInputError(f"bad type for sigma_target: {type(sigma_nb_target)}")
    if sigma_nb_target < 0.0:
        raise errors.SensInputError("value for sigma_target must be positive")

    list_values_results = ["dataframe", "string"]
    if results_format not in list_values_results:
        raise errors.SensInputError(f"value of results_format must be choosen in {list_values_results}")

    if tol_on_total_sensi_percent < 0.0 or tol_between_sensi_percent < 0.0:
        raise errors.SensInputError(f"value for tol must be >= 0.")

    case1_sens = case1.sensitivities.sort_values(by=["ISO", "REAC"]).reset_index(drop=True)
    case2_sens = case2.sensitivities.sort_values(by=["ISO", "REAC"]).reset_index(drop=True)

    case1_sens["key"] = case1_sens.apply(lambda row: (row["ISO"], row["REAC"]), axis=1)
    case2_sens["key"] = case2_sens.apply(lambda row: (row["ISO"], row["REAC"]), axis=1)
    commun_keys = set(case1_sens["key"]).intersection(set(case2_sens["key"]))
    # write_and_print(len(commun_keys))

    unique_keys_case1 = set(case1_sens["key"]) - set(case2_sens["key"])
    unique_keys_case2 = set(case2_sens["key"]) - set(case1_sens["key"])
    case1_sens_common = case1_sens[case1_sens["key"].isin(commun_keys)].drop(columns=["key"]).sort_values(by=["ISO", "REAC"]).reset_index(drop=True)
    case2_sens_common = case2_sens[case2_sens["key"].isin(commun_keys)].drop(columns=["key"]).sort_values(by=["ISO", "REAC"]).reset_index(drop=True)
    sum_abs_sensi_1 = case1_sens_common["SENSI_INTEGRAL_ABS"].sum()
    sum_abs_sensi_2 = case2_sens_common["SENSI_INTEGRAL_ABS"].sum()

    case_1_sens_unique = (
        case1_sens[case1_sens["key"].isin(unique_keys_case1)].drop(columns=["key"]).sort_values(by=["ISO", "REAC"]).reset_index(drop=True)
    )
    case_2_sens_unique = (
        case2_sens[case2_sens["key"].isin(unique_keys_case2)].drop(columns=["key"]).sort_values(by=["ISO", "REAC"]).reset_index(drop=True)
    )

    # write_and_print(case1_sens_common)
    # write_and_print(case2_sens_common)

    not_in_agreement = case1_sens_common.apply(
        lambda row: not compare_row(row, case2_sens_common.loc[row.name], sigma_nb_target=sigma_nb_target), axis=1
    )
    case1_not_in_agreement = case1_sens_common[not_in_agreement]
    case2_not_in_agreement = case2_sens_common[not_in_agreement]

    # write_and_print(case1_not_in_agreement)
    # write_and_print(case2_not_in_agreement)

    case_not_in_agreement = case1_not_in_agreement[["ISO", "REAC", "SENSI_INTEGRAL", "STD_INTEGRAL", "SENSI_INTEGRAL_ABS"]]
    case_not_in_agreement = case_not_in_agreement.rename(
        columns={
            "SENSI_INTEGRAL": "SENSI_INTEGRAL_1",
            "STD_INTEGRAL": "STD_INTEGRAL_1",
            "SENSI_INTEGRAL_ABS": "SENSI_INTEGRAL_ABS_1",
        }
    )
    case_not_in_agreement = pd.concat(
        [
            case_not_in_agreement,
            case2_not_in_agreement[["SENSI_INTEGRAL", "STD_INTEGRAL", "SENSI_INTEGRAL_ABS"]],
        ],
        axis=1,
    )
    case_not_in_agreement = case_not_in_agreement.rename(
        columns={
            "SENSI_INTEGRAL": "SENSI_INTEGRAL_2",
            "STD_INTEGRAL": "STD_INTEGRAL_2",
            "SENSI_INTEGRAL_ABS": "SENSI_INTEGRAL_ABS_2",
        }
    )
    case_not_in_agreement["%SUM_ABS_1"] = 100 * case_not_in_agreement["SENSI_INTEGRAL_ABS_1"] / float(sum_abs_sensi_1)
    case_not_in_agreement["%SUM_ABS_2"] = 100 * case_not_in_agreement["SENSI_INTEGRAL_ABS_2"] / float(sum_abs_sensi_2)
    case_not_in_agreement_reduce = case_not_in_agreement[
        (case_not_in_agreement["%SUM_ABS_1"] > tol_on_total_sensi_percent) | (case_not_in_agreement["%SUM_ABS_2"] > tol_on_total_sensi_percent)
    ]
    if tol_between_sensi_percent > 0.0:
        case_not_in_agreement_reduce = case_not_in_agreement_reduce[
            abs(1 - (case_not_in_agreement_reduce["SENSI_INTEGRAL_1"] / case_not_in_agreement_reduce["SENSI_INTEGRAL_2"]))
            >= (tol_between_sensi_percent / 100)
        ]

    # write_and_print("agreement")
    # write_and_print(case_not_in_agreement_reduce)

    if not case_not_in_agreement_reduce.empty:
        if results_format == "string":
            results = str(case_not_in_agreement_reduce)
        elif results_format == "dataframe":
            results = case_not_in_agreement_reduce
    else:
        results = None

    return results


def compare_row(row1, row2, sigma_nb_target):
    return is_values_in_agreement(
        row1["SENSI_INTEGRAL"],
        row2["SENSI_INTEGRAL"],
        row1["STD_INTEGRAL"],
        row2["STD_INTEGRAL"],
        sigma_nb_target=4,
    )


def is_values_in_agreement(val1, val2, sigma1, sigma2, sigma_nb_target):
    """test if val1 and val2 are in agreement within tol sigma"""

    if val1 < val2:
        val1, val2 = val2, val1

    if (val1 - sigma_nb_target * sigma1) <= (val2 + sigma_nb_target * sigma2):
        return True
    else:
        return False


# ====================================================================================================


def filter_cases(
    dict_data,
    to_keep=None,
    to_exclude=None,
):
    """
    Filter the list of cases based on keep and exclude patterns.

    Parameters
    -----------
    dict_data : dict
        [Required] Dictionnary with containing all cases
    to_keep : list of str, optional
        Patterns of cases to keep. Defaults to None, i.e. keeping all cases.
    to_exclude : list of str, optional
        Patterns of cases to exclude. Defaults to None, i.e. excluding no cases.

    Returns
    --------
    list
        List of case names.
    """
    cases = list(dict_data.keys())
    if to_keep:
        cases = [c for c in cases if any([k.lower() in c.lower() for k in to_keep])]
    if to_exclude:
        cases = [c for c in cases if not any([k.lower() in c.lower() for k in to_exclude])]
    return cases


def filter_common_data(dict_sensitivity):

    all_cases = [set(data["data"].keys()) for data in dict_sensitivity.values()]
    common_cases = set.intersection(*all_cases) if all_cases else set()

    filtered_dict = {}

    for combination, data in dict_sensitivity.items():
        filtered_data = {case: info for case, info in data["data"].items() if case in common_cases}

        if filtered_data:
            filtered_dict[combination] = filtered_data

    return filtered_dict


def find_sensitive_cases_for_combination(
    dict_data,
    cases,
    isotope,
    reaction,
    energy_region,
    sensitivity_threshold=None,
    fraction_threshold=None,
    integrate_data=False,
    normalize_lethargy=False,
):
    # Sub-function to process a single combination for find_sensitive_cases

    if sensitivity_threshold is None:
        sensitivity_threshold = 0

    if fraction_threshold is None:
        fraction_threshold = 0

    dict_sensitive = {
        "search_criteria": {
            "isotope": isotope,
            "reaction": reaction,
            "energy_region": energy_region,
            "sensitivity_threshold": sensitivity_threshold,
            "fraction_threshold": fraction_threshold,
            "integrate_data": integrate_data,
            "normalize_lethargy": normalize_lethargy,
        },
        "data": {},
    }

    if isinstance(isotope, str):
        isotope = convert_iso_string_to_id(isotope)

    if not isinstance(reaction, int):
        raise errors.UserInputError("Reactions must be specified as integers")

    for case in cases:
        case_data = dict_data[case]
        energy_indices, energy_bins = case_data.filter_energy_bins(energy_region)

        try:
            sensi_values = case_data.get_filtered_sensitivity_by_energy_region(
                isotope, reaction, energy_region, normalize_lethargy=normalize_lethargy
            )

            sensi_integral = case_data.get_integral_sensitivity_value(isotope, reaction, absolute_value=True)

            if integrate_data and not normalize_lethargy:
                sensi_values = [np.sum(sensi_values)]
                energy_bins = [energy_bins[0], energy_bins[-1]]
            if integrate_data and normalize_lethargy:
                raise errors.UserInputError("The integrate_data option cannot be enabled with data normalized to lethargy")
            sensitive_bins = [
                {"value": sensi_values[i], "e_bin": [energy_bins[i], energy_bins[i + 1]]}
                for i in range(len(sensi_values))
                if np.abs(sensi_values[i]) > sensitivity_threshold and np.abs(sensi_values[i]) / sensi_integral > fraction_threshold
            ]
            if sensitive_bins:
                max_index, max_value = max(enumerate(sensi_values), key=lambda x: np.abs(x[1]))
                max_sensitivity = {"value": max_value, "e_bin": [energy_bins[max_index], energy_bins[max_index + 1]]}
                dict_sensitive["data"][case] = {"sensitive_bins": sensitive_bins, "max_sensitivity": max_sensitivity}
        except errors.MissingDataError:
            # Skip the case if required data is missing
            continue

    return dict_sensitive


@log_exec()
def find_sensitive_cases(
    dict_data,
    cases,
    isotopes,
    reactions,
    energy_regions,
    sensitivity_threshold=None,
    fraction_threshold=None,
    integrate_data=True,
    normalize_lethargy=False,
):
    """
    Identifies cases with sensitivities above specified thresholds within selected isotopes, reactions, and energy regions.

    Parameters
    ----------
    dict_data : dict
        [Required] Dictionary containing sensitivity data for each case.
    cases : list
        [Required] List of case_id to analyze for sensitivity.
    isotopes : list of int or str
        [Required] Isotope(s) to analyze (e.g., 92235 or U235). Can be a single isotope or a list.
    reactions : list of int
        [Required] Nuclear reaction(s) to analyze (e.g., 18). Can be a single reaction or a list.
    energy_regions : list of str or list of lists
        [Required] Energy region(s) over which sensitivity is evaluated. Can be a string, a list of strings, or a list of sublists
        (each sublist should have length 2 energy limits [Upper limit, Lower limit]).
    sensitivity_threshold : float, optional
        Minimum sensitivity value for a bin to be considered significant. Default is None (interpreted as 0).
    fraction_threshold : float, optional
        Minimum fraction of total sensitivity for a bin to be considered significant. Default is None (interpreted as 0).
    integrate_data : bool, optional
        If True, integrates sensitivity values across the energy region. Cannot be used with `normalize_lethargy`.
        Default is True.
    normalize_lethargy : bool, optional
        If True, normalizes sensitivity values to lethargy. Cannot be used with `integrate_data`. Default is False.

    Returns
    -------
    dict
        A dictionary containing two keys:
        - "search_criteria": A dictionary summarizing the criteria used for the search.
        - "data": The filtered sensitivity data matching the specified criteria.
    """

    inputs = [isotopes, reactions, energy_regions]
    converted_inputs = [[arg] if isinstance(arg, (int, str)) else arg for arg in inputs]

    isotopes, reactions, energy_regions = converted_inputs

    if not all(isinstance(arg, list) for arg in [isotopes, reactions, energy_regions]):
        raise TypeError("All inputs (isotopes, reactions, energy_regions) must be of type list.")

    if not (len(isotopes) == len(reactions) == len(energy_regions)):
        raise ValueError("All input lists (isotopes, reactions, energy_regions) must have the same length.")

    for item in energy_regions:
        if isinstance(item, list):
            if len(item) != 2:
                raise ValueError("Each sublist in energy_regions must have a length of 2.")
        elif not isinstance(item, str):
            raise TypeError("Each element in energy_regions must be either a string or a list of length 2.")

    matched_cases = set(cases)
    dict_temp = {}

    for isotope, reaction, energy_region in zip(isotopes, reactions, energy_regions):

        sensitivity = find_sensitive_cases_for_combination(
            dict_data,
            cases,
            isotope,
            reaction,
            energy_region,
            sensitivity_threshold,
            fraction_threshold,
            integrate_data,
            normalize_lethargy,
        )

        if not isinstance(energy_region, (str)):
            energy_region_str = f"{energy_region[0]}_{energy_region[1]}"
        else:
            energy_region_str = energy_region

        dict_temp[f"{isotope}-{reaction}-{energy_region}"] = sensitivity

    filtered_data = filter_common_data(dict_temp)

    dict_sensitive = {
        "search_criteria": {
            "isotopes": isotopes,
            "reactions": reactions,
            "energy_regions": energy_regions,
            "sensitivity_threshold": sensitivity_threshold,
            "fraction_threshold": fraction_threshold,
            "integrate_data": integrate_data,
            "normalize_lethargy": normalize_lethargy,
        },
        "data": filtered_data,
    }

    return dict_sensitive


def build_sensitive_dataframe(dict_sensitive):
    """
    Converts a dictionary of sensitivity data into a DataFrame with each `case_id` as a row.
    Each unique `isotope-reaction-energy` combination creates two columns: one for `max_value` and one for `e_bin`.
    Also adds an `absolute_sum` column, which is the sum of absolute max sensitivity values for each `case_id`.

    Arugments :
    ---------
    dict_sensitive : dict
        [Required] Nested dictionary with sensitivity data organized by `isotope-energy` combinations and `case_id`s.

    Returns
    --------
    pd.DataFrame
        DataFrame with rows for each `case_id`, columns for `max_value`, `e_bin` per combination, and `absolute_sum` for total sensitivity. Sorted by `absolute_sum` in descending order.
        If no sensitive cases meet the criteria, the function warns and returns `None`.
    """
    records = []

    # Traverse the nested dictionary to extract values in a pivot-friendly format

    if not dict_sensitive.items():
        warn(f"No cases sensitive to the specified criteria were found")
        return None

    else:
        for combination, case_ids in dict_sensitive.items():

            for case_id, case_id_data in case_ids.items():
                # Create a dictionary for each case_id to store values under column names based on isotope-energy
                temp_record = {"case_id": case_id}

                # Extract max_sensitivity information
                max_sensitivity = case_id_data.get("max_sensitivity", {})
                max_sensitivity_value = max_sensitivity.get("value")
                max_sensitivity_e_bin = max_sensitivity.get("e_bin")

                # Add max sensitivity data to temp_record using isotope-energy as prefix
                temp_record[f"{combination}-max_value"] = max_sensitivity_value
                temp_record[f"{combination}-e_bin"] = max_sensitivity_e_bin

                records.append(temp_record)

        # Convert the list of dictionaries into a DataFrame
        df = pd.DataFrame(records)
        # Pivot the DataFrame to ensure that each `case_id` has a unique row
        df_pivoted = df.pivot_table(index="case_id", aggfunc="first").reset_index()

        # Calculate the sensitivity sum of absolute values
        # Select only the max_sensitivity_value columns for summing
        sensitivity_columns = [col for col in df_pivoted.columns if col.endswith("-max_value")]
        df_pivoted["absolute_sum"] = df_pivoted[sensitivity_columns].abs().sum(axis=1)

        # Sort by the sensitivity sum in descending order
        df_pivoted = df_pivoted.sort_values(by="absolute_sum", ascending=False).reset_index(drop=True)

        # Order columns as case_id, absolute_sum, alternating max_sensitivity and e_bin columns
        ordered_columns = ["case_id", "absolute_sum"]

        # Add max_sensitivity and e_bin columns in the desired alternating order
        for combination, case_ids in dict_sensitive.items():
            ordered_columns.append(f"{combination}-max_value")
            ordered_columns.append(f"{combination}-e_bin")

        # Reorder the columns
        df_pivoted = df_pivoted[ordered_columns]
        df_pivoted.index = range(1, len(df_pivoted) + 1)

        return df_pivoted


@log_exec()
def display_most_sensitive_cases(dict_sensitive, max_number_of_cases_displayed=50):
    """
    Display the most sensitive cases for each combination in the sensitivity dictionary, sorting sensitivities within each combination.
    This function prints a formatted table showing the rank, case name, maximum sensitivity value across bins,
    and the associated energy bin for each of the most sensitive cases, sorted within each combination.

    Parameters
    -----------
    dict_sensitive : dict
        [Required] Dictionary containing sensitivity data and search criteria.
    max_number_of_cases_displayed : int, optional
        Maximum number of cases to display per combination. Default is 50.
    """

    search_criteria = dict_sensitive["search_criteria"]
    isotopes = search_criteria["isotopes"]
    reactions = search_criteria["reactions"]
    energy_regions = search_criteria["energy_regions"]
    normalize_lethargy = search_criteria["normalize_lethargy"]
    integrate_data = search_criteria["integrate_data"]
    sensitivity_threshold = search_criteria["sensitivity_threshold"]
    fraction_threshold = search_criteria["fraction_threshold"]

    # Convert search criteria to strings for display
    sensitivity_threshold_str = f"sensitivity_threshold = {sensitivity_threshold}" if sensitivity_threshold else "no sensitivity_threshold"
    fraction_threshold_str = f"fraction_threshold = {fraction_threshold}" if fraction_threshold else "no fraction_threshold"
    data_type = "integral" if integrate_data else "per bin"
    normalization = "per unit lethargy" if normalize_lethargy else "raw"
    criteria_str = f"({sensitivity_threshold_str} ; {fraction_threshold_str})"
    settings_str = f"using {data_type}, {normalization} sensitivities:"

    if len(isotopes) > 1:
        message = (
            f"Displaying sensitive cases {criteria_str} for reactions {reactions} "
            f"of isotopes {isotopes} in energy regions {energy_regions}, {settings_str}\n"
            "\nInformation: only cases sensitive to all specified combinations will be displayed.\n"
        )
    else:
        message = (
            f"Displaying sensitive cases {criteria_str} for reaction {reactions[0]} "
            f"of isotope {isotopes[0]} in energy region {energy_regions[0]}, {settings_str}\n"
        )

    # Build dataframe
    df_sensitive = build_sensitive_dataframe(dict_sensitive["data"])

    if df_sensitive is not None:
        num_rows = df_sensitive.shape[0]

        if num_rows > max_number_of_cases_displayed:
            message += f"\n{num_rows} sensitive cases were found. Displaying the first {max_number_of_cases_displayed}:\n\n"
        else:
            message += f"\n{num_rows} sensitive cases were found:\n\n"

        # Set display options to show all rows and columns without truncation
        pd.options.display.max_rows = None
        pd.options.display.max_columns = None
        pd.options.display.width = None

        # Display
        message += f"{'=' * 40}\n{df_sensitive.head(max_number_of_cases_displayed)}\n{'=' * 40}\n"
        write_and_print(message)

        # Reset to default after displaying to avoid excessive outputs in future
        pd.reset_option("display.max_rows")
        pd.reset_option("display.max_columns")
        pd.reset_option("display.width")


@log_exec()
def plot_sensitivity_heatmap(
    dict_data,
    cases,
    isotopes,
    reactions,
    sensitivity_threshold=None,
    fraction_threshold=None,
    user_energy_bins=None,
    title=None,
    namefig=None,
    save_logfile=True,
):
    """
    Generates a sensitivity heatmap for given isotopes and reactions across specified energy bins.

    Parameters
    -----------
    dict_data : dict
        [Required] Dictionary containing the sensitivity data.
    cases : list
        [Required] List of cases to analyze.
    isotopes : list of int or str
        [Required] List of isotopes to be considered.
    reactions : list of int
        [Required] List of reactions to be considered.
    sensitivity_threshold : float, optional
        The minimum sensitivity value to consider a bin as significant. Default is None (interpreted as 0).
    fraction_threshold : float, optional
        The minimum fraction of the total sensitivity to consider a bin as significant. Default is None (interpreted as 0).
    user_energy_bins : list, optional
        Custom energy bins to use. If not provided, a default set of bins will be used. Default is None.
    title : str, optional
        Title of the heatmap. Default is None.
    namefig : str, optional
        Filename for the saved heatmap image. Default is None.
    save_logfile : bool, optional
        Whether to save a log file with details of sensitive cases. Default is True.

    Returns
    --------
        The function saves the heatmap as an image file and optionally a log file.
    """
    write_and_print(
        f"Generating a sensitivity heatmap for isotopes {isotopes} and reactions {reactions} (sensitivity_threshold = {sensitivity_threshold} ; fraction_threshold = {fraction_threshold}):"
    )

    if user_energy_bins is None:
        energy_bins = np.array([2e7, 1e6, 1e5, 1e4, 1e3, 1e2, 1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5])
    else:
        energy_bins = user_energy_bins

    num_bins = len(energy_bins) - 1
    num_reactions = len(isotopes) * len(reactions)

    data = np.full((num_reactions, num_bins), np.nan)

    reactions_str = []

    for reaction in reactions:
        if reaction == -2:
            reactions_str.append(101)
        else:
            reactions_str.append(reaction)

    iso_reac_list = [f"{isotope}-{reaction}" for isotope in isotopes for reaction in reactions_str]

    dict_sensitive_cases = {}
    total_sensitive_cases = []
    for i, (isotope, reaction) in enumerate([(iso, reac) for iso in isotopes for reac in reactions]):
        write_and_print(f"- Treating {isotope}-{reaction}")
        dict_sensitive_cases[f"{isotope}-{reaction}"] = {}
        dict_sensitive_cases[f"{isotope}-{reaction}"]["total"] = []

        for j in range(num_bins):
            current_bin = list(energy_bins[j : j + 2])

            sensitivity = find_sensitive_cases_for_combination(
                dict_data,
                cases,
                isotope,
                reaction,
                current_bin,
                sensitivity_threshold,
                fraction_threshold,
                integrate_data=True,
                normalize_lethargy=False,
            )

            sensitive_cases = sorted(
                ((case, np.round(details["max_sensitivity"]["value"], 6)) for case, details in sensitivity["data"].items()),
                key=lambda item: np.abs(item[1]),
                reverse=True,
            )

            dict_sensitive_cases[f"{isotope}-{reaction}"][f"{current_bin}"] = sensitive_cases

            for case in sensitive_cases:
                if case[0] not in dict_sensitive_cases[f"{isotope}-{reaction}"]["total"]:
                    dict_sensitive_cases[f"{isotope}-{reaction}"]["total"].append(case[0])

            num_sensitive = len(list(sensitivity["data"].keys()))

            if num_sensitive > 0:
                data[i, j] = num_sensitive

        total_sensitive_cases.append(len(dict_sensitive_cases[f"{isotope}-{reaction}"]["total"]))

    log_energy_bins = np.log10(energy_bins)

    fig, ax1 = plt.subplots(figsize=(12, 6))
    cax = ax1.pcolormesh(log_energy_bins, np.arange(len(iso_reac_list) + 1), data, cmap="jet")

    ax1.set_xlabel("Energy [eV]")
    ax1.set_ylabel("Isotope/Material and Interaction")

    ax1.set_yticks(np.arange(len(iso_reac_list)) + 0.5)
    ax1.set_yticklabels(iso_reac_list)

    energy_labels = [10**i for i in range(-5, 8)]
    ax1.set_xticks(np.log10(energy_labels))
    ax1.set_xticklabels([f"$10^{{{int(np.log10(e))}}}$" for e in energy_labels])

    # Add text annotations in each cell
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isnan(data[i, j]):
                ax1.text(
                    log_energy_bins[j] + (log_energy_bins[j + 1] - log_energy_bins[j]) / 2,
                    i + 0.5,
                    int(data[i, j]),
                    ha="center",
                    va="center",
                    color="white",
                )

    ax2 = ax1.twinx()
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_yticks(np.arange(len(iso_reac_list)) + 0.5)
    ax2.set_yticklabels([f"{int(total)}" for total in total_sensitive_cases])
    ax2.set_ylabel("Number of Unique Cases", labelpad=10)

    cbar = fig.colorbar(cax, ax=ax1, pad=0.10)
    cbar.set_label("Number of Cases")

    title = "" if title is None else title

    namefig = "heatmap-" + "_".join(map(str, isotopes)) + "-" + "_".join(map(str, reactions)) + ".png" if namefig is None else namefig

    ax1.set_title(title)

    plt.savefig(namefig, dpi=300, bbox_inches="tight")
    write_and_print(f"{namefig} successfully generated")

    if save_logfile is True:
        base, extension = os.path.splitext(namefig)
        namelog = f"log_{base}.txt"
        with open(namelog, "w") as logfile:
            logfile.write(
                f"Summary of sensitivity heatmap for isotopes {isotopes} and reactions {reactions} (sensitivity_threshold = {sensitivity_threshold} ; fraction_threshold = {fraction_threshold}):\n\n"
            )
            for idx, (iso_reac, details) in enumerate(dict_sensitive_cases.items()):
                logfile.write(f"{iso_reac}:\n")
                for bin_range, cases in details.items():
                    if bin_range != "total":
                        logfile.write(f"  Energy bin {bin_range} (sorted): {cases}\n")
                logfile.write(f"  Total sensitive cases (unsorted): {details['total']}")
                if idx < len(dict_sensitive_cases) - 1:
                    logfile.write("\n\n")
        write_and_print(f"{namelog} successfully generated")


@log_exec()
def sort_case_sensitivities(
    case,
    energy_region="all",
    integrate_data=True,
    normalize_lethargy=False,
    iso_list=None,
    reac_list=None,
):
    """
    Sorts case sensitivities by their maximum absolute sensitivity value.

    Parameters
    -----------
    case : str, Path, or Case instance
        [Required] The case object or path to the sensitivity data file.
    energy_region : str, optional, default="all"
        Energy region over which sensitivity is evaluated. Can be a string or a list with format [Upper limit, Lower limit].
    integrate_data : bool, optional, default=True
        If True, integrates the sensitivity data over the energy range.
    normalize_lethargy : bool, optional, default=False
        If True, normalizes data to lethargy (cannot be used with integrate_data).
    iso_list : list of int or str, optional
        List of isotopes to keep. If None, includes all isotopes from the case.
    reac_list : list of int, optional
        List of reactions to keep. If None, includes all reactions from the case.

    Returns
    --------
    pd.DataFrame
        A DataFrame sorted by maximum absolute sensitivity value, containing:
        - 'ISO': Isotope
        - 'REAC': Reaction
        - 'MAX_SENSI_VALUE': Maximum sensitivity value
        - 'ENERGY_BINS': Corresponding energy bins
    """

    if isinstance(case, (Path, str)):
        case = classes.Case(sdf_path=case)
    elif not isinstance(case, classes.Case):
        raise TypeError(f"Invalid sensitivity type for {case} - Choose Case, Path, or string")

    if not isinstance(iso_list, (list, type(None))):
        raise TypeError("iso_list must be None or a list")

    if not isinstance(reac_list, (list, type(None))):
        raise TypeError("reac_list must be None or a list")

    energy_indices, energy_bins = case.filter_energy_bins(energy_region)

    iso_reac_pairs = list(zip(case.sensitivities["ISO"], case.sensitivities["REAC"]))

    if not iso_list:
        iso_list = list({isotope for isotope, _ in iso_reac_pairs})
    else:
        for i, iso in enumerate(iso_list):
            if not isinstance(iso, int):
                iso_list[i] = convert_iso_string_to_id(iso)

    if not reac_list:
        reac_list = list({reaction for _, reaction in iso_reac_pairs})
    else:
        for reac in reac_list:
            if not isinstance(reac, int):
                raise errors.UserInputError("Reactions must be specified as integers")

    dict_data = {"ISO": [], "REAC": [], "MAX_SENSI_VALUE": [], "ENERGY_BINS": []}

    for isotope, reaction in iso_reac_pairs:
        try:
            if isotope in iso_list and reaction in reac_list:
                dict_data["ISO"].append(isotope)
                dict_data["REAC"].append(reaction)

                sensi_values = case.get_filtered_sensitivity_by_energy_region(isotope, reaction, energy_region, normalize_lethargy=normalize_lethargy)

                if integrate_data and not normalize_lethargy:
                    sensi_values = [np.sum(sensi_values)]
                    energy_bins = [energy_bins[0], energy_bins[-1]]
                elif integrate_data and normalize_lethargy:
                    raise errors.UserInputError("The integrate_data option cannot be enabled with data normalized to lethargy")

                max_index, max_value = max(enumerate(sensi_values), key=lambda x: np.abs(x[1]))
                dict_data["MAX_SENSI_VALUE"].append(max_value)
                dict_data["ENERGY_BINS"].append([energy_bins[max_index], energy_bins[max_index + 1]])

        except errors.MissingDataError:
            continue  # Skip missing data cases

    df = pd.DataFrame(dict_data)
    df.sort_values(by="MAX_SENSI_VALUE", ascending=False, inplace=True, key=abs)
    df.index = range(1, len(df) + 1)

    if df.empty:
        warn("No data matching the specified criteria was found")

    return df


@log_exec()
def display_dataframe(df, max_number_of_rows_displayed=30):
    """
    Displays a clean and formatted output of the first X rows of a DataFrame in the terminal.

    Parameters
    -----------
    df : pd.DataFrame
        [Required] The DataFrame to display.
    max_number_of_rows_displayed : int, optional (default=30)
        The maximum number of rows to display.
    """
    if df is None or df.empty:
        warn("Dataframe is empty.")
    else:
        num_rows = df.shape[0]

        if num_rows > max_number_of_rows_displayed:
            message = f"Showing the first {max_number_of_rows_displayed} rows out of {num_rows} total:\n\n"
        else:
            message = "Displaying the full dataframe:\n\n"

        # Set display options to show all rows and columns without truncation
        pd.options.display.max_rows = None
        pd.options.display.max_columns = None
        pd.options.display.width = None

        # Display
        message += f"{'=' * 40}\n{df.head(max_number_of_rows_displayed)}\n{'=' * 40}\n"
        write_and_print(message)

        # Reset to default after displaying to avoid excessive outputs in future
        pd.reset_option("display.max_rows")
        pd.reset_option("display.max_columns")
        pd.reset_option("display.width")
