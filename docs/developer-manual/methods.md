# Developer Manual - methods.py

Matrix multiplications are performed with the Python operator "A @ B", fully equivalent to numpy.matmul.

⚠️Note: operator precedence is very important when it comes to algorithmic matrix multiplication methods, to optimize calculation time.⚠️

---

## Acquisition and formatting of useful data (Pandas Dataframe)
Related functions: *format_ ... _to_dataframe (...)*

These functions read the various files containing covariance matrices and sensitivities, and format these data into Pandas DataFrames for easy manipulation.

Sensitivity DataFrames are built only from *.sdf* files as follows:

| ISO | REAC | SENSI |
| :-: | :-: | :-: |
|iso1|reac1[iso1] | $S^{1, 1}$
|iso1|reac2[iso1] | $S^{1, 2}$
|iso2|reac1[iso2] | $S^{2, 1}$
|iso2|reac2[iso2] | $S^{2, 2}$
|.|. | .
|.|. | .

Where:
- iso i: isotope identification number (Z x 1000 + A), as *string* or *int*;
- reac r[iso i]: reaction identification number (as per *reac_trad* dictionary), independent of isotope i;
- $S^{i, r}$: python list [. , . , ..] including the response's relative multi-group sensitivities to isotope i reaction r, in decreasing energy order;

***Note***: all isotopes and reactions encountered are stored in this dataframe, including the TOTAL reaction - ID: 1. However, a filter on the TOTAL reaction is applied when the sensitivity *numpy* vectors are being constructed.

Covariance matrix DataFrames are built as follows:

| ISO_H | REAC_H | ISO_V | REAC_V | STD |
| :-: | :-: | :-: | :-: | :-: |
|iso1_h|reac1_h[iso1_h] | iso1_v|reac1_v[iso1_v] |$Cov^{(1, 1), (1, 1)}$
|iso1_h|reac2_h[iso1_h] | iso1_v|reac2_v[iso1_v] |$Cov^{(1, 2), (1, 2)}$
|iso1_h|reac1_h[iso1_h] | iso1_v|reac2_v[iso1_v] |$Cov^{(1, 1), (1, 2)}$
|iso2_h|reac2_h[iso2_h] | iso2_v|reac2_v[iso2_v] |$Cov^{(2, 2), (2, 2)}$
|.|.|.|.|.
|.|.|.|.|.

Where:
- iso i_h: isotope identification number (Z x 1000 + A), as *string* or *int*, with energy groups decreasing horizontally in the $Cov^{(i_h, r_h), (i_v, r_v)}$ sub-matrix;
- reac r_h[iso i_h]: reaction identification number (as per *reac_trad* dictionary), with energy groups decreasing horizontally;
- iso i_v: isotope identification number (Z*1000 + A), as *string* or *int*, with energy groups decreasing vertically;
- reac r_v[iso i_v]: reaction identification number (as per *reac_trad* dictionary), with energy groups decreasing vertically;
- $Cov^{(i_h, r_h), (i_v, r_v)}$: 2D python list [[. , . , ..], [. , . , ..], ..] containing the multi-group covariance sub-matrix for reaction r_h of isotope i_h and reaction r_v of isotope i_v;

There is, among others, a format_..._to_dataframe(...) function for each type of covariance matrix. Currently readable matrices are:
- SCALE AMPX format (file form)
- SCALE binary format (file form)
- COMAC (folder of files)
- GENDF

**WARNING**: The **CALINS norm** is to format sensitivity and covariance data in **decreasing** energy group order from 20 MeV -> 0 MeV.

## Construction of matrices and vectors (Numpy Arrays)
Related functions: *make_ ... (...)*

These functions allow building vectors and matrices as *Numpy Arrays* while respecting the **alignment rule** (see chapter *Alignment rule for matrix calculations*).

To ensure that the alignment rule is respected between several sensitivity vectors and several covariance matrices, they are constructed simultaneously on the basis of a common isotope–reaction list. This list dictates the ordering of terms in the sensitivity vectors, as well as the placement of sub-matrices within the full covariance matrix.

This list is first built as either the intersection or the union of the isotope–reaction lists specific to each sensitivity vector involved (depending on the subsequent operations to be performed—by default the "union" option is applied, via the function get_common_iso_reac_list). If a covariance matrix is also involved, the isotope–reaction list is then reduced to its intersection with the isotope–reaction pairs present in the variance–covariance matrix data (this operation is performed in the function make_cov_matrix).

WARNING: If non-conventional isotopes such as "h-poly" (ID=1901) appear in the isotope–reaction list derived from sensitivity vectors but not in the variance–covariance DataFrame, the DataFrame is enriched so as to include variance–covariance data for this isotope, provided its base isotope (here h-1, ID=1001) is already present in the DataFrame. The data used for this non-conventional isotope are copied from those of the base isotope. The current method for detecting non-conventional isotopes is to check whether the hundreds digit of their ID is greater than or equal to 3 (example: h-poly / 1901 → 9 ≥ 3), according to our observations.
The same principle is applied for isotopes with a "natural" form (ID=xx000) in the variance–covariance data. If an isotope from sensitivity vectors is absent from the covariance DataFrame, the code checks whether data for its natural form are present. If so, the DataFrame is enriched with a copy of these data, associated with the originally missing isotope (and its reactions).

Note: At this stage, a filter excluding the TOTAL reaction is applied to the isotope–reaction list.

For sensitivity vectors, this list is iterated in order, and for each isotope–reaction pair the corresponding data are searched in its DataFrame. If present, the multigroup sensitivities are appended to the Numpy Array under construction. If not, zero multigroup sensitivities are inserted instead.

For the covariance matrix, initialization is done with a zero matrix of dimension 2D [(N_group × Nb_iso-reac), (N_group × Nb_iso-reac)]. The isotope–reaction list is then iterated through twice in nested *for* loops (for horizontal and vertical indices). The sub-matrix from the DataFrame corresponding to the two isotope–reaction pairs (i_h, r_h) and (i_v, r_v), ($Cov^{(i_h, r_h), (i_v, r_v)}$), is inserted into $Cov$ at the correct coordinates, if available. The code also checks whether the inverse pairing (i_v, r_v) and (i_h, r_h), i.e. ($Cov^{(i_v, r_v), (i_h, r_h)}$), exists in the DataFrame; if so, the transpose of the sub-matrix is inserted at the corresponding coordinates. Finally, if neither correspondence exists in the DataFrame, the sub-matrix remains zero at those coordinates.*

## Calculation functions
Functions of type: *calcul_... (...)*

These functions use the *make_...* functions to build sensitivity vectors and the covariance matrix based on a shared isotope-reaction list, to perform...

One function exists per similarity index. Implemented indices are:
- E
- C<sub>k</sub>
- G (CEA formula)
- Shared Sensitivities overlap index (Mariya BROVCHENKO formula)

These functions take as input at least two sensitivity vectors and a covariance matrix for C<sub>k</sub> calculations. Sensitivity data can either be a SDF-file path, or *Case* object (presented below), but not preconstructed *numpy* vectors.

A third function allows uncertainty calculation, with a sensitivity case and covariances matrix as input. Sensitivity data can either be a SDF-file path, or *Case* object (presented below), but not preconstructed *numpy* vectors. This function returns an *Uncertainty* object (described in chapter ***classes***).

## Utility functions
The *methods* module also contains various utility functions :
- to convert isotope identification numbers to strings and vice versa ;
- to search in a *Case* database for cases sensitive to specific criterion such as isotope, reaction, energy group and so on ;
- others ...