![](./images/Logo.PNG)
### **CAL**culations and **I**nvestigations on **N**uclear data uncertainties and **S**ensitivities

by French Authority for Nuclear Safety and Radiation Protection - ASNR

- [**INSTALLATION**](#installation)
- [**THEORY**](#theory)
  - [**Introduction**](#introduction)
  - [**Sensitivity vector**](#sensitivity-vector)
  - [**Covariance matrix**](#covariance-matrix)
  - [**Alignment rule for matrix calculations**](#alignment-rule-for-matrix-calculations)
  - [**Uncertainty calculation on the response: σrespND**](#uncertainty-calculation-on-the-response--σrespnd)
  - [**Assimilation of experimental data and bias calculation via GLLSM**](#assimilation-of-experimental-data-and-bias-calculation-via-gllsm)
    - [**Pre-sorting of benchmark cases to assimilate**](#pre-sorting-of-benchmark-cases-to-assimilate)
    - [**Assimilation formulas**](#assimilation-formulas)
  - [**Reflections on theory**](#reflections-on-theory)
- [**DEVELOPER MANUAL**](#developer-manual)
  - [**methods.py**](#methodspy)
    - [**Acquisition and formatting of useful data (Pandas Dataframe)**](#acquisition-and-formatting-of-useful-data-pandas-dataframe)
    - [**Construction of matrices and vectors (Numpy Arrays)**](#construction-of-matrices-and-vectors-numpy-arrays)
    - [**Calculation functions**](#calculation-functions)
    - [**Utility functions**](#utility-functions)
  - [**classes.py**](#classespy)
    - [**Case**](#case)
    - [**Uncertainty**](#uncertainty)
    - [**Bias**](#bias)
    - [**Assimilation**](#assimilation)
  - [**plots.py**](#plotspy)
  - [**errors.py**](#errorspy)


## **⚠️ Points of attention ⚠️**
(this chapter should be read with the theory chapter and the commands in USER_EXAMPLES.ipynb)
- Some particular isotopes present in sensitivity files are not present in the variance-covariance matrices, but TSURFER (SCALE) uses the variance-covariance data of the closest isotope. CALINS detects these cases if the hundreds digit of the ID is >= 3, examples:
  - be-9 bounded (ID=4309 or 4509) --> associated with covariances of be-9 (ID=4009)
  - h-1 bounded (ID=1901 or 1801) --> associated with covariances of h-1 (ID=1001)
- Isotopes with a 'natural' evaluated form, example:
  - c-12 (bounded?) (ID=6312) --> associated with covariances of c-0 (ID=6000)
- Some covariance matrices have data for the CAPTURE reaction (ID=101), others for the N,GAMMA reaction (ID=102) very close to capture. Be careful when associating these.
- *sdf* files from MCNP calculations can contain reactions with a negative ID:
  - CALINS associates reaction -2 to reaction 101 if the *mcnp* flag is enabled when creating the Case object. Other negative reactions simply trigger a warning.
- If the *sdf* file has multiple sensitivity profiles for the same isotope-reaction pair (same IDs), TSURFER (SCALE) seems to parse profiles differently:
  - *occurrences_rule = first*: the first occurrence is used as the sensitivity profile for that isotope-reaction pair.
  - *occurrences_rule = sum* (default): all occurrences are summed for the profile.
  - *occurrences_rule = last*: the last occurrence is used.

# **INSTALLATION**

Download the GitHub folder. Go to its root (if not compressed), then install the *calins* package:

```
pip install .
```
Or, if installing from the compressed archive directory:
```
pip install archive-name.zip
```

You can verify the package is installed by checking for auto-completion in a Python script.

# **THEORY**

⚠️The term "response" (or "resp" for "response") is used here for the integral (experimental or calculated) response in general.⚠️

## **Introduction**
When a response value is calculated, e.g., a k<sub>eff</sub>, there is a calculation bias Δresp (delta-response) between the simulated model's response value and the observed value.
The estimation method for Δ<sub>resp</sub> presented here is GLLSM (Generalized Linear Least Squares Method). It relies on the assumption that...
The method for calculating the uncertainty on the response calculation σ<sub>resp</sub><sup>ND</sup> due to nuclear data is a matrix computation called the "Sandwich Formula". It consists of...

Note: here "nuclear data" refers to microscopic cross-section data, as an example.
Finally, these formulas are valid for data expressed in relative values.

## **Sensitivity vector**
The sensitivity vector, associated with a unique modeling, consists of a series of sensitivity sub-vectors of resp<sup>calc</sup> to a multi-group cross-section for each isotope-reaction pair.
![](./images/Image_maker/Diapositive1.PNG)
The isotope-reaction pairs linked to the sub-vectors form a list called the "isotope-reaction list". This list and its order are important to respect the alignment rule...

## **Covariance matrix**
The covariance matrix (Cov) is built from the uncertainty data on multi-group cross-section values and the correlation data between these uncertainties.
![](./images/Image_maker/Diapositive2.PNG)
![](./images/Image_maker/Diapositive3.PNG)
The covariance matrix can also be interpreted as a composition of covariance sub-matrices, for each pair (i<sub>h</sub>, r<sub>h</sub>)-(i<sub>v</sub>, r<sub>v</sub>)...
![](./images/Image_maker/Diapositive4.PNG)

## **Alignment rule for matrix calculations**
To ensure each matrix multiplication is valid, the alignment of terms in the matrices/vectors involved must be respected. Thus, the linear application...

## **Uncertainty calculation on the response: σ<sub>resp</sub><sup>ND</sup>**
The sandwich formula allows the propagation of response sensitivities to nuclear data through the covariance matrix, to evaluate...

**Sandwich Formula:**

$$ {σ_{resp}}^{ND} = \sqrt{S \space · \space Cov \space · \space S^t} $$

$S$: sensitivity vector  $Cov$: covariance matrix

## **Assimilation of experimental data and bias calculation via GLLSM**
The GLLSM method consists in assimilating experimental simulation data (also called "benchmark" cases) whose response has been measured and calculated.

---
*SCALE-6.3.1 Manual*:
The GLLS approach accounts for potential variations in nuclear data and measured integral responses that minimize the differences between measured and calculated integral responses for a set of...

[^1]: R. N. Hwang. Topics in data adjustment theory and applications. In Proceedings of the Specialists' Meeting...

(...)

The GLLS procedure can modify the calculated value for the application if it is "similar" to certain experimental responses. In this case, the application's response...

---

This assimilation method yields a vector of global variations in nuclear data (from the assimilated benchmarks) $\Delta\mu_{XS}$, which can predict the bias...

The validity domain of this method is defined by several hypotheses:
- The total calculation bias is mostly made up of Δ<sub>resp</sub><sup>ND</sup>. This is verified especially for pointwise Monte Carlo calculations;
- The response variation caused by small variations in nuclear data is linear;
- The case under study is similar in sensitivity to the assimilated benchmarks (similarity is discussed).

### **Pre-sorting of benchmark cases to assimilate**
Assimilated benchmarks have a greater impact on the reduction of a case's a posteriori uncertainty σ<sub>resp</sub><sup>ND post</sup> the more similar they are to the case. This...

It is possible to pre-sort experiments to assimilate based on these physical criteria, i.e., similarities between benchmarks and the case. For example, classification...

Similarity indicators can be defined between two sensitivity vectors. They help pre-sort benchmarks, which, once assimilated, will best...

These indicators should be calculated for each benchmark to add to the assimilation list.

**Similarity index $E$ (value between 0 and 1):**
This index is a normalized dot product between a benchmark and the case.
$$ E = {{S_1} \space · \space {S_2}^t \over ||S_1|| \space · \space ||S_2||} $$

**Similarity index $C_k$ (value between 0 and 1):**
This index weights the dot product with covariance values.
$$ C_k = \sqrt{({{S_1} \space · \space Cov \space · \space {S_2}^t})^2 \over ({S_1} \space · \space Cov \space · \space {S_1}^t)({S_2} \space · \space Cov \space · \space {S_2}^t)} $$

**Overlap index G (value between 0 and 1):**
This index calculates the overlap rate of each sensitivity corresponding to an energy group, isotope, and reaction. The rate is calculated by taking a case as reference...

$$  SS = 1 - {\sum_{g}^{}\sum_{i, r}^{} \left\{\begin{matrix}
if \space \space {S_{ref}^{g, i, r}}\times {S_{comp}^{g, i, r}} > 0 \space and \space \left|{S_{ref}^{g, i, r}}\right|\geq \left|{S_{comp}^{g, i, r}}\right| : {S_{ref}^{g, i, r}} - {S_{comp}^{g, i, r}} \\
else \space : {S_{ref}^{g, i, r}}
\end{matrix}\right. \over 
\sum_{g}^{}\sum_{i, r}^{} {S_{ref}^{g, i, r}}  }$$

${S_{ref}}$: reference sensitivity vector  ${S_{comp}}$: comparison sensitivity vector

**Shared Sensitivity Overlap Index SS (formula by Mariya BROVCHENKO):**
This index calculates the overlap rate for each sensitivity for an energy group, isotope, and reaction, using the case as reference...

$$ SS =  {\sum_{g}^{}\sum_{i, r}^{} \left\{\begin{matrix}
if \space \space {S_{ref}^{g, i, r}}\times {S_{comp}^{g, i, r}} > 0 \space : min(\left|{S_{ref}^{g, i, r}}\right|, \left|{S_{comp}^{g, i, r}}\right|) \\
else \space : 0
\end{matrix}\right. \over 
\sum_{g}^{}\sum_{i, r}^{} \left| {S_{ref}^{g, i, r}} \right| }$$

${S_{ref}}$: reference sensitivity vector  ${S_{comp}}$: comparison sensitivity vector

**Sorting by the value of $\chi ^2$**

*To be written: see Chapter 6.8.4.2.1 of the SCALE 6.3.1 User Manual [^2]*
[^2]: https://scale-manual.ornl.gov/tsurfer.html#consistency-relations-and-chi-square-filtering

### **Assimilation formulas**
**Description of data used for GLLSM**:
- a "C-E benchmarks vector" $\Delta k_{C-E \space bench}$ containing the relative response differences ((C-E)/C) of selected benchmarks;
![](./images/Image_maker/Diapositive5.PNG)
- uncertainties σ<sub>resp</sub><sup>exp</sup>, to build the "experimental uncertainty matrix" $C_{bench}$. This matrix can include experimental correlations...
![](./images/Image_maker/Diapositive6.PNG)
- the "benchmark sensitivity matrix" $S_{bench} = ..|S_{bench}^j | S_{bench}^{j+1} | ..$ built from the sensitivity vectors $S_{bench}^j$ (j: benchmark index)
![](./images/Image_maker/Diapositive7.PNG)

⚠️Note that the $Cov$ matrix depends on the benchmarks used. It is defined by the union of isotope-reaction pairs present in the benchmarks and the case under study.

***

**Global nuclear data variation vector (from assimilated benchmarks) $\Delta\mu_{XS}$:**
$$ \Delta\mu_{XS} = - Cov \space · \space {S_{bench}}^t \space · \space  (C_{bench}+S_{bench} \space · \space Cov \space · \space {S_{bench}}^t)^{-1} \space · \space \Delta k_{C/E ~ bench} $$

**Adjusted covariance matrix $Cov'$:**
$$ Cov' = Cov - \Delta Cov_{assim} $$
$$ Cov' = Cov - Cov \space · \space {S_{bench}}^t \space · \space  (C_{bench}+S_{bench} \space · \space Cov \space · \space {S_{bench}}^t)^{-1} \space · \space  S_{bench} \space · \space Cov $$

The necessary inversion is not always mathematically possible. The impact of inversion methods leading to "pseudo-inverse" terms has not yet been studied. To verify...

***

**A posteriori bias Δresp<sup>post</sup>**, deviation $(k_{eff}^{expe}-k_{eff}^{calc}) \over k_{eff}^{calc}$ (unit: %, value relative to the calculated response):

$$ \Delta {resp}^{post} = {S_{cas}}  \space · \space  \Delta\mu_{XS} $$

**A priori uncertainty** (unit: %, value relative to the calculated response):

$$ {σ_{resp}}^{ND~ prior} = \sqrt{{S_{cas}} \space · \space Cov \space · \space {S_{cas}}^t} $$

**A posteriori uncertainty** (unit: %, value relative to the calculated response):

$$ {σ_{resp}}^{ND~ post} = \sqrt{{S_{cas}} \space · \space Cov' \space · \space {S_{cas}}^t} $$

***

## **Reflections on experimental correlations**
According to the above methodology, setting experimental correlations to 0 has several consequences. It has been observed that assimilation including several...

This raises questions about the similarity of assimilated benchmarks when using the GLLS method. If several experiments have intrinsically similar sensitivities and responses...

Studies on the impact of experimental correlations have already been conducted [^3].
[^3]: T. Nicol, C. Carmouze. Impact of experimental correlation on transposition method carry out with critical integral experiments. ICNC 2019...

# **DEVELOPER MANUAL**
Matrix multiplications are performed with the Python operator "A @ B", fully equivalent to numpy.matmul.

⚠️Note: operator precedence is very important when it comes to algorithmic matrix multiplication methods, to optimize calculation time.⚠️

---

## **methods.py**
### **Acquisition and formatting of useful data (Pandas Dataframe)**
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

***Note***: all isotopes and reactions encountered are stored in this dataframe, including the TOTAL reaction - ID: 1. However, a filter on the TOTAL reaction is applied...

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
- GENDF (*work in progress*)

**WARNING**: The **CALINS norm** is to format sensitivity and covariance data in **decreasing** energy group order from 20 MeV -> 0 MeV.

### **Construction of matrices and vectors (Numpy Arrays)**
Related functions: *make_ ... (...)*

These functions allow building vectors and matrices as *Numpy Arrays* while respecting the **alignment rule** (see chapter *Alignment rule for matrix calculations*).

To ensure alignment between sensitivity vectors and covariance matrices, they are built simultaneously based on an isotope-reaction list...

This list is first built as the intersection or union of each sensitivity vector's isotope-reaction lists, depending on the operations to be performed...

**NOTE**: if non-conventional isotopes like "h-poly" ID=1901 are present in the sensitivity vectors' isotope-reaction list and not in the covariance DataFrame...

The same principle applies for isotopes with a "natural" form (ID=xx000) in the variance-covariance data. If a sensitivity vector isotope isn't in the covariance matrix...

***Note:*** A filter excluding the TOTAL reaction is applied at this point on the isotope-reaction list.

For sensitivity vectors, this list is traversed in order, and for each isotope-reaction pair, data is searched for in its DataFrame...

For the matrix, it is initialized as a zero matrix, 2D dimensions *[(N_group x Nb_iso-reac), (N_group x Nb_iso-reac)]*. The isotope-reaction list is traversed twice...

### **Calculation functions**
Functions of type: *calcul_... (...)*

These functions use the *make_...* functions to build sensitivity vectors and the covariance matrix based on a shared isotope-reaction list, to perform...

One function exists per similarity index. Implemented indices are:
- E
- C<sub>k</sub>
- SS overlap index (CEA formula)
- SS overlap index (Mariya BROVCHENKO formula)

These functions take as input at least two sensitivity vectors and a covariance matrix for C<sub>k</sub> calculations. Sensitivity vectors can be paths to files...

A third function allows uncertainty calculation, with a sensitivity vector and covariance matrix as input. The sensitivity vector can be a file path...

### **Utility functions**
The *methods* module also contains utility functions to convert isotope identification numbers to strings and vice versa.

## **classes.py**

### **Case**
The *case* object enables building a study case (or benchmark case), from a *.sdf* sensitivity file, and storing relevant data on the case. Using a *case* object allows...

- *case.sdf_path*: path to the case;
- *case.casename*: case name (set as the base of the *.sdf* file path);
- *case.group_nb*: number of energy groups used;
- *case.e_bins*: bounds of the energy groups used;
- *case.iso_reac_list*: isotope-reaction list available in the sensitivity file, as isotope and reaction ID numbers;
- *case.sensitivities*: DataFrame of the case's sensitivities, obtained via *format_sensi_to_dataframe(...)* (*methods.py*)
- *case.resp_calc*: calculated/modelled response value;
- *case.sigma_resp_calc*: uncertainty on the calculated/modelled response;
- *case.resp_expe*: experimental response value (for benchmarks);
- *case.sigma_resp_expe*: uncertainty on the experimental response (for benchmarks);

Finally, the *case* object has a function to display and save as *.html* a histogram of integral sensitivities associated with the sensitivity profiles for each...
