# Developer Manual - classes.py

## Case
The *Case* object allows building a study case (or benchmark case) from a sensitivity file (*.sdf*) and storing the relevant data associated with the Case. Using an object makes this information easily accessible. The available attributes are as follows:
- *Case.sdf_path*: path to the case;
- *Case.casename*: name of the study case (initialized as the base name of the *.sdf* file path);
- *Case.group_nb*: number of energy groups in the used mesh;
- *Case.e_bins*: energy group boundaries of the used mesh;
- *Case.iso_reac_list*: list of isotope-reaction pairs for which data are available in the sensitivity file, in the form of isotope and reaction ID numbers;
- *Case.sensitivities*: DataFrame containing the sensitivities of the case, obtained using the function *format_sensi_to_dataframe(...)* (*methods.py*);
- *Case.resp_calc*: calculated/modeled response value;
- *Case.sigma_resp_calc*: uncertainty on the calculated/modeled response;
- *Case.resp_expe*: experimental response value (for benchmark cases);
- *Case.sigma_resp_expe*: uncertainty on the experimental response (for benchmark cases).

Finally, the *Case* object has a function to display and save, in *.html* format, a histogram of integral sensitivities along with the sensitivity profiles for each isotope and isotopic reaction.

## NDCovariances
The *NDCovariances* object provides a unified interface for loading and managing nuclear data covariance matrices from various file formats. This is the **recommended method** for handling covariance data in CALINS. 

The available attributes are:
- *NDCovariances.input_path*: path to the covariance file;
- *NDCovariances.format*: detected or specified format ('coverx', 'coverx_text', 'comac', 'gendf', 'xlsx', or 'auto');
- *NDCovariances.cov_dataf*: DataFrame containing the covariance data in standardized format;
- *NDCovariances.e_bins*: energy group boundaries;
- *NDCovariances.group_nb*: number of energy groups;
- *NDCovariances.header*: original file header (if available);
- *NDCovariances.iso_reac_list*: list of (isotope, reaction) pairs present in the covariance data, automatically extracted and filtered to exclude energy bin entries.

The object supports automatic format detection (using `format='auto'`) and provides `write_xlsx()` and `write_txt()` methods for exporting covariance data to Excel or COVERX text format, enabling full round-trip export/import cycles.

For detailed usage examples, see the [Loading Covariance Data](../usage-examples/loading-covariance-data.md) guide.

## Uncertainty
The *Uncertainty* object allows calculating an (absolute) uncertainty using the Sandwich formula and storing several properties related to this calculation. This object takes as input a sensitivity vector and a covariance matrix, already built as Numpy arrays, as well as the number of energy groups and the isotope-reaction list describing the construction of the vector and matrix. The available attributes are:
- *Uncertainty.resp_calc*: calculated response value of the study case;
- *Uncertainty.value*: calculated absolute uncertainty, in pcm;
- *Uncertainty.group_nb*: number of energy groups in the used mesh;
- *Uncertainty.e_bins*: energy group boundaries of the used mesh;
- *Uncertainty.iso_reac_list*: full isotope-reaction list, in the order of construction of the vector and matrix, given as isotope and reaction ID numbers;
- *Uncertainty.decomposition*: contribution to the relative uncertainty (%) squared, for each isotope-reaction pair, in the form of a DataFrame. This decomposition corresponds to the linear combination step in the final stage of the Sandwich formula squared calculation \((X \cdot S_{\text{cas}})\) as described below (taking covariances into account):

$$(\sigma_{\text{resp}}^{\text{ND}})^2 = S_{\text{cas}} \cdot Cov \cdot S_{\text{cas}}^t = X \cdot S_{\text{cas}}^t$$

with: $X = S_{\text{cas}} \cdot Cov$ (*DIM* : 1D -> I<sup>r</sup> x G)

Decomposition by iso-reac:

$$
d^1 = X^1 \cdot S_{\text{cas}}^1, \quad
d^2 = X^2 \cdot S_{\text{cas}}^2, \quad
d^3 = X^3 \cdot S_{\text{cas}}^3, \quad \ldots
$$

where:

$$
d^i = X^i \cdot S_{\text{cas}}^i
$$

\( S_{\text{cas}}^i \) is the sensitivity sub-vector for isotope-reaction pair *i*.

The *Uncertainty* object has a function to display and/or save several analysis elements in *.html* format (that can be activated with the parameter *output_html_path* in the function *calcul_uncertatiny()*):
- *Uncertainty.export_results()*:
  - the uncertainty;
  - the calculated response value and its calculation-scheme uncertainty;
  - the number of energy groups;
  - the sensitivity vector of the study case;
  - decomposition of uncertainty by isotope-reaction pair;
  - histograms of the covariance matrix integrals per isotope-reaction pair;
  - sub-matrices of covariances for user-selected isotope-reaction pairs.

Finally, the *Uncertainty* object is invoked in the *calcul_uncertainty(...)* function and twice as a property of the *Assimilation* object (for prior and posterior uncertainty).

## Bias
The *Bias* object is similar to the *Uncertainty* object. It allows computing the posterior bias \(\Delta resp^{post}\) using the bias formula (see chapter *THEORY*) and stores the same attributes as the *Uncertainty* object regarding the calculation. This object takes as input a sensitivity vector and the global nuclear data variation vector \(\Delta\mu_{XS}\) (after assimilation), already built as Numpy arrays, as well as the number of energy groups and isotope-reaction list describing the construction of the vectors. The available attributes are:
- *Bias.resp_calc*: calculated response value of the study case;
- *Bias.value*: calculated absolute bias, in pcm;
- *Bias.group_nb*: number of energy groups in the used mesh;
- *Bias.e_bins*: energy group boundaries of the used mesh;
- *Bias.iso_reac_list*: full isotope-reaction list, in the order of vector construction, given as isotope and reaction ID numbers;
- *Bias.decomposition*: contribution to the relative bias (%) for each isotope-reaction pair (as a DataFrame), described below:

$$
\Delta \text{resp}^{\text{post}} = S_{\text{cas}} \cdot \Delta\mu_{XS}
$$

Decomposition by iso-reac:

$$
d^1 = S_{\text{cas}}^1 \cdot \Delta\mu_{XS}^1, \quad
d^2 = S_{\text{cas}}^2 \cdot \Delta\mu_{XS}^2, \quad
d^3 = S_{\text{cas}}^3 \cdot \Delta\mu_{XS}^3, \quad \ldots
$$

where:

$$
d^i = S_{\text{cas}}^i \cdot \Delta\mu_{XS}^i
$$

\( S_{\text{cas}}^i \) is the sensitivity sub-vector for isotope-reaction pair *i*.

Finally, the *Bias* object is only invoked as a property of the *Assimilation* object.

---

## Assimilation
The *Assimilation* object performs GLLSM based on a list of benchmark cases, a study case, and a covariance matrix. Benchmark cases and the study case can be either paths to *SDF* files or *Case* objects. The covariance data must be provided as an *NDCovariances* object or an *Assimilation* object. This class uses functions from *methods*, *plots*, and *errors* to format the data, build vector/matrix objects, performs assimilation through GLLSM, and display useful data as plots.

Initialization of an *Assimilation* object consists of several main steps:
1. Formatting sensitivity data (benchmarks and study case) into a DataFrame;
2. Constructing the benchmark sensitivity matrix, the study case sensitivity vector, and the covariance matrix, following the same isotope-reaction order;
3. Constructing the benchmark response uncertainty matrix and the benchmark C/E vector, following the same benchmark case order;
4. *self.check_dimensions()* & *check_correspondences()*: verifying consistency of all vector/matrix dimensions and generating warnings for high sensitivities lacking variance-covariance data;
5. *self.calcul_prior_uncertainty()*: calculation and creation of an *Uncertainty* object;
6. *self.calcul_matrix_assimilation()*: calculation of \(\Delta \mu_{XS}\), \(\Delta Cov_{assim}\), application of a \(\chi^2\) filter if a target value is provided, and filter application on Ck if a target is provided;
7. *self.calcul_bias()*: calculation and creation of a *Bias* object;
8. *self.calcul_post_uncertainty()*: calculation and creation of an *Uncertainty* object.

The *Assimilation* object has two functions to display and/or save several analysis elements in *.html* format:
- the sensitivity vector of the study case;
- an output file including all assimilation parameters:
  - prior and posterior uncertainties, as well as bias;
  - initial and final \(\chi^2\) values;
  - number of removed benchmark cases;
  - number of energy groups;
  - list of benchmark cases with paths, calculated and experimental response values, ND prior and posterior uncertainty and bias, individual \(\chi^2\), similarity coefficients Ck with the study case, and a red background for excluded cases;
  - sensitivity vector of the study case;
  - decomposition by isotope-reaction pair for the two uncertainties and bias;
  - histograms of covariance matrix integrals per isotope-reaction pair, and of \(\Delta Cov_{assim}\) for comparison;
  - sub-matrices of covariances and \(\Delta Cov_{assim}\) for user-selected isotope-reaction pairs.
  - list of isotope-reaction pair included in benchmark cases, in covariances-data and in calculation for verification;