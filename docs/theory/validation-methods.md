# Validation Methods for the Upper Safety Limit (USL)

After performing data assimilation via GLLSM and obtaining a bias estimate, CALINS provides four independent methods to compute an **Upper Safety Limit (USL)** for the application case response. The USL are determined based on the assumption that the responses are the multiplication factor $k_{eff}$. The parametric and nonparametric methods are based on the formalism described by Kiedrowski et al. (NSE 181, 2015). The trending method implements the NUREG/CR-6698 single-sided lower tolerance band (also available in SCALE/VADER).

All three methods share the same general framework:

$$\text{USL} = 1 - \text{CM} - \text{MOS}$$

where the **Calculational Margin** (CM) accounts for bias, bias uncertainty, and an optional administrative margin of subcriticality (MOS) set to 0.05 by default.

---

## Common Definitions

For each non-removed benchmark $i$, the **scaled multiplication factor** is:

$$\tilde{k}_i = \frac{k_{\text{calc},i}}{k_{\text{bench},i}}$$

The **bias** of benchmark $i$ is $\beta_i = \tilde{k}_i - 1$, and the **combined uncertainty** is:

$$\sigma_i = \sqrt{\sigma_{\text{expe},i}^2 + \sigma_{\text{calc},i}^2}$$

The **nonconservative bias adjustment** common to all methods is:

$$\Delta_m = \max(0,\;\bar{\beta})$$

This ensures that if the mean bias is positive (calculation overestimates the response, i.e. is non-conservative), the adjustment increases the margin.

---

## 1. Parametric Method

### Hypotheses

- The scaled multiplication factors $\tilde{k}_i$ follow a **normal distribution**. A Shapiro-Wilk test is performed; if it fails ($p < 0.05$), the parametric method should not be trusted.
- All benchmarks are treated **equally** in weighting (weighted by inverse variance).
- The population of possible application responses is assumed Gaussian.

### Procedure

**Step 1 — Inverse-variance weights:**

$$w_i = \frac{1}{\sigma_i^2}, \quad W = \sum_i w_i$$

**Step 2 — Weighted mean:**

$$\bar{k} = \frac{\sum_i w_i \, \tilde{k}_i}{W}$$

**Step 3 — Bias:**

$$\beta = \bar{k} - 1$$

**Step 4 — Pooled standard deviation**:

The weighted sample variance is:

$$s_k^2 = \frac{1}{W} \cdot \frac{N}{N-1} \sum_i w_i (\tilde{k}_i - \bar{k})^2$$

The average variance of the mean is:

$$\bar{\sigma}_k^2 = \frac{N}{W}$$

The pooled standard deviation of the bias is:

$$\sigma_\beta = \sqrt{s_k^2 + \bar{\sigma}_k^2}$$

**Step 5 — Tolerance factor $\kappa$**:

Using the noncentral $t$-distribution with $\nu = N - 1$ degrees of freedom and noncentrality parameter $\delta = z_p \sqrt{N}$ (where $z_p$ is the $p$-quantile of the standard normal):

$$\kappa = \frac{t_{q,\nu,\delta}}{\sqrt{N}}$$

This factor guarantees that a fraction $p$ of the population is bounded at confidence level $q$ (defaults: $p = q = 0.99$).

**Step 6 — Calculational Margin and USL:**

$$\text{CM} = -\beta + \kappa \, \sigma_\beta + \Delta_m$$

$$\text{USL} = 1 - \text{CM} - \text{MOS}$$

### Summary

| Quantity | Formula |
|---|---|
| Bias | $\beta = \bar{k} - 1$ |
| Bias uncertainty | $\kappa \cdot \sigma_\beta$ |
| CM | $-\beta + \kappa \sigma_\beta + \Delta_m$ |
| Key hypothesis | **Normality** of $\tilde{k}_i$ |

---

## 2. Nonparametric Rank-Order Method

### Hypotheses

- **No distributional assumption** on $\tilde{k}_i$ (no normality required).
- The bias is estimated from the **worst case** (most negative $\tilde{k}_i$), making this method inherently conservative.
- A lookup table provides an additional margin $m_{NP}$ depending on the nonparametric confidence level $C_{NP}$.

### Procedure

**Step 1 — Worst-case bias**:

$$\beta = \min_i(\tilde{k}_i) - 1$$

This takes the most pessimistic benchmark as the bias estimate.

**Step 2 — Bias uncertainty:**

The uncertainty is based on the combined uncertainty of the worst-case benchmark, multiplied by a confidence factor $n_\sigma$ (default 2.6, corresponding to approximately 99% confidence for a normal distribution):

$$\sigma_\beta = n_\sigma \cdot \sigma_{\text{worst}}$$

**Step 3 — Nonparametric confidence level**:

$$C_{NP} = 1 - p_{\text{pop}}^N$$

This represents the probability that at least one of $N$ benchmarks exceeds the $p_{\text{pop}}$-quantile of the population. For example, with $N = 59$ benchmarks and $p_{\text{pop}} = 0.95$:

$$C_{NP} = 1 - 0.95^{59} \approx 0.953$$

**Step 4 — Nonparametric additional margin $m_{NP}$:**

Based on $C_{NP}$, a lookup table (Table I in reference) provides:

| $C_{NP}$ range | $m_{NP}$ |
|---|---|
| > 0.90 | 0.00 |
| 0.80 – 0.90 | 0.01 |
| 0.70 – 0.80 | 0.02 |
| 0.60 – 0.70 | 0.03 |
| 0.50 – 0.60 | 0.04 |
| 0.40 – 0.50 | 0.05 |
| ≤ 0.40 | Additional data needed |

If $C_{NP} \leq 0.4$, there are not enough benchmarks and the method cannot provide a result.

**Step 5 — Calculational Margin and USL:**

$$\text{CM} = -\beta + \sigma_\beta + \Delta_m + m_{NP}$$

$$\text{USL} = 1 - \text{CM} - \text{MOS}$$

### Summary

| Quantity | Formula |
|---|---|
| Bias | $\beta = \min(\tilde{k}_i) - 1$ (worst case) |
| Bias uncertainty | $n_\sigma \cdot \sigma_{\text{worst}}$ |
| CM | $-\beta + \sigma_\beta + \Delta_m + m_{NP}$ |
| Key hypothesis | **None** on distribution (conservative worst-case) |

---

---

## 3. Trending Method (Single-Sided Lower Tolerance Band)

### Hypotheses

- A **linear trend** of $\tilde{k}_i$ as a function of the trending parameter $x_i$ (here $x_i = Ck_i$, the similarity coefficient) is assumed.
- A **t-test for trend significance** is performed to determine whether the slope $\beta_1$ is statistically different from zero.
- Uses the NUREG/CR-6698 single-sided lower tolerance band formalism (also implemented in SCALE/VADER).
- All statistics use **inverse-variance weighting** (same as the parametric method).

### Procedure

**Step 1 — Inverse-variance weights and weighted means:**

$$w_i = \frac{1}{\sigma_i^2}, \quad W = \sum_i w_i$$

$$\bar{k} = \frac{\sum_i w_i \tilde{k}_i}{W}, \quad \bar{x} = \frac{\sum_i w_i x_i}{W}$$

**Step 2 — Weighted sums of squares:**

$$S_{xx} = \frac{N}{W} \sum_i w_i (x_i - \bar{x})^2$$

$$S_{kx} = \frac{N}{W} \sum_i w_i (x_i - \bar{x})(\tilde{k}_i - \bar{k})$$

**Step 3 — Linear fit coefficients:**

$$\beta_1 = \frac{S_{kx}}{S_{xx}}, \quad \beta_0 = \bar{k} - \beta_1 \bar{x}$$

The fitted value at any $x$ is: $k_{fit}(x) = \beta_0 + \beta_1 x$

**Step 4 — Variance of the fit and pooled standard deviation:**

$$\sigma_{fit}^2 = \frac{N}{(N-2) W} \sum_i w_i \left(\tilde{k}_i - k_{fit}(x_i)\right)^2$$

$$\bar{\sigma}^2 = \frac{N}{W}$$

$$S_p = \sqrt{\sigma_{fit}^2 + \bar{\sigma}^2}$$

**Step 5 — t-test for trend significance:**

$$t_{fit} = \frac{|\beta_1|}{\sigma_{fit} / \sqrt{S_{xx}}}$$

The null hypothesis $H_0: \beta_1 = 0$ is rejected (trend is significant) if:

$$t_{fit} > t_{1-\alpha/2,\; N-2}$$

where $\alpha = 1 - \text{confidence}$ (default confidence = 0.95).

**Step 6 — Evaluation point:**

The USL is evaluated at the maximum Ck value: $x_{eval} = \max(x_i)$, corresponding to the most similar benchmark.

**Step 7 — Bias at evaluation point:**

$$\beta(x_{eval}) = k_{fit}(x_{eval}) - 1$$

$$\Delta_m = \max(0, \beta(x_{eval}))$$

**Step 8 — Bias uncertainty (NUREG/CR-6698 tolerance band):**

$$\sigma_\beta(x) = S_p \left( \sqrt{2 F_{\alpha}(2, N-2) \left[\frac{1}{N} + \frac{(x - \bar{x}_{unw})^2}{S_{xx}}\right]} + \frac{z_p \sqrt{N-2}}{\sqrt{\chi^2_{1-\alpha, N-2}}} \right)$$

where:

- $F_{\alpha}(2, N-2)$ is the F-distribution quantile at confidence level $\alpha$
- $z_p$ is the standard normal quantile at proportion $p$ (default 0.95)
- $\chi^2_{1-\alpha, N-2}$ is the chi-squared quantile
- $\bar{x}_{unw}$ is the unweighted (arithmetic) mean of $x_i$

**Step 9 — Calculational Margin and USL:**

$$\text{CM} = -\beta(x_{eval}) + \sigma_\beta(x_{eval}) + \Delta_m$$

$$\text{USL} = 1 - \text{CM} - \text{MOS}$$

### Summary

| Quantity | Formula |
|---|---|
| Bias | $\beta(x) = k_{fit}(x) - 1$ (from linear regression) |
| Bias uncertainty | $\sigma_\beta(x)$ (NUREG/CR-6698 tolerance band) |
| CM | $-\beta(x) + \sigma_\beta(x) + \Delta_m$ |
| Key hypothesis | **Linear trend** of $\tilde{k}_i$ vs $Ck_i$ |
| Diagnostic | t-test for trend significance |

---

!!! note "Practical recommendations"
    - Use the **Parametric** method when Shapiro-Wilk passes and benchmarks are well-behaved (no significant trend).
    - Use the **Nonparametric** method as a conservative fallback when normality is not satisfied.
    - Use the **Trending** method when a significant linear trend is detected (t-test rejects $H_0$). It accounts for the dependence of bias on similarity.
    - If the trend is **not significant**, the parametric or nonparametric methods are more appropriate.

## References

- **B. C. Kiedrowski, F. B. Brown, A. C. Kahler, S. R. Bolding, M. E. Rising, and J. A. Favorite**, *Whisper: Sensitivity/Uncertainty-Based Computational Methods and Software for Determining Baseline Upper Subcritical Limits*, **Nuclear Science and Engineering, 181**, pp. 17-47 (2015).
- **J. C. Dean and R. W. Tayloe Jr.**, *Guide for Validation of Nuclear Criticality Safety Calculational Methodology*, **NUREG/CR-6698**, prepared for the US Nuclear Regulatory Commission by Science Applications International Corporation, Oak Ridge, TN, January 2001.
- **J. B. Clarity, W. J. Marshall, D. E. Mueller, S. S. Powers, B. T. Rearden, S. M. Bowman, and A. Barto**, *Determination of Bias and Bias Uncertainty for Criticality Safety Computational Methods*, **NUREG/CR-7311**, Oak Ridge National Laboratory, 2025.
