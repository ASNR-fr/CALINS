# Validation Methods for the Upper Safety Limit (USL)

After performing data assimilation via GLLSM and obtaining a bias estimate, CALINS provides three independent methods to compute an **Upper Safety Limit (USL)** for the application case response. The USL are determined based on the assumption that the responses are the multiplication factor $k_{eff}$. These methods are based on the formalism described by Kiedrowski et al. (NSE 181, 2015).

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

!!! note "Practical recommendations"
    - Use the **Parametric** method when Shapiro-Wilk passes and benchmarks are well-behaved.
    - Use the **Nonparametric** method as a conservative fallback when normality is not satisfied.

## References

- **B. C. Kiedrowski, F. B. Brown, A. C. Kahler, S. R. Bolding, M. E. Rising, and J. A. Favorite**, *Whisper: Sensitivity/Uncertainty-Based Computational Methods and Software for Determining Baseline Upper Subcritical Limits*, **Nuclear Science and Engineering, 181**, pp. 17-47 (2015).
