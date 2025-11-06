# Theory - Introduction

⚠️ The term *response* (or *resp*) is used here to denote the integral response (experimental or calculated). ⚠️

When a response value is calculated, e.g., a k<sub>eff</sub>, there is a calculation bias Δresp between the simulated model's response value and the experimental observed value. This bias is primarily composed of bias due to nuclear data Δ<sub>resp</sub><sup>ND</sup> and bias due to approximations in calculation schemes Δ<sub>resp</sub><sup>SC</sup>. The methods described here allow estimating the calculation bias and its associated uncertainty.

The estimation method Δ<sub>resp</sub> used is the GLLSM (Generalized Linear Least Squares Method). It assumes that Δ<sub>resp</sub><sup>ND</sup> dominates the total bias (valid e.g. for Monte Carlo). It also assumes that the response variation due to small nuclear data changes is linear.  

The uncertainty σ<sub>resp</sub><sup>ND</sup> due to nuclear data is computed using the propagation matrix product formula (i.e. "sandwich formula"), which propagates sensitivities (of the response to ND) through a variances-covariances of ND matrix.

Note: here "nuclear data" refers to microscopic cross-section data, as an example. Other data such as angular distributions could be included.
All formulas are valid for data expressed in relative values.