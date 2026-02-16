# Developer Manual - plots.py

The functions use the *Plotly* visualization package (v5.14), especially its *plotly.express* module. This package is ergonomic and explicit, allowing manipulation of each graph as a batch of HTML data, enabling the combination of multiple graphs in a single *.html* file if needed.

*Plots* functions are intended for use within objects but can also be used independently by users. Several functions are provided for different types of graphs:
- histogram of sensitivities for an application case, for each isotope and isotopic reaction;
- histogram of integral covariances of a matrix for each isotope and isotopic reaction;
- 2D colormap of a covariance sub-matrix for a given pair of (iso-reac)<sup>h</sup> - (iso-reac)<sup>v</sup>;
- HTML table generated from data in list format.