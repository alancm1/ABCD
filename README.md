# ABCD

This is a package for the change-point method ABCD (Adaptive Block-Based Change-point Detection), a flexible method for detecting spatially clustered change-points in high-dimensional time series.

## Description

ABCD is a method which uses an adaptive blocking strategy to capture sparse, spatially clustered change-points at multiple scales in a non-parametric fashion.  Distribution-free graph-based two-sample tests from Chu and Chen (2019) are used within this blocking framework to detect general distributional changes. ABCD controls Type 1 error via a permutation p-value. This package implements ABCD for high-dimensional and image time series, with options for automatic generation of blocking structures.

## Where to Start

The ABCD function can be used to detect sparse, spatially clustered change-points in an offline setting where dimensions are ordered on a grid in 1-dimensional Euclidean space. 

The ABCD_img function does the same for image data, that is, when dimensions are ordered on a grid in 2-dimensional Euclidean space. 

These are the main two functions implementing the ABCD method and should be your first port of call. Upon installing running ?ABCD or ?ABCD_img will provide more detailed documentation on appropriate arguments and basic examples of their usage. Further details on the ABCD method can be found in the corresponding paper.

### Dependencies

The ABCD R package requires the R package ade4 to be installed, which is used to construct k-MSTs to be used as similarity graphs within the method. 

### Installing

To install the program, install devtools and then run the following R code to install ABCD from Github:

```
install_github("alancm/ABCD")
```

### Reproducing Simulations

Within the folder sim_app_files, one may find .R scripts and R Markdown files which can be used to reproduce simulations found in sections 4 and 5 of the associated paper.

The file "abcd_sim_helpers_s.R" contains helper functions to run the simulations found in section 4.1. Additional dependencies for this R script include the R packages MASS, InspectChangepoint, factorcpt (found here: https://github.com/markov10000/factorcpt) and ggplot2.

The file "image_sim_helpers_s.R" contains helper functions to run the simulations found in section 4.2. Dependencies for this R script include the R packages MASS, RSatscan and L2HDChange, along with having the SaTScan software installed on your computer (a Windows machine is recommended for this).

The R Markdown file "simulation_testing_s.Rmd" includes code chunks to generate null Monte Carlo test statistics for each data type in section 4.1, as well as to run the full simulations found in Section 4.1 and 4.2 of the paper. This file requires the helper scripts listed above to be loaded.

The R Markdown file "data_app_clean.Rmd" includes code analyzing for change-points in the Natanz nuclear facility data application. Many external packages are required to actually run this code, along with the substantial remote sensing and construction polygon data.
