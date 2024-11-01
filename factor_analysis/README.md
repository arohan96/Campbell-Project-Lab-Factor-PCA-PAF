# Factor Analysis
This package contains MATLAB scripts that implement Principal Component Analysis and Principal Axis Factoring to a set of market returns in order to identify a subset of statistical factors that can be used for summarizing the returns of the given market. The package is broken down into the following modules:<br/>
### 1. pcaTestCase.m
This is the main entry point to the code. The script generates a synthetic return series and factor data and then calls **factorDecomposition.m** to apply either Principal Component Analysis or Principal Axis Factoring to the factor returns data depending on the parameter values provided in the script.
### 2. factorDecomposition.m
The driver script for dimensionality reduction of factor data. Computes correlation matrix and calls the **eigenValueDecomposition.m** applies either singular value decomposition or principal axis factoring to compute factor loadings. These factor loadings are, in turn, used to compute Portfolio Betas, Factor Returns, and Factor Volatilities. Also calls the **visualization.m** and **factorRotation.m** scripts (optional) for generating plots and applying factor rotation.
### 3. eigenValueDecomposition.m
Implements Singular Value Decomposition and Principal Axis Factoring to compute eigenvalues and eigenvectors
### 4. Visualization.m
Utility functions for creating plots
### 5. factorRotation.m
Utility functions for applying factor rotation
