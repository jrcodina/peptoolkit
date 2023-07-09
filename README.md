# PepToolkit

The `peptoolkit` R package provides a function to generate properties for peptide sequences for Principal Component Analysis (PCA).

## Installation

You can install the released version of peptoolkit from [CRAN](https://CRAN.R-project.org) with:

```r
install.packages("peptoolkit")
```

You can also install the development version from GitHub with:

```r
# install.packages("devtools") # Uncomment and run if you don't have the devtools package yet
devtools::install_github("jrcodina/peptoolkit")
```
# Example

This is a basic example which shows you how to use the main function:

```r
# Default usage
result <- peptoolkit::extract_features_QSAR(n = 3)

# Providing a custom peptide list
result <- peptoolkit::extract_features_QSAR(n = 3, custom.list = TRUE, PeList = c('ACA', 'ADE'))
```

Please refer to function documentation for more details on parameters and their usage.

## Citation

If you use `peptoolkit` in your research, please cite:

To cite package ‘peptoolkit’ in publications use:

  Codina J (2023). _peptoolkit: A Toolkit for Using Peptide
  Sequences in Machine Learning and Accelerate Virtual
  Screening_. R package version 0.0.1.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {peptoolkit: A Toolkit for Using Peptide Sequences in Machine Learning and
Accelerate Virtual Screening},
    author = {Josep-Ramon Codina},
    year = {2023},
    note = {R package version 0.0.1},
  }
