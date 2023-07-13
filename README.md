# PepToolkit

The `peptoolkit` R package is designed for the manipulation and analysis of peptides. The package provides a range of functionalities aimed at assisting researchers in peptide engineering and proteomics. The package allows users to manipulate peptides by adding amino acids at every position, count the occurrences of each amino acid at each position, and transform amino acid counts based on probabilities. Additionally, the package offers functionalities to select the best versus the worst peptides and further analyze these peptides. This includes counting specific residues, reducing peptide sequences, extracting features through One Hot Encoding (OHE), and utilizing Quantitative Structure-Activity Relationship (QSAR) properties. This package is intended for both researchers and bioinformatics enthusiasts who are working on peptide-based projects, specially for their use with machine learning.

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

  Codina J (2023). _peptoolkit: A Toolkit for Using Peptide
  Sequences in Machine Learning and Accelerate Virtual
  Screening_. R package version 0.0.1.

A BibTeX entry for LaTeX users is

```
  @Manual{,
    title = {peptoolkit: A Toolkit for Using Peptide Sequences in Machine Learning and Accelerate Virtual Screening},
    author = {Josep-Ramon Codina},
    year = {2023},
    note = {R package version 0.0.1},
  }
```
