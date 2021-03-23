# omicwas

Cell-Type-Specific Disease Association Testing in Bulk Omics Experiments

## Installation in R

In order to install the most recent version:

    install.packages("devtools")
    devtools::install_github("fumi-github/omicwas", build_vignettes = TRUE)

To install from CRAN archive (possibly a version older than github):

    install.packages("omicwas")

To uninstall package:

    remove.packages("omicwas")

## Usage

    library(omicwas)
    vignette("intro", package = "omicwas")

## Information

Please cite

[Takeuchi, F., Kato, N. Nonlinear ridge regression improves cell-type-specific differential expression analysis. BMC Bioinformatics 22, 141 (2021)](https://doi.org/10.1186/s12859-021-03982-3)
