# omicwas

Cell-Type-Specific Disease Association Testing in Bulk Omics Experiments

## Installation in R

In order to install the most recent version:

    install.packages("devtools")
    devtools::install_github("fumi-github/omicwas", build_vignettes = TRUE)

If you encounter dependency error for `sva` package,
install it from Bioconductor:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva")

To install from CRAN archive (possibly a version older than github):

    install.packages("omicwas")

To uninstall package:

    remove.packages("omicwas")

## Usage

    library(omicwas)
    vignette("intro", package = "omicwas")

## Information

[CSHL (14 Nov 2019) poster](http://103.253.147.127/PUBLICATIONS/191114cshl.pdf)
