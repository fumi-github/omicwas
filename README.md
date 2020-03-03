# omicwas

Cell-Type-Specific Disease Association Testing in Bulk Omics Experiments

## Installation in R

Run the following to install from CRAN archive:

    install.packages("omicwas")

In order to install the most recent version,
possibly not yet in CRAN:

    install.packages("devtools")
    devtools::install_github("fumi-github/omicwas")

If you encounter dependency error for `sva` package,
install it from Bioconductor:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva")

## Usage

    library(omicwas)
    ?ctassoc

## Information

[CSHL (14 Nov 2019) poster](http://103.253.147.127/PUBLICATIONS/191114cshl.pdf)
