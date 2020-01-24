# omicwas

Cell-Type-Specific Disease Association Testing in Bulk Omics Experiments

## Installation in R

First install `sva` which is in Bioconductor:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva")

Then install `omicwas`:

    install.packages("devtools")
    devtools::install_github("fumi-github/omicwas")

## Usage

    library(omicwas)
    ?ctassoc

## Information

[CSHL (14 Nov 2019) poster](http://103.253.147.127/PUBLICATIONS/191114cshl.pdf)

