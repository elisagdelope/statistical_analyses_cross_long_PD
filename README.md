# Table of contents
* [Introduction](#introduction)
* [Content](#content)
* [Data](#data)
* [Requirements](#requirements)
* [License](#license)

# Introduction
This repository contains the code for statistical analyses performed in Chapter 3 of my thesis "Cross-sectional and longitudinal profiling of PD transcriptomics and metabolomics".
The project consists on whole blood transcriptomics and blood plasma metabolomics cross-sectional and longitudinal profiling of Parkinson's disease patients and controls from the PPMI cohort and  the LuxPARK cohort respectively, to identify differential molecular and higher-level functional features in PD.

# Content
The code covers the following main tasks and analyses:

1. Loading, preparation and processing of datasets
2. Generate higher-level functional features (mean, median, sd, 1st principal component "pca", pathifier deregulation scores) for GOBP, GOCC, CORUM databases (PPMI), KEGG (LuxPARK)
3. Differential expression analysis for PD vs. control (PPMI) and differential abundance analysis for *de novo* PD vs. control and all PD vs. control (LuxPARK) for single molecules and higher-level functional representations
4. Longitudinal analyses: association with time, correlation, trend analysis over consecutive timepoints, for single molecules and higher-level functional representations
5. Pathway enrichment analysis using gsea, goana, mesh terms (PPMI)
6. Post-processing of results, visualizations (boxplots, trajectories)


There is a README.md inside each directory with corresponding explanation for each script:

*ppmi_analyses* contains code related to analysis on transcriptomics data from PPMI.

*luxpark_analyses* contains code related to analysis on metabolomics data from LuxPARK.


### Requirements
The code was mostly written and tested in R (R 4.0.3) on both current Mac (Ventura) and Linux operating systems (Ubuntu 23.04), relying on multiple R and BioConductor packages that are listed at the beginning of the corresponding R scripts. The code should be compatible with later versions of R installed on current Mac, Linux or Windows systems. R software packages loaded at the beginning of the R scripts must be installed before using the code. R packages available on CRAN can be installed with the command:
install.packages("BiocManager::install("PACKAGE_NAME")

R packages from Bioconductor can be installed with the following commands:

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
BiocManager::install("PACKAGE_NAME")

### License
TBD


