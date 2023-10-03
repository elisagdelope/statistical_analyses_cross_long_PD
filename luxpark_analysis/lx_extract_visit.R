# Title: lx_extract_visit.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script extracts clinical and log-transformed metabolomics data from a specific visit
# Usage: Rscript lx_extract_visit.R -v V0 -a "02-pred-V0-PD"
# Data: data from metabolites abundance, clinical data.

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(argparser, quietly = TRUE)



# I/O --------------------------------------------------------------------------
METAB.FILE <- "../data/00-cleansing/log_transformed_data_fte.tsv" # to be changed if de novo analysis!
PHENO.FILE <- "../data/00-cleansing/lx_pheno.tsv"

# Add command line arguments
p <- arg_parser("Filters expression data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--visit", help = "visit/timepoint name to perform deseq on", default = "V0", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-V0-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable
analysis_name <- argv$analysis # analysis name variable
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 

if (t %in% c("V0")) {
  OUT.PHENO.FILE <- file.path(OUT_DIR, paste0("pheno_", t, ".tsv"))
  OUT.METAB.FILE <- file.path(OUT_DIR, paste0("log_transformed_", t, ".tsv"))
} else {
  stop("Adequate arguments were not provided. Check R lx_generate_pathway_level.R --help for the right usage.")
} 



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# Data load --------------------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 
metab <- vroom(METAB.FILE, col_types = cols())



# filter visit 0 ---------------------------------------------------------------
pheno <- pheno %>%
  dplyr::filter(VISIT == t) 

metab <- metab %>%
  dplyr::filter(VISIT == t) 



# export outputs  --------------------------------------------------------------
readr::write_tsv(pheno, file = OUT.PHENO.FILE) 
readr::write_tsv(metab, file = OUT.METAB.FILE) 


# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


