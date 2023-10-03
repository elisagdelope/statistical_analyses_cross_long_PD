# Title: lx_denovo_filter.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script extracts clinical and log-transformed metabolomics data from de novo patients at v0
# Usage: Rscript lx_denovo_filter.R
# Data: data from metabolites de novo patients list, metabolites & pheno files.

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)


# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-V0-NOVO"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
DENOVO.FILE <- "../data/00-cleansing/lx_denovo_list.csv"
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv")
METAB.FILE <- file.path(OUT_DIR, "log_transformed_V0.tsv")
# ATTENTION: files will be rewritten, in-files = out-files



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# Data load --------------------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 
metab <- vroom(METAB.FILE, col_types = cols())
denovo_list <- vroom(DENOVO.FILE, delim = ",", col_names = FALSE, show_col_types = FALSE) %>% 
  pull()



# filter de novo patients ------------------------------------------------------
pheno <- pheno %>%
  dplyr::filter((PATIENT_ID %in% denovo_list) | (DIAGNOSIS == "HC"))    # keep HC & de novo

metab <- metab %>%
  dplyr::filter(PATIENT_ID %in% pheno$PATIENT_ID)



# export outputs  --------------------------------------------------------------
readr::write_tsv(pheno, file = PHENO.FILE) 
readr::write_tsv(metab, file = METAB.FILE) 



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


