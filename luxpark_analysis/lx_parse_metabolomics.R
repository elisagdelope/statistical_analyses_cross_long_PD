# Title: lx_parse_metabolomics.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script parses luxpark metabolomics data
# Usage: Rscript lx_parse_clinical.R
# Data: data from the original quantification file, chemical annotation & metadata.

# GC ----------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readxl)
library(readr)
library(dplyr)
library(vroom)


# I/O ------------------------------------------------------------------------
METABOLOMICS_PATH <- "../data/Metabolon/LUXE-01-20PHML+"  # Path to Salmon quantification files
OUT_DIR <- "../data/00-cleansing"
METABOLOMICS.FILE <- file.path(METABOLOMICS_PATH, "LUXE-01-20PHML+ DATA TABLES.XLSX")
PHENO.FILE <- file.path(OUT_DIR, "lx_pheno.tsv")



# generate metabolomics data ---------------------------------------------------
# add patient id and visit from pheno, restrictive to samples having both clinical and metabolon data
# rename colnames to "M"+chem_id for the sake of downstream processing and modelling

# extract chemical annotation
chem_data <- read_xlsx(METABOLOMICS.FILE, sheet = "Chemical Annotation")
chem_data <- chem_data %>%
  mutate(ANALYSIS_ID = paste("M", CHEM_ID, sep = ""))
chem_data[apply(chem_data, 1, function(i) any(grepl("\n|\r", i))),] %>%  # remove xlsx artifacts
  mutate_all(~str_replace(., "(\\,\\\r\\\n)", ",")) 
readr::write_delim(chem_data, file = file.path(OUT_DIR, "chemical_annotation.tsv"), delim = "\t")
rm(chem_data)


# load vars to merge from clinical data
pheno <- vroom(PHENO.FILE, show_col_types = FALSE)
pheno <- pheno %>% 
  select(all_of(c("PATIENT_ID", "VISIT", "SAMPLE_ID")))

# extract peak area
peak_area <- read_excel(METABOLOMICS.FILE, sheet = "Peak Area Data")
colnames(peak_area)[2:ncol(peak_area)] <- paste("M", colnames(peak_area)[2:ncol(peak_area)], sep = "")
colnames(peak_area)[1] <- "SAMPLE_ID"
peak_area <- pheno %>%
  inner_join(peak_area, by = "SAMPLE_ID")
readr::write_tsv(peak_area, file = file.path(OUT_DIR, "peak_area_data.tsv"))
rm(peak_area)

# extract Batch-normalized Data
batch_norm <- read_excel(METABOLOMICS.FILE, sheet = "Batch-normalized Data")
colnames(batch_norm)[2:ncol(batch_norm)] <- paste("M", colnames(batch_norm)[2:ncol(batch_norm)], sep = "")
colnames(batch_norm)[1] <- "SAMPLE_ID"
batch_norm <- pheno %>%
  inner_join(batch_norm, by = "SAMPLE_ID")
readr::write_tsv(batch_norm, file = file.path(OUT_DIR, "batch-norm_data.tsv"))
rm(batch_norm)

# extract Batch-norm Imputed Data
batch_norm_imputed <- read_excel(METABOLOMICS.FILE, sheet = "Batch-norm Imputed Data")
colnames(batch_norm_imputed)[2:ncol(batch_norm_imputed)] <- paste("M", colnames(batch_norm_imputed)[2:ncol(batch_norm_imputed)], sep = "")
colnames(batch_norm_imputed)[1] <- "SAMPLE_ID"
batch_norm_imputed <- pheno %>%
  inner_join(batch_norm_imputed, by = "SAMPLE_ID")
readr::write_tsv(batch_norm_imputed, file = file.path(OUT_DIR, "batch-norm_imputed_data.tsv"))
rm(batch_norm_imputed)

# extract Log Transformed Data
log_transformed <- read_excel(METABOLOMICS.FILE, sheet = "Log Transformed Data")
colnames(log_transformed)[2:ncol(log_transformed)] <- paste("M", colnames(log_transformed)[2:ncol(log_transformed)], sep = "")
colnames(log_transformed)[1] <- "SAMPLE_ID"
log_transformed <- pheno %>%
  inner_join(log_transformed, by = "SAMPLE_ID")
readr::write_tsv(log_transformed, file = file.path(OUT_DIR, "log_transformed_data.tsv"))
rm(log_transformed)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

