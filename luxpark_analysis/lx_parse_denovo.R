# Title: lx_parse_denovo.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script parses luxpark de novo file to retrieve de novo patients ID with matching metabolomics data
# Usage: Rscript lx_parse_denovo.R
# Data: data from the original de novo luxpark file.

# GC ----------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readxl)
library(readr)
library(dplyr)
library(stringr)



# I/O ------------------------------------------------------------------------
IN_DIR <- "../data/ADA-clinical"  # Path to Salmon quantification files
OUT_DIR <- "../data/00-cleansing"
DENOVO.FILE <- file.path(IN_DIR, "De novo_LuxPark_PD_4_2022.xlsm")
OUT.FILE <- file.path(OUT_DIR, "lx_denovo_list.csv")
METADATA.FILE <- file.path("../data/Metabolon/20210603 Sample metadata with ND.xlsx")



# generate de novo list ---------------------------------------------------
metadata <- read_excel(METADATA.FILE)
metadata <- metadata %>%
  dplyr::select(all_of(c("ND code", "Visit", "PARENT_SAMPLE_NAME"))) %>%
  dplyr::rename(PATIENT_ID = "ND code",
                SAMPLE_ID ="PARENT_SAMPLE_NAME") %>%
  rename_with(toupper) %>%
  mutate(VISIT = paste0(str_sub(VISIT, 1, 1), as.character(as.numeric(str_sub(VISIT, -1)) -1)))  # convert Visit 1, Visit 2,... -> V0, V1,...


# extract patients ID from file
denovo_data <- read_xlsx(DENOVO.FILE, sheet = "LuxPARK-DeNovoLP42022_DATA_LABE")
denovo_data <- denovo_data %>%
  dplyr::rename(PATIENT_ID = "Subject ID") %>%
  mutate(VISIT = paste0(str_sub(`Event Name`, 1, 1), as.character(as.numeric(str_sub(`Event Name`, -1)) -1)))  # convert Visit 1, Visit 2,... -> V0, V1,...


# add metabolomics sample ids, restrictive to samples having both clinical and metabolon data
denovo_data <- metadata %>%
  inner_join(denovo_data, by = c("PATIENT_ID", "VISIT"))
denovo_list <- denovo_data[["PATIENT_ID"]]



# export outputs  --------------------------------------------------------------
write(denovo_list, file = OUT.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

