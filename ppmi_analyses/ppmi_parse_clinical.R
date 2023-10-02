# Title: ppmi_parse_clinical.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script parses phenotype and clinical data from the original clinical dataset from PPMI
# Usage: R ppmi_parse_clinical.R
# Data: data from the original PPMI clinical file (curated data cuts)

# GC ----------------------------------------------------------------------
rm(list = ls())
gc(T)

# Packages ----------------------------------------------------------------
library(readr)
library(dplyr)

# I/O ---------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# setwd(args[1])  # Set root to RNAsequencing folder!

FILE.DATA <- "../data/Curated_Data_Cuts/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020_pat_extract.csv"
OUTPUT <- "../data/00-cleansing"

message(paste0("Set to ouput -> ", OUTPUT))

# Main --------------------------------------------------------------------
if (!dir.exists(OUTPUT)) {
  dir.create(OUTPUT, recursive = T)
}

# Read the PPMI curated data directly and save to disk
data <- readr::read_csv(FILE.DATA)
readr::write_tsv(data, file = file.path(OUTPUT, paste0("ppmi_raw_pheno.tsv")))

# Parse selected PPMI phenotypes so they fit with what we have from GEO datasets
pheno <- data.frame(
  GSEID = data$PATNO,
  visit = data$EVENT_ID,
  Diagnosis = data$APPRDX,
  Age = data$age,
  Gender = data$gen,
  Education = data$educ,
  is_IPD = (data$fampd_new %in% c(3)), # 3 = non family w/PD; 2 = non 1st degree family w/PD..
  Moca = data$moca,  # Already adjusted for education
  HoehnYahr = data$NHY,
  Updrs1 = data$updrs1_score,
  Updrs2 = data$updrs2_score,
  Updrs3 = data$updrs3_score,
  Updrs4 = data$updrs4_score,
  Medication = data$PD_MED_USE,
  AgeOnset = data$ageonset,
  MCI = data$MCI_testscores,
  Batch = NA
)

pheno <- pheno %>% 
  mutate(Diagnosis = case_when(
    as.character(Diagnosis) %in% c("2") ~ "HC",
    as.character(Diagnosis) %in% c("1") ~ "PD",
    TRUE ~ "Others"))

pheno <- pheno %>% 
  mutate(Gender = case_when(
    as.character(Gender) %in% c("1") ~ "M",
    as.character(Gender) %in% c("2") ~ "F",
    TRUE ~ as.character(Gender)))

if (!"diagnosis-others" %in% args) {
  pheno <- pheno %>% 
    dplyr::filter(Diagnosis %in% c("HC", "PD"))
}

readr::write_tsv(pheno, file = file.path(OUTPUT, paste0("ppmi_pheno.tsv")))

# Session info ------------------------------------------------------------
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()
