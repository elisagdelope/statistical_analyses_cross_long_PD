# Title: lx_parse_clinical.R
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
library(jsonlite)
library(stringr)



# I/O ------------------------------------------------------------------------
IN_DIR <- "../data/ADA-clinical"  # Path to Salmon quantification files
OUT_DIR <- "../data/00-cleansing"
JSON.FILE <- file.path(IN_DIR, "lux_park.clinical.json")
METADATA.FILE <- file.path("../data/Metabolon/20210603 Sample metadata with ND.xlsx")
OUT.FILE <- file.path(OUT_DIR, "lx_pheno.tsv")



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# load & transform data --------------------------------------------------------
datajson <- fromJSON(JSON.FILE)
clinical = do.call(cbind, datajson)
clinical = data.frame(clinical)

datajson <- fromJSON(JSON.FILE) # heavy load (~0.7GB)

# get clinical vars of interest
# subjet id, gender, visit, diagnosis, age, height, weight, bmi, hoehn and yahr, birth date, education level, high blood pressure, other diseases, visits dates, date dropout, date excluded.
clinical_vars <- c("cdisc_dm_usubjd", "cdisc_dm_sex", "redcap_event_name", "adiag_final_diagnosis", "diag_control", "control_q1", "diag_ipd", "meds_denovo", "question219_219",
  "sv_age", "status_weight", "status_height", "status_bmi", "hoehn_and_yahr_staging", "cdisc_dm_brthdtc", "cdisc_sc_sctestcd_edlevel", "scopa_q26c", "scopa_q26d", 
  "agecalc_visit_1_date_start", "agecalc_visit_2_date_start", "agecalc_visit_3_date_start", "agecalc_visit_4_date_start", "agecalc_visit_5_date_start", "agecalc_visit_6_date_start", "agecalc_visit_7_date_start",
  "dm_dropout_excluded2", "dm_dropout_excluded4")

clinical <- datajson %>%
  dplyr::select(all_of(clinical_vars))

# export the subset
readr::write_tsv(clinical, file = file.path(OUT_DIR, "lux_park.clinical_subset.tsv"))
rm(datajson)



# generate pheno file ----------------------------------------------------------
metadata <- read_excel(METADATA.FILE)
metadata <- metadata %>%
  dplyr::select(all_of(c("ND code", "Visit", "PARENT_SAMPLE_NAME"))) %>%
  dplyr::rename(PATIENT_ID = "ND code",
                SAMPLE_ID ="PARENT_SAMPLE_NAME") %>%
  rename_with(toupper) %>%
  mutate(VISIT = paste0(str_sub(VISIT, 1, 1), as.character(as.numeric(str_sub(VISIT, -1)) -1)))  # convert Visit 1, Visit 2,... -> V0, V1,...


pheno <- data.frame(
  PATIENT_ID = clinical$cdisc_dm_usubjd,
  VISIT = clinical$redcap_event_name,
  DIAGNOSIS = clinical$adiag_final_diagnosis,
  GENDER = clinical$cdisc_dm_sex,
  AGE = clinical$sv_age, 
  DENOVO = clinical$meds_denovo,
  FAM_HISTORY = clinical$question219_219,
  HOEHN_YAHR = clinical$hoehn_and_yahr_staging, 
  WEIGHT = clinical$status_weight, 
  HEIGHT = clinical$status_height, 
  BMI =  clinical$status_bmi, 
  EDUCATION = clinical$cdisc_sc_sctestcd_edlevel, 
  DATE_BIRTH = clinical$cdisc_dm_brthdtc, 
  HIGH_BLOOD_PRSR =  clinical$scopa_q26c, 
  OTHER_DISEASES = clinical$scopa_q26d
  ) 

pheno <- pheno %>%
  filter(DIAGNOSIS %in% c(1, 14)) %>%                   # 1 = PD, 14 = control
  mutate(
    DIAGNOSIS = case_when(                              
      as.character(DIAGNOSIS) %in% c("1") ~ "PD",        
      as.character(DIAGNOSIS) %in% c("14") ~ "HC"),
    GENDER = case_when(                                 # 1= male, 2= female
      as.character(GENDER) %in% c("1") ~ "M",
      as.character(GENDER) %in% c("2") ~ "F"),
    VISIT = paste0("V", VISIT),
    DENOVO = ifelse(DENOVO, 1, 0),                      # TRUE = 1, FALSE = 0
    FAM_HISTORY = ifelse(FAM_HISTORY, 1, 0),            # TRUE = 1, FALSE = 0
    HIGH_BLOOD_PRSR = ifelse(HIGH_BLOOD_PRSR, 1, 0),    # TRUE = 1, FALSE = 0
    OTHER_DISEASES = ifelse(OTHER_DISEASES, 1, 0)       # TRUE = 1, FALSE = 0
    )

# add metabolomics sample ids, restrictive to samples having both clinical and metabolon data
pheno <- metadata %>%
  inner_join(pheno, by = c("PATIENT_ID", "VISIT"))

  

# export outputs  --------------------------------------------------------------
readr::write_tsv(pheno, file = OUT.FILE) 



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


