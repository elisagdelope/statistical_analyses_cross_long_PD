# Title: ppmi_parse_clinical_targets.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script parses phenotype and clinical data from the original clinical dataset from PPMI
# Usage: R ppmi_parse_clinical_targets.R
# Data: data from the original PPMI clinical file (curated data cuts)

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(vroom)
library(dplyr)



# I/O --------------------------------------------------------------------------
OUT_DIR <- "../data/00-cleansing"
RAW.PHENO.FILE <- file.path(OUT_DIR, "ppmi_raw_pheno.tsv")
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
OUT.FILE <- file.path(OUT_DIR, "ppmi_clinical_targets_PD.tsv")



# Data load --------------------------------------------------------------------
pheno_raw <- vroom(RAW.PHENO.FILE, show_col_types = FALSE)
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddlf")) 



# Data transformation -----------------------------------------------------------
pheno_raw$patient_visit_id = paste0(pheno_raw$PATNO, ".", pheno_raw$EVENT_ID)

# select & filter variables 
clinical_vars <- c("patient_visit_id", "EVENT_ID", "PATNO", "age", "gen", "APPRDX",
                   "updrs3_score", 
                   "updrs3_score_on", 
                   "updrs2_score",
                   "updrs1_score",
                   "updrs_totscore",
                   "updrs_totscore_on",
                   "updrs4_score",
                   "upsit", 
                   "bjlot", 
                   "ess", # epworth
                   "ess_cat",  # epworth
                   "gds", # depression
                   "gds_cat", # depression
                   "MCI_testscores",
                   "cogstate",
                   "upsit_cat", 
                   "moca",  # moca
                   "quip",  #quip
                   "quip_any",     #quip         
                   "rem",     # rem
                   "rem_cat",  # rem
                   "rem_q6",  # rem
                   "scopa",
                   "tremor",
                   "tremor_on") # scopa

pheno_raw <- pheno_raw %>%
  dplyr::select(all_of(clinical_vars)) 

pheno_prog_PD <- pheno_raw %>% 
  inner_join(pheno %>% 
               dplyr::select(all_of(c("patient_visit_id"))),
             by = c("patient_visit_id")) %>%
  filter(APPRDX == "1") # filter for PD patients (remove HC)

pheno_prog_PD <- data.frame(
  PATIENT_VISIT_ID = pheno_prog_PD$patient_visit_id,
  PATIENT_ID = pheno_prog_PD$PATNO,
  VISIT = pheno_prog_PD$EVENT_ID,
  AGE = pheno_prog_PD$age,
  GENDER = pheno_prog_PD$gen,
  DIAGNOSIS = pheno_prog_PD$APPRDX,
  UPDRS1 = pheno_prog_PD$updrs1_score,
  UPDRS2 = pheno_prog_PD$updrs2_score,
  UPDRS3 = pheno_prog_PD$updrs3_score,
  UPDRS4 = pheno_prog_PD$updrs4_score,
  UPDRS3_ON = pheno_prog_PD$updrs3_score_on,
  UPDRS_TOTAL = pheno_prog_PD$updrs_totscore,
  UPDRS_TOTAL_ON = pheno_prog_PD$updrs_totscore_on,
  UPSIT = pheno_prog_PD$upsit, 
  UPSIT_CAT = pheno_prog_PD$upsit_cat, 
  BJLOT = pheno_prog_PD$bjlot, 
  ESS = pheno_prog_PD$ess, # epworth
  ESS_CAT = pheno_prog_PD$ess_cat,  # epworth
  GDS = pheno_prog_PD$gds, # depression
  GDS_CAT = pheno_prog_PD$gds_cat, # depression
  MCI = pheno_prog_PD$MCI_testscores,
  COGSTATE = pheno_prog_PD$cogstate,
  MOCA = pheno_prog_PD$moca,  # moca
  QUIP = pheno_prog_PD$quip,  #quip
  QUIP_ANY = pheno_prog_PD$quip_any,     #quip         
  REM = pheno_prog_PD$rem,     # rem
  REM_CAT = pheno_prog_PD$rem_cat,  # rem
  REM_Q6 = pheno_prog_PD$rem_q6,  # rem
  SCOPA = pheno_prog_PD$scopa,
  TREMOR = pheno_prog_PD$tremor,
  TREMOR_ON = pheno_prog_PD$tremor_on # scopa
  ) %>%
  rename_with(toupper)

pheno_prog_PD <- pheno_prog_PD %>% 
  mutate(
    DIAGNOSIS = case_when(
    as.character(DIAGNOSIS) %in% c("2") ~ "HC",
    as.character(DIAGNOSIS) %in% c("1") ~ "PD"),
    GENDER = case_when(
      as.character(GENDER) %in% c("1") ~ "M",
      as.character(GENDER) %in% c("2") ~ "F")
    )



# output results ---------------------------------------------------------------
readr::write_tsv(pheno_prog_PD, file = OUT.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


