# Title: lx_generate_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates pathway-level aggregated metabolomics data (mean, median & st) across metabolites of pathways (defined in annotation file)
# Usage: Rscript lx_generate_pathway_level.R
# Usage: Rscript lx_generate_pathway_level.R -v "V0" -a "01-dea-V0-PD" # for timepoint-specific analysis!
# Data: data from metabolomics at metabolite level and annotation file containings pathway annotation.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(argparser)
library(pathifier)
library(caret)



# I/O --------------------------------------------------------------------------
# Add command line arguments
p <- arg_parser("Generates aggregated metabolomics statistics", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--visit", help = "visit/timepoint name", default = "ALL", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-TS-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable
analysis_name <- argv$analysis # analysis name variable

ANNOTATION.PATH <- "../data/00-cleansing/chemical_annotation.tsv"
METAB.PATH <- "../data/00-cleansing/log_transformed_data_fte.tsv" # to be changed if de novo analysis!
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level")

if (t %in% c("V0")) {
  METAB.PATH <- file.path(OUT_DIR, paste0("log_transformed_", t, ".tsv"))
} else if (t != "ALL") {
  stop("Adequate arguments were not provided. Check R lx_generate_pathway_level.R --help for the right usage.")
} 

source("func_summarize_pathway_level.R")



# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------------
annotation <- vroom(ANNOTATION.PATH, col_types = cols())
metab <- vroom(METAB.PATH, col_types = cols()) 



# Generate pathway-level metabolomics ------------------------------------------

# flip the dataset as rows should correspond to metabolites and columns to samples
metab.t = metab %>%
  dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>% 
  pivot_longer(!SAMPLE_ID, names_to = "METABOLITES", values_to = "COUNT") %>% 
  pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
  column_to_rownames("METABOLITES")

# generate list of lists where list name is pw and list components are ids:
PW = list()
for (p in unique(annotation$SUB_PATHWAY[!is.na(annotation$SUB_PATHWAY)])) {
  PW[[p]] <- annotation %>% 
    filter(SUB_PATHWAY == p) %>%
    pull(ANALYSIS_ID)
}

PW_AGG <- list()
for (st in c("mean", "median", "sd", "pca")) { # add pathifier?
  PW_AGG[[st]] <- summarize_pathway_level(metab.t, PW, type = st, minsize = 6)
  
  # flip the matrix (samples, features) to be consistent with saved metabolite-level format
  PW_AGG[[st]] <- as.data.frame(t(PW_AGG[[st]])) %>%
    rownames_to_column(var = "SAMPLE_ID") %>%
    inner_join(metab[, c("PATIENT_ID", "VISIT", "SAMPLE_ID")], 
               by = "SAMPLE_ID") %>% 
    relocate(c("PATIENT_ID", "VISIT"))
  
  # output results
  if (t != "ALL") {
    readr::write_tsv(PW_AGG[[st]], file = file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, "_", t, ".tsv")))
  } else {
    readr::write_tsv(PW_AGG[[st]], file = file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv")))
  }
}


# generate pathifier scores
st = "pathifier"

# remove near zero variance metabolites (strict "near zero" => "almost zero" variance)
nzv <- nearZeroVar(metab %>% dplyr::select(-any_of(c("SAMPLE_ID", "PATIENT_ID", "VISIT"))), 
                   freqCut = 100, # default = 20
                   uniqueCut = 5)  # default = 10
#ATTENTION: to generate pw level for de novo, the default nearZeroVar thresholds were used.
metab.t <- metab.t[-nzv, ]
pathifier_agg <- quantify_pathways_deregulation(data = as.matrix(metab.t),
                                                allgenes = rownames(metab.t),
                                                syms = PW, 
                                                pathwaynames = names(PW),
                                                normals = NULL, 
                                                logfile = file.path(OUT_DIR_PATHWAY, paste("lx", st, "log", sep = ".")), 
                                                attempts = 5,
                                                min_exp = -10) # = remove effect if min_exp 

pathifier_scores <- data.frame(Reduce(rbind, pathifier_agg$scores))
colnames(pathifier_scores) <- colnames(metab.t) 
rownames(pathifier_scores) <- names(pathifier_agg$scores)
PW_AGG[[st]] <- pathifier_scores

# flip the matrix (samples, features) to be consistent with saved metabolite-level format
PW_AGG[[st]] <- as.data.frame(t(PW_AGG[[st]])) %>%
  rownames_to_column(var = "SAMPLE_ID") %>%
  inner_join(metab[, c("PATIENT_ID", "VISIT", "SAMPLE_ID")], 
             by = "SAMPLE_ID") %>% 
  relocate(c("PATIENT_ID", "VISIT"))

# output results
if (t != "ALL") {
  readr::write_tsv(PW_AGG[[st]], file = file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, "_", t, ".tsv")))
} else {
  readr::write_tsv(PW_AGG[[st]], file = file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv")))
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


