# Title: lx_filter_treateffect.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script aims to remove L-dopa associated treatment effects. Correlation filter with 3-methoxytyrosine; 
# Tyrosine and Tryptophan pathways are removed from metabolomics abundance matrix.
# Usage: Rscript lx_filter_treateffect.R
# Data: data from metabolites abundance

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
METAB.INFILE <- "../data/00-cleansing/log_transformed_data.tsv"
METAB.OUTFILE <- "../data/00-cleansing/log_transformed_data_fte.tsv"
ANNOTATION.FILE <- "../data/00-cleansing/chemical_annotation.tsv"
M1342.FILE <- "../data/00-cleansing/M1342.tsv"


# Data processing -----------------------------------------------------
data <- vroom(METAB.INFILE, col_types = cols())
annotation <- vroom(ANNOTATION.FILE, col_types = cols())

# extract 3-methoxytyrosine
M1342 <- data %>%
  dplyr::select(c("SAMPLE_ID", "M1342"))

# remove metabolites correlated with "M1342" (3-methoxytyrosine)
data.cor = cor(data[,-c(1:3)]) 
mcor <- data.cor[, "M1342"][abs(data.cor[, "M1342"]) > 0.2] # remove any potential correlated metabolite
mcor <-mcor[!is.na(mcor)]
data <- data %>%
  dplyr::select(-any_of(unique(names(mcor))))

# remove metabolites in Tyrosine and Tryptophan metabolism pathways 
pw_filter <- annotation %>%
  filter(SUB_PATHWAY %in% c("Tyrosine Metabolism", "Tryptophan Metabolism")) %>%
  pull(ANALYSIS_ID)
data <- data %>%
  dplyr::select(-any_of(pw_filter))



# export outputs  --------------------------------------------------------------
readr::write_tsv(data, file = METAB.OUTFILE)
readr::write_tsv(M1342, file = M1342.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


