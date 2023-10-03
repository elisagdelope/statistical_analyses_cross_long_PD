# Title: lx_annotate_DEA.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script adds annotation data to input file, typically metabolite-level DA results.
# Usage: Rscript lx_annotate_DEA.R -f "lt_limma_V0.tsv" -a "01-dea-V0-PD"
# Data: data resulting from metabolites-level DA & annotation data.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(readr)
library(plyr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(stringr)
library(argparser)



# I/O --------------------------------------------------------------------------
IN_DIR <- "../data/00-cleansing/"
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")

# Add command line arguments
p <- arg_parser("Filters expression data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--file", help = "file name of DA to annotate", default = "lt_limma_V0.tsv", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-V0-PD", type = "character", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
f <- argv$file # file name variable
analysis_name <- argv$analysis # analysis name variable

OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
IN.DEA.FILE <- file.path(OUT_DIR, f)
OUT.DEA.FILE <- paste0(str_extract(IN.DEA.FILE,  ".+?(?=\\.[a-z])"), "_annotated.tsv")



# Data load & transformation ---------------------------------------------------
annotation <- vroom(ANNOTATION.FILE, col_types = cols())
dea_res <- vroom(IN.DEA.FILE, col_types = cols())

annotation_vars <- c("ANALYSIS_ID", "SUPER_PATHWAY",	"SUB_PATHWAY", "CHEMICAL_NAME",	"KEGG",	"PUBCHEM", "CHEM_ID",	"CAS", "HMDB", "TYPE")
annotation <- annotation %>%
  dplyr::select(all_of(annotation_vars))

# join with metabolite-level annotations with limma t-test results
dea_res <- dea_res %>% 
  dplyr::select(-any_of(c("AveExpr", "t", "B", "SUM_TREND",	"IS_CONSISTENT"))) %>%
  left_join(annotation, by = c("METABOLITES_ID" = "ANALYSIS_ID"))



# export outputs  --------------------------------------------------------------
readr::write_tsv(dea_res, file = OUT.DEA.FILE) 



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

