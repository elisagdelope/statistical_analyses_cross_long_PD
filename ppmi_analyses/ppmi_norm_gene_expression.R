# Title: ppmi_norm_gene_expression.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs normalization of data for downstream visualization AND aggregation. 
# Usage: R ppmi_norm_gene_expression.R
# Usage: R ppmi_norm_gene_expression.R -v "BL" -a "01-dea-BL-PD" # for timepoint-specific analysis!
# Data: data from expression at gene level (STAR) -already filtered for low expression genes- and clinical (pheno) data.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)




# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(DESeq2)
library(argparser)



# I/O --------------------------------------------------------------------------
# Add command line arguments
p <- arg_parser("Normalizes expression data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--visit", help = "visit/timepoint name to perform deseq on", default = "ALL", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-TS-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable
analysis_name <- argv$analysis # analysis name variable

OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
if (t %in% c("BL", "V04", "V06", "V08")) {
  EXPRESSION.FILE <- file.path(OUT_DIR, paste0("flt_star_", t, ".tsv"))
  PHENO.FILE <- file.path(OUT_DIR, paste("ppmi_pheno", t, "dsq.tsv", sep = "_")) 
  if (!file.exists(PHENO.FILE)){
    PHENO.FILE <- file.path(OUT_DIR, paste0("pheno_", t, ".tsv")) 
  }
  OUT.EXPRESSION.FILE <- file.path(OUT_DIR, paste0("flt_norm_star_", t, ".tsv"))
} else if (t == "ALL") {
  EXPRESSION.FILE <- file.path(OUT_DIR, "flt_star_all_TS.tsv")
  OUT.EXPRESSION.FILE <- file.path(OUT_DIR, "flt_norm_star_all_TS.tsv")
} else if (t != "ALL") {
  stop("Adequate arguments were not provided. Check R ppmi_filter_gene_expression.R --help for the right usage.")
}



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# Data load --------------------------------------------------------------------

pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) %>%
  rename_with(toupper)
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>% 
  rename_with(toupper) %>%
  column_to_rownames(var = 'GENEID')

# 1. Match the metadata and counts data
# make sure that we have sample names that match between the two files, and that the samples are in the right order
all(colnames(expression) == pheno$PATIENT_VISIT_ID)

# 2. Create DESEq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = pheno, design = ~ AGE + GENDER + DIAGNOSIS)

# 3. perform median of ratios method
dds <- DESeq2::estimateSizeFactors(dds)

# 4. retrieve the normalized counts matrix from dds
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts) %>% 
  rownames_to_column(var = 'GENEID') 



# output results ---------------------------------------------------------------
readr::write_tsv(normalized_counts, file = OUT.EXPRESSION.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

