# Title: ppmi_deseq_visit.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs DESEQ analysis on specific timepoint data for diagnosis with Age & Gender as covariates.
# Usage: R ppmi_deseq_visit.R "BL"
# Data: data from gene expression + pheno (ready to deseq files) of STAR+FC 

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)

# Packages ---------------------------------------------------------------------
library(readr)
library(dplyr)
library(vroom)
library(tidyr)
library(stringr)
library(tibble)
library(DESeq2)
library(edgeR)
library(argparser, quietly = TRUE)

# I/O --------------------------------------------------------------------------

# Add command line arguments
p <- arg_parser("Performs DESEQ DEA analysis on data at specific timepoint", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "visit", help = "visit/timepoint name to perform deseq on", default = "BL", type = "string", short = "v")
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable

analysis_name <- paste0("01-dea-", t, "-PD")
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")

EXPRESSION.FILE <- file.path(OUT_DIR, paste0("flt_star_", t, ".tsv")) # filtered expression data
PHENO.FILE <- file.path(OUT_DIR, paste("ppmi_pheno", t, "dsq.tsv", sep = "_"))  



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS)) | (!dir.exists(OUT_DIR_ENRICHMENT))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
  dir.create(OUT_DIR_ENRICHMENT, recursive = T)
}



# Data load --------------------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) %>%
  rename_with(toupper)
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>%
  column_to_rownames(var = "geneid")



# DEA of PD/HC for expression counts -------------------------------------------

E.STAR <- DESeq2::DESeqDataSetFromMatrix(expression, 
                                         colData = pheno, 
                                         design = ~ AGE + GENDER + DIAGNOSIS)

# pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
keep <- rowSums(counts(E.STAR)) >= 10
E.STAR <- E.STAR[keep,]

dds <- DESeq(E.STAR) # complete DESEQ process
res <- results(dds, contrast = c("DIAGNOSIS", "PD", "HC"))

res <- as.data.frame(res) %>% 
  rownames_to_column(var = "geneid") %>%
  arrange(pvalue)

# remove genes with padj = NA and pvalue < 0.05 #### CHANGED THRESHOLD PVALUE<0.05 ==> PADJ<0.05
res.flt <- res %>%
            filter(padj < 0.05) %>% ##### CHANGED!
            filter(!is.na(padj)) 

message("DESeq analysis completed")


# output results ---------------------------------------------------------------
readr::write_tsv(res, file = file.path(OUT_DIR, paste0("DESeq_res_", t, ".tsv")))
readr::write_tsv(res.flt, file = file.path(OUT_DIR, paste0("genes_deseq_", t, ".tsv")))



# Session info ------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

