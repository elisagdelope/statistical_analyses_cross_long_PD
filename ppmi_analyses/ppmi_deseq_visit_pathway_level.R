# Title: ppmi_deseq_visit_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs DESEQ analysis on specific timepoint data at aggregated level (mean, median & st): GO BP, GO CC, CORUM. DESEQ for diagnosis with Age & Gender as covariates.
# Usage: R ppmi_deseq_visit_pathway_level.R "BL"
# Data: data from the gene expression + pheno (ready to deseq files) of patients(PD) & controls(HC)

# GC ----------------------------------------------------------------------
rm(list=ls())
gc(T)


# Packages ----------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(DESeq2)
library(edgeR)
library(argparser, quietly = TRUE)



# I/O ---------------------------------------------------------------------

# Add command line arguments
p <- arg_parser("Performs DESEQ DEA analysis on data at specific timepoint", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "visit", help = "visit/timepoint name to perform deseq on", default = "BL", type = "string", short = "v")
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable

analysis_name <- paste0("01-dea-", t, "-PD")
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- file.path(OUT_DIR, paste("ppmi_pheno", t, "dsq.tsv", sep = "_"))  
AGG.EXPRESSION.FILE <- list()
for (db in c("GOBP", "GOCC", "CORUM")) {
  AGG.EXPRESSION.FILE[[db]] <- list()
  for (st in c("mean", "median", "sd", "pathifier")) {
    AGG.EXPRESSION.FILE[[db]][[st]] <- file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_"))
  }
}



# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------

pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) %>%
  rename_with(toupper)


for (db in c("GOBP", "GOCC", "CORUM")) {
  var_id <- paste0(db, "_name")
  
  for (st in names(AGG.EXPRESSION.FILE[[db]])) { 
    expression <- vroom(AGG.EXPRESSION.FILE[[db]][[st]] , col_types = cols(), delim = "\t") %>%  # load data
      column_to_rownames(var = var_id) 
    
    if (nrow(expression) > 2) {
      
      if ((st == "pathifier") & (db != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
        E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression * 100), 
                                                  colData = pheno, 
                                                  design = ~ AGE + GENDER + DIAGNOSIS)
      } else {
        E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression), 
                                                  colData = pheno, 
                                                  design = ~ AGE + GENDER + DIAGNOSIS)
      }
      
      # pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
      keep <- rowSums(counts(E.STAR)) >= 10
      E.STAR <- E.STAR[keep,]
      
      # deseq ATTENTION: DATA ALREADY NORMALIZED
      #    DDS <- DESeq(E.STAR)
      #    DDS <- estimateSizeFactors(E.STAR)
      sizeFactors_nonnorm <- rep(1, ncol(expression)) # create a vector of 1s to avoid normalization
      names(sizeFactors_nonnorm) <- colnames(expression)
      sizeFactors(E.STAR) <- sizeFactors_nonnorm
      DDS <- estimateDispersions(E.STAR)
      DDS <- nbinomWaldTest(DDS, maxit=500)
      
      res <- results(DDS, c("DIAGNOSIS", "PD", "HC")) #, independentFiltering = FALSE)??
      res <- as.data.frame(res) %>% 
        rownames_to_column(var = var_id) %>%
        arrange(pvalue)
      
      # remove genes with padj = NA and pvalue > 0.05
      res.flt <- res %>%
        filter(pvalue < 0.05) %>%
        filter(!is.na(padj)) 
      print(paste("DESeq analysis for", db, st, t, "completed"))

      # output results
      readr::write_tsv(as.data.frame(res), file = file.path(OUT_DIR_PATHWAY, paste0(paste("DESeq_res", db, st, t, sep = "_"), ".tsv")))
      readr::write_tsv(res.flt, file = file.path(OUT_DIR_PATHWAY, paste0(paste(db, st, "deseq", t, sep = "_"), ".tsv")))
      
    } else {
      message(paste(db, st, "failed: Need at least two features to fit a model", sep = " "))
      next
    }
    
  } 
}




# Session info ------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


