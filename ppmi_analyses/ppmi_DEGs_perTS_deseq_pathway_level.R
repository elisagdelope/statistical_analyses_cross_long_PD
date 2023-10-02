# Title: ppmi_DEGs_perTS_deseq_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script computes deseq - DEAs PD-HC for each timepoint at aggregated levels, then filters those common to all timepoints
# Usage: R ppmi_DEGs_perTS_deseq_pathway_level.R 
# Data: data from expression at aggregated levels and pheno.


# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(DESeq2)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
AGG.EXPRESSION.FILE <- list()
for (db in c("GOBP", "GOCC", "CORUM")) {
  AGG.EXPRESSION.FILE[[db]] <- list()
  for (st in c("mean", "median", "sd", "pathifier")) {
    AGG.EXPRESSION.FILE[[db]][[st]] <- file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_"))
  }
}



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------------
pheno_dsq <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) 
pheno_dsq <- pheno_dsq %>% 
  rename_with(toupper)



# Subset data per timepoint ---------------------------------------
visits <- as.character(unique(pheno_dsq$VISIT))
classes <- c("PD", "HC") # please mind the order! it will be used for deseq -> reference/control group goes second

AGG.EXPRESSION <- list()
GENES <- list()
DEA_visits <- list()
for (db in c("GOBP", "GOCC", "CORUM")) {
  AGG.EXPRESSION[[db]] <- list()
  DEA_visits[[db]] <- list()
  
  for (st in names(AGG.EXPRESSION.FILE[[db]])) { 
    AGG.EXPRESSION[[db]][[st]] <- vroom(AGG.EXPRESSION.FILE[[db]][[st]] , col_types = cols(), delim = "\t") # load data
    DEA_visits[[db]][[st]] <- list()
    var_id = paste0(db, "_name")
    E.STAR <- list()
    DDS <- list()
    RES <- list()
    for (v in visits) {
      # subset expression data for each visit
      PHENO <- pheno_dsq %>% 
        dplyr::filter(VISIT == v)
      EXPR <- AGG.EXPRESSION[[db]][[st]] %>% 
        column_to_rownames(var = var_id) %>%
        dplyr::select(matches(v)) 
      
      # create deseq object
      if ((st == "pathifier") & (db != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
        E.STAR[[v]] <- DESeq2::DESeqDataSetFromMatrix(round(EXPR * 100), 
                                                      colData = PHENO, 
                                                      design = ~ AGE + GENDER + DIAGNOSIS)
      } else {
        E.STAR[[v]] <- DESeq2::DESeqDataSetFromMatrix(round(EXPR), 
                                                      colData = PHENO, 
                                                      design = ~ AGE + GENDER + DIAGNOSIS)
      }
      
      # pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
      keep <- rowSums(counts(E.STAR[[v]])) >= 10
      E.STAR[[v]] <- E.STAR[[v]][keep,]
      
      # deseq ATTENTION: DATA ALREADY NORMALIZED
      #    DDS[[v]] <- DESeq(E.STAR[[v]])
      #    DDS <- estimateSizeFactors(E.STAR)
      sizeFactors_nonnorm <- rep(1, ncol(EXPR)) # create a vector of 1s to avoid normalization
      names(sizeFactors_nonnorm) <- colnames(EXPR)
      sizeFactors(E.STAR[[v]]) <- sizeFactors_nonnorm
      DDS[[v]] <- estimateDispersions(E.STAR[[v]])
      DDS[[v]] <- nbinomWaldTest(DDS[[v]], maxit=500)
      
      res.flt <- results(DDS[[v]], contrast = c("DIAGNOSIS", classes)) #, independentFiltering = FALSE)??
      res.flt <- as.data.frame(res.flt) %>% 
        rownames_to_column(var = var_id)
      RES[[v]] <- res.flt[order(res.flt$pvalue),]
      
      # output results           
      readr::write_tsv(as.data.frame(RES[[v]]), file = file.path(OUT_DIR_PATHWAY, paste("DESeq_res", db, st, paste0(classes, collapse = ""), paste0(v, ".tsv"), sep = "_")))
      print(paste("DESeq analysis for", paste0(classes, collapse = "-"), v, "completed"))
    
        
      ###### 
      # filter resulting DEGs
      
      # remove NAs in padj
      RES[[v]] <- RES[[v]][!is.na(RES[[v]]$padj),]
      # cut off nominal pvalue < 0.05
      RES[[v]] <- RES[[v]][RES[[v]]$pvalue < 0.05,]
      # output filtered result - NOT NECESSARY
      # readr::write_tsv(RES[[v]], file = file.path(OUT_DIR_PATHWAY, paste("DESeq_fltranked", db, st, paste0(classes, collapse = ""), paste0(v, ".tsv"), sep = "_")))
      
      GENES[[v]] <- RES[[v]][[var_id]]
    }

    # merge resulting DEGs on all visits
    common_genes <- Reduce(intersect, GENES)
    for (v in names(RES)) {
      RES[[v]] <- RES[[v]] %>%
        filter(get(var_id) %in% common_genes) %>%
        mutate(!!paste0('rank.',v) := rank(pvalue)) %>%
        dplyr::select(matches(paste0(var_id,'|p|rank')))
    }
    # concatenate datasets
    DEA_visits[[db]][[st]] <- as.data.frame(Reduce(function(...) inner_join(..., by = var_id), RES))
    
    # output results ---------------------------------------------------------------
    readr::write_tsv(DEA_visits[[db]][[st]], file = file.path(OUT_DIR_PATHWAY, paste(db, st, "consistent_DEGsPDHC_TS.tsv", sep = "_")))
    
  }
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


