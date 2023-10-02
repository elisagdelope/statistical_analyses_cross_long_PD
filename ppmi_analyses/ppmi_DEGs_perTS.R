# Title: ppmi_DEGs_perTS.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script computes DESEQ DEAs PD-HC for each timepoint at gene level, then filters those common to all timepoints
# Usage: R ppmi_DEGs_perTS.R 
# Data: data from expression at gene level and pheno.


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
library(edgeR)
library(SummarizedExperiment)


# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"

EXPRESSION.FILE <- file.path(OUT_DIR, "flt_star_all_TS.tsv")
OUT.FILE.STARTWITH <- file.path(OUT_DIR, "DESeq_res_")




# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------
e_level = "GENE"
target = "DIAGNOSIS"
if (e_level == "GENE") { 
  var_id = "geneid"
} else { 
  var_id = paste0(e_level, "_name")
}

pheno_dsq <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl"))
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>% 
  column_to_rownames(var = var_id) 



# Subset data per timepoint ---------------------------------------
visits <- str_extract(sort(colnames(expression)), '[A-Z0-9]*$')
visits <- unique(visits[(!visits == "") & (!is.na(visits))])
classes <- c("PD", "HC") # please mind the order! it will be used for deseq -> reference/control group goes second

pheno_dsq <- pheno_dsq %>% 
  rename_with(toupper)


RES <- list()
GENES <- list()
for (v in visits) {
  # subset expression data for each visit
  PHENO <- pheno_dsq %>% 
    dplyr::filter(VISIT == v)
  STAR <- expression %>%
    dplyr::select(matches(v)) 
  
  E.STAR <- DESeq2::DESeqDataSetFromMatrix(STAR, 
                                                colData = PHENO, 
                                                design = ~ AGE + GENDER + DIAGNOSIS)
  
  # pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
  keep <- rowSums(counts(E.STAR)) >= 10 
  E.STAR <- E.STAR[keep,]
                        
  # deseq
  DDS <- estimateSizeFactors(E.STAR)
  DDS <- estimateDispersions(DDS)
  DDS <- nbinomWaldTest(DDS, maxit=500)
#  DDS <- DESeq(E.STAR)
  res.flt <- results(DDS, contrast = c(target, classes[1], classes[2])) # change order if necessary, pass by argument
  res.flt <- as.data.frame(res.flt) %>% 
    rownames_to_column(var = var_id)
  RES[[v]] <- res.flt[order(res.flt$pvalue),]
  
  rm(DDS, E.STAR, STAR, PHENO)
  # output results
  readr::write_tsv(as.data.frame(RES[[v]]), file = file.path(OUT_DIR, paste0("DESeq_res_", classes[1], classes[2], "_", v, ".tsv")))
  print(as.data.frame(RES[[v]][1:10,]))
  print(paste("DESeq analysis for", v, "completed"))
  
  
  
  
  ###### 
  # filter resulting DEGs
  
  # remove NAs in padj
  RES[[v]] <- as.data.frame(RES[[v]])
  RES[[v]] <- RES[[v]][!is.na(RES[[v]]$padj),]
  # cut off nominal pvalue < 0.05
  RES[[v]] <- RES[[v]][RES[[v]]$pvalue < 0.05,]
  # output filtered result
  readr::write_tsv(as.data.frame(RES[[v]]), file = file.path(OUT_DIR, paste0("DESeq_fltranked_", classes[1], classes[2], "_", v, ".tsv")))
  
  GENES[[v]] <- RES[[v]] %>% 
    pull(var_id)
}

# merge resulting DEGs on all visits
common_genes <- Reduce(intersect, GENES)
for (v in names(RES)) {
  RES[[v]] <- RES[[v]] %>%
    filter(geneid %in% common_genes) %>%
    mutate(!!paste0('rank.',v) := rank(pvalue)) %>%
    dplyr::select(matches(paste0(var_id,'|p|rank')))
}
 
# concatenate datasets
DEA_visits <- Reduce(function(...) inner_join(..., by = 'geneid'), RES)


# output results ---------------------------------------------------------------
readr::write_tsv(DEA_visits, file = file.path(OUT_DIR, "genes_consistent_DEGsPDHC_TS.tsv"))
