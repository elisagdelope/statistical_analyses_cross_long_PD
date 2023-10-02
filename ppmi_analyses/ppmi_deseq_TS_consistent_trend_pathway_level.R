# Title: ppmi_deseq_TS_consistent_trend_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters pathways/protein complexes/cellular locations whose aggregated expression has consistent sign of foldchange across all the Bayesian Moderated t-statistics on consecutive timepoints.
# Usage: R ppmi_deseq_TS_consistent_trend_pathway_level.R
# Data: data from ebayes t-statistic at aggregated level: GOBP, GOCC and CORUM on all aggregated stats.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)


# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)


# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

DEA.FILES <- list()
labels <- c("l1","l2")
DEA.FILES[[labels[1]]] <- list.files(path = OUT_DIR_PATHWAY, pattern = "^.+deseq_PD.tsv$", full.names = TRUE) # please note it's important that l1 refers to the patient/symptom-specific class
DEA.FILES[[labels[2]]] <- list.files(path = OUT_DIR_PATHWAY, pattern = "^.+deseq_HC.tsv$", full.names = TRUE) 




# Data load and transformation -------------------------------------------------

for (db in c("GOBP", "GOCC", "CORUM")) {
  var_id <- paste0(db, "_name")
  for (st in c("mean", "median", "sd", "pathifier")) {
    DEA.FULL.C <- list()
    G.LIST <- list()
    for (class in labels) {
      FILE <- str_subset(DEA.FILES[[class]], paste(db, st, sep = "_"))
      
      if ((length(FILE) != 0) && (file.exists(FILE))) {
        DEA_full <- vroom(FILE, col_types = cols())
        
        # remove genes with padj = NA for all DEAs
        DEA_full <- DEA_full %>% filter_at(vars(starts_with('padj_')), all_vars (!is.na(.)))
        
        # filter genes with consistent sign of foldchange across the 3 DEAs
        DEA.FULL.C[[class]] <- DEA_full %>%
          mutate(sum_trend = rowSums(sign(DEA_full %>% dplyr::select(matches('^log2FoldChange_[A-Z]+.*\\d$')))),
                 up = (sum_trend == -3), # note log2foldchange(t1,t2) < 0 ==> t1 < t2 ==> trend is upwards ; log2foldchange(t1,t2) > 0 ==> t1 > t2 ==> trend is downwards
                 down = (sum_trend == 3),
                 is_consistent = (sum_trend %in% c(3, -3))) %>%
          filter(is_consistent == TRUE) %>% 
          rowwise() %>%
          mutate(pvalue_median = median(c_across(matches("pvalue"))),
                 padj_median = median(c_across(matches("adj"))),
                 logFC_median = median(c_across(matches("log2FoldChange")))) %>%
          dplyr::arrange(pvalue_median, desc(abs(logFC_median)))
        
        # store list of consistent genes
        G.LIST[[class]] <- DEA.FULL.C[[class]][[var_id]]
        
      } else {
        message(paste("File not found for", db, st, sep = " "))
        next
      }
    }
    
    
    # Class-specific consistent trend list
    # for each db-st combination having DEAs for 2 classes, compare lists in both classes and their sign of trend
    
    if (length(DEA.FULL.C) == 2) {
      SPEC <- list()
      common_genes <- Reduce(intersect, G.LIST) 
      for (c in seq(1,length(labels))) { 
        if (c < length(labels)) { # if class-1
          gene_spec <- Reduce(setdiff, G.LIST[[c+1]], G.LIST[[c]])  # elements in class l1 but not in class l2 
          gene_diftrend <- DEA.FULL.C[[labels[c]]] %>%      # elements in common between classes with different trends. 
            filter(get(var_id) %in% common_genes)  %>%
            dplyr::select(all_of(c(var_id, "up", "down"))) %>%
            inner_join(DEA.FULL.C[[labels[c+1]]] %>%
                         filter(get(var_id) %in% common_genes) %>%
                         dplyr::select(all_of(c(var_id, "up", "down"))),
                       by = var_id,
                       suffix = paste0("_", labels)) %>%
            filter(((get(paste0("up_", labels[c])) == TRUE) & (get(paste0("up_", labels[c+1])) == FALSE)) | 
                     ((get(paste0("up_", labels[c])) == FALSE) & (get(paste0("up_", labels[c+1])) == TRUE))) %>%
            pull(var_id)
          out_name <- "consistent_trend_TS.tsv"
        } else { # class-2
          gene_spec <- Reduce(setdiff, G.LIST[[c-1]], G.LIST[[c]]) # elements in class l2 but not in class l1
          out_name <- "HC_consistent_trend_TS.tsv"
        }
        
        gene_spec <- c(gene_spec, gene_diftrend)
        
        # subset identified class-specific genes from dataset
        SPEC[[labels[c]]] <- DEA.FULL.C[[labels[c]]] %>%
          filter(get(var_id) %in% gene_spec) %>%
          dplyr::select(-c(is_consistent, sum_trend))
        
        # output results       
        readr::write_tsv(SPEC[[labels[c]]], file = file.path(OUT_DIR_PATHWAY, paste(db, st, out_name, sep = "_")))
        
      }
    }
  }
}
      
