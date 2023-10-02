# Title: ppmi_deseq_TS_consistent_trend.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters genes with consistent sign of foldchange across all the consecutive DEAs for each class separately, and subsets identified class 1-specific (e.g. PD-specific) genes with consistent sign of foldchange across all the consecutive DEAs.
# Usage: R ppmi_deseq_TS_consistent_trend.R
# Data: data resulting from DESEQ at consecutive timepoints (DESeq_res_t1t2)

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

DEA.FILES <- list()
labels <- c("l1","l2")
DEA.FILES[[labels[1]]] <- list.files(path = OUT_DIR, pattern = "^DESeq_res_.+_PD.tsv$", full.names = TRUE) # please note it's important that l1 refers to the patient/symptom-specific class
DEA.FILES[[labels[2]]] <- list.files(path = OUT_DIR, pattern = "^DESeq_res_.+_HC.tsv$", full.names = TRUE) 



# Data load --------------------------------------------------------------------

DEAS <- list()
for (class in labels) {
  DEAS[[class]] <- list()
  for (f in DEA.FILES[[class]]) {
    vs <- str_split(f, pattern = "_")[[1]][3]
    label <- str_split(str_split(f, pattern = "_")[[1]][4], pattern = "\\.")[[1]][1]
    f.name <- paste(vs, label, sep = "_")
    DEAS[[class]][[f.name]] <- vroom(f, col_types = cols())
  }
}




# Data joins & transformations -------------------------------------------------
# Identify genes with consistent sign of trend across consecutive timepoints

DEA_full <- list()
G.LIST <- list()
DEA_full_consistent <- list()

for  (class in labels) {
  # identify common genes in all DEAs
  GENES <- list()
  for (t in names(DEAS[[class]])) {
    GENES[[t]] <- DEAS[[class]][[t]] %>% 
      pull(geneid)
  }
  common_genes <- Reduce(intersect, GENES)
  
  # filter common genes in all DEAs & remove NAs
  flt <- c()
  for (t in names(DEAS[[class]])) {
    DEAS[[class]][[t]] <- DEAS[[class]][[t]] %>% 
      dplyr::filter(geneid %in% common_genes) %>% 
      dplyr::select(geneid, log2FoldChange, pvalue, padj) %>%
      as.data.frame()
    flt <- c(flt, DEAS[[class]][[t]][is.na(DEAS[[class]][[t]]$padj),]$geneid) # append all geneids with NAs in padj variable
  }
  flt <- unique(flt)
  DEAS[[class]][[1]] <- DEAS[[class]][[1]] %>% 
    dplyr::filter(!geneid %in% flt)
  
  # join all DEAs' results
  DEA_full[[class]] <- Reduce(function(...) inner_join(..., by = 'geneid', suffix = c("_1", "_2")), DEAS[[class]]) %>%
    dplyr::rename(log2FoldChange_3 = log2FoldChange,
                  pvalue_3 = pvalue,
                  padj_3 = padj) 
  
  # filter genes with consistent sign of foldchange across the 3 DEAs
  DEA_full_consistent[[class]] <- DEA_full[[class]] %>%
    mutate(sum_trend = rowSums(sign(DEA_full[[class]] %>% dplyr::select(matches('^log2FoldChange_\\d$')))),
           up = (sum_trend == -3), # note log2foldchange(t1,t2) < 0 ==> t1 < t2 ==> trend is upwards ; log2foldchange(t1,t2) > 0 ==> t1 > t2 ==> trend is downwards
           down = (sum_trend == 3),
           is_consistent = (sum_trend %in% c(3, -3))) %>%
    filter(is_consistent == TRUE) %>% 
    rowwise() %>%
    mutate(pvalue_median = median(c_across(matches("pvalue"))),
           padj_median = median(c_across(matches("adj"))),
           log2FoldChange_median = median(c_across(matches("log2FoldChange_")))) %>%
    dplyr::arrange(pvalue_median, desc(abs(log2FoldChange_median)))
  # store list of consistent genes
  G.LIST[[class]] <- DEA_full_consistent[[class]]$geneid
}



# compare lists in both classes and their sign of the trend --------------------

# elements in class l1 but not in class l2. 
gene_spec <- Reduce(setdiff, G.LIST) 

# elements in common between classes with different trends. 
common_genes <- Reduce(intersect, G.LIST)
gene_diftrend <- DEA_full_consistent[[labels[1]]] %>%
  filter(geneid %in% common_genes) %>%
  dplyr::select(c(geneid, up, down)) %>%
  inner_join(DEA_full_consistent[[labels[2]]] %>%
               filter(geneid %in% common_genes) %>%
               dplyr::select(c(geneid, up, down)),
             by = 'geneid',
             suffix = paste0("_", labels)) %>%
  filter(((get(paste0("up_", labels[1])) == TRUE) & (get(paste0("up_", labels[2])) == FALSE)) | ((get(paste0("up_", labels[1])) == FALSE) & (get(paste0("up_", labels[2])) == TRUE))) %>%
  pull(geneid)

gene_spec <- c(gene_spec, gene_diftrend)

# subset identified class 1-specific genes from dataset
DEA_full_consistent[[labels[1]]] <- DEA_full_consistent[[labels[1]]] %>%
  filter(geneid %in% gene_spec) %>%
  dplyr::select(-c(is_consistent, sum_trend)) %>% 
  filter(pvalue_median < 0.05) # significancy filter: pvalue < 0.05



# output results ---------------------------------------------------------------
readr::write_tsv(DEA_full_consistent[[labels[1]]], file = file.path(OUT_DIR, "genes_consistent_trend_TS.tsv"))



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

