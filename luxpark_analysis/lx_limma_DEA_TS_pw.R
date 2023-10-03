# Title: lx_limma_DEA_TS_pw.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs DA on metabolomics data at pathway level via limma.
# Usage: Rscript lx_limma_DEA_TS_pw.R  
# Data: data from metabolites abundance, clinical data.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(edgeR)
library(limma)
library(argparser)
library(caret)
library(VIM)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
IN_DIR <- "../data/00-cleansing/"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level")
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")
PHENO.FILE <- file.path(IN_DIR, "lx_pheno.tsv")
METAB.FILE <- file.path(IN_DIR, "log_transformed_data_fte.tsv") # to be change if de novo!
M1342.FILE <- file.path(IN_DIR, "M1342.tsv") # M1342 aka 3-methoxytyrosine: confounder
target = "DIAGNOSIS" # cmd line arg?

AGG.METAB.FILE <- list()
for (st in c("mean", "median", "sd", "pathifier", "pca")) {
  AGG.METAB.FILE[[st]] <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv"))
}



# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------------
annotation <- vroom(ANNOTATION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = cols())
metabolites <- vroom(METAB.FILE, col_types = cols())
M1342_df <- vroom(M1342.FILE, col_types = cols()) 



# Processing ---------------------------------------------------------------
var_id = "PATHWAY_NAME"

# kNN imputation for BMI variable
pheno <- VIM::kNN(pheno, variable = "BMI", k= 5, imp_var = F)

# define timepoints
visits <- sort(unique(pheno$VISIT))
visits <- visits[-length(visits)] # remove V4 as it only has 2 observations! maybe v3 should also be removed (only 10 observations!)
v1 = visits[1:length(visits) -1]
v2 = visits[2:length(visits)]

METAB_class <- list()
PHENO_class <- list()
DEAS <- list()
class = "PD"  # in luxpark there is only PD samples for longitudinal data
for (st in names(AGG.METAB.FILE)) { 
  DEAS[[st]] <- list()
  metab <- vroom(AGG.METAB.FILE[[st]], col_types = cols(), delim = "\t") 
  
  for (v in seq(1:length(v1))) {
    # identify individuals by class and visit
    class.flt <- pheno %>% 
      dplyr::filter(VISIT %in% c(v1[v], v2[v])) %>% 
      dplyr::filter(get(target) == class) %>% 
      pull("SAMPLE_ID") %>% 
      unique()
    PHENO_class[[paste0(v1[v],v2[v])]] <- pheno %>% 
      dplyr::filter(VISIT %in% c(v1[v], v2[v])) %>% 
      filter(get(target) == class) 
    METAB_class[[paste0(v1[v],v2[v])]] <- metab %>%
      dplyr::filter(VISIT %in% c(v1[v], v2[v])) %>%
      dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>%
      dplyr::filter(get("SAMPLE_ID") %in% class.flt) 

    if (nrow(METAB_class[[paste0(v1[v],v2[v])]]) > 2) {
      # Pre-processing ---------------------------------------------------------------
      
      # remove near Zero variance PW
      nzv = nearZeroVar(METAB_class[[paste0(v1[v],v2[v])]], names = TRUE)
      if (length(nzv) > 0) {
        METAB_class[[paste0(v1[v],v2[v])]] <- METAB_class[[paste0(v1[v],v2[v])]] %>%
          dplyr::select(-any_of(nzv))
      }
      
      # confounding factor 3-methoxytyrosine
      M1342 <- M1342_df %>%
        filter(SAMPLE_ID %in% METAB_class[[paste0(v1[v],v2[v])]]$SAMPLE_ID) %>%
        pull(M1342)
      
      # flip the dataset as rows should correspond to metabolites and columns to samples
      METAB_class[[paste0(v1[v],v2[v])]] <- METAB_class[[paste0(v1[v],v2[v])]] %>% 
        pivot_longer(!SAMPLE_ID, names_to = var_id, values_to = "COUNT") %>% 
        pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
        column_to_rownames(var_id)
      
      # Apply Bayesian Moderated t-statistic -----------------------------------------
      # empirical bayes t-test: visit + confounding factors
      design <- model.matrix(~PHENO_class[[paste0(v1[v],v2[v])]]$VISIT 
                             + PHENO_class[[paste0(v1[v],v2[v])]]$AGE 
                             + PHENO_class[[paste0(v1[v],v2[v])]]$GENDER 
                             + PHENO_class[[paste0(v1[v],v2[v])]]$BMI 
                             + M1342) 
      colnames(design) <- c(v1[v], v2[v], "AGE", "GENDER", "BMI", "M1342")
      contrast.matrix <- limma::makeContrasts(paste0(v1[v], '-', v2[v]), levels = design) # comparison for visits
      
      xfit <- limma::lmFit(METAB_class[[paste0(v1[v],v2[v])]], design)
      xfit <- limma::contrasts.fit(xfit, contrast.matrix)
      ebayes <- limma::eBayes(xfit)
      lm_summary <- summary(decideTests(ebayes))
      print(lm_summary)
      DEAS[[st]][[paste0(v1[v],v2[v])]] <- topTable(ebayes, coef=1, number = nrow(ebayes), sort.by = "P")
      DEAS[[st]][[paste0(v1[v],v2[v])]] <- as.data.frame(DEAS[[st]][[paste0(v1[v],v2[v])]]) %>%
        rownames_to_column(var = var_id)
      print(head(DEAS[[st]][[paste0(v1[v],v2[v])]]))
      
      
      # output results
      readr::write_tsv(DEAS[[st]][[paste0(v1[v],v2[v])]], file = file.path(OUT_DIR_PATHWAY, paste0("lt_PW_limma_", st, "_", v1[v], v2[v], "_", class, ".tsv")))
      print(paste("limma analysis for", v1[v], "-", v2[v], "-", class, "completed"))
      
    } else {
      message(paste("limma analysis for", v1[v], "-", v2[v], "-", class, "failed: Need at least two features to fit a model"))
      next
    }
  }
  
}
      



# Data joins & transformations -------------------------------------------------
# Identify genes with consistent sign of trend across consecutive timepoints

for (st in names(AGG.METAB.FILE)) { 
    
  # identify common genes to all DE analyses & remove NAs
  GENES <- list()
  for (t in names(DEAS[[st]])) {
    GENES[[t]] <- DEAS[[st]][[t]] %>% 
      pull(var_id)
  }
  common_genes <- Reduce(intersect, GENES)  
  flt <- c()
  for (t in names(DEAS[[st]])) {
    DEAS[[st]][[t]] <- DEAS[[st]][[t]] %>% 
      dplyr::filter(get(var_id) %in% common_genes) %>% 
      dplyr::select(!!var_id, logFC, P.Value, adj.P.Val) %>%
      as.data.frame()
    flt <- c(flt, DEAS[[st]][[t]][is.na(DEAS[[st]][[t]]$adj.P.Val),][[var_id]]) # append all metabolites with NAs in padj variable
  }
  flt <- unique(flt)
  DEAS[[st]][[1]] <- DEAS[[st]][[1]] %>% 
    dplyr::filter(!var_id %in% flt)
  
  # join all DEAs' results
  DEA_full <- Reduce(function(...) inner_join(..., by = var_id, suffix = c("_1", "_2")), DEAS[[st]]) %>%
    dplyr::rename(logFC_3 = logFC,
                  P.Value_3 = P.Value,
                  adj.P.Val_3 = adj.P.Val) 
  
  # filter genes with consistent sign of foldchange across the 3 DEAs
  DEA_consistent <- DEA_full %>%
    mutate(sum_trend = rowSums(sign(DEA_full %>% dplyr::select(matches('^logFC_\\d$')))),
           up = (sum_trend == -3), # note logFC(t1,t2) < 0 ==> t1 < t2 ==> trend is upwards ; logFC(t1,t2) > 0 ==> t1 > t2 ==> trend is downwards
           down = (sum_trend == 3),
           is_consistent = (sum_trend %in% c(3, -3))) %>%
    filter(is_consistent == TRUE) %>% 
    rowwise() %>%
    dplyr::mutate(pvalue_median = median(c_across(matches("P.Value"))),
           padj_median = median(c_across(matches("adj"))),
           logFC_median = median(c_across(matches("logFC_")))) %>%
    dplyr::arrange(pvalue_median, desc(abs(logFC_median))) %>%
    rename_with(toupper)
  
  
  
  # output results ---------------------------------------------------------------
  readr::write_tsv(DEA_consistent, file = file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_consistent_trend_TS.tsv")))

}




# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()




  
  
  

    


