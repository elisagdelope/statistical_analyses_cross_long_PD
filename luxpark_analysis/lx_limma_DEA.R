# Title: lx_limma_DEA.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs DA on metabolomics data via limma.
# Usage: Rscript lx_limma_DEA.R -a "01-dea-V0-PD"
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
# Add command line arguments
p <- arg_parser("Generates aggregated metabolomics statistics", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-V0-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
analysis_name <- argv$analysis # analysis name variable

IN_DIR <- "../data/00-cleansing/"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")
M1342.FILE <- file.path(IN_DIR, "M1342.tsv") # M1342 aka 3-methoxytyrosine: confounder
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv")
METAB.FILE <- file.path(OUT_DIR, "log_transformed_V0.tsv")
OUT.DEA.FILE <- file.path(OUT_DIR, "lt_limma_V0.tsv")
 


# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# Data load --------------------------------------------------------------------
annotation <- vroom(ANNOTATION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 
metab <- vroom(METAB.FILE, col_types = cols())
metab <- metab %>%
  dplyr::select(-any_of(c("PATIENT_ID", "VISIT")))



# Pre-processing ---------------------------------------------------------------
var_id = "METABOLITES_ID"

# remove near Zero variance metabolites
nzv = nearZeroVar(metab, names = TRUE)
if (length(nzv) > 0) {
  metab <- metab %>%
    dplyr::select(-any_of(nzv))
}

if (!str_detect(toupper(analysis_name), "NOVO")) {
  # confounding factor 3-methoxytyrosine
  M1342 <- vroom(M1342.FILE, col_types = cols()) 
  M1342 <- M1342 %>%
    filter(SAMPLE_ID %in% metab$SAMPLE_ID) %>%
    pull(M1342)
}

# flip the dataset as rows should correspond to metabolites and columns to samples
metab.t = metab %>% 
  pivot_longer(!SAMPLE_ID, names_to = var_id, values_to = "COUNT") %>% 
  pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
  column_to_rownames(var_id)

# kNN imputation for BMI variable
pheno <- VIM::kNN(pheno, variable = "BMI", k= 5, imp_var = F)



# Apply Bayesian Moderated t-statistic -----------------------------------------
# empirical bayes t-test
if (!str_detect(toupper(analysis_name), "NOVO"))  {
  design <- model.matrix(~0 + pheno$DIAGNOSIS + pheno$AGE + pheno$GENDER + pheno$BMI + M1342) # diagnosis + confounding factors
  colnames(design) <- c(levels(pheno$DIAGNOSIS), "AGE", "GENDER", "BMI", "M1342")
} else {
  design <- model.matrix(~0 + pheno$DIAGNOSIS + pheno$AGE + pheno$GENDER + pheno$BMI) # diagnosis + confounding factors
  colnames(design) <- c(levels(pheno$DIAGNOSIS), "AGE", "GENDER", "BMI")
}
contrast.matrix <- limma::makeContrasts('PD-HC', levels = design) # comparison for diagnosis, important to specify the right order PD/HC

xfit <- limma::lmFit(metab.t, design)
xfit <- limma::contrasts.fit(xfit, contrast.matrix)
ebayes <- limma::eBayes(xfit)
lm_summary <- summary(decideTests(ebayes))
print(lm_summary)
lm_res <- topTable(ebayes, coef=1, number = nrow(ebayes), sort.by = "P")
lm_res <- lm_res[((lm_res$adj.P.Val < 0.05) & (!is.na(lm_res$adj.P.Val))),]  # remove genes with padj = NA / pvalue > 0.05
lm_res <- lm_res %>%
  rownames_to_column(var = var_id)
print(head(lm_res))

# export outputs  --------------------------------------------------------------
readr::write_tsv(lm_res, file = OUT.DEA.FILE) 



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

