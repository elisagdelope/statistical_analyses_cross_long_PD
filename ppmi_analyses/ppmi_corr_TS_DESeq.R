# Title: ppmi_corr_TS_DESeq.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs a linear model of expression & time per diagnosis class adjusting for confounders (via deseq/limma), filters by significance of variables, and subsets identified class 1-specific (e.g. PD-specific) genes with significant association with time.
# Usage: Rscript ppmi_corr_TS_DESeq.R -l gobp -s sd
# Data: data from gene expression
# Note: A similar version of ppmi_corr_TS.R using DESeq2 for the linear model.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(tibble)
library(argparser)
library(DESeq2)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"

# Add command line arguments
p <- arg_parser("Plots expression of class 1-specific features", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level

if ((e_level == "GENE") | (e_level == "G")) { # gene level alone / lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR, "flt_norm_star_all_TS.tsv") # normalized?
  var_id = "geneid"
  OUT.FILE <- file.path(OUT_DIR, paste(e_level, "timecorr_PD_TS_DESeq.tsv", sep = "_"))
  
} else { # not gene level & not lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_")) 
  var_id = paste0(e_level, "_name")
  OUT.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "timecorr_PD_TS_DESeq.tsv", sep = "_"))
  
  if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier", paste(c("mean", "median", "sd", "pathifier"), "lasso", sep = "_"))) | (!file.exists(EXPRESSION.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_corr_TS.R --help for the right usage.")
  }
}



# Data load --------------------------------------------------------------------
expression <- vroom(EXPRESSION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("cfffdfildddddddddl"))



# Data filtering and transformations -------------------------------------------
# for each feature, generate linear regression model versus time

expr <- expression %>%
  pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "counts") %>%
  left_join(pheno, by = c("sample" = "patient_visit_id")) %>%
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) %>%
  dplyr::select(all_of(c(var_id, "sample", "counts", "timepoint", "Diagnosis", "Age", "Gender"))) %>%
  rename_with(toupper)

pheno <- pheno %>%
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) %>%
  rename_with(toupper)
 
labels <- c("PD","HC")
G.LIST <- list()
RES <- list()
for (class in labels){
  
  pheno_class <- pheno %>%
    filter(DIAGNOSIS == class)
  expression_class <- expression %>%
    column_to_rownames(var_id) %>%
    dplyr::select(matches(paste(unique(pheno_class$PATIENT_VISIT_ID), collapse = "|")))
  
  pheno_class <- pheno_class %>%
    group_by(GSEID) %>%
    mutate(SUBJECT_ID = cur_group_id()) # individual ID as int not factor (otherwise DESeq2 doesn't like it)
  
  if ((st == "pathifier") & (e_level != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
    E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression_class * 100), 
                                             colData = pheno_class, 
                                             design = ~ AGE + GENDER + TIMEPOINT) #  + SUBJECT_ID
  } else {
    E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression_class), 
                                             colData = pheno_class, 
                                             design = ~ AGE + GENDER + TIMEPOINT) #  + SUBJECT_ID
  }
  
  
  # pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
  keep <- rowSums(counts(E.STAR)) >= 10
  E.STAR <- E.STAR[keep,]
  
  # deseq ATTENTION: DATA ALREADY NORMALIZED
  #    DDS <- DESeq(E.STAR)
  #    DDS <- estimateSizeFactors(E.STAR)
  sizeFactors_nonnorm <- rep(1, ncol(expression_class)) # create a vector of 1s to avoid normalization
  names(sizeFactors_nonnorm) <- colnames(expression_class)
  sizeFactors(E.STAR) <- sizeFactors_nonnorm
  DDS <- estimateDispersions(E.STAR)
  DDS <- nbinomWaldTest(DDS, maxit=500)
  
  res <- results(DDS) #, c("DIAGNOSIS", "PD", "HC")) #, independentFiltering = FALSE)??
  RES[[class]] <- as.data.frame(res) %>% 
    rownames_to_column(var = var_id) %>%
    filter(!is.na(padj)) %>%
    arrange(pvalue)
  
  G.LIST[[class]] <- RES[[class]][which(RES[[class]]$padj < 0.1),] %>%
    pull(var_id)
}


# filter corr elements in class PD but not in class HC
PD_spec <- setdiff(G.LIST[["PD"]], G.LIST[["HC"]])  
RES[["PD"]] <- RES[["PD"]] %>%
  filter(!!as.symbol(var_id) %in% PD_spec)



# output results ---------------------------------------------------------------
readr::write_tsv(RES[["PD"]], file = OUT.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()



