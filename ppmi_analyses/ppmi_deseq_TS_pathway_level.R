# Title: ppmi_deseq_TS_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script Prepares data and performs DESEQ analysis on consecutive timepoints (t1 vs t2; t2 vs t3; t3 vs t4) for visit with Age & Gender as covariates separately for binary classes.
# Usage: R ppmi_deseq_TS_pathway_level.R diagnosis -l gobp -s sd 
# Data: data from the gene expression + pheno (ready to deseq files) of patients(PD) & controls(HC)

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
library(argparser, quietly = TRUE)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"

# Add command line arguments
p <- arg_parser("Prepares data and performs DEAs on consecutive timepoints for binary classes", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "target", help = "target variable (diagnosis, updrs, ...)", default = "diagnosis", type = "string", short = "t")
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
target = toupper(argv$target) # target variable
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level


EXPRESSION.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
OUT.FILE.STARTWITH <- file.path(OUT_DIR_PATHWAY, paste("DESeq_res", e_level, st, sep = "_"))
if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier", paste(c("mean", "median", "sd", "pathifier"), "lasso", sep = "_"))) | (!file.exists(EXPRESSION.FILE))) { 
  stop("Adequate arguments were not provided. Check R ppmi_deseq_TS_pathway_level.R --help for the right usage.")
}


# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------

var_id = paste0(e_level, "_name")

pheno_dsq <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl"))
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>% 
  column_to_rownames(var = var_id) 



# Subset data per consecutive timepoints ---------------------------------------
visits <- str_extract(sort(colnames(expression)), '[A-Z0-9]*$')
visits <- unique(visits[(!visits == "") & (!is.na(visits))])
v1 = visits[1:length(visits) -1]
v2 = visits[2:length(visits)]

STAR <- list()
PHENO <- list()
for (v in seq(1:length(v1))) {
  # subset expression data for consecutive visits
  PHENO[[paste0(v1[v],v2[v])]] <- pheno_dsq %>% 
    dplyr::filter(visit %in% c(v1[v], v2[v]))
  STAR[[paste0(v1[v],v2[v])]] <- expression %>%
    dplyr::select(matches(paste0(v1[v],'|', v2[v]))) 
}
rm(expression)



# DEA of consecutive visits for PD/HC ------------------------------------------

pheno_dsq <- pheno_dsq %>% 
  rename_with(toupper)

STAR_class <- list()
PHENO_class <- list()
E.STAR <- list()
for (class in as.character(unique(pheno_dsq[[target]]))) {
  # subset consecutive visits by class and DESEQ
  DDS <- list()
  RES <- list()
  
  for (v in seq(1:length(v1))) {
    # identify individuals by class
    class.flt <- PHENO[[paste0(v1[v],v2[v])]] %>% 
      rename_with(toupper) %>% 
      filter(get(target) == class) %>% 
      pull(GSEID) %>% 
      unique()
    PHENO_class[[paste0(v1[v],v2[v])]] <- PHENO[[paste0(v1[v],v2[v])]] %>% 
      rename_with(toupper) %>% 
      filter(get(target) == class) 
    STAR_class[[paste0(v1[v],v2[v])]] <- STAR[[paste0(v1[v],v2[v])]] %>%
      dplyr::select(matches(class.flt)) 
    
    if ((st == "pathifier") & (e_level != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
      E.STAR[[paste0(v1[v],v2[v])]] <- DESeq2::DESeqDataSetFromMatrix(round(STAR_class[[paste0(v1[v],v2[v])]] * 100), 
                                                                      colData = PHENO_class[[paste0(v1[v],v2[v])]], 
                                                                      design = ~ AGE + GENDER + VISIT)
      
    } else {
      E.STAR[[paste0(v1[v],v2[v])]] <- DESeq2::DESeqDataSetFromMatrix(round(STAR_class[[paste0(v1[v],v2[v])]]), 
                                                                      colData = PHENO_class[[paste0(v1[v],v2[v])]], 
                                                                      design = ~ AGE + GENDER + VISIT)
      
    }
    
    # pre-filtering has been already applied with edgeR::filterbyExpr but just in case:
    keep <- rowSums(counts(E.STAR[[paste0(v1[v],v2[v])]])) >= 10
    E.STAR[[paste0(v1[v],v2[v])]] <- E.STAR[[paste0(v1[v],v2[v])]][keep,]
    
    # deseq ATTENTION: DATA ALREADY NORMALIZED
#    DDS[[paste0(v1[v],v2[v])]] <- DESeq(E.STAR[[paste0(v1[v],v2[v])]])
#    DDS <- estimateSizeFactors(E.STAR)
    sizeFactors_nonnorm <- rep(1, ncol(STAR_class[[paste0(v1[v],v2[v])]])) # create a vector of 1s to avoid normalization
    names(sizeFactors_nonnorm) <- colnames(STAR_class[[paste0(v1[v],v2[v])]])
    sizeFactors(E.STAR[[paste0(v1[v],v2[v])]]) <- sizeFactors_nonnorm
    DDS[[paste0(v1[v],v2[v])]] <- estimateDispersions(E.STAR[[paste0(v1[v],v2[v])]])
    DDS[[paste0(v1[v],v2[v])]] <- nbinomWaldTest(DDS[[paste0(v1[v],v2[v])]], maxit=500)
    
    res.flt <- results(DDS[[paste0(v1[v],v2[v])]], contrast = c("VISIT", v1[v], v2[v])) #, independentFiltering = FALSE)?? -> en ppio no 
    res.flt <- as.data.frame(res.flt) %>% 
      rownames_to_column(var = var_id)
    RES[[paste0(v1[v],v2[v])]] <- res.flt[order(res.flt$pvalue),]
    
    # output results
    readr::write_tsv(as.data.frame(RES[[paste0(v1[v],v2[v])]]), file = paste0(OUT.FILE.STARTWITH, "_", v1[v], v2[v], "_", class, ".tsv"))
    print(paste("DESeq analysis for", v1[v], "-", v2[v], "-", class, "completed"))
  }
  
  # merge consecutive DEAs of each class in the same object
  e.deseq <- data.frame()
  e.deseq <- list(RES[[paste0(v1[1],v2[1])]] %>%
                     rename_at(vars(-any_of(var_id)),function(x) paste0(x,"_",paste0(v1[1],v2[1]))),
                   RES[[paste0(v1[2],v2[2])]] %>%
                     rename_at(vars(-any_of(var_id)),function(x) paste0(x,"_",paste0(v1[2],v2[2]))), 
                   RES[[paste0(v1[3],v2[3])]] %>%
                     rename_at(vars(-any_of(var_id)),function(x) paste0(x,"_",paste0(v1[3],v2[3])))) %>%
    purrr::reduce(left_join, by = var_id)
  
  # output merged results 
  readr::write_tsv(e.deseq, file = file.path(OUT_DIR_PATHWAY, paste(e_level, st, "deseq", paste0(class, ".tsv"), sep = "_")))
  print(paste("deseq analysis for", class, "in", e_level, st, "completed"))
  
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


