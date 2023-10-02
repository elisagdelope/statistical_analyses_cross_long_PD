# Title: ppmi_filter_gene_expression.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script filters out low expression genes from raw gene expression counts. 
# Usage: R ppmi_filter_gene_expression.R 
# Usage: R ppmi_filter_gene_expression.R -v "BL" -a "01-dea-BL-PD" -p "yes" # for timepoint-specific analysis!
# Data: data from expression at gene level (STAR) and clinical (pheno) data.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(DESeq2)
library(edgeR)
library(argparser)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT.EXPRESSION.FILE <- file.path(OUT_DIR,"flt_star_all_TS.tsv")
STAR.FILE <- "../data/01-dea-PDHC/02-outfiles/fc_counts_dsq.tsv"
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"

# Add command line arguments
p <- arg_parser("Filters expression data", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--visit", help = "visit/timepoint name to perform deseq on", default = "ALL", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-TS-PD", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--print", help = "boolean string (yes/no) whether to print time-visit extract from original data (no filtered)", default = "NO", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable
analysis_name <- argv$analysis # analysis name variable
print_extract <- toupper(argv$print) # print visit extract

if (t %in% c("BL", "V04", "V06", "V08")) {
  OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
  OUT.EXPRESSION.FILE <- file.path(OUT_DIR, paste0("flt_star_", t, ".tsv"))
  OUT.PHENO.FILE <- paste("ppmi_pheno", t, "dsq.tsv", sep = "_") 
} else if (t == "ALL") {
  OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
  OUT.EXPRESSION.FILE <- file.path(OUT_DIR,"flt_star_all_TS.tsv")
} else if (t != "ALL") {
  stop("Adequate arguments were not provided. Check R ppmi_filter_gene_expression.R --help for the right usage.")
} 





# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = T)
}



# Data load --------------------------------------------------------------------

expression <- vroom(STAR.FILE, col_types = cols())  %>%
  mutate(geneid = str_replace(geneid, '\\.[0-9]*$', '')) %>%
  column_to_rownames(var = "geneid")
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) 




# When timepoint-specific analysis ---------------------------------------------
if (t != "ALL") {
  expression <- expression %>%
    dplyr::select(matches(t))
  pheno <- pheno %>% 
    dplyr::filter(visit == t)
  
  if (print_extract == "YES") {
    readr::write_tsv(expression %>% 
                       rownames_to_column(var = "geneid"), file = file.path(OUT_DIR, paste0("fc_counts_", t, "_dsq.tsv"))) 
    readr::write_tsv(pheno, file = file.path(OUT_DIR, OUT.PHENO.FILE)) 
  }
  
  message(paste("Data filtered by timepoint", t))
}



# Filter expression data -------------------------------------------------------
# filtering rows with no or nearly no information about the amount of gene expression with edgeR::filterByExpr()

expression.deseq <- DESeq2::DESeqDataSetFromMatrix(expression, 
                                                   colData = pheno, 
                                                   design = ~ Age + Gender) # only account for confounders
expression.dgelist <- DGEList(counts(expression.deseq), group = colData(expression.deseq)$Diagnosis) #  The function accesses the group factor to compute the minimum group size, but the filtering is performed independently of which sample belongs to which group so that no bias is introduced
gene.flt <- edgeR::filterByExpr(expression.dgelist, min.count = 10)
expression.flt <- expression.deseq[gene.flt, ]

expression <- as.data.frame(assay(expression.flt)) %>%
  rownames_to_column(var = "geneid")



# output results ---------------------------------------------------------------
readr::write_tsv(expression, file = OUT.EXPRESSION.FILE)



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

