# Title: ppmi_deseq_TS_consistent_enrichment.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs enrichment analysis at gene level.
# Usage: Rscript ppmi_deseq_TS_consistent_enrichment.R -a "01-dea-TS-PD" -i "../data/01-dea-TS-PD/02-outfiles/genes_consistent_trend_TS.tsv" -o "consistent_trend"
# Data: data from expression at gene level, list of genes, entrezID/ensemblID correspondance.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)


# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(meshr)
library(DESeq2)
library(edgeR)
library(argparser)
library(AnnotationHub)
library(MeSHDbi)



# I/O --------------------------------------------------------------------------
source("func_download_ensembldb.R")
source("func_enrichment.R")

# Add command line arguments
p <- arg_parser("Enrichment analysis at gene level", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--input", help = "full name of input file", default = "../data/01-dea-TS-PD/02-outfiles/genes_consistent_trend_TS.tsv", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--output", help = "pattern name of output file", default = "consistent_trend", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "name of analysis (directory name)", default = "01-dea-TS-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
IN.FILE = argv$input # input file with genes of interest
OUT.FILE = argv$output # string to tag output filename
analysis_name <- argv$analysis

OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
REFERENCES_DIR <- "../references"
FLT_EXPRESSION.FILE <- list.files(path = OUT_DIR, pattern = "^flt_star_(.{2}|.{6})\\.tsv$", full.names = TRUE) # takes flt_star_BL.tsv or flt_star_all_TS.tsv accordingly
ENSEMBL.FILE <- sort(list.files("../references/", pattern = "^ids2ensembldb_.+\\.tsv"), decreasing = TRUE)[1]  # take newest version



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS)) | (!dir.exists(OUT_DIR_ENRICHMENT))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
  dir.create(OUT_DIR_ENRICHMENT, recursive = T)
}



# Data load --------------------------------------------------------------------

c_genes <- vroom(IN.FILE, col_types = cols())
expression <- vroom(FLT_EXPRESSION.FILE, col_types = cols()) 

# Entrez - ensembl conversion db
if (file.exists(file.path(REFERENCES_DIR, ENSEMBL.FILE))) {
  ids2ensembldb <- vroom(paste0("../references/", ENSEMBL.FILE), col_types = cols())
} else {
  ids2ensembldb <- download_ids2ensembldb()
}
message("Ensembl ID names loaded")

# mesh data
ah <- AnnotationHub()
dbfile1 <- query(ah, c("MeSHDb", "MeSH.db", "v002"))[[1]]
MeSH.db <- MeSHDbi::MeSHDb(dbfile1)
dbfile2 <- query(ah, c("MeSHDb", "Homo sapiens", "v002"))[[1]]
MeSH.Hsa.eg.db <- MeSHDbi::MeSHDb(dbfile2)



# Standard Enrichment analysis -------------------------------------------------
# Test for over-representation of gene ontology (GO) terms and KEGG pathways in DE genes

# Translate Ensembl IDs to EntrezGene IDs
ix <- match(c_genes$geneid, ids2ensembldb$Ensembl.ID)
c_genes$entrezid <- ids2ensembldb$Entrez.ID[ix]
c_genes_entrez <- c_genes$entrezid #### ALERT: THERE CAN BE NAs - i.e. genes WITH NO ENTREZID
print(paste0("Selected consistent genes (n = ", length(c_genes_entrez), "; entrezID NAs = ", nrow(c_genes[is.na(c_genes$entrezid),]), ")"))
# remove NAs from list of genes
c_genes_entrez <- c_genes_entrez[!is.na(c_genes_entrez)]

ix <- match(expression$geneid, ids2ensembldb$Ensembl.ID)
expression$entrezid <- ids2ensembldb$Entrez.ID[ix]
all_genes <- expression$entrezid
rm(expression)
# mesh data
ah <- AnnotationHub()
dbfile1 <- query(ah, c("MeSHDb", "MeSH.db", "v002"))[[1]]
MeSH.db <- MeSHDbi::MeSHDb(dbfile1)
dbfile2 <- query(ah, c("MeSHDb", "Homo sapiens", "v002"))[[1]]
MeSH.Hsa.eg.db <- MeSHDbi::MeSHDb(dbfile2)
target_meshterms = c("Parkinson Disease", "Parkinson Disease, Secondary", "Parkinsonian Disorders", "Neurodegenerative Diseases")
meshParams <- new("MeSHHyperGParams", 
                  geneIds = c_genes_entrez, 
                  universeGeneIds = unique(all_genes),
                  annotation = "MeSH.Hsa.eg.db",
                  meshdb = "MeSH.db", # new requirement for BioC 3.14+
                  category = "C", 
                  database = "gendoo", 
                  pvalueCutoff = 1,
                  pAdjust = "none")

# call enrichment analysis functions
ENRICH <- list()
if (length(c_genes_entrez) > 0) {
  ENRICH[["KEGG"]] <- gsea_kegg(c_genes_entrez)
  ENRICH[["GO"]] <- gsea_go(c_genes_entrez)
  ENRICH[["MeSH"]] <- mesh_enrich(c_genes_entrez, all_genes, target_meshterms, meshParams)
}


# output results ---------------------------------------------------------------
for (i in names(ENRICH)) {
  readr::write_tsv(ENRICH[[i]], file = file.path(OUT_DIR_ENRICHMENT, paste0("enrichment_", OUT.FILE, "_", i, ".tsv")))
  print(knitr::kable(ENRICH[[i]][1:30,], row.names = F))
}



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

