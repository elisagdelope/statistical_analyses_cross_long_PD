# Title: ppmi_generate_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates expression data at aggregated level (mean, median & st) across gene members of aggregations: GO BP, GO CC, CORUM levels.
# Usage: R ppmi_generate_pathway_level.R
# Usage: R ppmi_generate_pathway_level.R -v "BL" -a "01-dea-BL-PD" # for timepoint-specific analysis!
# Data: data from expression at gene level (STAR), entrezID/ensemblID correspondance, CORUMdb file, links from GOBP, GOCC and CORUM databases.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)



# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(vroom)
library(stringr)
library(tibble)
library(pathifier)
library(argparser)



# I/O --------------------------------------------------------------------------
# Add command line arguments
p <- arg_parser("Generates aggregated expression statistics", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--visit", help = "visit/timepoint name to perform deseq on", default = "ALL", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-TS-PD", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
t <- toupper(argv$visit) # visit name variable
analysis_name <- argv$analysis # analysis name variable

REFERENCES_DIR <- "../references"
ENSEMBL.FILE <- sort(list.files(REFERENCES_DIR, pattern = "^ids2ensembldb_.+\\.tsv"), decreasing = TRUE)[1]  # take newest version
CORUM.DB.FILE <- file.path(REFERENCES_DIR, "CORUMdb.txt.zip")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level")
FLT_NORM_EXPRESSION.FILE <- file.path(OUT_DIR, "flt_norm_star_all_TS.tsv")

if (t %in% c("BL", "V04", "V06", "V08")) {
  FLT_NORM_EXPRESSION.FILE <- file.path(OUT_DIR, paste0("flt_norm_star_", t, ".tsv"))
} else if (t != "ALL") {
  stop("Adequate arguments were not provided. Check R ppmi_filter_gene_expression.R --help for the right usage.")
} 

source("func_download_ensembldb.R")
source("func_summarize_pathway_level.R")

LINKDB <- list()
LINKDB[["GOBP"]] <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/c5.go.bp.v7.4.entrez.gmt"
LINKDB[["GOCC"]] <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.4/c5.go.cc.v7.4.entrez.gmt"
LINKDB[["CORUM"]] <- "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip"



# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}



# Data load --------------------------------------------------------------------

expression <- vroom(FLT_NORM_EXPRESSION.FILE, col_types = cols()) %>%
  rename_with(toupper)

# Entrez - ensembl conversion db
if (file.exists(paste0("../references/", ENSEMBL.FILE))) {
  ids2ensembldb <- vroom(paste0("../references/", ENSEMBL.FILE), col_types = cols())
} else {
  ids2ensembldb <- download_ids2ensembldb()
}
message("Ensembl ID names loaded")

if (file.exists(CORUM.DB.FILE)) {
  corum_db <- vroom(CORUM.DB.FILE, col_types = cols())
} else {
  download.file(LINKDB[["CORUM"]], CORUM.DB.FILE, method = "curl")
  corum_db <- vroom(CORUM.DB.FILE, col_types = cols())
}
message("CORUM db loaded")



# Aggregated statistics of expression across gene members of cellular pathways: GO BP, GO CC, CORUM -------------

# get entrez ids and set them as rownames; removes entries with no entrezID
expression <- expression %>%
#  mutate(gene_id = str_replace(geneid, '\\.[0-9]*$', '')) %>% 
  mutate(ENTREZID = ids2ensembldb$Entrez.ID[match(GENEID, ids2ensembldb$Ensembl.ID)]) %>%
  filter(!is.na(ENTREZID)) %>%
  column_to_rownames(var = "ENTREZID") %>%
  dplyr::select(matches('^[0-9]{4}\\.[A-Z]+'))


# Compute aggregated statistics for each db: mean, median & st
DB <- list()
AGG.EXPRESSION <- list()
for (db in c("GOCC", "CORUM", "GOBP")) {
  if (db == "CORUM") {
    downdb <- as.data.frame(corum_db)
    downdb <- as.data.frame(downdb) %>%
      dplyr::select(ComplexName, `subunits(Entrez IDs)`) %>%  # select(matches('subunits|^Complex[A-Z]+')) %>%
      group_by(ComplexName) %>%
      summarise(entrez_ids = paste(`subunits(Entrez IDs)`, collapse = ";")) # collapse all subunits (genes) from complex having the same complex name
    DB[[db]] <- list()
    for (i in seq(1:length(downdb$ComplexName))) {
      DB[[db]][i] <- downdb[i,]$entrez_ids
      DB[[db]][i] <- strsplit(DB[[db]][[i]],";")
      DB[[db]][[i]] <- unique(DB[[db]][[i]])
    }    
    names(DB[[db]]) <- downdb$ComplexName
    
  } else {
    downdb <- sapply(readLines(LINKDB[[db]]), function(x) strsplit(x, "\t")[[1]])
    DB[[db]] <- sapply(as.matrix(downdb), function(x) x[3:length(x)])
    names(DB[[db]]) = sapply(as.matrix(downdb), function(x) x[1])
  }
  
  AGG.EXPRESSION[[db]] <- list()
  for (st in c("mean", "median", "sd", "pathifier", "pca")) {
    var_id <- paste0(db, "_NAME")
    
    if (st == "pathifier") {
      pathifier_agg <- quantify_pathways_deregulation(data = as.matrix(expression),
                                                      allgenes = rownames(expression),
                                                      syms = DB[[db]], 
                                                      pathwaynames = names(DB[[db]]),
                                                      normals = NULL, 
                                                      logfile = file.path(OUT_DIR_PATHWAY, paste("ppmi", st, db, "log", sep = ".")), 
                                                      attempts = 5)
      pathifier_scores <- data.frame(Reduce(rbind, pathifier_agg$scores))
      colnames(pathifier_scores) <- colnames(expression) 
      rownames(pathifier_scores) <- names(pathifier_agg$scores)
      AGG.EXPRESSION[[db]][[st]] <- pathifier_scores
      
     } else {
      AGG.EXPRESSION[[db]][[st]] <- summarize_pathway_level(expression, DB[[db]], type = st)
     }
    
    # output results
    readr::write_tsv(as.data.frame(AGG.EXPRESSION[[db]][[st]]) %>%
                       rownames_to_column(var = var_id), file = file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_")))
    
  }
  
}




# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

