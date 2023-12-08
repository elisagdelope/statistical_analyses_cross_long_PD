# Title: ppmi_deseq_salmon_star.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script processes data and performs DESEQ analysis for STAR+featureCounts and Salmon quantification pipelines on all samples (PD/HC with Age & Gender as covariates)
# Usage: R ppmi_deseq_salmon_star.R
# Data: data from the quantification file at gene level, both STAR+FC & Salmon + pheno

# GC ----------------------------------------------------------------------
rm(list=ls())
gc(T)


# Packages ----------------------------------------------------------------
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(DESeq2)
library(edgeR)
library(SummarizedExperiment)


# I/O ---------------------------------------------------------------------
analysis_name <- "01-dea-PDHC"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles") 
OUT_DIR_ENRICHMENT <- paste0("../data/", analysis_name, "/03-enrichment_tables")
PHENO.FILE <- "../data/00-cleansing/ppmi_pheno.tsv"
SALMON.FILE <- "../data/00-cleansing/ppmi_rnaseq_salmon_genequant_counts.tsv"
STAR.FILE <- "../data/00-cleansing/ppmi_rnaseq_fc_counts.tsv"
# ENSEMBL.FILE <- sort(list.files("../references/", pattern = "^ids2ensembldb_.+\\.tsv"), decreasing = TRUE)[1]  # take newest version

# Main --------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS)) | (!dir.exists(OUT_DIR_ENRICHMENT))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
  dir.create(OUT_DIR_ENRICHMENT, recursive = T)
}

# Data load --------------------------------------------------------------
pheno <- vroom(PHENO.FILE, col_types = "cffdfildddddddddl")
salmon_genequant <- vroom(SALMON.FILE, col_types = cols())
fc_counts <- vroom(STAR.FILE, col_types = cols())


# prepare data for deseq analysis ----------------------------------------

# match pheno, quant and count data -> all same ids should be in all datasets but unfortunately they are not.
# assign rownames to gene ids, remove extra cols in fc counts and set salmon quants as integer

pheno <- pheno %>%
  unite("patient_visit_id", c(GSEID, visit), sep = ".", remove = FALSE)

ids_phenoquant <- sort(colnames(salmon_genequant)[match(pheno$patient_visit_id, colnames(salmon_genequant))])
ids_phenocounts <- sort(colnames(fc_counts)[match(pheno$patient_visit_id, colnames(fc_counts))])

if ((length(ids_phenoquant) != length(pheno$patient_visit_id)) || (length(ids_phenocounts) != length(pheno$patient_visit_id))) {
  message("###--------- ALERT: IDs in phenotype are not matching those in fc_counts or salmon quants-----------###")
  if (!identical(ids_phenoquant, ids_phenocounts)) {
    message("###--------- ALERT: IDs in fc_counts and salmon quants are not matching-----------###")
    message("###--------- ALERT: Execution will stop -----------###")
    stop()
  } else {
    message("###--------- ALERT: Non matching IDs will be removed from the corresponding datasets -----------###")
    
    pheno_dsq <- pheno %>%
      filter(patient_visit_id %in% ids_phenocounts)
    
    salmon_genequant_dsq <- salmon_genequant %>%
      dplyr::select(id, ids_phenoquant) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "id") %>%
      mutate(across(where(is.double), as.integer))
    
    fc_counts_dsq <- fc_counts %>%
      dplyr::select(Geneid, ids_phenocounts) %>% 
      remove_rownames %>% 
      column_to_rownames(var = "Geneid") %>%
      mutate(across(where(is.double), as.integer))
  }
} else {
  message("Data in phenotype, fc_counts and salmon quants is matching")
  
  salmon_genequant_dsq <- salmon_genequant %>% 
    remove_rownames %>% 
    column_to_rownames(var = "id") %>%
    mutate(across(where(is.double), as.integer))
  
  fc_counts_dsq <- fc_counts %>% 
    dplyr::select(Geneid, 7:ncol(.)) %>% 
    remove_rownames %>% 
    column_to_rownames(var = "Geneid") %>%
    mutate(across(where(is.double), as.integer))
}

# export dsq data
readr::write_tsv(pheno_dsq, file = file.path(OUT_DIR, "ppmi_pheno_dsq.tsv"))
readr::write_tsv(salmon_genequant_dsq %>% 
                   rownames_to_column(var = "geneid"), file = file.path(OUT_DIR, "salmon_genequant_dsq.tsv"))
readr::write_tsv(fc_counts_dsq %>% 
                   rownames_to_column(var = "geneid") , file = file.path(OUT_DIR, "fc_counts_dsq.tsv"))

message("Data prepared for DESeq analysis")
rm(pheno, salmon_genequant, fc_counts)


# filtering & mean variance plots ----------------------------------------

# filtering of rows with no or nearly no information about the amount of gene expression with edgeR::filterByExpr() 
# By default, the function keeps genes with about 10 read counts or more in a minimum number of samples, where the number of samples is chosen according to the minimum group sample size.
# mean - variance plot of filtered/unfiltered counts in salmon and star

# build DESeq2 object
E = list()
E[["salmon"]] <- DESeq2::DESeqDataSetFromMatrix(salmon_genequant_dsq, colData = pheno_dsq, design = ~ Diagnosis + Age + Gender)
E[["STAR"]] <- DESeq2::DESeqDataSetFromMatrix(fc_counts_dsq, colData = pheno_dsq, design = ~ Diagnosis + Age + Gender)

E.flt <- list()
E.flt.vst <- list()
plot_name <- "mean-variance.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 12, height = 8)
par(mfrow = c(1,2))
for (pipeline in c("STAR", "salmon")) {
  exprs <- DGEList(counts(E[[pipeline]]), group = colData(E[[pipeline]])$Diagnosis)
  v <- voom(exprs, plot = T)
  
  gene.flt <- edgeR::filterByExpr(exprs)
  E.flt[[pipeline]] <- E[[pipeline]][gene.flt, ]
  E.flt.vst[[pipeline]] <- vst(E.flt[[pipeline]], blind = FALSE)
  exprs.flt <- DGEList(assay(E.flt.vst[[pipeline]]), group = colData(E.flt.vst[[pipeline]])$Diagnosis)
  v.flt <- voom(exprs.flt, plot = T)
  
#  readr::write_tsv(as.data.frame(assay(E.flt.vst[[pipeline]])) %>%
#                     rownames_to_column(var = "geneid"), file = file.path(OUT_DIR, paste0("flt_transf_", pipeline, ".tsv")))
}
dev.off()
rm(E.flt.vst, E, salmon_genequant_dsq, fc_counts_dsq, pheno_dsq)

# correlation of counts in salmon and star + fcounts ------------------

ix <- match(rownames(E.flt[["STAR"]]), rownames(E.flt[["salmon"]]))

df <- data.frame(
  gene = rep(rownames(E.flt[["STAR"]]), ncol(E.flt[["STAR"]])),
  STAR = c(counts(E.flt[["STAR"]])),
  salmon = c(counts(E.flt[["salmon"]])[ix, ])
)

plot_name <- "correlation_STAR-salmon.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, plot_name), width = 10, height = 8)
ggplot(data = df, aes(x = STAR, y = salmon)) +
  geom_point(shape = 16, alpha = 0.1) +
  theme_bw() +
  xlab("STAR + featureCounts count") +
  ylab("salmon count") +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  geom_smooth(method = "lm", formula = 'y ~ x',  orientation = "y")
dev.off()



# DEA of PD/controls for star and salmon counts ------------------------------

DDS <- list()
RES <- list()

for (pipeline in c("STAR", "salmon")) {
  DDS[[pipeline]] <- DESeq(E.flt[[pipeline]])
  res.flt <- results(DDS[[pipeline]], contrast = c("Diagnosis", "PD", "HC"))
  res.flt$gene_id <- str_replace(rownames(res.flt), '\\.[0-9]*$', '')
#  res.flt$gene_name <- ensembl.db$hgnc_symbol[match(res.flt$gene_id, ensembl.db$ensembl_gene_id)]
  RES[[pipeline]] <- res.flt[order(res.flt$pvalue),]
  readr::write_tsv(as.data.frame(RES[[pipeline]]) %>% rownames_to_column( var = "geneid") , file = file.path(OUT_DIR, paste0("DESeq_res_", pipeline, ".tsv")))
  message(paste0("DESeq analysis for ", pipeline, " completed"))
}


# Session info ------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


