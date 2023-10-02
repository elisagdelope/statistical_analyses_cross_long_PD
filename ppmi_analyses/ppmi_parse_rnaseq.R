# Title: ppmi_parse_rnaseq.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script parses both STAR+FeatureCounts and Salmon transcript counts into format for DEG discovery & does summarization from transcript level to gene level
# Usage: R ppmi_parse_rnaseq.R
# Data: data from the original quantification files (salmon) + counts files (STAR+FC) + genecode annotation

# GC ----------------------------------------------------------------------
rm(list=ls())
gc(T)


# Packages ----------------------------------------------------------------
library(GenomicFeatures)
library(readr)
library(tximport)
library(rlang)
library(ggplot2)

# # I/O ---------------------------------------------------------------------
# args = commandArgs(trailingOnly=TRUE)
# setwd(args[1])  # Set root to RNAsequencing folder!

QUANT.PATH <- "../data/PPMI_RNAseq_IR3_pat_extract/quant"  # Path to Salmon quantification files
COUNT.PATH <- "../data/PPMI_RNAseq_IR3_pat_extract/counts"  # Path to FeatureCounts files
ANNOTATION.FILE <- "../references/gencode.v29.ALLannotation21042021.gff3.gz" 
# PPMITX2TX <- "Data/00-cleansing/ppmi/ppmitx2tx.tsv"  # This file maps tx ID used by PPMI with tx ID from GENECODE
OUTPUT <- "../data/00-cleansing"

message(paste0("Set to ouput -> ", OUTPUT))

# Main --------------------------------------------------------------------
if (!dir.exists(OUTPUT)) {
  dir.create(OUTPUT, recursive = T)
}

# Get gene metadata -------------------------------------------------------
# Read gene annotations from GENCODE
if (file.exists(file.path(OUTPUT, "tx2gene.tsv"))) {
  tx2gene <- readr::read_tsv(file.path(OUTPUT, "tx2gene.tsv"))
} else {
  txdb <- GenomicFeatures::makeTxDbFromGFF(file=ANNOTATION.FILE)
  k <- biomaRt::keys(txdb, keytype = "TXNAME")
  tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")
  tx2gene <- subset(tx2gene, GENEID %in% tx2gene$GENEID[!grepl('_PAR_',  tx2gene$GENEID)])
  message("tx names were filtered to be unique")
  readr::write_tsv(tx2gene, file.path(OUTPUT, "tx2gene.tsv"))
  rm(k)
  gc()
}
message("GENCODE tx names loaded")

# # Match the tx name they used in PPMI with GENECODE tx ID
# if (file.exists(file.path(OUTPUT, "ppmitx2gene.tsv"))) {
#    ppmitx2gene <- readr::read_tsv(file.path(OUTPUT, "ppmitx2gene.tsv"))
# } else {
#   ppmitx2tx <- readr::read_tsv(PPMITX2TX)
#   ppmitx2tx2gene <- merge(tx2gene, ppmitx2tx, by="TXNAME")
#   ppmitx2tx2gene$TXNAME <- ppmitx2tx2gene$PPMI_TXNAME
#   ppmitx2tx2gene$PPMI_TXNAME <- NULL
#   ppmitx2gene <- ppmitx2tx2gene
#   rm(ppmitx2tx2gene, tx2gene, ppmitx2tx, txdb)
#   gc()
#   readr::write_tsv(ppmitx2gene, file.path(OUTPUT, "ppmitx2gene.tsv"))
# }
# message("PPMI tx names matched to GENCODE tx names")


# Extract meta data -----------------------------------------------------
# Parse file names
files <- c(list.files(QUANT.PATH), list.files(COUNT.PATH))
files <- grep(pattern = 'transcripts|featureCounts', x = files, value = TRUE) # salmon genes files will not be processed
files_split <- strsplit(files, "[.]")
data_meta <- data.frame(
  "fp"=rep(NA, length(files)),
  "type"=rep(NA, length(files)),
  "patient_id"=rep(NA, length(files)),
  "visit"=rep(NA, length(files)),
  "sample_id"=rep(NA, length(files)),
  "seq_facility"=rep(NA, length(files)),
  "pvid"=rep(NA, length(files))
)

for (i in 1:length(files)) {
  fp = files_split[[i]]
  ext <- paste0(fp[(length(fp)-2):length(fp)], collapse=".")
  if (ext=="salmon.transcripts.sf") { # Salmon transcript quant file
    data_meta$type[i] <- "salmon"
    data_meta$fp[i] <- file.path(QUANT.PATH, paste0(fp, collapse="."))
  } else if (ext=="featureCounts.GencodeV29.txt") {  # featureCounts file
    data_meta$type[i] <- "featureCounts"
    data_meta$fp[i] <- file.path(COUNT.PATH, paste0(fp, collapse="."))
  }

  if (!is.na(data_meta$type[i])) {
    data_meta$patient_id[i] <- fp[2]
    data_meta$visit[i] <- fp[3]
    data_meta$sample_id[i] <- fp[4]
    data_meta$seq_facility[i] <- fp[5]
    data_meta$pvid[i] <- paste(fp[2], fp[3], sep=".")
  }
}

rm(ext, files, fp, files_split)
gc()

data_meta <- data_meta[complete.cases(data_meta[3:7]), ]
message("Metadata extracted")


# Extract Salmon data -----------------------------------------------------
# Select only salmon data and remove pooled data
salmon_txquant_meta <- data_meta[(data_meta$type == "salmon") & (data_meta$visit!="POOL"),]

# Load RNA-seq data
salmon_txquant <- tximport::tximport(salmon_txquant_meta$fp, type = "salmon", txOut=T)  # We will simply return the standard Counts/TPM
salmon_txquant_bkp <- salmon_txquant

# Compute length scaled TPMs <- raw counts scaled by library and transcript size
lsTPM <- tximport::makeCountsFromAbundance(salmon_txquant$counts, salmon_txquant$abundance, salmon_txquant$length, countsFromAbundance = "lengthScaledTPM")

# Add row and column names
salmon_txquant$abundance <- as.data.frame(salmon_txquant$abundance) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

salmon_txquant$counts <- as.data.frame(salmon_txquant$counts) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

salmon_txquant$length <- as.data.frame(salmon_txquant$length) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

lsTPM <- as.data.frame(lsTPM) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

# Save to disk
readr::write_tsv(salmon_txquant_meta, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_meta.tsv")))
readr::write_tsv(salmon_txquant$abundance, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_txquant_abundance.tsv")))
readr::write_tsv(salmon_txquant$counts, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_txquant_counts.tsv")))
readr::write_tsv(salmon_txquant$length, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_txquant_length.tsv")))
readr::write_tsv(lsTPM, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_txquant_lsTPM.tsv")))

# Restore backup
salmon_txquant <- salmon_txquant_bkp
rm(salmon_txquant_bkp, lsTPM)
gc()

# Summarize to gene-level quant from raw transcript data (not scaled)
salmon_genequant <- tximport::summarizeToGene(salmon_txquant, tx2gene = tx2gene)  # We will simply collapse transcript into genes and return the standard Counts/TPM
rm(salmon_txquant)
gc()

# Compute length scaled TPMs at gene-level
lsTPM <- tximport::makeCountsFromAbundance(salmon_genequant$counts, salmon_genequant$abundance, salmon_genequant$length, countsFromAbundance = "lengthScaledTPM")

# test to scale first and collapse scaled transcripts to gene
# salmon_txquant_genetest <- tximport::tximport(salmon_txquant_meta$fp, type = "salmon", txOut=T, countsFromAbundance = "lengthScaledTPM")
# lsTPM_genetest <- tximport::summarizeToGene(salmon_txquant_genetest, tx2gene = tx2gene)  # We will simply collapse scaled transcript into genes 
# test worked out: same results as first summarize to gene-level + scale at gene level

# Add row and column names
salmon_genequant$abundance <- as.data.frame(salmon_genequant$abundance) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

salmon_genequant$counts <- as.data.frame(salmon_genequant$counts) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

salmon_genequant$length <- as.data.frame(salmon_genequant$length) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

lsTPM <- as.data.frame(lsTPM) %>%
  rename_with(~ salmon_txquant_meta$pvid) %>%
  rownames_to_column( var = "id") 

# Save to disk
readr::write_tsv(salmon_genequant$abundance, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_genequant_abundance.tsv")))
readr::write_tsv(salmon_genequant$counts, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_genequant_counts.tsv")))
readr::write_tsv(salmon_genequant$length, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_genequant_length.tsv")))
readr::write_tsv(lsTPM, file = file.path(OUTPUT, paste0("ppmi_rnaseq_salmon_genequant_lsTPM.tsv")))

rm(salmon_genequant, salmon_txquant_meta, lsTPM)
gc()
message("Salmon quant data parsed successfully")


# Extract featureCounts data ----------------------------------------------

# Select only featureCounts data and remove pooled data
fc_counts_meta <- data_meta[(data_meta$type == "featureCounts") & (data_meta$visit != "POOL"),]

for (i in 1:nrow(fc_counts_meta)) {
  d <- readr::read_tsv(fc_counts_meta$fp[i], skip = 1)
  colnames(d)[length(d)] <- fc_counts_meta$pvid[i]
  if (i == 1) {
    fc_counts <- d
  } else {
    d <- d[, c(1,7)]
    fc_counts <- merge(fc_counts, d, by = "Geneid")
  }
}
fc_counts <- fc_counts[!grepl("_PAR_", fc_counts$Geneid),] # removing those gene ids containing _PAR_Y

# Save to disk
readr::write_tsv(fc_counts_meta, file = file.path(OUTPUT, paste0("ppmi_rnaseq_fc_meta.tsv")))
readr::write_tsv(as.data.frame(fc_counts), file = file.path(OUTPUT, paste0("ppmi_rnaseq_fc_counts.tsv")))

message("featureCounts data parsed successfully")


# Session info ------------------------------------------------------------
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

