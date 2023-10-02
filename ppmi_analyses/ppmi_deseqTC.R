# For a control and treatment time series, design formula containing the condition factor, the time factor, and the interaction of the two. 
# using the likelihood ratio test with a reduced model which does not contain the interaction terms 
# will test whether the condition induces a change in gene expression at any time point after the reference level time point (time 0). 
# We use a design formula that models the strain difference at time 0, the difference over time, and any strain-specific differences over time (the interaction term strain:minute).

# likelihood ratio test, where we remove the strain-specific differences over time. Genes with small p values from this test are those which at one or more time points after time 0 showed a strain-specific effect. 
# Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both strains.

# Title: ppmi_deseqTC.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script deseq2 time course analysis from raw counts.
# Usage: R ppmi_deseqTC.R
# Data: data from expression at gene level (STAR) -already filtered for low expression genes- and clinical (pheno) data.

# GC ---------------------------------------------------------------------------
rm(list=ls())
gc(T)





# Packages ----------------------------------------------------------------
library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(DESeq2)
library(pheatmap)


# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 

EXPRESSION.FILE <- file.path(OUT_DIR, "flt_star_all_TS.tsv")
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
OUT.EXPRESSION.FILE <- file.path(OUT_DIR, "genes_deseqTC.tsv")



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) | (!dir.exists(OUT_DIR_PLOTS))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------

pheno_dsq <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) %>% 
  column_to_rownames(var = 'patient_visit_id') %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3)))
expression <- vroom(EXPRESSION.FILE, col_types = cols()) %>% 
  column_to_rownames(var = 'geneid') 

dds <- DESeq2::DESeqDataSetFromMatrix(countData = expression, colData = pheno_dsq, design = ~ Age + Gender + Diagnosis + timepoint + Diagnosis:timepoint)
dds <- DESeq(dds, test="LRT", reduced = ~ Age + Gender + Diagnosis + timepoint)

res <- results(dds)
head(res[order(res$padj),], 4)
genes_TC <- rownames(as.data.frame(res) %>%
  filter(padj < 0.05) %>%
  arrange(padj) )

res_TC <- as.data.frame(res) %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  rownames_to_column('geneid')


# output results ---------------------------------------------------------------
readr::write_tsv(res_TC, file = OUT.EXPRESSION.FILE)


# plot gene expression with plotCounts -----------------------------------------
# for genes with significant adjusted p value. 

NAME.PLOT = "GENE_plotCounts_deseqTC.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT), width = 14, height = 8)
for (gene in genes_TC) {
  fiss <- plotCounts(dds, gene, # which.min(resTC$padj),
                     intgroup = c("timepoint","Diagnosis"), returnData = TRUE)
  fiss$timepoint <- as.numeric(as.character(fiss$timepoint))
  p <- ggplot(fiss,
         aes(x = timepoint, y = count, color = Diagnosis, group = Diagnosis)) + 
    geom_point() + 
    stat_summary(fun=median, geom="line") +
    stat_summary(fun=mean, geom="line", linetype="dotted") +
    scale_y_log10() +
    labs(title = paste("plotCounts with median (-) and mean (···) of ", gene),
         x = "Time points", 
         y = "Gene counts") 
  
  print(p)
  print(paste0("Plotting ", which(gene==genes_TC), "/", length(genes_TC)))
}
dev.off()


# Heatmap of log2 fold changes  -----------------------------------------------------
# for genes with significant adjusted p value. 

betas <- coef(dds)
colnames(betas)
mat <- betas[genes_TC, -c(1)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

NAME.PLOT2 = "GENE_heatmapFC_deseqTC.pdf"
pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT2), width = 14, height = 20)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)
dev.off()


# Another option is to model the counts as a smooth function of time, and to include an interaction term of the condition with the smooth function. 
# It is possible to build such a model using spline basis functions within R, and another, more modern approach is using Gaussian processes (Tonner et al. 2017).

# We can plot the counts for the groups over time using ggplot2, for the gene with the smallest adjusted p value, testing for condition-dependent time profile and accounting for differences at time 0 (figure below). Keep in mind that the interaction terms are the difference between the two groups at a given time after accounting for the difference at time 0.



