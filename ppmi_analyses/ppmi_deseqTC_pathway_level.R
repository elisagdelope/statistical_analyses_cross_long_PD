# Title: ppmi_deseqTC_pathway_level.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script deseq2 time course analysis from raw counts.
# Usage: R ppmi_deseqTC_pathway_level.R
# Data: data from expression at aggregated level (GOBP, GOCC, CORUM) -already normalized for sequencing depth and RNA composition- and clinical (pheno) data.

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
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 

PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"
AGG.EXPRESSION.FILE <- list()
for (db in c("GOBP", "GOCC", "CORUM")) {
  AGG.EXPRESSION.FILE[[db]] <- list()
  for (st in c("mean", "median", "sd", "pathifier")) {
    AGG.EXPRESSION.FILE[[db]][[st]] <- file.path(OUT_DIR_PATHWAY, paste(db, st, "expression.tsv", sep = "_"))
  }
}



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR_PATHWAY))) {
  dir.create(OUT_DIR_PATHWAY, recursive = T)
}




# Data load --------------------------------------------------------------------

pheno_dsq <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl")) %>% 
  column_to_rownames(var = 'patient_visit_id') %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3)))



AGG.EXPRESSION <- list()
for (db in c("GOBP", "GOCC", "CORUM")) {
  AGG.EXPRESSION[[db]] <- list()
  var_id <- paste0(db, "_name")
  
  for (st in names(AGG.EXPRESSION.FILE[[db]])) { 
    expression <- vroom(AGG.EXPRESSION.FILE[[db]][[st]] , col_types = cols(), delim = "\t") %>% 
      column_to_rownames(var = var_id) 
    
    if (nrow(expression) > 2) {
      
      if ((st == "pathifier") & (db != "GENE")) {  # for pathifier scores [0,1], they are multiplied by 100 to avoid a matrix of 1s/0s
        E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression * 100), 
                                                 colData = pheno_dsq, 
                                                 design = ~ Age + Gender + Diagnosis + timepoint + Diagnosis:timepoint)
      } else {
        E.STAR <- DESeq2::DESeqDataSetFromMatrix(round(expression), 
                                                 colData = pheno_dsq, 
                                                 design = ~ Age + Gender + Diagnosis + timepoint + Diagnosis:timepoint)
      }
      
      sizeFactors_nonnorm <- rep(1, ncol(expression)) # create a vector of 1s to avoid normalization
      names(sizeFactors_nonnorm) <- colnames(expression)
      sizeFactors(E.STAR) <- sizeFactors_nonnorm
      dds <- estimateDispersions(E.STAR)
      dds <- nbinomLRT(dds, reduced = ~ Age + Gender + Diagnosis + timepoint, maxit=500)
      
      res <- results(dds)
      head(res[order(res$padj),], 4)
      genes_TC <- rownames(as.data.frame(res) %>%
                             filter(pvalue < 0.05) %>% # CHANGED BECAUSE PADJ < 0.05 WAS TOO RESTRICTIVE!
                             arrange(padj) )
      
      res_TC <- as.data.frame(res) %>%
        filter(pvalue < 0.05) %>% # CHANGED BECAUSE PADJ < 0.05 WAS TOO RESTRICTIVE!
        arrange(padj) %>%
        rownames_to_column(var_id)
      
      
      # output results ---------------------------------------------------------
      readr::write_tsv(res_TC, file = file.path(OUT_DIR_PATHWAY, paste(db, st, "deseqTC.tsv", sep = "_")))
      print(paste("TC deseq analysis for", db, st, "completed"))
      
      
      
      
      # plot gene expression with plotCounts -----------------------------------
      # for genes with significant p value. 
      
      if (nrow(res_TC) > 0) {
        NAME.PLOT = paste(db, st, "plotCounts_deseqTC.pdf", sep = "_") 
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
        # for genes with significant p value. 
        if (nrow(res_TC) > 2) { 
          
          betas <- coef(dds)
          colnames(betas)
          mat <- betas[genes_TC, -c(1)]
          thr <- 3 
          mat[mat < -thr] <- -thr
          mat[mat > thr] <- thr
          
          NAME.PLOT2 = paste(db, st, "plotCounts_heatmapFC_deseqTC.pdf", sep = "_") 
          pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT2), width = 14, height = 20)
          pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
                   cluster_col=FALSE)
          dev.off()
        } else {
          message("Heatmap was not generated, need at least 3 features to compare.")
        }
      } else {
        message("Plots were not generated due to inexistent significant features.")
      }
        
      
      
    } else {
      message(paste(db, st, "failed: Need at least two features to fit a model", sep = " "))
      next
    }
  }
}





# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()


      
      
      