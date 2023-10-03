# Title: lx_boxplots.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates plots of distribution of metabolomics data at v0.
# Usage: Rscript lx_boxplots.R -a "01-dea-V0-NOVO" -l pw -s sd 
# Data: data from expression at different levels, pheno and features to plot.

# GC ---------------------------------------------------------------------------
rm(list = ls())
gc(T)


# Packages ---------------------------------------------------------------------
#library(tidyverse)
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(FactoMineR)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(viridis)
library(sparcl)
library(MASS)
library(argparser)



# I/O --------------------------------------------------------------------------
# Add command line arguments
p <- arg_parser("Boxplots", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--analysis", help = "analysis name (directory name)", default = "01-dea-V0-PD", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--level", help = "level of features (metab / gobp / gocc / corum)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd / pca / pathifier)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
analysis_name <- argv$analysis # analysis name variable
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level

OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- file.path(OUT_DIR, "pheno_V0.tsv")
target = "DIAGNOSIS"

if (e_level == "METAB") {
  METAB.FILE <- file.path(OUT_DIR, "log_transformed_V0.tsv")
  DEA.FILE <- file.path(OUT_DIR, "lt_limma_V0.tsv")
  NAME.PLOT1 <- "lt_V0_mplots.pdf" 
  NAME.PLOT2 <- "lt_V0_iplots.pdf"
  var_id = "METABOLITES_ID"
} else { # aggregations
  METAB.FILE <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, "_V0.tsv"))
  DEA.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_limma_", st, "_V0.tsv"))
  NAME.PLOT1 <- paste("lt", e_level, st, "V0_mplots.pdf", sep = "_") 
  NAME.PLOT2 <- paste("lt", e_level, st, "V0_iplots.pdf", sep = "_")
  var_id = "PATHWAY_NAME"
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(METAB.FILE)) | (!file.exists(DEA.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
  }
}



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR_PLOTS)) {
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------
metab <- vroom(METAB.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("cccffdiiiddddidii")) 
features_dea <- vroom(DEA.FILE, col_types = cols())



# multiple BOX plot ------------------------------------------------------------
if (nrow(features_dea) > 12) {
  features_dea_top <- features_dea[[var_id]][1:12]
} else {
  features_dea_top <- features_dea[[var_id]]
}

genes_box <- metab %>%
  pivot_longer(cols = -any_of(c("PATIENT_ID", "VISIT", "SAMPLE_ID")), names_to = var_id, values_to = "COUNTS") %>%
  left_join(pheno, by = c("SAMPLE_ID", "PATIENT_ID", "VISIT")) %>% 
  mutate(TIMEPOINT = as.integer(case_when(
    as.character(VISIT) == "V0" ~ 0,
    as.character(VISIT) == "V1" ~ 1,
    as.character(VISIT) == "V2" ~ 2,
    as.character(VISIT) == "V3" ~ 3))) 
genes_box$DIAGNOSIS <- factor(genes_box$DIAGNOSIS , levels=c("HC", "PD")) # set specific order for box plots categories

pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT1), width = 14, height = 8)
genes_box %>%
  filter(get(var_id) %in% features_dea_top) %>% 
  ggplot(aes(x=DIAGNOSIS, 
             y=COUNTS,
             fill=DIAGNOSIS)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.8) +
  geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = DIAGNOSIS), size = 0.5) + # , aes(color = diagnosis)
  labs(title = "Distribution of feature abundance by diagnosis and timepoints",
       x = "Time points", 
       y = "log-transformed abundance") + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))

genes_box %>%
  filter(get(var_id) %in% features_dea_top) %>% 
  ggplot(aes(x=DIAGNOSIS, 
             y=COUNTS,
             fill=DIAGNOSIS)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
  labs(title = "Distribution of feature abundance by diagnosis and timepoints",
       x = "Time points", 
       y = "log-transformed abundance")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))
dev.off()



# individual BOX plots ---------------------------------------------------------

pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT2), width = 14, height = 8)
for (feature in features_dea[[var_id]]) {
  p <- genes_box %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x=DIAGNOSIS, 
               y=COUNTS,
               fill=DIAGNOSIS)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    labs(title = feature,
         x = "Time points", # should be Diagnosis
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18)
  
  
  # remove outliers and zoom y axis
  sts <- boxplot.stats(genes_box[genes_box[[var_id]] == feature,]$COUNTS)$stats  # Compute lower and upper whisker limits
  #p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p)
  
  p <- genes_box %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x=DIAGNOSIS, 
               y=COUNTS,
               fill=DIAGNOSIS)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = DIAGNOSIS), size = 0.5) + # , aes(color = diagnosis)
    labs(title = feature,
         x = "Time points", 
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18)
  print(p)
}
dev.off()



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()
