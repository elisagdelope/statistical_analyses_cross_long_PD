# Title: lx_boxplots_TS.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates plots of distribution of temporal metabolomics data.
# Usage: Rscript lx_boxplots_TS.R -l pw -s sd -o "correlated"
# Data: data from expression at different levels, pheno and features to plot.

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
library(argparser)



# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
IN_DIR <- "../data/00-cleansing/"
PHENO.FILE <- file.path(IN_DIR, "lx_pheno.tsv")

target = "DIAGNOSIS"

# Add command line arguments
p <- arg_parser("Nested CV", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (metab / gobp / gocc / corum)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd / pathifier / pca)", default = "mean", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--output", help = "pattern name of output file", default = "consistent_trend", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level
TYPE.FILE = tolower(argv$output)

if (e_level == "METAB") {
  METAB.FILE <- file.path(IN_DIR, "log_transformed_data_fte.tsv")
  FEAT.FILE <- file.path(OUT_DIR, paste0("lt_", TYPE.FILE, "_TS.tsv"))
  NAME.PLOT1 <- paste("lt_TS", TYPE.FILE, "mplots.pdf", sep = "_") ########P <-> C
  NAME.PLOT2 <- paste("lt_TS", TYPE.FILE, "iboxplots.pdf", sep = "_") 
  NAME.PLOT3 <- paste("lt_TS", TYPE.FILE, "lineplots.pdf", sep = "_") 
  NAME.PLOT4 <- paste("lt_TS", TYPE.FILE, "iboxplotsdots.pdf", sep = "_") 
  var_id = "METABOLITES_ID"
} else { # aggregations
  METAB.FILE <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv"))
  FEAT.FILE <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_", TYPE.FILE, "_TS.tsv"))
  NAME.PLOT1 <- paste("lt", e_level, st, TYPE.FILE, "TS_mplots.pdf", sep = "_") 
  NAME.PLOT2 <- paste("lt", e_level, st, TYPE.FILE, "TS_iboxplots.pdf", sep = "_")
  NAME.PLOT3 <- paste("lt", e_level, st, TYPE.FILE, "TS_lineplots.pdf", sep = "_")
  NAME.PLOT4 <- paste("lt", e_level, st, TYPE.FILE, "TS_iboxplotsdots.pdf", sep = "_")
  var_id = "PATHWAY_NAME"
  if ((!e_level %in% c("PW", "METAB")) | (!st %in% c("mean", "median", "sd", "pathifier", "pca")) | (!file.exists(METAB.FILE)) | (!file.exists(FEAT.FILE))) { 
    stop("Adequate arguments were not provided. Check R ppmi_rnaseq_binaryclass.R --help for the right usage.")
  }
}



# Main -------------------------------------------------------------------------
if (!dir.exists(OUT_DIR_PLOTS)) {
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------
metab <- vroom(METAB.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = cols())
features_dea <- vroom(FEAT.FILE, col_types = cols())



# Data filtering and transformations -------------------------------------------
if (nrow(features_dea) > 12) {
  features_dea_top <- features_dea[[var_id]][1:12]
} else {
  features_dea_top <- features_dea[[var_id]]
}

#c_genes_sub <- c_genes[[var_id]][1:20]

metab_pheno <- metab %>%
  pivot_longer(cols = -any_of(c("PATIENT_ID", "VISIT", "SAMPLE_ID")), names_to = var_id, values_to = "COUNTS") %>%
  left_join(pheno, by = c("SAMPLE_ID", "PATIENT_ID", "VISIT")) %>%
  filter(VISIT != "V4") # few data in V4, wasn't used for the TS analysis
metab_pheno$DIAGNOSIS <- factor(metab_pheno$DIAGNOSIS , levels=c("HC", "PD")) # set specific order for box plots categories

metab_traj <- metab_pheno %>%
  group_by_at(c(var_id, "VISIT", "DIAGNOSIS")) %>% 
  dplyr::summarise(MEDIAN_COUNTS = median(COUNTS)) %>% 
  mutate(TIMEPOINT = as.integer(case_when(
    as.character(VISIT) == "V0" ~ 0,
    as.character(VISIT) == "V1" ~ 1,
    as.character(VISIT) == "V2" ~ 2,
    as.character(VISIT) == "V3" ~ 3))) %>%
  drop_na() # remove metab that have NAs at one or more timepoints

metab_traj_mean <- metab_pheno %>% 
  group_by_at(c(var_id, "VISIT", "DIAGNOSIS")) %>% 
  dplyr::summarise(MEAN_COUNTS = mean(COUNTS)) %>% 
  mutate(TIMEPOINT = as.integer(case_when(
    as.character(VISIT) == "V0" ~ 0,
    as.character(VISIT) == "V1" ~ 1,
    as.character(VISIT) == "V2" ~ 2,
    as.character(VISIT) == "V3" ~ 3))) %>%
  drop_na() # remove metab that have NAs at one or more timepoints

metab_box <- metab_pheno %>% 
  mutate(TIMEPOINT = as.integer(case_when(
    as.character(VISIT) == "V0" ~ 0,
    as.character(VISIT) == "V1" ~ 1,
    as.character(VISIT) == "V2" ~ 2,
    as.character(VISIT) == "V3" ~ 3))) 


# PLOT DISTRIBUTIONS -----------------------------------------------------------

pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT1), width = 14, height = 8)

# multiple line plot
metab_traj %>%
  filter(get(var_id) %in% features_dea_top) %>% 
  ggplot(aes(x = TIMEPOINT,
             y = MEDIAN_COUNTS,
             fill = DIAGNOSIS,
             group = interaction(var_id, DIAGNOSIS))) +
  geom_line(aes(color = DIAGNOSIS), size = 0.5) + #, linetype=geneid
  labs(title = "Trajectories of metabolites abundance profiles (median) - PD patients",
       x = "Time points", 
       y = "Log-transformed abundance")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id))) + 
  scale_x_continuous(labels= unique(as.character(metab_pheno$VISIT)))

## time series plot(for each gene: median exprs per visit)
#palette <- brewer.pal(6, "Set1") 
#metab_traj2 <- metab_pheno %>% 
#  filter(get(var_id) %in% features_dea_top) %>% 
#  group_by_at(c(var_id, "VISIT", "DIAGNOSIS")) %>% 
#  dplyr::summarise(MEDIAN_COUNTS = median(COUNTS)) %>%
#  pivot_wider(names_from=VISIT, values_from=MEDIAN_COUNTS) %>%
#  unite("id_diag", !!var_id:DIAGNOSIS) %>%
#  column_to_rownames(var = "id_diag") %>% 
#  drop_na() # remove metab that have NAs at one or more timepoints
#
#ispd <- rep("blue", nrow(metab_traj2))
#parcoord(metab_traj2, 
#         col=ispd,
#         main = "Trajectories of 20 metabolite profiles (median) - PD patients")


# multiple BOX plot
metab_box %>%
  filter(get(var_id) %in% features_dea_top) %>% 
  ggplot(aes(x=VISIT, 
             y=COUNTS,
             fill=DIAGNOSIS)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.8) +
  geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = DIAGNOSIS), size = 0.5) + # , aes(color = diagnosis)
  labs(title = "Distribution of feature abundance by diagnosis and timepoints",
       x = "Time points", 
       y = "log-transformed abundance")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))

metab_box %>%
  filter(get(var_id) %in% features_dea_top) %>% 
  ggplot(aes(x=VISIT, 
             y=COUNTS,
             fill=DIAGNOSIS)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
  labs(title = "Distribution of feature abundance by diagnosis and timepoints",
       x = "Time points", 
       y = "log-transformed abundance")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))

dev.off()


# plot patients counts for each gene level as points + boxplots colored based on diagnosis 
pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT2), width = 14, height = 8)
for (feature in features_dea[[var_id]]) {
  
  p <- metab_box %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x=VISIT, 
               y=COUNTS,
               fill=DIAGNOSIS)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    labs(title = feature,
         x = "Time points", 
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18)
  
  # remove outliers and zoom y axis
  sts <- boxplot.stats(metab_box[metab_box[[var_id]] == feature,]$COUNTS)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p)
  
}
dev.off()

# trajectory lines (median, mean)
pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT3), width = 14, height = 8)
for (feature in features_dea[[var_id]]) {
  
  p <- metab_traj %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x = TIMEPOINT,
               y = MEDIAN_COUNTS,
               fill = DIAGNOSIS,
               group = interaction(var_id, DIAGNOSIS))) +
    geom_line(aes(color = DIAGNOSIS), size = 0.5) + #, linetype=geneid
    labs(title = paste("Trajectory (median) -", feature),
         x = "Time points", 
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18) + 
    scale_x_continuous(labels= unique(as.character(metab_pheno$VISIT)))
  print(p)
  
  pp <- metab_traj_mean %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x = TIMEPOINT,
               y = MEAN_COUNTS,
               fill = DIAGNOSIS,
               group = interaction(var_id, DIAGNOSIS))) +
    geom_line(aes(color = DIAGNOSIS), size = 0.5) + #, linetype=geneid
    labs(title = paste("Trajectory (mean) -", feature),
         x = "Time points", 
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18) + 
    scale_x_continuous(labels= unique(as.character(metab_pheno$VISIT)))
  print(pp)
}

dev.off()




pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT4), width = 14, height = 8)
for (feature in features_dea[[var_id]]) {
  p <- metab_box %>%
    filter(get(var_id) == feature) %>% 
    ggplot(aes(x=VISIT, 
               y=COUNTS,
               fill=DIAGNOSIS)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = DIAGNOSIS), size = 0.5) + # , aes(color = diagnosis)
    labs(title = feature,
         x = "Time points", 
         y = "log-transformed abundance")  + 
    theme_minimal(base_size = 18)
  
  # remove outliers and zoom y axis
  sts <- boxplot.stats(metab_box[metab_box[[var_id]] == feature,]$COUNTS)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p)
  
}
dev.off()



# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()
