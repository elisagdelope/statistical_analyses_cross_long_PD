# Title: ppmi_boxplots_TS.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates plots of distribution of expression in classes.
# Usage: R ppmi_boxplots_TS.R diagnosis -l gobp -s sd -i "../data/01-dea-TS-PD/02-outfiles/genes_consistent_trend_TS.tsv" -o "consistent_trend"
# Data: data from expression at different levels, pheno and class-specific features with consistent trend.

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
analysis_name <- "01-dea-TS-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- "../data/01-dea-PDHC/02-outfiles/ppmi_pheno_dsq.tsv"

# Add command line arguments
p <- arg_parser("Plots expression of class 1-specific features", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "target", help = "target variable (diagnosis, updrs, ...)", default = "diagnosis", type = "string", short = "t")
p <- add_argument(parser = p, arg = "--input", help = "full name of input file", default = "../data/01-dea-TS-PD/02-outfiles/genes_consistent_trend_TS.tsv", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--output", help = "pattern name of output file", default = "consistent_trend", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
target = toupper(argv$target) # target variable
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level
GENES.FILE = argv$input # input file with genes of interest
OUT.FILE = argv$output # string to tag output filename


if ((e_level == "GENE") | (e_level == "G")) { # gene level alone / lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR, "flt_norm_star_all_TS.tsv") 
#  CONSISTENT.FILE <- file.path(OUT_DIR, "genes_consistent_trend_TS.tsv")
  NAME.PLOT1 <- paste(e_level, OUT.FILE, "DEAs_TS_plots.pdf", sep = "_") ########P <-> C
  NAME.PLOT2 <- paste(e_level, OUT.FILE, "DEAs_TS_boxplots.pdf", sep = "_") 
  NAME.PLOT3 <- paste(e_level, OUT.FILE, "DEAs_TS_lineplots.pdf", sep = "_") 
  NAME.PLOT4 <- paste(e_level, OUT.FILE, "DEAs_TS_boxplotsdots.pdf", sep = "_") 

} else { # not gene level & not lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
#  CONSISTENT.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "consistent_trend_TS.tsv", sep = "_"))
  NAME.PLOT1 <- paste(e_level, st, OUT.FILE, "DEAs_TS_plots.pdf", sep = "_") 
  NAME.PLOT2 <- paste(e_level, st, OUT.FILE, "DEAs_TS_boxplots.pdf", sep = "_") 
  NAME.PLOT3 <- paste(e_level, st, OUT.FILE, "DEAs_TS_lineplots.pdf", sep = "_") 
  NAME.PLOT4 <- paste(e_level, st, OUT.FILE, "DEAs_TS_boxplotsdots.pdf", sep = "_") 
    
  if ((!e_level %in% c("GOBP", "GOCC", "CORUM")) | (!st %in% c("mean", "median", "sd", "pathifier", paste(c("mean", "median", "sd", "pathifier"), "lasso", sep = "_"))) | (!file.exists(EXPRESSION.FILE))) { 
    stop("Adequate arguments were not provided. Check R test_boxplots.R --help for the right usage.")
  }
}



# Main -------------------------------------------------------------------------
if ((!dir.exists(OUT_DIR)) || (!dir.exists(OUT_DIR_PLOTS))) {
  dir.create(OUT_DIR, recursive = T)
  dir.create(OUT_DIR_PLOTS, recursive = T)
}



# Data load --------------------------------------------------------------------
expression <- vroom(EXPRESSION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = c("ccffdfildddddddddl"))
c_genes <- vroom(GENES.FILE, col_types = cols())


# Data filtering and transformations -------------------------------------------

if (e_level == "GENE") { 
  var_id = "geneid"
} else { 
  var_id = paste0(e_level, "_name")
}

expression <- expression %>% 
                filter(get(var_id) %in% c_genes[[var_id]])
pheno <- pheno %>% 
  rename_with(toupper)

c_genes_sub <- c_genes[[var_id]][1:20]
exprs_pheno <- expression %>%
  pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "counts") %>%
  left_join(pheno, by = c("sample" = "PATIENT_VISIT_ID")) %>%
  rename_with(tolower)
var_id = tolower(var_id)
  
genes_traj <- exprs_pheno %>% 
  group_by_at(c(var_id, "visit", "diagnosis")) %>% 
  summarise(median_counts = median(counts)) %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) %>%
  drop_na() # remove genes that have NAs at one or more timepoints

genes_traj_mean <- exprs_pheno %>% 
  group_by_at(c(var_id, "visit", "diagnosis")) %>% 
  summarise(mean_counts = mean(counts)) %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) %>%
  drop_na() # remove genes that have NAs at one or more timepoints

genes_box <- exprs_pheno %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) 

c_genes <- c_genes %>%
  rename_with(tolower)


# PLOT DISTRIBUTIONS -----------------------------------------------------------

pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT1), width = 14, height = 8)

# multiple line plot
# aggregate patients counts for each gene level (mean and median) + draw 2 trajectory lines colored based on diagnosis
genes_traj %>%
  filter(get(var_id) %in% c_genes_sub) %>% 
  ggplot(aes(x = timepoint,
             y = median_counts,
             fill = diagnosis,
             group = interaction(var_id, diagnosis))) +
  geom_line(aes(color = diagnosis), size = 0.5) + #, linetype=geneid
  labs(title = "Trajectories of gene expression profiles (median) by diagnosis",
       x = "Time points", 
       y = "Gene counts")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id))) + 
  scale_x_continuous(labels= unique(as.character(exprs_pheno$visit)))

# time series plot(for each gene: median exprs per visit)
palette <- brewer.pal(6, "Set1") 
genes_traj2 <- exprs_pheno %>% 
  filter(get(var_id) %in% c_genes_sub) %>% 
  group_by_at(c(var_id, "visit", "diagnosis")) %>% 
  summarise(median_counts = median(counts)) %>%
  pivot_wider(names_from=visit, values_from=median_counts) %>%
  unite("id_diag", !!var_id:diagnosis) %>%
  column_to_rownames(var = "id_diag") %>% 
  drop_na() # remove genes that have NAs at one or more timepoints

ispd <- rep(c("red", "blue"), nrow(genes_traj2))
parcoord(genes_traj2, 
         col=ispd,
         main = "Trajectories of 20 gene expression profiles (median) by diagnosis")


# multiple BOX plot
# plot patients counts for each gene level as points + boxplots colored based on diagnosis

genes_box %>%
  filter(get(var_id) %in% c_genes_sub[1:12]) %>% 
  ggplot(aes(x=visit, 
             y=counts,
             fill=diagnosis)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.8) +
  geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = diagnosis), size = 0.5) + # , aes(color = diagnosis)
  labs(title = "Distribution of gene counts by diagnosis and timepoints",
       x = "Time points", 
       y = "Gene counts")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))

genes_box %>%
  filter(get(var_id) %in% c_genes_sub[1:12]) %>% 
  ggplot(aes(x=visit, 
             y=counts,
             fill=diagnosis)) +
  geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
  labs(title = "Distribution of gene counts by diagnosis and timepoints",
       x = "Time points", 
       y = "Gene counts")  + 
  theme_minimal() + 
  facet_wrap(as.formula(paste("~", var_id)))

dev.off()

pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT2), width = 14, height = 8)
for (gene in c_genes[[var_id]]) {
  
  p <- genes_box %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x=visit, 
               y=counts,
               fill=diagnosis)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    labs(title = paste("Distribution of gene counts by diagnosis and timepoints -", gene),
         x = "Time points", 
         y = "Gene counts")  + 
    theme_minimal() 

  # remove outliers and zoom y axis
  sts <- boxplot.stats(genes_box[genes_box[[var_id]] == gene,]$counts)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
#  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p1)
  
}
dev.off()

# trajectory lines (median, mean)
pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT3), width = 14, height = 8)
for (gene in c_genes[[var_id]]) {
  
  p <- genes_traj %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x = timepoint,
               y = median_counts,
               fill = diagnosis,
               group = interaction(var_id, diagnosis))) +
    geom_line(aes(color = diagnosis), size = 0.5) + #, linetype=geneid
    labs(title = paste("Trajectories of gene expression profiles (median) by diagnosis -", gene),
         x = "Time points", 
         y = "Gene counts")  + 
    theme_minimal() + 
    scale_x_continuous(labels= unique(as.character(exprs_pheno$visit)))
  print(p)
  
  pp <- genes_traj_mean %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x = timepoint,
               y = mean_counts,
               fill = diagnosis,
               group = interaction(var_id, diagnosis))) +
    geom_line(aes(color = diagnosis), size = 0.5) + #, linetype=geneid
    labs(title = paste("Trajectories of gene expression profiles (mean) by diagnosis -", gene),
         x = "Time points", 
         y = "Gene counts")  + 
    theme_minimal() + 
    scale_x_continuous(labels= unique(as.character(exprs_pheno$visit)))
  print(pp)
}

dev.off()




pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT4), width = 14, height = 8)
for (gene in c_genes[[var_id]]) {
  p <- genes_box %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x=visit, 
               y=counts,
               fill=diagnosis)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = diagnosis), size = 0.5) + # , aes(color = diagnosis)
    labs(title = paste("Distribution of gene counts by diagnosis and timepoints -", gene),
         x = "Time points", 
         y = "Gene counts")  + 
    theme_minimal() 

  # remove outliers and zoom y axis
  sts <- boxplot.stats(genes_box[genes_box[[var_id]] == gene,]$counts)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p1)
  
}
dev.off()




# Session info -----------------------------------------------------------------
rm(list = ls())
gc(T)
cat("\n================\n  SESSION INFO\n================\n")
sessionInfo()

