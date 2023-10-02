# Title: ppmi_boxplots_BL.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script generates plots of distribution of expression in classes.
# Usage: R ppmi_boxplots_BL.R diagnosis -l gobp -s sd 
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
analysis_name <- "01-dea-BL-PD"
OUT_DIR_PLOTS <- paste0("../data/", analysis_name, "/01-plots")
OUT_DIR <- paste0("../data/", analysis_name, "/02-outfiles")
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level") 
PHENO.FILE <- file.path(OUT_DIR, "ppmi_pheno_BL_dsq.tsv")

# Add command line arguments
p <- arg_parser("Plots expression of class 1-specific features", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (gene / gobp / gocc / corum)", default = "gene", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd)", default = "mean", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "target", help = "target variable (diagnosis, updrs, ...)", default = "diagnosis", type = "string", short = "t")
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
target = toupper(argv$target) # target variable
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level


if ((e_level == "GENE") | (e_level == "G")) { # gene level alone / lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR, "flt_norm_star_BL.tsv") 
  NAME.PLOT <- paste(e_level, "BL_boxplots.pdf", sep = "_") 
  DEA.FILE <- file.path(OUT_DIR, "genes_deseq_BL.tsv") 
  
} else { # not gene level & not lasso ft selection
  EXPRESSION.FILE <- file.path(OUT_DIR_PATHWAY, paste(e_level, st, "expression.tsv", sep = "_"))
  NAME.PLOT <- paste(e_level, st, "BL_boxplots.pdf", sep = "_") 
  DEA.FILE <- file.path(OUT_DIR_PATHWAY, paste("DESeq_res", e_level, st, "BL.tsv", sep = "_")) 
  
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
dea_genes <- vroom(DEA.FILE, col_types = cols())



# Data filtering and transformations -------------------------------------------

if (e_level == "GENE") { 
  var_id = "geneid"
} else { 
  var_id = paste0(e_level, "_name")
}

expression <- expression %>% 
  filter(get(var_id) %in% dea_genes[[var_id]])
pheno <- pheno %>% 
  rename_with(toupper)

dea_genes_sub <- dea_genes[[var_id]][1:20]
exprs_pheno <- expression %>%
  pivot_longer(cols = stringr::str_subset(colnames(expression), "[0-9]{4}\\.[A-Z]"), names_to = "sample", values_to = "counts") %>%
  left_join(pheno, by = c("sample" = "PATIENT_VISIT_ID")) %>%
  rename_with(tolower)
var_id = tolower(var_id)


genes_box <- exprs_pheno %>% 
  mutate(timepoint = as.integer(case_when(
    as.character(visit) == "BL" ~ 0,
    as.character(visit) == "V04" ~ 1,
    as.character(visit) == "V06" ~ 2,
    as.character(visit) == "V08" ~ 3))) 

dea_genes <- dea_genes %>%
  rename_with(tolower)



pdf(file = file.path(OUT_DIR_PLOTS, NAME.PLOT), width = 14, height = 8)
for (gene in dea_genes[[var_id]]) {
  # boxplot
  p <- genes_box %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x=diagnosis, 
               y=counts,
               fill=diagnosis)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    labs(title = gene,
         x = "Diagnosis", 
         y = "Gene counts")  + 
    theme_minimal(base_size = 18)
  
  # remove outliers and zoom y axis
  sts <- boxplot.stats(genes_box[genes_box[[var_id]] == gene,]$counts)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p1)
  
  # boxplot with dots
  p <- genes_box %>%
    filter(get(var_id) == gene) %>% 
    ggplot(aes(x=diagnosis, 
               y=counts,
               fill=diagnosis)) +
    geom_boxplot(width = 0.5, color = "grey", alpha = 0.9) +
    geom_point(position=position_jitterdodge(dodge.width = 0.65), aes(color = diagnosis), size = 0.5) + # , aes(color = diagnosis)
    labs(title = gene,
         x = "Diagnosis", 
         y = "Gene counts")  + 
    theme_minimal(base_size = 18) 
  
  # remove outliers and zoom y axis
  sts <- boxplot.stats(genes_box[genes_box[[var_id]] == gene,]$counts)$stats  # Compute lower and upper whisker limits
  p1 = p + coord_cartesian(ylim = sts[c(1,5)] * 1.05)
  #  p1 = p + coord_cartesian(ylim = c(sts[2]/2,max(sts)*1.05))
  print(p1)
  
}
dev.off()
