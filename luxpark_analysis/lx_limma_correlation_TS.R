# Title: lx_limma_correlation_TS.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script performs linear model of time and metabolomics abundance at metabolite and pathway level via limma.
# Usage: Rscript lx_limma_correlation_TS.R -l pw -s mean  
# Data: data from metabolites abundance, clinical data.

# Packages ---------------------------------------------------------------------
library(tidyr)
library(dplyr)
library(stringr)
library(vroom)
library(tibble)
library(argparser)
library(edgeR)
library(limma)
library(statmod)
library(caret)


# I/O --------------------------------------------------------------------------
analysis_name <- "01-dea-TS-PD"
IN_DIR <- "../data/00-cleansing"
OUT_DIR <- paste0("../data/", analysis_name , "/02-outfiles") 
OUT_DIR_PATHWAY <- paste0("../data/", analysis_name , "/04-pathway_level")
ANNOTATION.FILE <- file.path(IN_DIR, "chemical_annotation.tsv")
PHENO.FILE <- file.path(IN_DIR, "lx_pheno.tsv")
METAB.FILE <- file.path(IN_DIR, "log_transformed_data_fte.tsv")
M1342.FILE <- file.path(IN_DIR, "M1342.tsv") # M1342 aka 3-methoxytyrosine: confounder
FILE.OUT <- file.path(OUT_DIR, paste0("lt_correlated_TS.tsv"))
var_id = "METABOLITES_ID"


# Add command line arguments
p <- arg_parser("Performs linear model of time vs metabolomics abundance at metabolite and pathway level via limma", hide.opts = FALSE)
p <- add_argument(parser = p, arg = "--level", help = "level of features (metab / pw)", default = "metab", type = "string", nargs = '+')
p <- add_argument(parser = p, arg = "--stat", help = "aggregated statistic (mean / median / sd / pathifier / pca)", default = "mean", type = "string", nargs = '+')
argv <- parse_args(p, argv = commandArgs(trailingOnly = TRUE))
e_level = toupper(argv$level) # gene (g) / aggregations
st = tolower(argv$stat) # stat at aggregation level

if (e_level == "PW") {
  METAB.FILE <- file.path(OUT_DIR_PATHWAY, paste0("log_transformed_PW_", st, ".tsv"))
  FILE.OUT <- file.path(OUT_DIR_PATHWAY, paste0("lt_PW_", st, "_correlated_TS.tsv"))
  var_id = "PATHWAY_NAME"
}


# Data load --------------------------------------------------------------------
annotation <- vroom(ANNOTATION.FILE, col_types = cols())
pheno <- vroom(PHENO.FILE, col_types = cols()) 
metab <- vroom(METAB.FILE, col_types = cols()) 
metab = metab[!duplicated(metab[["SAMPLE_ID"]]), ]
M1342 <- vroom(M1342.FILE, col_types = cols()) 



# Processing ---------------------------------------------------------------

# kNN imputation for BMI variable
pheno <- VIM::kNN(pheno, variable = "BMI", k= 5, imp_var = F)

# remove patients that DON'T have data for at least 2 visits 
n_tp = 2
patients_list <- metab %>%
  group_by(PATIENT_ID) %>% 
  filter(n() >= n_tp) %>%   
  ungroup() %>% 
  dplyr::select(PATIENT_ID) %>%
  distinct() %>%
  pull()

metab <- metab %>% 
  filter(PATIENT_ID %in% patients_list)
pheno <- pheno %>% 
  filter(PATIENT_ID %in% patients_list)

# timepoint as a visit number
#patient_timepoints <- split(pheno$VISIT, pheno$PATIENT_ID)
#timepoints_list <- list()
#for (n in names(patient_timepoints)) {
#  timepoints_list[[n]] <- seq_along(patient_timepoints[[n]]) - 1
#}
#pheno$TIMEPOINT <- unlist(timepoints_list)

# timepoint as time elapsed (years)
pheno <- pheno %>%
  mutate(TIMEPOINT = as.integer(case_when(
    as.character(VISIT) == "V0" ~ 0,
    as.character(VISIT) == "V1" ~ 1,
    as.character(VISIT) == "V2" ~ 2,
    as.character(VISIT) == "V3" ~ 3,
    as.character(VISIT) == "V4" ~ 4)))
RES <- list()
class = "PD"  # in luxpark there is only PD samples for longitudinal data
pheno_class <- pheno %>%
  filter(DIAGNOSIS == class) # remove HC (existing at V0)
metab_class <- metab %>%
  dplyr::filter(SAMPLE_ID %in% pheno_class$SAMPLE_ID)



# Pre-processing ---------------------------------------------------------------
# remove near zero variance features 
nzv = nearZeroVar(metab_class, names = TRUE)
if (length(nzv) > 0) {
  metab_class <- metab_class %>%
    dplyr::select(-any_of(nzv))
}

# Add confounding factor 3-methoxytyrosine (M1342) to pheno df
pheno_class <- pheno_class %>%
  left_join(M1342, 
            by = "SAMPLE_ID")

# generate patient_ids in numeric format (1:nbpatients)
patient_ids <- metab_class %>%
  group_by(PATIENT_ID) %>%
  mutate(SUBJECT_ID = cur_group_id()) %>% # individual ID as int not factor (otherwise DESeq2 doesn't like it)
  pull(SUBJECT_ID)

# flip the dataset as rows should correspond to metabolites and columns to samples
metab_class <- metab_class %>% 
  dplyr::select(-any_of(c("PATIENT_ID", "VISIT"))) %>% 
  pivot_longer(!SAMPLE_ID, names_to = var_id, values_to = "COUNT") %>% 
  pivot_wider(names_from = "SAMPLE_ID", values_from = "COUNT") %>%
  column_to_rownames(var_id)


# Apply Bayesian Moderated t-statistic -----------------------------------------
# empirical bayes t-test: visit + confounding factors

design <- model.matrix(~ pheno_class$AGE 
                       + pheno_class$GENDER 
                       + pheno_class$BMI 
                       + pheno_class$M1342
                       + pheno_class$TIMEPOINT
                       )  
colnames(design) <- c("intercept", "AGE", "GENDERM", "BMI", "M1342", "TIMEPOINT")
# accounting for within-patient correlations via duplicateCorrelations (since there is not the same number of samples per timepoint here)
corfit <- limma::duplicateCorrelation(metab_class, design, block = patient_ids)
xfit <- limma::lmFit(metab_class, design, block = patient_ids, correlation = corfit$consensus)
ebayes <- limma::eBayes(xfit)
res <- topTable(ebayes, coef = "TIMEPOINT", number = nrow(ebayes), adjust = "fdr", sort.by = "P") # coef indicates the column number of interest (account for intercept)
# ###
# head(res)
# plot(y=metab_class["M100015723",], x = pheno_class$TIMEPOINT)
# ###

RES[[class]] <- as.data.frame(res) %>% 
  rownames_to_column(var = var_id) %>%
  filter(!is.na(adj.P.Val)) %>%
  arrange(P.Value)




# output results ---------------------------------------------------------------
readr::write_tsv(RES[[class]], file = FILE.OUT)


