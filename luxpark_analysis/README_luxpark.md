Scripts on data processing and analysis of transcriptomics data from the PPMI cohort. 

### Parsing data

##### lx_parse_clinical.R
Extract and parse clinical data from json original file.

##### lx_parse metabolomics.R
Extract and parse metabolomics and annotation data from xlsx file(s).

##### lx_parse_denovo.R
Parses luxpark de novo file to retrieve de novo patients ID with matching metabolomics data

##### lx_filter_treateffect.R
Addresses L-dopa associated treatment effects using a correlation filter with 3-methoxytyrosine and removing metabolites of the Tyrosine and Tryptophan metabolism, generates M1342 data (3-methoxytyrosine profiling) for further use as confounder.


### DEA baseline time (T0) PD/HC
DEA on T=0 to extract features with differential metabolomics abundance in PD/HC. *The DA is performed on log-transformed peak area data.

##### lx_extract_visit.R 
Extracts clinical and metabolomics data at a specific timepoint.
OPTIONS -v "V0" for timepoint-specific analysis

##### lx_denovo_filter.R
Extracts clinical and log-transformed metabolomics data from de novo patients at v0

##### lx_generate_pathway_level.R 
Generates expression data at aggregated level (mean, median, sd, 1st principal component "pca", pathifier deregulation scores") from log-transformed peak area data of KEGG pathways. Outputs aggregated expression for all aggregated stats.

##### lx_limma_DEA.R, lx_limma_DEA_pw.R
Perform limma Bayes moderated t-test between PD/HC for single metabolites (lx_limma_DEA.R) and pathway level (lx_limma_DEA_pw.R), using Age, gender, BMI and 3-Methoxytyrosine as confounders.

##### lx_annotate_DEA.R
Adds annotation data to input file, typically metabolite-level DA results.

##### lx_descriptive_plots_V0.R
Generates exploratory and descriptive plots of distribution of metabolomics data at v0.

##### lx_boxplots.R
Generates boxplots of PD/HC metabolomics data at v0.



### Longitudinal analysis PD/HC
Association with time and trend analyses for the study of longitudinal metabolomics samples (short series).

##### lx_limma_DEA_TS.R, lx_limma_DEA_TS_pw.R
Perform moderated Bayes t-test from limma package between consecutive timepoints (t1 vs t2; t2 vs t3; t3 vs t4) performed separately for binary classes; uses Age, gender, BMI and 3-Methoxytyrosine as confounders. Selects PD-specific features with a sustained trend for single metabolites and pathway aggregations. *The DA is performed on log-transformed peak area data.

##### lx_limma_correlation_TS.R -l pw -s {mean, median,...}
Fits a (limma) linear model to test for association between expression & time per diagnosis class, filters by significance. Test on median of counts, and subsets identified class 1-specific (e.g. PD-specific) genes with significant association with time.

##### lx_descriptive_plots_V0.R
Generates exploratory and descriptive plots of distribution of metabolomics data at all timepoints.

##### lx_boxplots_TS.R
Generates plots (boxplots, trajectories) of distribution of temporal metabolomics data.


