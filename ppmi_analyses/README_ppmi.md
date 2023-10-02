Scripts on data processing and analysis of transcriptomics data from the PPMI cohort. 

### Parsing data

##### ppmi_parse_clinical.R
Parses phenotype and clinical data from the original clinical dataset from PPMI. Outputs ppmi_raw_pheno.tsv, which contains the whole original dataset, and ppmi_pheno.tsv, which contains a reduced set of relevant clinical features for PD.


##### ppmi_parse_rnaseq.R
Parses both STAR+FeatureCounts and Salmon transcript counts into format that can be readily used in R for differentially expressed gene discovery. Also does a summarization from transcript level to gene level (which is what was actually used for analysis). Outputs a file with transcript to gene correspondance, several files with quantification info (raw counts, abundance, length, length-scaled TPM) at transcript and gene level (salmon) and gene level quantification for STAR+FC.


### DEA Baseline PD/HC

##### ppmi_filter_gene_expression.R
Filters low reads from raw gene expression counts. 
OPTIONS -v "BL" for timepoint-specific analysis

##### ppmi_norm_gene_expression.R
Normalizes filtered gene expression counts. 
OPTIONS -v "BL" for timepoint-specific analysis

##### ppmi_deseq_visit.R 
Performs DESEQ analysis on specific timepoint data for diagnosis with Age & Gender as covariates.
OPTIONS -v "BL" for timepoint-specific analysis

##### ppmi_generate_pathway_level.R 
Generates expression data at aggregated level (mean, median, sd, pca) from normalized counts across gene members of cellular pathways: GO BP, GO CC, CORUM + Pathway deregulation scores (pathifier). Outputs aggregated expression for all databases and all aggregated stats.

##### ppmi_deseq_visit_pathway_level.R 
Performs differential analysis for higher-level functional features on specific timepoint data for diagnosis with Age & Gender as covariates.
OPTIONS -v "BL" for timepoint-specific analysis

