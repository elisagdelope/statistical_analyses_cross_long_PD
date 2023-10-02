## ppmi_parse_clinical.R
Parses phenotype and clinical data from the original clinical dataset from PPMI. Outputs ppmi_raw_pheno.tsv, which contains the whole original dataset, and ppmi_pheno.tsv, which contains a reduced set of relevant clinical features for PD.


## ppmi_parse_rnaseq.R
Parses both STAR+FeatureCounts and Salmon transcript counts into format that can be readily used in R for differentially expressed gene discovery. Also does a summarization from transcript level to gene level (which is what was actually used for analysis). Outputs a file with transcript to gene correspondance, several files with quantification info (raw counts, abundance, length, length-scaled TPM) at transcript and gene level (salmon) and gene level quantification for STAR+FC.
