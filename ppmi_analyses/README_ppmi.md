Scripts on data processing and analysis of transcriptomics data from the PPMI cohort. 

### Parsing data

##### ppmi_parse_clinical.R
Parses phenotype and clinical data from the original clinical dataset from PPMI. Outputs ppmi_raw_pheno.tsv, which contains the whole original dataset, and ppmi_pheno.tsv, which contains a reduced set of relevant clinical features for PD.

##### ppmi_parse_rnaseq.R
Parses both STAR+FeatureCounts and Salmon transcript counts into format that can be readily used in R for differentially expressed gene discovery. Also does a summarization from transcript level to gene level (which is what was actually used for analysis). Outputs a file with transcript to gene correspondance, several files with quantification info (raw counts, abundance, length, length-scaled TPM) at transcript and gene level (salmon) and gene level quantification for STAR+FC.

##### ppmitx2tx.sh
Parse transcript names from PPMI data to be compatible with GENECODE transcript IDs.

##### ppmi_deseq_salmon_star.R
Process data and compare the results of performing DESEQ analysis for STAR+featureCounts and Salmon quantification pipelines on all samples (PD/HC with Age & Gender as covariates). 

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


### Longitudinal analysis PD/HC

##### ppmi_deseq_TS.R (-l gene)
Performs DESEQ analysis on consecutive timepoints (t1 vs t2; t2 vs t3; t3 vs t4) with design = ~ Age + Gender + visit (visit number with Age & Gender as covariates) separately for binary classes. Outputs DESEQ results for each class. 

##### ppmi_deseq_TS_consistent_trend.R
Filters genes with consistent sign of foldchange across all the consecutive DEAs for each class separately, and subsets identified class 1-specific (e.g. PD-specific) genes with consistent sign of foldchange across all the consecutive DEAs. Outputs results of DEAs of those genes with consistent trend across all DEAs on consecutive timepoints that are specific of class 1 (e.g. PD).

##### ppmi_deseq_TS_consistent_enrichment.R
Performs enrichment analysis on KEGG, GO and MESH databases of genes.

##### ppmi_corr_TS_deseq.R
Fits a (deseq) linear model to test for association between expression & time per diagnosis class, filters by significance. Test on median of counts, and subsets identified class 1-specific (e.g. PD-specific) genes with significant association with time.

##### ppmi_boxplots_BL.R
Generates some visuals (boxplots) from file of genes from DEAs anaysis at baseline time.

##### ppmi_boxplots_TS.R
Generates some visuals (boxplots, trajectories) from file of genes from consecutive DEAs anaysis or from association with time.

##### ppmi_deseqTC.R
Performs deseq2 time course analysis (tests for differential trajectories, whichever the trajectory is) from raw counts and generates some visuals (plotCounts, heatmap of log2FC).

##### ppmi_DEGs_perTS.R
Computes DESEQ DEAs PD-HC for each timepoint at gene level, then filters those common to all timepoints.

##### ppmi_deseq_TS_pathway_level.R
Performs DESEQ analysis on expression at aggregated levels (CORUM, GOCC, GOBP) from consecutive timepoints (t1 vs t2; t2 vs t3; t3 vs t4) with design = ~ Age + Gender + visit (visit number with Age & Gender as covariates) separately for binary classes. Outputs DESEQ results for each class.

##### ppmi_deseq_TS_consistent_trend_pathway_level.R
Filters pathways/protein complexes/cellular locations whose aggregated expression has consistent sign of foldchange across all the Bayesian Moderated t-statistics on consecutive timepoints. Outputs results of the Bayesian Moderated t-statistics of those aggregations with consistent trend across all Bayesian Moderated t-statistics on consecutive timepoints.

##### ppmi_deseqTC_pathway_level.R
Performs time course deseq2-normalization DEA (time + age + gender as covariates) on expression at aggregated levels (CORUM, GOCC, GOBP) 

##### ppmi_DEGs_perTS_deseq_pathway_level.R
Computes deseq2-normalization - DEAs PD-HC for each timepoint at aggregated levels, then filters those common to all timepoints

