# Capstone2024ccRCC

## Setup 
First, clone this repository on your local machine. 

### Data Procurement
The analyses in this project rely on data from The Cancer Genome Atlas (TCGA). You will need to download data from the following links, unzip them, and save them in the "Input" folder. DO NOT change any of the original file names. If you download multiple versions of the same file, ensure that the file name does *not* include a trailing copy number (ie. "(1)"). 

1. TCGA & Clinical Data - [TCGA_&_Clinical_Data.tsv] (on box)
2. miRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.mirna.tsv.gz
3. mRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz
4. DNA 450k methylation Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.methylation450.tsv.gz
5. Mutation Data - https://cbioportal-datahub.s3.amazonaws.com/kirc_tcga_pan_can_atlas_2018.tar.gz
6. RPPA Data - "data_rppa.txt" within the folder downloaded on step 5
7. mtDNA Data - Download "Supplementary file 1" at https://elifesciences.org/articles/10769/figures#SD4-data. You can also find the link in this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4775221/.

## RUN
Hooray! You have downloaded all the requisite data, now let's generate results. The folders for generating each figure are labeled with numbers. **The folders should be run in order they are labeled in ** If a figure requires multiple scripts, the scripts will be labeled in order as well (ie 1.datapreprocessing should be run first). **Each script will only have one variable that needs to be modified, called "base_dir". This will be the path to wherever you are storing this repo.** Be sure to include "/" at the end of the base_dir variable (Each path is currently formatted for mac/linux paths). 

NOTE: some scripts may include "TIME_INTENSE" in their name. If so, they have a command that takes >15 minutes. If you wish to bypass this command, these scripts will have this line of code annotated with "TIME INTENSE", and optional command immediately afterwards to pre-load in results. 

Outputs, including plots and data, will either be printed out or saved in the "output" folder of this repo. 

### Figure 1 (mRNA)
1. Install and open Jupyter (as step 2 will need to be run in Jupyter notebooks)
2. Run Figure1_preprocessing.ipynb with "base_dir"
3. Run Figure1.rmd with "base_dir"

### Figure 2 (Kaplan-Meier Clinical, BMI distribution)
1. Run Figure2.Rmd with "base_dir". *Note that the plots will all output into the output/plots folder except for the kaplan meier because it is composed of a graph and text and is not compatible with ggsave. When run in R-studio you can still see the plot just fine and Dr. Payne said that would be okay.
2. Make sure to install the forestplot and survminer packages if they do not already exist

### miRNA
1. Run Supplemental_Figure_miNA_DESeq2.R with "base_dir"
2. Run Supplemental_Figure_miRNA_limma.R with "base_dir". Note that this file does not have an output, just a print statement with results from the limma test.

### mtDNA
1. Run MTDNA_CopyNumberAnalysis.Rmd with "base_dir" after downloading the Supplementary file 1.

### Methylation
1. Run TIME_INTENSE_MethylGSA.Rmd with "base_dir"

### Figure 3 (iClusterBayes)
1. Run Figure3.Rmd with "base_dir"

### RPPA
1. Run "RPPA.R" with "base_dir"

### Mutation
Navigate to the scripts/4.mutations directory and run each script, updating "base_dir" each time:
1. Mutation_Methylation_Table_Creation.R	
	* Compare your output with mut_meth_table_comparison.tsv (found in Output/Data)
2. Mutation_Freq_BMI_Analysis.R
	* Compare output with mut_freq_plot_comparison.pdf (found in Output/Plots)
3. Mutation_3p_Status_BMI_Analysis.R
	* Compare output with mut_3pstatus_plot_comparison.pdf (found in Output/Plots)
4. Cancer_Kegg_Hallmark_Pathway_Analyses.R
	* Compare output with cancer_mut_BMI_analysis_comparison.csv, hallmark_mut_BMI_analysis_comparison.csv, and kegg_mut_BMI_analysis_comparison.csv (found in Output/Data)
5. GSEA_Adj_Mut_Pathway_Tables_Creation.R
	* Compare output with Cancer_Analyses_Table_Comparison.csv, Hallmark_Analyses_Table_Comparison.csv, and Kegg_Analyses_Table_Comparison.csv (found in Output/Plots)
6. Minimax_Multiomics_Analysis.R
	* Compare output with Cancer_Minimax_comparison.csv, Kegg_Minimax_Comparison.csv, and Cancer_Minimax_comparison.csv (found in Output/Data)
7. Indel_Proportion_T-test_Analysis.R
	* Compare output with Indel_	Proportion_Boxplot_comparison.pdf (found in Output/Plots)

