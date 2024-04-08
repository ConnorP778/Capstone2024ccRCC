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


## RUN
Hooray! You have downloaded all the requisite data, now let's generate results. The folders for generating each figure are labeled with numbers. **The folders should be run in this order: 1, 2, 5, 6, 4, 3.** If a figure requires multiple scripts, the scripts will be labeled in order as well (ie 1.datapreprocessing should be run first). **Each script will only have one variable that needs to be modified, called "base_dir". This will be the path to wherever you are storing this repo.** Be sure to include "/" at the end of the base_dir variable. 

NOTE: some scripts may include "TIME_INTENSE" in their name. If so, they have a command that takes >15 minutes. If you wish to bypass this command, these scripts will have this line of code annotated with "TIME INTENSE", and optional command immediately afterwards to pre-load in results. 

Outputs, including plots and data, will either be printed out or saved in the "output" folder of this repo. 

### Figure 1 (mRNA)
1. Install and open Jupyter (as step 2 will need to be run in Jupyter notebooks)
2. Run Figure1_preprocessing.ipynb with "base_dir"
3. Run Figure1.rmd with "base_dir"

## Figure 2 (Kaplan-Meier Clinical)
1. Run Figure2.Rmd with "base_dir"

### miRNA
1. Run Supplemental_Figure_miNA_DESeq2.R with "base_dir"
2. Run Supplemental_Figure_miRNA_limma.R with "base_dir". Note that this file does not have an output, just a print statement with results from the limma test.

## Methylation
1. Run MethylGSA.Rmd with "base_dir"

### Mutation
Steps for downloading the TCGA mutation data
1. Navigate to cbioportal.org
2. Select 'Kidney' in the column on the left side of the website
3. Choose 'Kidney Renal Clear Cell Carcinoma (TCGA, PanCancer Atlas)' in the center column
4. Click the 'Explore Selected Studies' button located at the bottom of the website
5. Once the page loads, locate and click the download button in the top left portion of the webiste, which is situated right next to the title
6. After downloading, unzip the folder containing the data
7. Move the unzipped folder to the 'Input' folder within the base directory of your project. The path to the destination folder will be:
   path_to_TCGA_data <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/')

### Figure 3 (iClusterBayes)
1. Run Figure3.Rmd with "base_dir"
