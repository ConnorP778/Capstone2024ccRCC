# Capstone2024ccRCC

## Setup 
First, clone this repository on your local machine. 

### Data Procurement
The analyses in this project rely on data from The Cancer Genome Atlas (TCGA). You will need to download data from the following links, unzip them, and save them in the "Input" folder. DO NOT change any of the original file names. If you download multiple versions of the same file, ensure that the file name does *not* include a trailing copy number (ie. "(1)"). 

1. TCGA & Clinical Data - [TCGA_&_Clinical_Data.tsv]
2. miRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.mirna.tsv.gz
3. mRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz
4. DNA 450k methylation Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.methylation450.tsv.gz


## RUN
Hooray! You have downloaded all the requisite data, now let's generate results. The folders for generating each figure are labeled in the order they need to be run (ie 1.Figure1 should be run first). If a figure requires multiple scripts, the scripts will be labeled in order as well (ie 1.datapreprocessing should be run first). **Each script will only have one variable that needs to be modified, called "base_dir". This will be the path to wherever you are storing this repo.** Be sure to include "/" at the end of the base_dir variable. 

NOTE: some scripts may include "TIME_INTENSE" in their name. If so, they have a command that takes >15 minutes. If you wish to bypass this command, these scripts will have this line of code annotated with "TIME INTENSE", and optional command immediately afterwards to pre-load in results. 

All output, including plots and data will be saved in the "output" folder of this repo. 

### Figure 1 (mRNA)
1. Run Figure1_preprocessing.ipynb with "base_dir"
2. Run Figure1.rmd with "base_dir"


### miRNA data scripts
Download the TCGA miRNA data found at this link: https://xenabrowser.net/datapages/?dataset=TCGA-KIRC.mirna.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443. Click the first link after the "download" section of the page. Unzip this data. Put this data, together with the TCGA_&_Clinical_Data.tsv file, into a folder that is set as the working directory for R studio. Then run the DESeq2 script found in the miRNA file, making sure the file names match yours. To confirm the DESeq2 results, run the limma script.


### Mutation data scripts
Steps for downloading the TCGA mutation data
1. Navigate to cbioportal.org
2. Select 'Kidney' in the column on the left side of the website
3. Choose 'Kidney Renal Clear Cell Carcinoma (TCGA, PanCancer Atlas)' in the center column
4. Click the 'Explore Selected Studies' button located at the bottom of the website
5. Once the page loads, locate and click the download button in the top left portion of the webiste, which is situated right next to the title
6. After downloading, unzip the folder containing the data
7. Move the unzipped folder to the 'Input' folder within the base directory of your project. The path to the destination folder will be:
   path_to_TCGA_data <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/')
