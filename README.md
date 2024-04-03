# Capstone2024ccRCC

## Setup 
First, clone this repository on your local machine. The path to this repository will be used as the "base_dir" variable in each script. This will be the only variable you need to modify.

### Data Procurement
The analyses in this project rely on data from The Cancer Genome Atlas (TCGA). You will need to download data from the following links, unzip them, and save them in the "Input" folder. DO NOT change any of the original file names. If you download multiple versions of the same file, ensure that the file name does *not* include a trailing copy number (ie. "(1)"). 

1. TCGA & Clinical Data - [TCGA_&_Clinical_Data.tsv]
2. miRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.mirna.tsv.gz
3. mRNA Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz
4. DNA 450k methylation Data - https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.methylation450.tsv.gz


## miRNA data scripts
Download the TCGA miRNA data found at this link: https://xenabrowser.net/datapages/?dataset=TCGA-KIRC.mirna.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443. Click the first link after the "download" section of the page. Unzip this data. Put this data, together with the TCGA_&_Clinical_Data.tsv file, into a folder that is set as the working directory for R studio. Then run the DESeq2 script found in the miRNA file, making sure the file names match yours. To confirm the DESeq2 results, run the limma script.

## Figure 1 (mRNA)
1. First, download the *Figure1* folder.
2. Download the TCGA mRNA data found at this link: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-KIRC.htseq_counts.tsv.gz, unzip and place in *Figure1/xena/* subfolder
3. Open *Figure1/Figure1_preprocessing.ipynb* in your preferred IDE, and set the base_dir variable as the path to *Figure1* folder on your computer. Don't forget "/" at the end of base_dir! Run the script.
4. Open *Figure1/Figure1.RMD* in your preferred IDE, and set the base_dir variable as the path to *Figure1* folder on your computer. Run the script.
5. All Figure 1 plots will now be visible in the *Figure1/Plots* subfolder, and driving genes are visible at *Figure1/Data/leading_genes.csv*. 
