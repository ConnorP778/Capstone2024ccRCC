# script to create GSEA, Adjusted, and Mutation Pathway Tables
library(dplyr)
library(tidyr)


# FIXME: ensure Carter's files are in the Output/Data folder and with the correct names
# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'


# read in data
mut_kegg <- read.csv(file.path(base_dir, 'Output/Plots/kegg_mut_BMI_analysis.csv'))
mut_hallmark <- read.csv(file.path(base_dir, 'Output/Plots/hallmark_mut_BMI_analysis.csv'))
mut_cancer <- read.csv(file.path(base_dir, 'Output/Plots/cancer_mut_BMI_analysis.csv'))
gsea_kegg <- read.csv(file.path(base_dir, 'Output/Data/KEGG_Pathways_GSEA.csv'))
gsea_hallmark <-read.csv(file.path(base_dir, 'Output/Data/Hallmark_Pathways_GSEA.csv'))
gsea_cancer <- read.csv(file.path(base_dir, 'Output/Data/Cancer_Pathways_GSEA.csv'))


# Badi's script needs to generate these before I run this analysis line 226-230 MethylGSA
adj_kegg <- read.csv(file.path(base_dir, 'Output/Data/Kegg_pathways_adj_results.csv'))
adj_hallmark <- read.csv(file.path(base_dir, 'Output/Data/Hallmark_pathways_adj_results.csv'))
adj_cancer <- read.csv(file.path(base_dir, 'Output/Data/Cancer_pathways_adj_results.csv'))


# Define the file paths for saving the CSV of the analyses
kegg_pathway_table <- file.path(base_dir, "Output/Plots/Kegg_Analyses_Table.csv")
cancer_pathway_table <- file.path(base_dir, "Output/Plots/Cancer_Analyses_Table.csv")
hallmark_pathway_table <- file.path(base_dir, "Output/Plots/Hallmark_Analyses_Table.csv")


# Rename 'pval' column to designate analysis and data types, select 'pathway' and 'pval' columns
gsea_kegg <- gsea_kegg %>%
  rename(pval_gsea = pval) %>%
  select(pathway, pval_gsea)

gsea_hallmark <- gsea_hallmark %>%
  rename(pval_gsea = pval) %>%
  select(pathway, pval_gsea)

gsea_cancer <- gsea_cancer %>%
  rename(pval_gsea = pval) %>%
  select(pathway, pval_gsea)

mut_kegg <- mut_kegg %>%
  rename(pval_mut = pval) %>%
  select(pathway, pval_mut)

mut_hallmark <- mut_hallmark %>%
  rename(pval_mut = pval) %>%
  select(pathway, pval_mut)

mut_cancer <- mut_cancer %>%
  rename(pval_mut = pval) %>%
  select(pathway, pval_mut)

adj_kegg <- adj_kegg %>%
  rename(pval_adj = pvalue, pathway = ID) %>%
  select(pathway, pval_adj)

adj_hallmark <- adj_hallmark %>%
  rename(pval_adj = pvalue, pathway = ID) %>%
  select(pathway, pval_adj)

adj_cancer <- adj_cancer %>%
  rename(pval_adj = pvalue, pathway = ID) %>%
  select(pathway, pval_adj)


# Merge mut_kegg, gsea_kegg, adj_kegg based on the pathway column
merged_kegg <- merge(mut_kegg, gsea_kegg, by = "pathway", all = TRUE)
merged_kegg <- merge(merged_kegg, adj_kegg, by = "pathway", all = TRUE)


# Merge mut_hallmark, gsea_hallmark, adj_hallmark based on the pathway column
merged_hallmark <- merge(mut_hallmark, gsea_hallmark, by = "pathway", all = TRUE)
merged_hallmark <- merge(merged_hallmark, adj_hallmark, by = "pathway", all = TRUE)


# Merge mut_cancer, gsea_cancer, adj_cancer based on the pathway column
merged_cancer <- merge(mut_cancer, gsea_cancer, by = "pathway", all = TRUE)
merged_cancer <- merge(merged_cancer, adj_cancer, by = "pathway", all = TRUE)


# Write the merged data to a new CSV file
write.csv(merged_kegg, file = kegg_pathway_table, row.names = FALSE)
write.csv(merged_cancer, file = cancer_pathway_table, row.names = FALSE)
write.csv(merged_hallmark, file = hallmark_pathway_table, row.names = FALSE)

