# script to run MiniMax multiomics analysis
library(tidyverse)
library(pathwayMultiomics)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'


# Define the file path for saving the CSV of the analysis
kegg_minimax_pathway <- file.path(base_dir, "Output/Data/Kegg_Minimax.csv")
cancer_minimax_pathway <- file.path(base_dir, "Output/Data/Cancer_Minimax.csv")
hallmark_minimax_pathway <- file.path(base_dir, "Output/Data/Hallmark_Minimax.csv")


# read in data
cancer_df <- read.csv(file.path(base_dir, "Output/Plots/Cancer_Analyses_Table.csv"))
hallmark_df <- read.csv(file.path(base_dir, "Output/Plots/Hallmark_Analyses_Table.csv"))
kegg_df <- read.csv(file.path(base_dir, "Output/Plots/Kegg_Analyses_Table.csv"))


# Replace NA values in pval_mut or pval_gsea columns with 1
cancer_df$pval_mut <- ifelse(is.na(cancer_df$pval_mut), 1, cancer_df$pval_mut)
cancer_df$pval_gsea <- ifelse(is.na(cancer_df$pval_gsea), 1, cancer_df$pval_gsea)
cancer_df$pval_adj <- ifelse(is.na(cancer_df$pval_adj), 1, cancer_df$pval_adj)

hallmark_df$pval_mut <- ifelse(is.na(hallmark_df$pval_mut), 1, hallmark_df$pval_mut)
hallmark_df$pval_gsea <- ifelse(is.na(hallmark_df$pval_gsea), 1, hallmark_df$pval_gsea)
hallmark_df$pval_adj <- ifelse(is.na(hallmark_df$pval_adj), 1, hallmark_df$pval_adj)

kegg_df$pval_mut <- ifelse(is.na(kegg_df$pval_mut), 1, kegg_df$pval_mut)
kegg_df$pval_gsea <- ifelse(is.na(kegg_df$pval_gsea), 1, kegg_df$pval_gsea)
kegg_df$pval_adj <- ifelse(is.na(kegg_df$pval_adj), 1, kegg_df$pval_adj)


# calculate minimax statistic
cancer_minimax_df <- cancer_df %>%
  select(pathway, pval_mut, pval_gsea, pval_adj) %>%
  rename(mut = pval_mut, gsea = pval_gsea, adj = pval_adj) %>%
  MiniMax()

hallmark_minimax_df <- hallmark_df %>%
  select(pathway, pval_mut, pval_gsea, pval_adj) %>%
  rename(mut = pval_mut, gsea = pval_gsea, adj = pval_adj) %>%
  MiniMax()

kegg_minimax_df <- kegg_df %>%
  select(pathway, pval_mut, pval_gsea, pval_adj) %>%
  rename(mut = pval_mut, gsea = pval_gsea, adj = pval_adj) %>%
  MiniMax()


# adjust pathway p-values and filter to the most significant
cancer_res_df <- cancer_minimax_df %>%
  mutate(MiniMaxFDR = p.adjust(MiniMaxP, method = "fdr"))

hallmark_res_df <- hallmark_minimax_df %>%
  mutate(MiniMaxFDR = p.adjust(MiniMaxP, method = "fdr"))

kegg_res_df <- kegg_minimax_df %>%
  mutate(MiniMaxFDR = p.adjust(MiniMaxP, method = "fdr"))


# Save the data frames as CSV files
write.csv(cancer_res_df, file = cancer_minimax_pathway, row.names = FALSE)
write.csv(hallmark_res_df, file = hallmark_minimax_pathway, row.names = FALSE)
write.csv(kegg_res_df, file = kegg_minimax_pathway, row.names = FALSE)
