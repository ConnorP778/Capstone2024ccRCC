# Comprehensive script for Mutation Analyses


# use appropriate packages
library(dplyr)
library(tidyr)



# Download appropriate data from TCGA

# Step 1: navigate to cbioportal.org
# Step 2: check 'Kidney' in the column on the left side of the website
# step 3: check 'Kidney Renal Clear Cell Carcinoma (TCGA, PanCancer Atlas) in center column
# Step 4: click 'Explore Selected Studies' button at bottom of the website
# Step 5: click the download button in the top left portion of the website (right next to the title)
# Step 6: unzip the downloaded folder
# Step 7: enter the address to this TCGA data as the 'path_to_TCGA_data' below
path_to_TCGA_data <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Reproducibility/kirc_tcga_pan_can_atlas_2018'

mutations_file <- 'data_mutations.txt'
path_to_mutations_file <- file.path(path_to_TCGA_data, mutations_file)
# Step 8: download the clinical data FIXME: WHERE TO GET THIS FROM
# Step 9: enter the address to this clinical data as the 'path_to_clinical_data' below
path_to_clinical_data_file <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/TCGA_&_Clinical_Data.tsv'



# ANALYSIS 1: Mutation Methylation Table

#FIXME: find this path from Badi's script
mut_meth_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/kirc_mut_meth.tsv'



# ANALYSIS 2: Mutation BMI frequency Analysis

# read mutation data
data_mutation <- as.data.frame(read.delim(path_to_mutations_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(Hugo_Symbol)
selected_columns <- data_mutation[c("Hugo_Symbol", "Tumor_Sample_Barcode")]
# Remove "-01" from the end of every entry in the "Tumor_Sample_Barcode" column
data_mutation$Tumor_Sample_Barcode <- gsub("-01$", "", data_mutation$Tumor_Sample_Barcode)


# count unique entries of "Tumor_Sample_Barcode"
barcode_counts <- table(data_mutation$Tumor_Sample_Barcode)

# Convert the table to a data frame, format patient IDs appropriately
barcode_counts_df <- data.frame(Patient_ID = names(barcode_counts), Mutation_Count = as.vector(barcode_counts))
barcode_counts_df$'Patient_ID' <- sub("-01$", "", barcode_counts_df$'Patient_ID')

# read clinical data
data_clinical <- as.data.frame(read.delim(path_to_clinical_data_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(X.Patient.Identifier) %>%
  rename(PATIENT_ID = X.Patient.Identifier, SUBTYPE = Subtype)


# Filter rows with "KIRC" value in the "CANCER_TYPE_ACRYONYM" column and remove rows with NA values
kirc_data <- data_clinical[!is.na(data_clinical$SUBTYPE) & data_clinical$SUBTYPE == "KIRC", ]
clinical_filtered_columns <- kirc_data[c("PATIENT_ID", "BMI", "BMI.Classification")]

# Create a new column "Simplified_Classifications"
clinical_filtered_columns <- clinical_filtered_columns %>%
  mutate(Simplified_Classifications = case_when(
    BMI.Classification %in% c("OB I", "OB II", "OB III") ~ "OB",  # Combine OB I, OB II, and OB III into OB
    TRUE ~ BMI.Classification  # For all other cases, keep the original classification
  ))

# Merge the data frames based on the specified columns
merged_data <- merge(clinical_filtered_columns, barcode_counts_df, by.x = "PATIENT_ID", by.y = "Patient_ID", all = TRUE)

# Convert BMI classification to a factor
merged_data$Simplified_Classifications <- factor(merged_data$Simplified_Classifications, levels = c("NW", "OW", "OB"))

# Fit logistic regression model
model <- glm(Simplified_Classifications ~ Mutation_Count, data = merged_data, family = binomial)

# Summary of the model
summary_model <- summary(model)

# Extract the p-value
p_value <- summary_model$coefficients["Mutation_Count", "Pr(>|z|)"]

# Define the file path for saving the PDF of the analysis
mutation_frequency_file <- 'frequency_analysis.pdf'
path_to_mutations_file <- file.path(path_to_TCGA_data, mutation_frequency_file)

# Open a PDF device to start recording the plot
pdf(path_to_mutations_file)

# Visualize the association with custom x-axis range
plot(merged_data$Mutation_Count, merged_data$Simplified_Classifications, 
     xlab = "Mutation Frequency", ylab = "BMI Classification", col = "blue",
     xlim = c(0, 200), yaxt = "n")

# Add a red line indicating the linear trend
abline(lm(merged_data$Simplified_Classifications ~ merged_data$Mutation_Count), col = "red")

# Define the y-axis tick labels
axis(2, at = 1:length(levels(merged_data$Simplified_Classifications)), 
     labels = levels(merged_data$Simplified_Classifications))

# Create the title including the p-value
title <- paste("\nAssociation between Mutation Frequency and BMI Classification\n(p-value =", round(p_value, 4), ")")

# Add title to the plot
title(title)  # Adjust the line parameter as needed to position the title

# Close the PDF device to save the plot
dev.off()



# ANALYSIS 3: 3p Loss Analysis

# file paths
arm_level_file <- 'data_armlevel_cna.txt'
path_to_arm_level_file <- file.path(path_to_TCGA_data, arm_level_file)

# read in cna arm data
data_cna_armlevel <- as.data.frame(read.delim(path_to_arm_level_file, header = TRUE, sep = "\t", dec = ".")) %>%
  filter(NAME == '3p') %>%
  select(-ENTITY_STABLE_ID, -NAME, -DESCRIPTION) %>%
  pivot_longer(cols = everything(), names_to = "PatientID", values_to = "Status.3p")

# replace . with - and remove the -01 so patient IDs align
data_cna_armlevel$PatientID <- gsub('\\.', '-', data_cna_armlevel$PatientID)
data_cna_armlevel$PatientID <- sub("-01$", "", data_cna_armlevel$PatientID)

# Filter data_cna_armlevel to keep only values of PatientID that are in kirc_patient_ids
data_cna_armlevel_filtered <- data_cna_armlevel[data_cna_armlevel$'PatientID' %in% kirc_patient_ids, ]

# calculate proportion of 3p that are loss
loss_count <- data_cna_armlevel_filtered %>%
  filter(Status.3p == "Loss") %>%
  summarise(loss_counts = n()) %>%
  pull(loss_counts)  # Extracting the value as a vector

total_count <- data_cna_armlevel_filtered %>%
  summarise(total_entries = n()) %>% 
  pull(total_entries)  # Extracting the value as a vector

proportion_loss <- loss_count / total_count

# find IDs of the patients with 3p loss, ignoring NA values
threeP_loss_ids <- data_cna_armlevel_filtered$PatientID[!is.na(data_cna_armlevel_filtered$Status.3p) & data_cna_armlevel_filtered$Status.3p == "Loss"]

# open mutation/methylation table
mut_meth_table <- read.delim(mut_meth_path, header = TRUE, sep = "\t")

# merge the mut_meth_table with the dna_cna_armlevel_filtered table
mut_meth_3p_table <- left_join(mut_meth_table, data_cna_armlevel_filtered, by = c("TCGA.Number" = "PatientID"))

# calculate proportion of lost that are VHL mutated
loss_VHL_count <- mut_meth_3p_table %>%
  filter(Status.3p == "Loss", VHL.Mutation.Status == "Mutated") %>%
  summarise(loss_VHL_count = n()) %>%
  pull(loss_VHL_count)
proportion_loss_VHL <- loss_VHL_count / loss_count

# calculate proportion of lost that are hypomethylated
loss_hypometh_count <- mut_meth_3p_table %>%
  filter(Status.3p == "Loss", VHL.Methylation.Status == "low") %>%
  summarise(loss_hypometh_count = n()) %>%
  pull(loss_hypometh_count)
proportion_loss_hypometh <- loss_hypometh_count / loss_count

# Create the contingency table
contingency_table <- table(mut_meth_3p_table$Status.3p, mut_meth_3p_table$VHL.Mutation.Status)

print(contingency_table)
# Perform the chi-square test, ignoring NA values
chi_sq_test <- chisq.test(contingency_table)

# Print the results
print(chi_sq_test)

# Select specific columns from clinical_data, add the bmi data to the table
clinical_subset <- clinical_data[c("PATIENT_ID", "BMI", "BMI.Classification")]

# Merge the two datasets based on the common columns
merged_data <- merge(clinical_subset, mut_meth_3p_table, by.x = "PATIENT_ID", by.y = "TCGA.Number")

# Check if the variables are factors
merged_data$Status.3p <- factor(merged_data$Status.3p)

# Perform logistic regression
logit_model <- glm(Status.3p ~ BMI, data = merged_data, family = "binomial")

# Check the summary of the logistic regression model
summary(logit_model)

# Extract the p-value for the coefficient of BMI
p_value <- summary(logit_model)$coefficients["BMI", "Pr(>|z|)"]

# Define the file path for saving the PDF of the analysis
chr_3p_file <- 'chr_3p_analysis.pdf'
path_to_chr_3p_file <- file.path(path_to_TCGA_data, chr_3p_file)

# Open a PDF device to start recording the plot
pdf(path_to_chr_3p_file)

# Visualize the association with custom x-axis range
plot(merged_data$Status.3p, merged_data$BMI, 
     xlab = "Chr 3p Status", ylab = "BMI", col = "blue")

# Add a red line indicating the linear trend
abline(lm(merged_data$Status.3p ~ merged_data$BMI), col = "red")

# Create the title including the p-value
title <- paste("Association between Chr 3p status and BMI (p-value =", round(p_value, 4), ")")

# Add title to the plot
title(title, line = 2.5)  # Adjust the line parameter as needed to position the title

# Close the PDF device to save the plot
dev.off()



# ANALYSIS 4: Mutation, KEGG, Hallmark Pathway Analyses

# download file with genes in each pathway
# FIXME: these should be in the box, or in Carter's code
# file paths
# FIXME: these should align with the downloaded files that were just downloaded
cancer_pathway_genes_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/PathwayGenes.csv'
hallmark_pathway_genes_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/hallmark.csv'
kegg_pathway_genes_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/kegg.csv'

# Read pathway genes
cancer_pathway_genes <- read.csv(cancer_pathway_genes_path, fileEncoding = "UTF-8")
kegg_pathway_genes <- read.csv(kegg_pathway_genes_path, fileEncoding = "UTF-8")
hallmark_pathway_genes <- read.csv(hallmark_pathway_genes_path, fileEncoding = "UTF-8")

# Filter obese and normal weight patient IDs
ob_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification %in% c("OB I", "OB II", "OB III")]
nw_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification == "NW" & complete.cases(data_clinical$BMI.Classification)]

# Initialize list to store cancer pathway results
cancer_results_list <- list()

# Loop through each cancer pathway
for (pathway_name in names(cancer_pathway_genes)) {
  pathway_genes_current <- cancer_pathway_genes[[pathway_name]]
  
  # Filter mutation data for NW patients
  nw_mutations <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% nw_patient_ids, Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Filter mutation data for OB patients
  ob_mutations <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% ob_patient_ids, Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Calculate counts
  nw_count_mutated <- sum(nw_mutations$has_mutation == "Mutated", na.rm = TRUE)
  ob_count_mutated <- sum(ob_mutations$has_mutation == "Mutated", na.rm = TRUE)
  nw_count_unmutated <- length(nw_patient_ids) - nw_count_mutated
  ob_count_unmutated <- length(ob_patient_ids) - ob_count_mutated
  
  # Create counts matrix
  counts_matrix <- matrix(c(nw_count_mutated, ob_count_mutated, nw_count_unmutated, ob_count_unmutated), nrow = 2, byrow = TRUE,
                          dimnames = list(c("Mutated", "Unmutated"), c("Normal Weight", "Obese")))
  
  # Perform chi-squared test, p-value is 0.05
  p_value <- perform_chi_squared_test(counts_matrix)
  significant_association <- p_value < 0.05
  
  # Store results
  cancer_results_list[[pathway_name]] <- data.frame(pathway = pathway_name, Significant_Association = significant_association, pval = p_value, stringsAsFactors = FALSE)
}

# Combine results into a single dataframe
cancer_results_pathway_chisq <- bind_rows(cancer_results_list)

# Print or use results_pathway_chisq as needed
print(cancer_results_pathway_chisq)

# Define the file path for saving the CSV of the analysis
cancer_mut_file <- 'cancer_pathway_analysis.csv'
cancer_pathway_mutation_df <- file.path(path_to_TCGA_data, cancer_mut_file)

# Save the dataframe as a CSV file
write.table(cancer_results_pathway_chisq, file = cancer_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)


# Initialize list to store results
hallmark_results_list <- list()

# Loop through each KEGG pathway
for (i in 1:nrow(hallmark_pathway_genes)) {
  pathway_name <- hallmark_pathway_genes$name[i]
  pathway_genes_current <- unlist(strsplit(gsub("\"", "", pathway_genes$value_str[i]), ", "))
  
  # Check if each patient has mutations in the specified genes for NW patients
  nw_mutation_status <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% nw_patient_ids & Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Check if each patient has mutations in the specified genes for OB patients
  ob_mutation_status <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% ob_patient_ids & Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Calculate counts for NW patients
  nw_count_mutated <- nrow(subset(nw_mutation_status, !is.na(has_mutation)))
  nw_patient_ids_count <- length(nw_patient_ids)
  nw_count_unmutated <- nw_patient_ids_count - nw_count_mutated
  
  # Calculate counts for OB patients
  ob_count_mutated <- nrow(subset(ob_mutation_status, !is.na(has_mutation)))
  ob_patient_ids_count <- length(ob_patient_ids)
  ob_count_unmutated <- ob_patient_ids_count - ob_count_mutated
  
  # Create a counts matrix
  counts_matrix <- t(matrix(c(nw_count_mutated, ob_count_mutated, nw_count_unmutated, ob_count_unmutated), nrow = 2))
  colnames(counts_matrix) <- c("Normal Weight", "Obese")
  
  # Perform chi-squared test and store the result in results_df
  p_value <- perform_chi_squared_test(counts_matrix)
  significant_association <- p_value < 0.05
  
  # Store results
  hallmark_results_list[[pathway_name]] <- data.frame(pathway = pathway_name, Significant_Association = significant_association, pval = p_value, stringsAsFactors = FALSE)
}

# Combine results into a single dataframe
results_hallmark_pathways_chisq <- bind_rows(hallmark_results_list)

# Print or use results_kegg_pathways_chisq as needed
print(results_hallmark_pathways_chisq)

# Define the file path for saving the CSV of the analysis
hallmark_mut_file <- 'hallmark_pathway_analysis.csv'
hallmark_pathway_mutation_df <- file.path(path_to_TCGA_data, hallmark_mut_file)

# Save the dataframe as a CSV file
write.table(results_hallmark_pathways_chisq, file = hallmark_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)


# Initialize list to store results
kegg_results_list <- list()

# Loop through each KEGG pathway
for (i in 1:nrow(kegg_pathway_genes)) {
  pathway_name <- kegg_pathway_genes$name[i]
  pathway_genes_current <- unlist(strsplit(gsub("\"", "", pathway_genes$value_str[i]), ", "))
  
  # Check if each patient has mutations in the specified genes for NW patients
  nw_mutation_status <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% nw_patient_ids & Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Check if each patient has mutations in the specified genes for OB patients
  ob_mutation_status <- data_mutation %>%
    filter(Tumor_Sample_Barcode %in% ob_patient_ids & Hugo_Symbol %in% pathway_genes_current) %>%
    group_by(Tumor_Sample_Barcode) %>%
    summarise(has_mutation = ifelse(any(Hugo_Symbol %in% pathway_genes_current), "Mutated", "Not Mutated"))
  
  # Calculate counts for NW patients
  nw_count_mutated <- nrow(subset(nw_mutation_status, !is.na(has_mutation)))
  nw_patient_ids_count <- length(nw_patient_ids)
  nw_count_unmutated <- nw_patient_ids_count - nw_count_mutated
  
  # Calculate counts for OB patients
  ob_count_mutated <- nrow(subset(ob_mutation_status, !is.na(has_mutation)))
  ob_patient_ids_count <- length(ob_patient_ids)
  ob_count_unmutated <- ob_patient_ids_count - ob_count_mutated
  
  # Create a counts matrix
  counts_matrix <- t(matrix(c(nw_count_mutated, ob_count_mutated, nw_count_unmutated, ob_count_unmutated), nrow = 2))
  colnames(counts_matrix) <- c("Normal Weight", "Obese")
  
  # Perform chi-squared test and store the result in results_df
  p_value <- perform_chi_squared_test(counts_matrix)
  significant_association <- p_value < 0.05
  
  # Store results
  kegg_results_list[[pathway_name]] <- data.frame(pathway = pathway_name, Significant_Association = significant_association, pval = p_value, stringsAsFactors = FALSE)
}

# Combine results into a single dataframe
results_kegg_pathways_chisq <- bind_rows(kegg_results_list)

# Print or use results_kegg_pathways_chisq as needed
print(results_kegg_pathways_chisq)

# Define the file path for saving the CSV of the analysis
kegg_mut_file <- 'kegg_pathway_analysis.csv'
kegg_pathway_mutation_df <- file.path(path_to_TCGA_data, kegg_mut_file)

# Save the dataframe as a CSV file
write.table(results_kegg_pathways_chisq, file = kegg_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)



# ANALYSIS 5: GSEA Mutation Table

# FIXME: link to Carter's GSEA and Adj data
# load the individual csv files
mut_kegg <- read.csv(kegg_pathway_mutation_df)
mut_hallmark <- read.csv(hallmark_pathway_mutation_df)
mut_cancer <- read.csv(cancer_pathway_mutation_df)
gsea_kegg <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/KEGG_Pathways_GSEA.csv')
gsea_hallmark <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Hallmark_Pathways_GSEA.csv')
gsea_cancer <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Cancer_Pathways_GSEA.csv')
adj_kegg <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Kegg_pathways_adj_results.csv')
adj_hallmark <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Hallmark_pathways_adj_results.csv')
adj_cancer <- read.csv('/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Cancer_pathways_adj_results.csv')

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

# Define the file path for saving the CSV of the analysis
kegg_table_file <- 'kegg_analysis_table.csv'
kegg_pathway_table <- file.path(path_to_TCGA_data, kegg_table_file)
cancer_table_file <- 'cancer_analysis_table.csv'
cancer_pathway_table <- file.path(path_to_TCGA_data, cancer_table_file)
hallmark_table_file <- 'hallmark_analysis_table.csv'
hallmark_pathway_table <- file.path(path_to_TCGA_data, hallmark_table_file)

# Write the merged data to a new CSV file
write.csv(merged_kegg, file = kegg_pathway_table, row.names = FALSE)
write.csv(merged_cancer, file = cancer_pathway_table, row.names = FALSE)
write.csv(merged_hallmark, file = hallmark_pathway_table, row.names = FALSE)



# ANALYSIS 6: Multiomics

library(tidyverse)
library(pathwayMultiomics)

# read in data
cancer_df <- read.csv(cancer_pathway_table)
hallmark_df <- read.csv(hallmark_pathway_table)
kegg_df <- read.csv(kegg_pathway_table)

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

# Define the file path for saving the CSV of the analysis
kegg_minimax_file <- 'kegg_minimax.csv'
kegg_minimax_pathway <- file.path(path_to_TCGA_data, kegg_minimax_file)
cancer_minimax_file <- 'cancer_minimax.csv'
cancer_minimax_pathway <- file.path(path_to_TCGA_data, cancer_minimax_file)
hallmark_minimax_file <- 'hallmark_minimax.csv'
hallmark_minimax_pathway <- file.path(path_to_TCGA_data, hallmark_minimax_file)

# Save the data frames as CSV files
write.csv(cancer_res_df, file = cancer_minimax_pathway, row.names = FALSE)
write.csv(hallmark_res_df, file = hallmark_minimax_pathway, row.names = FALSE)
write.csv(kegg_res_df, file = kegg_minimax_pathway, row.names = FALSE)



# ANALYSIS 7: Indel SNV Analysis

library(ggplot2)

# Add INDEL_COUNT and SNP_COUNT columns with zero in each entry
filtered_data <- kirc_data %>%
  select(PATIENT_ID) %>%
  mutate(INDEL_COUNT = 0, SNP_COUNT = 0)
# Loop through each row in data_mutation
for (i in 1:nrow(data_mutation)) {
  # Get the Tumor_Sample_Barcode and Variant_Type for the current row
  tumor_sample <- data_mutation$Tumor_Sample_Barcode[i]
  variant_type <- data_mutation$Variant_Type[i]
  
  # Check if the Tumor_Sample_Barcode is in the PATIENT_ID column of filtered_data
  if (tumor_sample %in% filtered_data$PATIENT_ID) {
    # Get the corresponding row in filtered_data
    patient_row <- which(filtered_data$PATIENT_ID == tumor_sample)
    
    # Update counts based on Variant_Type
    if (variant_type == "SNP") {
      filtered_data$SNP_COUNT[patient_row] <- filtered_data$SNP_COUNT[patient_row] + 1
    } else if (variant_type == "INS" || variant_type == "DEL") {
      filtered_data$INDEL_COUNT[patient_row] <- filtered_data$INDEL_COUNT[patient_row] + 1
    }
  }
}

# Add a new column called BMI_CLASS to filtered_data
filtered_data$BMI_CLASS <- NA  # Initialize the column with NA

# Update BMI_CLASS based on patient IDs
filtered_data$BMI_CLASS[filtered_data$PATIENT_ID %in% ob_patient_ids] <- "OB"
filtered_data$BMI_CLASS[filtered_data$PATIENT_ID %in% nw_patient_ids] <- "NW"


# Create a new column called PROPORTION
filtered_data <- filtered_data %>%
  mutate(PROPORTION = INDEL_COUNT / (INDEL_COUNT + SNP_COUNT))


boxplot(PROPORTION ~ BMI_CLASS, data = filtered_data,
        main = "Proportion of Indels by BMI Classification",
        xlab = "BMI Classification", ylab = "Proportion of Indels",
        col = c("orange", "blue"),
        notch = TRUE)

# Perform a t-test to compare proportions
t_test <- t.test(PROPORTION ~ BMI_CLASS, data = filtered_data)

# Add p-value to the plot
p_value <- format.pval(t_test$p.value)
text(1.5, max(filtered_data$PROPORTION) - 0.05, paste("p-value =", p_value), adj = 0.5)

# Define the file path for saving the PDF of the analysis
indel_snv_file <- 'indel_prop_boxplot.pdf'
indel_snv_pdf_pathway <- file.path(path_to_TCGA_data, indel_snv_file)

# Save the plot as a file
dev.print(pdf, file = indel_snv_pdf_pathway, width = 8, height = 6)

