# script to perform chi square analysis on mutation counts for cancer, kegg, and hallmark pathways
library(dplyr)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'
path_to_mutations_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_mutations.txt')
path_to_clinical_data_file <- file.path(base_dir, 'Input/TCGA_&_Clinical_Data.tsv')
cancer_pathway_genes_path <- file.path(base_dir, 'Input/pathway_annotations/cancer.csv')
hallmark_pathway_genes_path <- file.path(base_dir, 'Input/pathway_annotations/hallmark.csv')
kegg_pathway_genes_path <- file.path(base_dir, 'Input/pathway_annotations/kegg.csv')
cancer_pathway_mutation_df <- file.path(base_dir, 'Output/Plots/cancer_mut_BMI_analysis.csv')
hallmark_pathway_mutation_df <- file.path(base_dir, 'Output/Plots/hallmark_mut_BMI_analysis.csv')
kegg_pathway_mutation_df <- file.path(base_dir, 'Output/Plots/kegg_mut_BMI_analysis.csv')


# read mutation data
data_mutation <- as.data.frame(read.delim(path_to_mutations_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(Hugo_Symbol)
selected_columns <- data_mutation[c("Hugo_Symbol", "Tumor_Sample_Barcode")]
# Remove "-01" from the end of every entry in the "Tumor_Sample_Barcode" column
data_mutation$Tumor_Sample_Barcode <- gsub("-01$", "", data_mutation$Tumor_Sample_Barcode)


# read pathway genes
cancer_pathway_genes <- read.csv(cancer_pathway_genes_path, fileEncoding = "UTF-8")
kegg_pathway_genes <- read.csv(kegg_pathway_genes_path, fileEncoding = "UTF-8")
hallmark_pathway_genes <- read.csv(hallmark_pathway_genes_path, fileEncoding = "UTF-8")


# read clinical data
data_clinical <- read.delim(path_to_clinical_data_file, header = TRUE, sep = "\t") %>%
  arrange(X.Patient.Identifier) %>%
  rename(PATIENT_ID = X.Patient.Identifier, SUBTYPE = Subtype)


# Filter obese and normal weight patient IDs
ob_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification %in% c("OB I", "OB II", "OB III")]
nw_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification == "NW" & complete.cases(data_clinical$BMI.Classification)]


# function to run chi_square analysis
perform_chi_squared_test <- function(counts_matrix) {
  chi_square_test <- chisq.test(counts_matrix)
  p_value <- chi_square_test$p.value
  return(p_value)
}


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


# Save the dataframe as a CSV file
write.table(cancer_results_pathway_chisq, file = cancer_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)


# Initialize list to store results
hallmark_results_list <- list()


# Loop through each KEGG pathway
for (i in 1:nrow(hallmark_pathway_genes)) {
  pathway_name <- hallmark_pathway_genes$name[i]
  pathway_genes_current <- unlist(strsplit(gsub("\"", "", hallmark_pathway_genes$value_str[i]), ", "))
  
  
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


# Save the dataframe as a CSV file
write.table(results_hallmark_pathways_chisq, file = hallmark_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)


# Initialize list to store results
kegg_results_list <- list()


# Loop through each KEGG pathway
for (i in 1:nrow(kegg_pathway_genes)) {
  pathway_name <- kegg_pathway_genes$name[i]
  pathway_genes_current <- unlist(strsplit(gsub("\"", "", kegg_pathway_genes$value_str[i]), ", "))
  
  
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


# Save the dataframe as a CSV file
write.table(results_kegg_pathways_chisq, file = kegg_pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)
