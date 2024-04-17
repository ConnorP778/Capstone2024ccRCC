library(dplyr)

# file addresses
mutation_file_address <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/kirc_tcga_pan_can_atlas_2018/data_mutations.txt'
clinical_data_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/PatientDataWithBMI.tsv'
pathway_genes_path <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/PathwayGenes.csv'
pathway_mutation_df <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/pathway_mutations.csv'


# Read pathway genes
pathway_genes <- read.csv(pathway_genes_path, fileEncoding = "UTF-8")

# Read clinical data
data_clinical <- read.delim(clinical_data_path, header = TRUE, sep = "\t", dec = ".") %>%
  arrange(PATIENT_ID)

# Filter obese and normal weight patient IDs
ob_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification %in% c("OB I", "OB II", "OB III")]
nw_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification == "NW" & complete.cases(data_clinical$BMI.Classification)]

# Initialize list to store results
results_list <- list()

# Loop through each pathway
for (pathway_name in names(pathway_genes)) {
  pathway_genes_current <- pathway_genes[[pathway_name]]
  
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
  
  # Perform chi-squared test
  p_value <- perform_chi_squared_test(counts_matrix)
  significant_association <- p_value < 0.05
  
  # Store results
  results_list[[pathway_name]] <- data.frame(Pathway = pathway_name, Significant_Association = significant_association, p_value = p_value, stringsAsFactors = FALSE)
}

# Combine results into a single dataframe
results_pathway_chisq <- bind_rows(results_list)

# Print or use results_pathway_chisq as needed
print(results_pathway_chisq)

# Save the dataframe as a CSV file
write.table(results_pathway_chisq, file = pathway_mutation_df, sep = ",", quote = FALSE, row.names = FALSE)
