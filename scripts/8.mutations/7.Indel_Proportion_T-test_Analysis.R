# script to perform Indel and SNV proportion t-test analysis
library(ggplot2)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'
path_to_mutations_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_mutations.txt')
path_to_clinical_data_file <- file.path(base_dir, 'Input/TCGA_&_Clinical_Data.tsv')


# Define the file path for saving the PDF of the analysis
indel_snv_pdf_pathway <- file.path(base_dir, 'Output/Plots/Indel_Proportion_Boxplot.pdf')


# read clinical data
data_clinical <- as.data.frame(read.delim(path_to_clinical_data_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(X.Patient.Identifier) %>%
  rename(PATIENT_ID = X.Patient.Identifier, SUBTYPE = Subtype)


# Filter rows with "KIRC" value in the "CANCER_TYPE_ACRYONYM" column and remove rows with NA values
kirc_data <- data_clinical[!is.na(data_clinical$SUBTYPE) & data_clinical$SUBTYPE == "KIRC", ]
clinical_filtered_columns <- kirc_data[c("PATIENT_ID", "BMI", "BMI.Classification")]


# Add INDEL_COUNT and SNP_COUNT columns with zero in each entry
filtered_data <- kirc_data %>%
  select(PATIENT_ID) %>%
  mutate(INDEL_COUNT = 0, SNP_COUNT = 0)


# Filter obese and normal weight patient IDs
ob_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification %in% c("OB I", "OB II", "OB III")]
nw_patient_ids <- data_clinical$PATIENT_ID[data_clinical$BMI.Classification == "NW" & complete.cases(data_clinical$BMI.Classification)]


# read mutation data
data_mutation <- as.data.frame(read.delim(path_to_mutations_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(Hugo_Symbol)
selected_columns <- data_mutation[c("Hugo_Symbol", "Tumor_Sample_Barcode")]


# Remove "-01" from the end of every entry in the "Tumor_Sample_Barcode" column
data_mutation$Tumor_Sample_Barcode <- gsub("-01$", "", data_mutation$Tumor_Sample_Barcode)


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


# Save the plot as a file
dev.print(pdf, file = indel_snv_pdf_pathway, width = 8, height = 6)
