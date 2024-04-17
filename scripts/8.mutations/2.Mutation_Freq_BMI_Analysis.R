# script for mutation frequency and BMI classification logistic regression analysis (matches original script)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'
path_to_mutations_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_mutations.txt')
path_to_clinical_data_file <- file.path(base_dir, 'Input/TCGA_&_Clinical_Data.tsv')
output_plot_dir <- file.path(base_dir, 'Output/Plots/mut_freq_plot.pdf')


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


# Open a PDF device to start recording the plot
pdf(output_plot_dir)


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
