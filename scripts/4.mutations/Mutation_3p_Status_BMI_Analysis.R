# script for 3p status and BMI classification logistic regression analysis (matches original script)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'
path_to_arm_level_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_armlevel_cna.txt')
path_to_clinical_data_file <- file.path(base_dir, 'Input/TCGA_&_Clinical_Data.tsv')
output_plot_dir <- file.path(base_dir, 'Output/Plots/mut_3pstatus_plot.pdf')
mut_meth_path <- file.path(base_dir, 'Output/Plots/mut_meth_table.tsv')


# read clinical data
clinical_data <- read.delim(path_to_clinical_data_file, header = TRUE, sep = "\t") %>%
  arrange(X.Patient.Identifier) %>%
  rename(PATIENT_ID = X.Patient.Identifier, SUBTYPE = Subtype)

# Filter the data frame to select rows where SUBTYPE is "KIRC" and extract PATIENT_ID
kirc_patient_ids <- clinical_data$PATIENT_ID[clinical_data$SUBTYPE == "KIRC"]


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


# Open a PDF device to start recording the plot
pdf(output_plot_dir)

# Visualize the association with custom x-axis range
plot(merged_data$Status.3p, merged_data$BMI, 
     xlab = "Chr 3p Status", ylab = "BMI", col = "blue")

# Add a red line indicating the linear trend
abline(lm(merged_data$Status.3p ~ merged_data$BMI), col = "red")

# Create the title including the p-value
title <- paste("Association between Chr 3p status and BMI (p-value =", round(p_value, 4), ")")

# Add title to the plot
title(title, line = 2.5)

# Close the PDF device to save the plot
dev.off()
