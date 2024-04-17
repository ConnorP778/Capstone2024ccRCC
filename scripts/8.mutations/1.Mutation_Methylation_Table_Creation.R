# script to create a table with VHL mutation and methylation status for each patient (matches original script)


# install packages, if you don't already have them
install.packages("dplyr")
install.packages("tidyr")

library(dplyr)
library(tidyr)


# Update the base_dir to match where the github repo was downloaded
base_dir <- '/Users/tufts/OneDrive/Winter 2024/Bioinformatic Capstone/Mutation Project/Capstone2024ccRCC'
path_to_mutations_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_mutations.txt')
path_to_methylation_file <- file.path(base_dir, 'Input/kirc_tcga_pan_can_atlas_2018/data_methylation_hm27_hm450_merged.txt')
path_to_clinical_data_file <- file.path(base_dir, 'Input/TCGA_&_Clinical_Data.tsv')
output_table_dir <- file.path(base_dir, 'Output/Data/mut_meth_table.tsv')


# read clinical data
clinical_data <- read.delim(path_to_clinical_data_file, header = TRUE, sep = "\t") %>%
  arrange(X.Patient.Identifier) %>%
  rename(PATIENT_ID = X.Patient.Identifier, SUBTYPE = Subtype)

# Filter the data frame to select rows where SUBTYPE is "KIRC" and extract PATIENT_ID
kirc_patient_ids <- clinical_data$PATIENT_ID[clinical_data$SUBTYPE == "KIRC"]

# read methylation data
data_methylation <- as.data.frame(read.delim(path_to_methylation_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(NAME) %>%
  na.omit()

# process methylatyon data
VHL_data = data_methylation %>%
  filter(NAME == 'VHL') %>%
  select(-ENTITY_STABLE_ID,-TRANSCRIPT_ID) %>%
  pivot_longer(TCGA.3Z.A93Z.01:TCGA.T7.A92I.01, names_to = 'Patient', values_to = "Meth_Level")  %>%
  group_by(DESCRIPTION,Patient) %>%
  summarise(average = mean(Meth_Level)) %>%
  mutate(type_meth = case_when(average < 0.2 ~ 'low',
                               average <= 0.5 ~ 'med',
                               average > 0.5 ~ 'high'))


# replace . with - so patient IDs align
VHL_data$Patient <- gsub('\\.', '-', VHL_data$Patient)


# read mutation data
data_mutation <- as.data.frame(read.delim(path_to_mutations_file, header = TRUE, sep = "\t", dec = ".")) %>%
  arrange(Hugo_Symbol)

# process mutation data
VHL_mut_data <- data_mutation %>%
  filter(Hugo_Symbol == 'VHL') %>%
  select(Hugo_Symbol, Tumor_Sample_Barcode)  

VHL_mut_data <- distinct(VHL_mut_data, Tumor_Sample_Barcode, .keep_all = TRUE)
VHL_mut_data$Hugo_Symbol <- gsub('VHL', 'Mutated', VHL_mut_data$Hugo_Symbol)


# merge mutation/methylation dataframes based on patient IDs
joined_df <- left_join(VHL_data, VHL_mut_data, by = c("Patient" = "Tumor_Sample_Barcode")) %>%
  mutate(Hugo_Symbol = coalesce(Hugo_Symbol, "Not Mutated")) %>%
  rename('TCGA Number' = Patient, 
         'VHL Methylation Status' = type_meth,
         'VHL Mutation Status' = Hugo_Symbol)

# remove the average column and the -01 from each patient ID, keep only the KIRC patients
joined_df <- joined_df[, -1] %>%
  select(-average)
joined_df$'TCGA Number' <- sub("-01$", "", joined_df$'TCGA Number')
kirc_df <- joined_df[joined_df$'TCGA Number' %in% kirc_patient_ids, ]

# Save the dataframe as a TSV file
write.table(kirc_df, file = output_table_dir, sep = "\t", quote = FALSE, row.names = FALSE)
