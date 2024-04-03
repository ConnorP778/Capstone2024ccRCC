library(tidyverse)
# library(readr)
# 
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 2) {
#   stop("Please provide exactly two file paths as command-line arguments.")
# }
# 
# patient_file <- args[1]
# clinical_file <- args[2]
# 
# patientData <- read_tsv(patient_file)
# clinicalData <- read_tsv(clinical_file)

#TODO: Place your filepaths for the TCGA and Clinical data here
patientData = read_tsv("Downloads/kirc_tcga_pan_can_atlas_2018/data_clinical_patient.txt")
clinicalData = read_tsv("Downloads/011724_Clinical TCGA data.xlsx - Merged Complete.tsv")

# Define the classify_bmi function
classify_bmi = function(bmi) {
  if (is.na(bmi)) {
    return('NA')  # Handle missing values
  } else if (bmi < 18.5) {
    return('Underweight')
  } else if (bmi >= 18.5 & bmi < 25) {
    return('NW')
  } else if (bmi >= 25 & bmi < 30) {
    return('OW')
  } else if (bmi >= 30 & bmi < 35) {
    return('OB I')
  } else if (bmi >= 35 & bmi < 40) {
    return('OB II')
  } else if (bmi >= 40) {
    return('OB III')
  } else {
    return('NA')
  }
}

mergedData = full_join(patientData, clinicalData, by = c("#Patient Identifier" = "TCGA Sample Code")) %>%
  filter(!`#Patient Identifier` %in% c("#Identifier to uniquely specify a patient.", "#STRING", "#1", "PATIENT_ID")) %>%
  #Combine Age at Diagnosis Values
  #This renames the Age columns from the TCGA and Clinical
  #It then creates a new column, "Age at Diagnosis", which is a combination of the two
  #If they are different from each other, it defers to the TCGA data. If TCGA data isn't present, it uses the Clinical
  rename(`Age at Diagnosis (Clinical)` = `Age at Diagnosis`) %>%
  rename(`Age at Diagnosis (TCGA)` = `Diagnosis Age`) %>%
  mutate(`Age at Diagnosis (TCGA)` = as.integer(as.numeric(`Age at Diagnosis (TCGA)`)),
         `Age at Diagnosis (Clinical)` = as.integer(as.numeric(`Age at Diagnosis (Clinical)`))) %>%
  mutate(`Age at Diagnosis` = coalesce(`Age at Diagnosis (TCGA)`, `Age at Diagnosis (Clinical)`)) %>%
  #Combine Sex values
  rename(`Sex (TCGA)` = `Sex`) %>%
  rename(`Sex (Clinical)` = `TC`) %>%
  mutate(`Sex` = str_to_title(coalesce(`Sex (TCGA)`, `Sex (Clinical)`))) %>%
  #Combine Race Values
  rename(`Race (TCGA)` = `Race Category`) %>%
  rename(`Race (Clinical)` = `Race`) %>%
  mutate(`Race` = str_to_title(coalesce(`Race (TCGA)`, `Race (Clinical)`))) %>%
  #rename Ethnicity, keeping only the TCGA ethnicity
  rename(`Ethnicity (Clinical)` = `Ethnicity`) %>% #All are NA
  mutate(`Ethnicity (TCGA)` = `Ethnicity Category`) %>%
  rename(`Ethnicity` = `Ethnicity Category`) %>%
  #Combine Cancer Stage
  rename(`Cancer Stage (TCGA)` = `Neoplasm Disease Stage American Joint Committee on Cancer Code`) %>%
  rename(`Cancer Stage (Clinical)` = `Tumor Stage`) %>%
  mutate(`Cancer Stage` = toupper(coalesce(`Cancer Stage (TCGA)`, `Cancer Stage (Clinical)`))) %>%
  #Combine Tumor Stage
  rename(`Tumor Stage (TCGA)` = `American Joint Committee on Cancer Tumor Stage Code`) %>%
  rename(`Tumor Stage (Clinical)` = `T Stage`) %>%
  mutate(`Tumor Stage` = toupper(coalesce(`Tumor Stage (TCGA)`, `Tumor Stage (Clinical)`))) %>%
  #Combine Days to Last Followup
  rename(`Days to Last Followup (TCGA)` = `Last Communication Contact from Initial Pathologic Diagnosis Date`) %>%
  rename(`Days to Last Followup (Clinical)` = `Days to Last followup`) %>%
  mutate(`Overall Survival (Months)` = as.double(as.numeric(`Overall Survival (Months)`))) %>%
  mutate(`Days to Last Followup (TCGA)` = as.integer(as.numeric(`Days to Last Followup (TCGA)`))) %>%
  mutate(`Days to Last Followup (Clinical)` = as.integer(as.numeric(`Days to Last Followup (Clinical)`))) %>%
  mutate(`Days to Last Followup` = ceiling(coalesce(`Days to Last Followup (TCGA)`, 
                                                    `Overall Survival (Months)`*30.4, `Days to Last Followup (Clinical)`))) %>%
  mutate(`Overall Survival (Months)` = coalesce(`Overall Survival (Months)`, round(`Days to Last Followup`/30.4, 8))) %>%
  #Combine Metastasis
  rename(`Metastasis (TCGA)` = `American Joint Committee on Cancer Metastasis Stage Code`) %>%
  rename(`Metastasis (Clinical)` = `Pathologic Distant Metastasis`) %>%
  mutate(`Metastasis` = coalesce(`Metastasis (TCGA)`, `Metastasis (Clinical)`)) %>%
  #Combine Lymph Node Metastasis
  rename(`Lymph Node Metastasis (TCGA)` = `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code`) %>%
  rename(`Lymph Node Metastasis (Clinical)` = `Pathologic Lymph Nodes`) %>%
  mutate(`Lymph Node Metastasis` = coalesce(`Lymph Node Metastasis (TCGA)`, `Lymph Node Metastasis (Clinical)`)) %>%
  #Rename Tumor Status
  rename(`Tumor Status-Tissue Collection (Clinical)` = `Tumor Status-Tissue Collection`) %>%
  rename(`Tumor Status-Last Followup (Clinical)` = `Tumor Status-Last Followup`) %>%
  rename(`Person Neoplasm Cancer Status (TCGA)` = `Person Neoplasm Cancer Status`) %>%
  #Rename Vital Status
  rename(`Vital Status-Enrollment (Clinical)` = `Vital Status-Enrollment`) %>%
  rename(`Vital Status Followup (Clinical)` = `Vital Status Followup`) %>%
  #Add BMI Classification and Group
  mutate(`BMI Classification` = sapply(BMI, classify_bmi), .after = BMI) %>%
  mutate(`BMI Group` = case_when(`BMI Classification` %in% c("OB I", "OB II", "OB III") 
                                 ~ "OB", TRUE ~ `BMI Classification`), .after = `BMI Classification`) %>%
  #Split Overall Living Status
  separate(`Overall Survival Status`, into = c("Overall Survival Status Number", "Overall Survival Status Classification"), sep = ":") %>%
  mutate(`Overall Survival Status Classification` = coalesce(`Overall Survival Status Classification`, `Vital Status Followup (Clinical)`)) %>%
  #TODO: Make another sapply function for Deceeased vs Living
  mutate(`Overall Survival Status Number` = case_when(
    is.na(`Overall Survival Status Number`) & `Vital Status Followup (Clinical)` == "DECEASED" ~ 1,
    is.na(`Overall Survival Status Number`) & `Vital Status Followup (Clinical)` == "LIVING" ~ 0,
    TRUE ~ as.numeric(`Overall Survival Status Number`))) %>%  #Split Disease-specific Survival status
  separate(`Disease-specific Survival status`, into = c("Disease-specific Survival status Number", "Disease-specific Survival status Classification"), sep = ":") %>%
  separate(`Disease Free Status`, into = c("Disease Free Status Number", "Disease Free Status Classification"), sep = ":") %>%
  separate(`Progression Free Status`, into = c("Progression Free Status Number", "Progression Free Status Classification"), sep = ":") %>%
  rename(`Survival Since Diagnosis (Months)` = `Overall Survival (Months)`) #%>%
  #Reorder columns
  # select(`#Patient Identifier`, Subtype, `Age at Diagnosis`, Sex, Race, Ethnicity, `Days to Last Followup`, `Cancer Stage`, `Tumor Stage`,
  #        everything(), 
  #        `Age at Diagnosis (TCGA)`, `Age at Diagnosis (Clinical)`, `Sex (TCGA)`, `Sex (Clinical)`, `Race (TCGA)`, `Race (Clinical)`, `Ethnicity (Clinical)`,
  #        `Days to Last Followup (TCGA)`, `Days to Last Followup (Clinical)`, `Cancer Stage (TCGA)`, `Cancer Stage (Clinical)`, `Tumor Stage (TCGA)`, `Tumor Stage (Clinical)`,
  #        `TCGA PanCanAtlas Cancer Type Acronym`, `Other Patient ID`, `American Joint Committee on Cancer Publication Version Type`, )

#TODO: Add age quartile
# Print and save the results

#print(mergedData)



write_tsv(mergedData, "Downloads/TCGA_&_Clinical_Data.tsv")

#Move columns we'll never use to the end
