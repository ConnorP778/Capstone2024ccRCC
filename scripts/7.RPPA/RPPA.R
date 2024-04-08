library(tidyverse)

# Change this to match your base directory
base_dir <- "/Users/ConnorP778/repo/GitHub/Capstone2024ccRCC"


rppa_data = read.table(file = paste0(base_dir, "data_rppa.txt"), header=TRUE, row.names=1, na.strings="NA", sep="\t")

col_names <- colnames(rppa_data)

# Apply transformations to each column name
new_col_names <- gsub("\\.", "-", substr(col_names, 1, nchar(col_names) - 3))

# Assign the new column names to the RPPA data table
colnames(rppa_data) <- new_col_names

clinical_data = read_tsv("TCGA_&_Clinical_Data_TEMP.txt")

clinical_data = rename(clinical_data, BMI_Group = 'BMI Group')
clinical_data = rename(clinical_data, Patient_ID = '#Patient Identifier')

filter(clinical_data, BMI_Group == "OB") %>%
  pull(Patient_ID) -> ob_ids

filter(clinical_data, BMI_Group == "NW") %>%
  pull(Patient_ID) -> nw_ids


ob_data <- rppa_data[, colnames(rppa_data) %in% ob_ids]
nw_data <- rppa_data[, colnames(rppa_data) %in% nw_ids]

merged_data <- merge(ob_data, nw_data, by = "row.names", all = FALSE)

merged_data <- na.omit(merged_data)

protein_names <- merged_data[, 1]

# Separate the merged data into two matrices for the comparison
mat1 <- as.matrix(merged_data[, colnames(merged_data) %in% colnames(ob_data)])
mat2 <- as.matrix(merged_data[, colnames(merged_data) %in% colnames(nw_data)])


row.names(mat1) <- protein_names
row.names(mat2) <- protein_names

all_data <- cbind(mat1, mat2)

design <- model.matrix(~ 0 + factor(c(rep("Group1", ncol(mat1)), rep("Group2", ncol(mat2)))))

# Fit linear model
fit <- lmFit(all_data, design)

# Performs empirical Bayes moderation
fit <- eBayes(fit)

# Extract differentially expressed proteins
de_proteins <- topTable(fit, coef = 1, n = Inf)

de_proteins <- transform(de_proteins, Significance = ifelse(abs(logFC) > 1 & adj.P.Val < 0.05, "Significant", "Not Significant"))

# Volcano plot
volcano = ggplot(de_proteins, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(size = 2, aes(color = Significance)) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
  theme_bw()







                



