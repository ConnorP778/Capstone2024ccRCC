# Package names
packages <- c("dplyr", "tidyverse", "BiocManager", "DESeq2")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

library(DESeq2)
library(dplyr)
library(tidyverse)
library(BiocManager)

#Change this to match your base directory 
base_dir <- "/Users/elisamcrae/repo/GitHub/Capstone2024ccRCC/"

# Read data
data <- read.table(file = paste0(base_dir, "Input/TCGA-KIRC.mirna.tsv"), header = TRUE)
miRNA_ID = data$miRNA_ID

# Clean data
data <- select(data, -contains(".11"), -miRNA_ID)

# Read clinical data and preprocess
clinical <- read_tsv(file = paste0(base_dir, "Input/TCGA_&_Clinical_Data.tsv"))
#If you run into a problem on this line, sometimes this line works one way or the other. Change it to be "bmiGroup" = "BMI Group".
clinical <- rename(clinical, c("bmiGroup" = "BMI Group"))
clinical <- rename(clinical, c("TCGA.Sample.Code" = "#Patient Identifier"))
c2 = clinical

# Mutate and preprocess clinical data
clinical <- mutate(clinical, bmiGroup = fct_recode(bmiGroup, "Normal" = "NW", "Obese" = "OB", "Overweight" = "OW")) %>%
  select(TCGA.Sample.Code, bmiGroup) %>%
  filter(!is.na(bmiGroup)) %>% filter(bmiGroup == "Obese" | bmiGroup == "Normal")

# Process exprBT
exprBT <- data
colnames(exprBT) <- substr(colnames(exprBT), 1, 12) %>% str_replace_all("\\.", "-")
rownames(exprBT) <- miRNA_ID
exprBT <- round((2^exprBT-1),0)
exprBT <- exprBT[, !duplicated(names(exprBT))]

# Match columns and rename
for (col_name in colnames(exprBT)) {
  if (col_name %in% clinical$TCGA.Sample.Code) {
    obesity <- clinical[clinical$TCGA.Sample.Code == col_name, "bmiGroup"]
    new_col_name <- paste(col_name, obesity, sep = "_")
    colnames(exprBT)[colnames(exprBT) == col_name] <- new_col_name
  } else {
    exprBT <- select(exprBT, -col_name)
  }
}

# Create metadata
meta <- as.tibble(t(exprBT))
meta$rowNames <- row.names(t(exprBT))

# Separate rowNames into ID and bmiGroup columns
meta <- separate(meta, rowNames, into = c("ID", "bmiGroup"), sep = "_", remove = FALSE) %>%
  select(-rowNames)

## Create DESeq object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = exprBT, colData = meta, design = ~ bmiGroup)
dds <- DESeq2::DESeq(dds)
dds <- dds[rowSums(counts(dds)) > 100, ]

# Extract results from DESeq analysis
res1 <- as.data.frame(DESeq2::results(dds))

#### NOW FOR NORMAL AND OVERWEIGHT
c2 <- mutate(c2, bmiGroup = fct_recode(bmiGroup, "Normal" = "NW", "Obese" = "OB", "Overweight" = "OW")) %>%
  select(TCGA.Sample.Code, bmiGroup) %>%
  filter(!is.na(bmiGroup)) %>% filter(bmiGroup == "Overweight" | bmiGroup == "Normal")

# Process exprBT
exprBT <- data
colnames(exprBT) <- substr(colnames(exprBT), 1, 12) %>% str_replace_all("\\.", "-")
rownames(exprBT) <- miRNA_ID
exprBT <- round((2^exprBT-1),0)
exprBT <- exprBT[, !duplicated(names(exprBT))]

# Match columns and rename
for (col_name in colnames(exprBT)) {
  if (col_name %in% c2$TCGA.Sample.Code) {
    obesity <- c2[c2$TCGA.Sample.Code == col_name, "bmiGroup"]
    new_col_name <- paste(col_name, obesity, sep = "_")
    colnames(exprBT)[colnames(exprBT) == col_name] <- new_col_name
  } else {
    exprBT <- select(exprBT, -col_name)
  }
}

# Create metadata
meta1 <- as.tibble(t(exprBT))
meta1$rowNames <- row.names(t(exprBT))

# Separate rowNames into ID and bmiGroup columns
meta1 <- separate(meta1, rowNames, into = c("ID", "bmiGroup"), sep = "_", remove = FALSE) %>%
  select(-rowNames)

## Create DESeq object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = exprBT, colData = meta1, design = ~ bmiGroup)
dds <- DESeq2::DESeq(dds)
dds <- dds[rowSums(counts(dds)) > 100, ]

# Extract results from DESeq analysis
res2 <- as.data.frame(DESeq2::results(dds))
deseq2.sig.res2 <- (res2[!is.na(res2$padj) & res2$padj <= 0.1, ])

print(deseq2.sig.res2)

res1$BMI = "Obese"
res2$BMI = "Overweight"
summary(res2)
res = rbind(res1, res2)

# Create a volcano plot
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1.5, aes(color = ifelse(padj < 0.1 & abs(log2FoldChange) > 1, "Significant", "Non-significant"))) +
  scale_color_manual(values = c("black", "red"), breaks = c("Non-significant", "Significant"), name = "Significance") +  
  facet_wrap(~BMI) +
  theme_bw() +
  labs(x = ~Log[2]*" Fold Change", y = ~-Log[10]*"(p-value)", title = NULL) #+

print(volcano_plot)

ggsave("miRNA_Volcano_Plot.png", plot = volcano_plot, width = 6, height = 4, dpi = 300)

