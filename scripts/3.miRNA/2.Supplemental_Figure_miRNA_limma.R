#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0267291#sec006
#https://blog.devgenius.io/differential-gene-expression-analysis-using-limma-step-by-step-358da9d41c4e

# Package names
packages <- c("dplyr", "limma", "tidyverse", "BiocManager", "DESeq2", "edgeR")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Load necessary libraries
library(dplyr)
library(limma)
library(edgeR)

#Change this to match your base directory 
base_dir <- "/Users/elisamcrae/repo/GitHub/Capstone2024ccRCC/"

# Read data
data <- read.table(file = paste0(base_dir, "scripts/inputs/TCGA-KIRC.mirna.tsv"), header = TRUE) 

# Remove patient data from normal samples
data <- data %>% select(-contains(".11"), -miRNA_ID)

# Convert expression data to appropriate format
exprBT <- 2^data - 1
rownames(exprBT) <- data$miRNA_ID
colnames(exprBT) <- gsub("\\.", "-", substr(colnames(exprBT), 1, 12))

# Read clinical data
clinical <- read_tsv(paste0(base_dir, "Input/TCGA_&_Clinical_Data.tsv"))
clinical <- rename(clinical, c("bmiGroup" = "BMI Group"))
clinical <- rename(clinical, c("TCGA.Sample.Code" = "#Patient Identifier"))
clinical <- mutate(clinical, bmiGroup = fct_recode(bmiGroup, "Normal" = "NW", "Obese" = "OB", "Overweight" = "OW")) %>%
  select(TCGA.Sample.Code, bmiGroup) %>%
  filter(!is.na(bmiGroup)) %>% filter(bmiGroup == "Obese" | bmiGroup == "Normal")

# Remove duplicated columns
exprBT <- exprBT[, !duplicated(names(exprBT))]

# Update column names based on clinical data
for (col_name in colnames(exprBT)) {
  if (col_name %in% clinical$TCGA.Sample.Code) {
    obesity <- clinical[clinical$TCGA.Sample.Code == col_name, "bmiGroup"]
    new_col_name <- paste(col_name, obesity, sep = "_")
    colnames(exprBT)[colnames(exprBT) == col_name] <- new_col_name
  } else {
    exprBT <- select(exprBT, -col_name)
  }
}

# Create DGEList object
expLIMMA <- exprBT
x <- DGEList(expLIMMA)

# Set up level information for samples
snames <- colnames(x)
group <- substr(snames, 14, nchar(snames))
x$samples$group <- group

# Filter lowly expressed genes
keep.exprs <- filterByExpr(x, group = group)
x <- x[keep.exprs, , keep.lib.sizes = FALSE]

# Calculate library size factors
x <- calcNormFactors(x)

# Create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- gsub("group", "", colnames(design))
colnames(design) <- c("Normal", "Obese")

# Create contrasts
if (length(unique(group)) >= 2) {
  contr.matrix <- makeContrasts(ObesevsNormal = Normal - Obese,
                                levels = colnames(design))
} else {
  print("Error: Not enough levels in the group variable to create contrasts.")
}

# Fit linear model
v <- x$counts
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)

# Empirical Bayes
efit <- eBayes(vfit)

# Treat differentially expressed genes
tfit <- treat(vfit, lfc = 1)
print(summary(decideTests(tfit)))
DEGsTreat <- topTreat(tfit, n = Inf)
