library(tidyverse)
library(dplyr)
library(readr)

#TODO: Change to your own base directory
base_dir = ""

data = read_tsv(paste0(base_dir, "Input/TCGA_&_Clinical_Data.tsv")) %>%
  filter(`BMI Classification` != "Underweight")

mtData = read_csv(paste0(base_dir, "Input/Supplementary_file_1.csv")) %>%
  filter(Study=="KIRC") %>%
  mutate(`#Patient Identifier` = paste0("TCGA-", `Sample ID`)) %>%
  rename(`Mitochondria Log2 Fold Change` = `Log2 Fold Change`) %>%
  select(`#Patient Identifier`, `Mitochondria Log2 Fold Change`)

mergedData = full_join(data, mtData)

mergedData <- mergedData %>%
  drop_na(`BMI Group`, `Mitochondria Log2 Fold Change`)

ggplot(mergedData, aes(x = `BMI Group`, y = `Mitochondria Log2 Fold Change`)) +
  geom_boxplot() +
  stat_summary(fun=mean, color ="#FDAE6B", shape = 23, fill ="#FDAE6B") +
  theme_bw()

ggsave(paste0(base_dir, "Output/Plots/ccRCC_mtDNA_CopyNumberAnalysis.png"), device = "png", height = 8)
