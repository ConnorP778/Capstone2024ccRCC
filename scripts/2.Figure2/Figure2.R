install.packages('survminer')
install.packages('forestplot')

library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(readr)
library(forestplot)

#Add Base Directory Here
base_dir = "C:\\Users\\biqe2\\Downloads\\Capstone2024ccRCC-main\\Capstone2024ccRCC-main\\"

data = read_tsv(paste0(base_dir,"Input/TCGA_&_Clinical_Data.tsv"))

quartiles <- quantile(data$`Age at Diagnosis`, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)

# Add new columns that contain BMI Group, Age group, and Cancer Stage Group
data$`BMI Classification` = factor(data$`BMI Classification`, 
                                   levels = c('Underweight', 'NW', 'OW', 'OB I', 'OB II', 'OB III'))
data$Sex = as.factor(toupper(data$Sex))
data$`Cancer Stage` = factor(toupper(data$`Cancer Stage`), levels = c("STAGE I", "STAGE II", "STAGE III", "STAGE IV"))
data = data %>% filter(!is.na(`BMI Classification`))
data = mutate(data, `Max Tumor Dimension` = as.numeric(`Max Tumor Dimension`)) %>%
  filter(`BMI Classification` != "Underweight") %>%
  rename(`PatientID` = `#Patient Identifier`) %>%
  mutate(`BMI Group` = case_when(`BMI Classification` %in% c("OB I", "OB II", "OB III") ~ "OB", TRUE ~ `BMI Classification`)) %>%
  mutate(`BMI Group` = factor(toupper(`BMI Group`), levels = c("NW", "OW", "OB"))) %>%
  filter(`Survival Since Diagnosis (Months)` != 0) %>%
  #mutate(`Age Quartile2` = ntile(`Age at Diagnosis`, 4)) %>%
  mutate(`Above or Below 60` = case_when(`Age at Diagnosis` >= 60 ~ "Above 60", `Age at Diagnosis` < 60 ~ "Below 60")) %>%
  mutate(`Above or Below 60` = factor(`Above or Below 60`, levels = c("Below 60", "Above 60"))) %>%
  mutate(`Cancer Stage Group` = case_when(
      `Cancer Stage` %in% c("STAGE I", "STAGE II") ~ "Stage I/II",
      `Cancer Stage` %in% c("STAGE III", "STAGE IV") ~ "Stage III/IV", TRUE ~ `Cancer Stage`)) %>%
  mutate(`Cancer Stage Group` = factor(`Cancer Stage Group`, levels = c("Stage I/II", "Stage III/IV"))) #%>%
  #write_tsv("Downloads/TCGA_Clinical_Combined_With_BMI.tsv")

#Colors for survival curve
nw_color <- "#56B4E9"  # Pastel Green
ow_color <- "#009E73"  # Pastel Blue
ob_colors <- c("#f7c8c9", "#D55E00", "#ff3a3a")  # Gradient of pastel reds
colorBlindFriendly = c("#009E73", "#56B4E9", "#D55E00")


#Get BMI Counts
data_counts <- data %>%
  group_by(`BMI Group`) %>%
  summarise(Count = n())
data <- merge(data, data_counts, by = "BMI Group", all.x = TRUE)

ggplot(data, aes(x=`BMI Group`)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = ..count.., y = ..count..), vjust = -0.5) +
  theme_bw() +
  labs(y = "Count")
ggsave(paste0(base_dir, "Output/Plots/ccRCC_BMI_Histogram.png"), device = "png", height = 8)

summary_data <- data %>%
  group_by(`BMI Group`) %>%
  summarise(mean_BMI = mean(`BMI`), count = n())

ggplot(data, aes(x = `BMI Group`, y = `BMI`)) +
  #geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.5) +
  geom_violin(alpha = 0.3) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill="orange") + 
  geom_text(data = summary_data, aes(label = count, y = mean_BMI), vjust = 2) +
  theme_bw()
ggsave(paste0(base_dir, "Output/Plots/ccRCC_BMI_ViolinPlot.png"), device = "png", height = 8)


set.seed(123)
data2 <- data.frame(
  `BMI Classification` = rep(c("NW", "OW", "OB I", "OB II", "OB III"), each = 20),
  `Survival Since Diagnosis (Months)` = rnorm(100),
  `Overall Survival Status Classification` = rep(c("Living", "Deceased"), each = 50))

nw_color <- "#56B4E9"  # Pastel Green
ow_color <- "#009E73"  # Pastel Blue
ob_colors <- c("#f7c8c9", "#D55E00", "#ff3a3a")  # Gradient of pastel reds

colorBlindFriendly = c("#009E73", "#56B4E9", "#D55E00")


#Kaplan meier (4 Kaplan meier NW to OW, OB1,2,3)
data <- data %>%
  rename(Months = `Survival Since Diagnosis (Months)`) %>%
  rename(Survival = `Overall Survival Status Number`) %>%
  rename(BMI_Class = `BMI Classification`) %>%
  rename(BMI_Group = `BMI Group`)


group_sizes <- data %>%
  group_by(BMI_Group) %>%
  summarise(Count = n())
totalObservations = sum(group_sizes$Count)
print(totalObservations)
group_labels <- setNames(group_sizes$Count, group_sizes$BMI_Group)

fit <- survfit(Surv(Months, Survival) ~ BMI_Group, data = data)
custom_legend_labs <- sapply(levels(data$BMI_Group), function(x) paste(x, " (", group_labels[x], ")", sep=""))

survPlot = ggsurvplot(fit, data = data, pval = TRUE, conf.int = FALSE,
           title = paste("ccRCC Survival by BMI Group (", totalObservations, ")", sep = ""),
           #title = NULL,
           break.time.by = 12,
           xlab = "Time (months)", ylab = "Survival probability",
           legend.title = "BMI Classification:",
           legend.labs = custom_legend_labs,
           risk.table = TRUE,
           risk.table.title = "Number at Risk",
           risk.table.col = "strata",
           tables.theme = theme_cleantable(),
           palette = c(nw_color, ow_color, ob_colors[2]))#c("black", "black", "black"))

ggsave(paste0(base_dir, "Output/Plots/ccRCC_KaplanMeier_Table.png"), device = "png", width = 8, height = 3)

survPlot$plot = survPlot$plot + theme(legend.title = element_text(size = 12))
ggsave(paste0(base_dir, "Output/Plots/ccRCC_KaplanMeier_Curve.png"), plot = survPlot$plot, device = "png", width = 8, height = 6, dpi = 300)

print(survPlot)

cox_model = coxph(Surv(Months, Survival) ~ BMI_Group + `Above or Below 60` + Sex + `Cancer Stage Group`, data = data) # + `Tumor Grade`
summary(cox_model)

hazard_ratios <- exp(coef(cox_model))
conf_intervals <- exp(confint(cox_model))

# Create a data frame for forest plot
forest_data <- data.frame(
  variables = names(hazard_ratios),
  HR = hazard_ratios,
  lower = conf_intervals[, 1],
  upper = conf_intervals[, 2]
)

data2 = filter(data, !is.na(`Tumor Grade`)) %>%
  mutate(`Tumor Grade` = factor(`Tumor Grade`, levels = c("G1", "G2", "G3", "G4"))) %>%
  mutate(`Tumor Grade Group` = case_when(
    `Tumor Grade` %in% c("G1", "G2") ~ "G1/G2", #Less Aggressive Cancers
    `Tumor Grade` %in% c("G3", "G4") ~ "G3/G4", TRUE ~ `Tumor Grade`)) %>% #More Aggressive Cancers
  mutate(`Tumor Grade Group` = factor(`Tumor Grade Group`, levels = c("G1/G2", "G3/G4")))


cox_model2 <- coxph(Surv(Months, Survival) ~ BMI_Group + `Above or Below 60` + Sex + `Cancer Stage Group` + `Tumor Grade Group`, data = data2)

# summary(cox_model2)

hazard_ratios2 <- exp(coef(cox_model2))
conf_intervals2 <- exp(confint(cox_model2))

# Create a data frame for forest plot
forest_data2 <- data.frame(
  variables = names(hazard_ratios2),
  HR = hazard_ratios2,
  lower = conf_intervals2[, 1],
  upper = conf_intervals2[, 2]
)

p_values <- summary(cox_model2)$coefficients[, "Pr(>|z|)"]

significance_levels <- c("***", "**", "*")  # You can customize this based on your preference

forest_data2$significance <- ifelse(p_values < 0.001, significance_levels[1],
                                    ifelse(p_values < 0.01, significance_levels[2],
                                           ifelse(p_values < 0.05, significance_levels[3], "")))
forest_data2$pValue <- p_values

ggplot(data = forest_data2, aes(x = variables, y = HR, label = significance)) +
  geom_point(size = 3) +  # This adds the dots for Hazard Ratios
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # This adds the whiskers for confidence intervals
  coord_flip() +  # Flips the coordinates to make the plot horizontal
  theme_bw() +
  geom_text(aes(label = significance), vjust = -0.5, color = "black", size = 4) +  # Adjust vjust to position text
  labs(x = "Variables", y = "Hazard Ratio") +
  geom_hline(yintercept = 1, linetype = "solid", color = "grey35") +
  scale_y_continuous(breaks = seq(floor(min(forest_data2$HR)), ceiling(max(forest_data2$upper)), by = 0.5)) +
  scale_x_discrete(labels = rev(c("Male\n(Ref: Female)", "OW\n(Ref: NW)", "OB\n(Ref: NW)",
                                  "Tumor Grades G3/G4\n(Ref: Grades G1/G2)", "Cancer Stages III/IV\n(Ref: Stages I/II)",
                                  "Age Above 60\n(Ref: Age Below 60)"))) + #, "Grade 2", "Grade 3", "Grade 4"
  annotate("text", x = 5.5, y = 3.55, label = "* = p < 0.05  \n*** = p < 0.001", hjust = 1, size = 3)

ggsave(paste0(base_dir, "Output/Plots/ccRCC_coxPropModel.png"), device = "png", height = 8)

print(forest_data2)
