# Load required libraries
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(readr)
library(pheatmap)
library(forcats)

### Multivariate COX MODEL

# 1. Clean cluster_df — keep one sample per patient
cluster_df_clean <- cluster_df %>%
  mutate(patient = str_extract(sample_id, "[0-9]{4}") %>% paste0("MMRF_", .)) %>%
  distinct(patient, .keep_all = TRUE)  # one row per patient

# 2. Clean annotation_col
annotation_annot <- annotation_col %>%
  rownames_to_column("sample_id") %>%
  mutate(patient = str_extract(sample_id, "[0-9]{4}") %>% paste0("MMRF_", .)) %>%
  distinct(patient, .keep_all = TRUE)

# First join: survival + annotation
merged_surv_annot <- survival_clean %>%
  left_join(annotation_annot, by = "patient")

# Then add Cluster safely
cox_data <- merged_surv_annot %>%
  left_join(cluster_df_clean %>% select(patient, Cluster), by = "patient")

names(cox_data)
# Make sure "Cluster" is there

cox_data <- cox_data %>%
  mutate(
    Cluster = Cluster.y  # keep the one from cluster_df_clean
  ) %>%
  select(-Cluster.x, -Cluster.y) %>%  # drop the extras
  drop_na(time, status, Cluster)

cox_data <- cox_data %>%
  mutate(
    Cluster = factor(Cluster),
    relapsed = factor(relapsed),
    NDS2 = factor(NDS2),
    MYC = factor(MYC),
    del1p22 = factor(del1p22),
    MAF = factor(MAF),
    MAFB = factor(MAFB),
    del17p13 = factor(del17p13),
    Gain1q21 = factor(Gain1q21),
  )

cox_model <- coxph(
  Surv(time, status) ~ Cluster + relapsed + NDS2 + MYC +
    MAF + MAFB + del1p22 + del17p13 + Gain1q21,
  data = cox_data
)

summary(cox_model)

ggforest(cox_model,
         data = cox_data,
         main = "Hazard Ratios for NK Ligand Clusters and Cytogenetics",
         fontsize = 1.0,
         cpositions = c(0.05, 0.22, 0.4),  # adjust label + CI positioning
         refLabel = "Reference",           # label for reference level
         noDigits = 3                      # round estimates
)

#––– assume you’ve already fit:
# cox_model <- coxph(Surv(time, status) ~ Cluster + relapsed + NDS2 + MYC + del17p13 + Gain1q21,
#                    data = cox_data)

# load the plotting helper
library(forestmodel)

# build a forest plot
forest_model(cox_model) +
  ggtitle("Multivariate Cox Model\nNK Ligand Clusters + Cytogenetic Risks") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 10)
  )
