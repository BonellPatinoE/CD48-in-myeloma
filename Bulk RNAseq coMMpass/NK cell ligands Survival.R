##############
## Survival ##
###############

library(tidyverse)
library(survival)
library(survminer)

# Load survival file
survival <- read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv")
survival <- data.frame(survival)

# Clean survival data
survival_clean <- survival %>%
  mutate(
    patient = PUBLIC_ID,
    status = case_when(
      !is.na(deathdy) & deathdy > 1 ~ 1,  # died
      TRUE ~ 0                           # alive
    ),
    time = coalesce(deathdy, lstalive, lvisitdy)  # fallback to other columns
  ) %>%
  select(patient, time, status)

# Map sample_id → patient (e.g., "MMRF_1052")
cluster_df <- cluster_df %>%
  mutate(patient = str_extract(sample_id, "[0-9]{4}") %>% paste0("MMRF_", .))

survival_merged <- cluster_df %>%
  select(patient, Cluster) %>%
  left_join(survival_clean, by = "patient") %>%
  drop_na(time, status, Cluster)

# Build survival object
surv_object <- Surv(time = survival_merged$time/30.44, event = survival_merged$status)

# Fit survival curves
fit <- survfit(surv_object ~ Cluster, data = survival_merged)

# Plot
ggsurvplot(fit,
           data = survival_merged,
           pval = TRUE,
           risk.table = TRUE,
           palette = RColorBrewer::brewer.pal(3, "Dark2"),
           title = "Overall Survival by NK-Ligand Cluster",
           xlab = "Months",
           ylab = "Survival Probability")


# Use the survival_merged data from your Kaplan-Meier plot
pairwise_survdiff(Surv(time, status) ~ Cluster, data = survival_merged, p.adjust.method = "BH")
