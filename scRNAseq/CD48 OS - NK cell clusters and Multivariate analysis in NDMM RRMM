# =============================================================================
# KAPLAN–MEIER + MULTIVARIATE COX ONLY
# NK Ligand Clustering + CD48 — MMRF CoMMpass IA19
#
# ─────────────────────────────────────────────────────────────────────────────
# 0. LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────
library(DESeq2)
library(gprofiler2)
library(readr)
library(pheatmap)
library(forcats)
library(RColorBrewer)
library(survival)
library(survminer)
library(forestmodel)
library(ggpubr)
library(rstatix)
library(patchwork)
library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS: p-value formatting + pairwise label for KM plots
# ─────────────────────────────────────────────────────────────────────────────
format_p <- function(p) {
  if (is.na(p))   return("NA")
  if (p < 0.001)  return("<0.001")
  if (p < 0.01)   return(paste0("= ", round(p, 3)))
  return(paste0("= ", round(p, 2)))
}
make_pw_label <- function(pw_obj) {
  m <- pw_obj$p.value
  paste0(
    "A vs B: p ", format_p(m["B", "A"]), "\n",
    "A vs C: p ", format_p(m["C", "A"]), "\n",
    "B vs C: p ", format_p(m["C", "B"])
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD RAW COUNTS AND BUILD VST MATRIX
# ─────────────────────────────────────────────────────────────────────────────
raw <- read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
counts <- as.matrix(raw[, -1])
rownames(counts) <- raw$Gene
ids <- data.frame(sample_id = colnames(counts))
rownames(ids) <- ids$sample_id
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ids, design = ~1)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# ─────────────────────────────────────────────────────────────────────────────
# 2. CONVERT ENSEMBL IDs TO HGNC SYMBOLS
# ─────────────────────────────────────────────────────────────────────────────
conversion <- gconvert(query = rownames(vsd_mat), organism = "hsapiens",
                       target = "HGNC", mthreshold = 1, filter_na = FALSE)
vsd_df <- as.data.frame(vsd_mat)
vsd_df$gene_symbol <- conversion$target[match(rownames(vsd_df), conversion$input)]

# ─────────────────────────────────────────────────────────────────────────────
# 3. BUILD NK LIGAND MATRIX AND Z-SCORE SCALE
# ─────────────────────────────────────────────────────────────────────────────
nk_ligands <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E",
                "PCNA", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
                "NECTIN2", "PVR", "CD48")
nk_matrix <- vsd_df %>%
  filter(gene_symbol %in% nk_ligands) %>%
  select(all_of(colnames(vsd_mat)))
rownames(nk_matrix) <- vsd_df$gene_symbol[vsd_df$gene_symbol %in% nk_ligands]
nk_matrix_scaled <- t(scale(t(nk_matrix)))
nk_matrix_scaled[nk_matrix_scaled >  4] <-  4
nk_matrix_scaled[nk_matrix_scaled < -4] <- -4

# ─────────────────────────────────────────────────────────────────────────────
# 4. LOAD CYTOGENETIC DATA
# ─────────────────────────────────────────────────────────────────────────────
translocations_df <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_tumor_only_mm_igtx_pairoscope.tsv")
translocations_df <- data.frame(translocations_df)
translocations_df2 <- data.frame(
  sample_id = translocations_df$SAMPLE,
  patient   = str_extract(translocations_df$SAMPLE, "[0-9]{4}"),
  visit     = str_extract(translocations_df$SAMPLE, "_[0-9]_") %>% str_extract("[0-9]"),
  type      = str_extract(translocations_df$SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
  NSD2      = translocations_df$NSD2_CALL,
  MYC       = translocations_df$MYC_CALL,
  MAF       = translocations_df$MAF_CALL,
  CCND1     = translocations_df$CCND1_CALL,
  MAFB      = translocations_df$MAFB_CALL
) %>%
  filter(type == "BM")

seqFISH_raw <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_gatk_cna_seqFISH.tsv")
seqFISH_raw <- data.frame(seqFISH_raw)
seqFISH <- seqFISH_raw %>%
  transmute(
    sample_id = SAMPLE,
    type      = str_extract(SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
    del17p13  = SeqWGS_Cp_17p13_20percent,
    Gain1q21  = SeqWGS_Cp_1q21_20percent,
    del1p22   = SeqWGS_Cp_1p22_20percent
  ) %>%
  filter(type == "BM")

cyto_combined <- translocations_df2 %>%
  left_join(seqFISH, by = "sample_id", suffix = c("_trans", "_cna"))

# ─────────────────────────────────────────────────────────────────────────────
# 5. RELAPSE STATUS
# ─────────────────────────────────────────────────────────────────────────────
ids <- ids %>%
  mutate(
    visit    = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"),
    relapsed = case_when(
      visit == "1"                 ~ "Newly Diagnosed",
      visit %in% as.character(2:6) ~ "Relapsed",
      TRUE                         ~ NA_character_
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6. BUILD ANNOTATION COLUMN
# ─────────────────────────────────────────────────────────────────────────────
annotation_col <- cyto_combined %>%
  left_join(ids %>% select(sample_id, relapsed), by = "sample_id") %>%
  filter(sample_id %in% colnames(nk_matrix_scaled)) %>%
  select(sample_id, relapsed, NSD2, MYC, MAF, CCND1, MAFB, del17p13, Gain1q21, del1p22) %>%
  column_to_rownames("sample_id")

common_samples   <- intersect(colnames(nk_matrix_scaled), rownames(annotation_col))
nk_matrix_scaled <- nk_matrix_scaled[, common_samples]
annotation_col   <- annotation_col[common_samples, ]
annotation_col[] <- lapply(annotation_col, as.factor)

annotation_colors <- list(
  relapsed  = c("Newly Diagnosed" = "white", "Relapsed" = "black"),
  NSD2      = c("0" = "white", "1" = "#e41a1c"),
  MYC       = c("0" = "white", "1" = "#377eb8"),
  MAF       = c("0" = "white", "1" = "#4daf4a"),
  CCND1     = c("0" = "white", "1" = "#984ea3"),
  MAFB      = c("0" = "white", "1" = "#ff7f00"),
  del17p13  = c("0" = "white", "1" = "#a65628"),
  Gain1q21  = c("0" = "white", "1" = "#f781bf"),
  del1p22   = c("0" = "white", "1" = "#999999")
)
custom_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# ─────────────────────────────────────────────────────────────────────────────
# 7. CLUSTERING — silent heatmap, extract cluster assignments (REQUIRED)
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)
heat <- pheatmap(nk_matrix_scaled,
                 annotation_col    = annotation_col,
                 annotation_colors = annotation_colors,
                 cluster_rows      = TRUE,
                 cluster_cols      = TRUE,
                 cutree_cols       = 3,
                 show_colnames     = FALSE,
                 fontsize_row      = 9,
                 color             = custom_palette,
                 silent            = TRUE)
sample_clusters <- cutree(heat$tree_col, k = 3)
cluster_df <- data.frame(
  sample_id = names(sample_clusters),
  Cluster   = LETTERS[sample_clusters]
) %>%
  mutate(Cluster = factor(Cluster, levels = c("A", "B", "C")))

# Panel C ASSIGNMENTS ONLY (heatmap plot removed). Needed by the Cox join below.
annotation_col$Cluster <- cluster_df$Cluster[match(rownames(annotation_col),
                                                   cluster_df$sample_id)]
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = c("A", "B", "C"))

# ─────────────────────────────────────────────────────────────────────────────
# 8. SURVIVAL DATA  — coalesce(deathdy, lstalive, lvisitdy)
# ─────────────────────────────────────────────────────────────────────────────
Sur_df <- read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv")
Sur_df <- data.frame(Sur_df)
survival_clean <- Sur_df %>%
  mutate(
    patient = PUBLIC_ID,
    status  = case_when(!is.na(deathdy) & deathdy > 1 ~ 1, TRUE ~ 0),
    time    = coalesce(deathdy, lstalive, lvisitdy)
  ) %>%
  select(patient, time, status) %>%
  distinct(patient, .keep_all = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 9. KM SURVIVAL DATA — cluster joined directly to survival (no cyto filter)
# ─────────────────────────────────────────────────────────────────────────────
cluster_df_surv <- cluster_df %>%
  mutate(patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}")))
ids_relapse <- ids %>%
  mutate(patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}")))

survival_merged <- cluster_df_surv %>%
  select(patient, Cluster) %>%
  left_join(survival_clean, by = "patient") %>%
  left_join(ids_relapse %>% select(patient, relapsed) %>%
              distinct(patient, .keep_all = TRUE),
            by = "patient") %>%
  drop_na(time, status, Cluster) %>%
  mutate(
    time_months = time / 30.44,
    Cluster     = factor(Cluster, levels = c("A", "B", "C")),
    relapsed    = factor(relapsed, levels = c("Newly Diagnosed", "Relapsed"))
  )

cat("\n--- KM survival_merged summary ---\n")
cat("Total patients:", nrow(survival_merged), "\n")
print(table(survival_merged$Cluster))
print(table(survival_merged$relapsed))

survival_NDMM <- survival_merged %>%
  filter(as.character(relapsed) == "Newly Diagnosed") %>% droplevels()
survival_RRMM <- survival_merged %>%
  filter(as.character(relapsed) == "Relapsed") %>% droplevels()
cat("NDMM KM n =", nrow(survival_NDMM), "\n")
cat("RRMM KM n =", nrow(survival_RRMM), "\n")

# CD48 RAW COUNTS — all visits (for top/bottom 20% KM)
cd48_raw_all <- raw %>%
  filter(Gene == "ENSG00000117091") %>%
  pivot_longer(-Gene, names_to = "sample_id", values_to = "raw_count") %>%
  select(-Gene) %>%
  mutate(
    patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}")),
    visit   = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"),
    type    = str_extract(sample_id, "_[A-Z]{2}_") %>% str_extract("[A-Z]+")
  ) %>%
  filter(type == "BM")

# ═════════════════════════════════════════════════════════════════════════════
# KAPLAN–MEIER
# ═════════════════════════════════════════════════════════════════════════════

# ── CD48 top vs bottom 20% — ALL PATIENTS (baseline visit) ───────────────────
cd48_baseline <- cd48_raw_all %>%
  filter(visit == "1") %>%
  arrange(patient, visit) %>%
  distinct(patient, .keep_all = TRUE)
n_group <- floor(nrow(cd48_baseline) * 0.20)
surv_20pct <- bind_rows(
  cd48_baseline %>% slice_max(raw_count, n = n_group) %>% mutate(CD48_group = "Top 20%"),
  cd48_baseline %>% slice_min(raw_count, n = n_group) %>% mutate(CD48_group = "Bottom 20%")
) %>%
  left_join(survival_clean, by = "patient") %>%
  drop_na(time, status) %>%
  mutate(time_months = time / 30.44,
         CD48_group  = factor(CD48_group, levels = c("Bottom 20%", "Top 20%")))
cat("\n--- CD48 top vs bottom 20% (all patients) log-rank ---\n")
print(survdiff(Surv(time_months, status) ~ CD48_group, data = surv_20pct))
p_km_cd48 <- ggsurvplot(
  survfit(Surv(time_months, status) ~ CD48_group, data = surv_20pct),
  data = surv_20pct, pval = TRUE, pval.method = TRUE, conf.int = FALSE,
  risk.table = TRUE, risk.table.height = 0.25, break.time.by = 12,
  palette = c("#F8766D", "#00BFC4"),
  legend.labs = c("Bottom 20%", "Top 20%"),
  legend.title = "CD48 mRNA (CoMMpass)", title = "",
  ggtheme = theme_classic(), xlab = "Time (months)")
print(p_km_cd48)

# ── NK-ligand cluster KM — ALL PATIENTS ──────────────────────────────────────
km_all <- survfit(Surv(time_months, status) ~ Cluster, data = survival_merged)
cat("\n--- Cluster KM (all patients) global log-rank ---\n")
print(survdiff(Surv(time_months, status) ~ Cluster, data = survival_merged))
pw_all <- pairwise_survdiff(Surv(time_months, status) ~ Cluster,
                            data = survival_merged, p.adjust.method = "BH")
print(pw_all)
p_all <- ggsurvplot(km_all, data = survival_merged, pval = TRUE, pval.method = TRUE,
                    conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.28,
                    break.time.by = 30, palette = brewer.pal(3, "Dark2"),
                    legend.labs = c("Cluster=A", "Cluster=B", "Cluster=C"),
                    title = "Overall Survival by NK-Ligand Cluster",
                    ggtheme = theme_classic(), xlab = "Months", ylab = "Survival Probability")
p_all$plot <- p_all$plot +
  annotate("label", x = max(survival_merged$time_months, na.rm = TRUE) * 0.55, y = 0.25,
           label = make_pw_label(pw_all), size = 3.5, hjust = 0,
           fill = "white", label.size = 0.3)
print(p_all)

# ── NK-ligand cluster KM — NDMM ──────────────────────────────────────────────
km_NDMM <- survfit(Surv(time_months, status) ~ Cluster, data = survival_NDMM)
cat("\n--- Cluster KM NDMM global log-rank ---\n")
print(survdiff(Surv(time_months, status) ~ Cluster, data = survival_NDMM))
pw_NDMM <- pairwise_survdiff(Surv(time_months, status) ~ Cluster,
                             data = survival_NDMM, p.adjust.method = "BH")
print(pw_NDMM)
p_NDMM <- ggsurvplot(km_NDMM, data = survival_NDMM, pval = TRUE, pval.method = TRUE,
                     conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.28,
                     break.time.by = 30, palette = brewer.pal(3, "Dark2"),
                     legend.labs = paste0("Cluster ", levels(survival_NDMM$Cluster)),
                     title = "Overall Survival by NK-Ligand Cluster — NDMM",
                     ggtheme = theme_classic(), xlab = "Months", ylab = "Survival Probability")
p_NDMM$plot <- p_NDMM$plot +
  annotate("label", x = max(survival_NDMM$time_months, na.rm = TRUE) * 0.55, y = 0.25,
           label = make_pw_label(pw_NDMM), size = 3.5, hjust = 0,
           fill = "white", label.size = 0.3)
print(p_NDMM)

# ── NK-ligand cluster KM — RRMM ──────────────────────────────────────────────
km_RRMM <- survfit(Surv(time_months, status) ~ Cluster, data = survival_RRMM)
cat("\n--- Cluster KM RRMM global log-rank ---\n")
print(survdiff(Surv(time_months, status) ~ Cluster, data = survival_RRMM))
pw_RRMM <- pairwise_survdiff(Surv(time_months, status) ~ Cluster,
                             data = survival_RRMM, p.adjust.method = "BH")
print(pw_RRMM)
p_RRMM <- ggsurvplot(km_RRMM, data = survival_RRMM, pval = TRUE, pval.method = TRUE,
                     conf.int = FALSE, risk.table = TRUE, risk.table.height = 0.28,
                     break.time.by = 30, palette = brewer.pal(3, "Dark2"),
                     legend.labs = paste0("Cluster ", levels(survival_RRMM$Cluster)),
                     title = "Overall Survival by NK-Ligand Cluster — RRMM",
                     ggtheme = theme_classic(), xlab = "Months", ylab = "Survival Probability")
p_RRMM$plot <- p_RRMM$plot +
  annotate("label", x = max(survival_RRMM$time_months, na.rm = TRUE) * 0.55, y = 0.25,
           label = make_pw_label(pw_RRMM), size = 3.5, hjust = 0,
           fill = "white", label.size = 0.3)
print(p_RRMM)

# ── CD48 top/bottom 20% KM — NDMM and RRMM separately ────────────────────────
run_cd48_km <- function(visit_filter, label) {
  cohort <- cd48_raw_all %>%
    filter(visit %in% visit_filter) %>%
    arrange(patient, visit) %>%
    distinct(patient, .keep_all = TRUE)
  cat("\n---", label, ": patients before survival join:", nrow(cohort), "---\n")
  n <- floor(nrow(cohort) * 0.20)
  dat <- bind_rows(
    cohort %>% slice_max(raw_count, n = n) %>% mutate(CD48_group = "Top 20%"),
    cohort %>% slice_min(raw_count, n = n) %>% mutate(CD48_group = "Bottom 20%")
  ) %>%
    left_join(survival_clean, by = "patient") %>%
    drop_na(time, status) %>%
    mutate(time_months = time / 30.44,
           CD48_group  = factor(CD48_group, levels = c("Bottom 20%", "Top 20%")))
  cat(label, ": patients after survival join:", nrow(dat), "\n")
  if (nrow(dat) == 0) stop("No patients remain for ", label)
  cat("\n---", label, ": CD48 top vs bottom 20% log-rank ---\n")
  print(survdiff(Surv(time_months, status) ~ CD48_group, data = dat))
  p <- ggsurvplot(survfit(Surv(time_months, status) ~ CD48_group, data = dat),
                  data = dat, pval = TRUE, pval.method = TRUE, conf.int = FALSE,
                  risk.table = TRUE, risk.table.height = 0.25, break.time.by = 12,
                  palette = c("#F8766D", "#00BFC4"),
                  legend.labs = c("Bottom 20%", "Top 20%"),
                  title = paste0("CD48 Expression (Top vs Bottom 20%) — ", label),
                  ggtheme = theme_classic(), xlab = "Time (months)")
  print(p)
}
run_cd48_km("1",               "NDMM")
run_cd48_km(as.character(2:6), "RRMM")

# ═════════════════════════════════════════════════════════════════════════════
# MULTIVARIATE COX  (CD48 + Clusters + Cytogenetics)  — cytogenetic-filtered path
# ═════════════════════════════════════════════════════════════════════════════
survival_cox <- Sur_df %>%
  mutate(
    patient = paste0("MMRF_", str_extract(PUBLIC_ID, "[0-9]{4}")),
    status  = case_when(!is.na(deathdy) & deathdy > 1 ~ 1, TRUE ~ 0),
    time    = coalesce(deathdy, lstalive, lvisitdy) / 30.44
  ) %>%
  select(patient, time, status) %>%
  distinct(patient, .keep_all = TRUE)

cluster_df_clean <- cluster_df %>%
  mutate(patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}"))) %>%
  distinct(patient, .keep_all = TRUE)

annotation_annot <- annotation_col %>%
  rownames_to_column("sample_id") %>%
  mutate(patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}"))) %>%
  distinct(patient, .keep_all = TRUE)

# Strip any pre-existing Cluster from annotation_annot, then bring in the
# authoritative Cluster from cluster_df_clean -> single Cluster col, no suffixes.
cox_data <- survival_cox %>%
  left_join(annotation_annot %>% select(-any_of("Cluster")), by = "patient") %>%
  left_join(cluster_df_clean %>% select(patient, Cluster), by = "patient") %>%
  drop_na(time, status, Cluster) %>%
  mutate(
    Cluster  = factor(Cluster,  levels = c("A", "B", "C")),
    relapsed = factor(relapsed, levels = c("Newly Diagnosed", "Relapsed")),
    NSD2     = factor(NSD2),  MYC   = factor(MYC),   MAF     = factor(MAF),
    CCND1    = factor(CCND1), MAFB  = factor(MAFB),  del1p22 = factor(del1p22),
    del17p13 = factor(del17p13), Gain1q21 = factor(Gain1q21)
  )

cd48_vst_patient <- vsd_df %>%
  filter(gene_symbol == "CD48") %>%
  select(-gene_symbol) %>%
  pivot_longer(everything(), names_to = "sample_id", values_to = "CD48_vst") %>%
  mutate(
    patient = paste0("MMRF_", str_extract(sample_id, "[0-9]{4}")),
    visit   = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"),
    type    = str_extract(sample_id, "_[A-Z]{2}_") %>% str_extract("[A-Z]+")
  ) %>%
  filter(type == "BM") %>%
  arrange(patient, visit) %>%
  distinct(patient, .keep_all = TRUE) %>%
  select(patient, CD48_vst)

cox_data_cd48 <- cox_data %>%
  left_join(cd48_vst_patient, by = "patient") %>%
  drop_na(CD48_vst)
cat("\n--- Cox n =", nrow(cox_data_cd48), "---\n")
print(table(cox_data_cd48$Cluster))

# ── Full multivariate Cox ────────────────────────────────────────────────────
cox_full <- coxph(
  Surv(time, status) ~ CD48_vst + Cluster + relapsed +
    NSD2 + MYC + MAF + MAFB + del1p22 + del17p13 + Gain1q21,
  data = cox_data_cd48)
cat("\n--- Full Cox summary ---\n")
print(summary(cox_full))
print(forest_model(cox_full) +
        ggtitle("Multivariate Cox Model\nCD48 + Clusters + Cytogenetics") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5)))

# ── Multivariate Cox — NDMM ──────────────────────────────────────────────────
cox_data_NDMM <- cox_data_cd48 %>%
  filter(as.character(relapsed) == "Newly Diagnosed") %>% droplevels()
cox_data_RRMM <- cox_data_cd48 %>%
  filter(as.character(relapsed) == "Relapsed") %>% droplevels()
cat("\nCox NDMM n =", nrow(cox_data_NDMM), "  Cox RRMM n =", nrow(cox_data_RRMM), "\n")

cox_NDMM <- coxph(
  Surv(time, status) ~ CD48_vst + Cluster + NSD2 + MYC + MAF + MAFB +
    del1p22 + del17p13 + Gain1q21, data = cox_data_NDMM)
cat("\n--- Multivariate Cox — NDMM ---\n")
print(summary(cox_NDMM))
print(forest_model(cox_NDMM) +
        ggtitle("Multivariate Cox — NDMM\nCD48 + NK Clusters + Cytogenetics") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5)))

# ── Multivariate Cox — RRMM ──────────────────────────────────────────────────
cox_RRMM <- coxph(
  Surv(time, status) ~ CD48_vst + Cluster + NSD2 + MYC + MAF + MAFB +
    del1p22 + del17p13 + Gain1q21, data = cox_data_RRMM)
cat("\n--- Multivariate Cox — RRMM ---\n")
print(summary(cox_RRMM))
print(forest_model(cox_RRMM) +
        ggtitle("Multivariate Cox — RRMM\nCD48 + NK Clusters + Cytogenetics") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5)))

# =============================================================================
