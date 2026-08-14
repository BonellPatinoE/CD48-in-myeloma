# ============================================================
# CD48 EXPRESSION IN vk*MYC MYELOMA PLASMA CELLS 
# ============================================================

set.seed(1234)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(scales)
})

# ============================================================
# CONSTANTS
# ============================================================
# Control is excluded from myeloma comparisons —
# no MYC+ plasma cells exist in control mice by definition
comparisons_myeloma <- list(
  c("ND", "MM"),
  c("ND", "AMG"),
  c("MM", "AMG")
)

myeloma_levels <- c("ND", "MM", "AMG")
myeloma_colors <- c("ND" = "#FC8D62", "MM" = "#8DA0CB", "AMG" = "#E78AC3")
all_colors     <- c("Control" = "#66C2A5", "ND" = "#FC8D62",
                    "MM"      = "#8DA0CB", "AMG" = "#E78AC3")

# Seurat v5 helpers
seurat_v5  <- packageVersion("Seurat") >= "5.0.0"
get_layer  <- function(obj, layer = "data",   assay = "RNA") {
  if (seurat_v5) GetAssayData(obj, assay = assay, layer = layer)
  else           GetAssayData(obj, assay = assay, slot  = layer)
}
get_counts <- function(obj, assay = "RNA") {
  if (seurat_v5) GetAssayData(obj, assay = assay, layer = "counts")
  else           GetAssayData(obj, assay = assay, slot  = "counts")
}

# ============================================================
# 1. LOAD & GATE CELLS
# ============================================================
cat("Loading object...\n")
seurat_obj <- readRDS("merged_harmony_annotated_filtered.rds")
if (seurat_v5) seurat_obj <- JoinLayers(seurat_obj)

# Subset B cells
b_cells        <- subset(seurat_obj, subset = broad_celltype == "B cells")
if (seurat_v5) b_cells <- JoinLayers(b_cells)
b_cells$group  <- factor(b_cells$group, levels = c("Control","ND","MM","AMG"))

# Gate on human MYC raw counts
myc_raw        <- as.numeric(get_counts(b_cells)["MYC", ])
names(myc_raw) <- colnames(b_cells)

myeloma_cells  <- names(myc_raw)[myc_raw > 0]
normal_b_cells <- names(myc_raw)[myc_raw == 0]

cat("MYC+ myeloma PCs: ", length(myeloma_cells), "\n")
cat("MYC- normal B:    ", length(normal_b_cells), "\n\n")

# Subset objects
myeloma_pc       <- subset(b_cells, cells = myeloma_cells)
myeloma_pc$group <- factor(myeloma_pc$group, levels = myeloma_levels)
myeloma_pc$group <- droplevels(myeloma_pc$group)

normal_b         <- subset(b_cells, cells = normal_b_cells)
normal_b$group   <- factor(normal_b$group,
                           levels = c("Control","ND","MM","AMG"))

cat("MYC+ cells per group:\n")
print(table(myeloma_pc$group))

# ============================================================
# 2. CHECK Cd48
# ============================================================
cd48_gene <- "Cd48"
if (!cd48_gene %in% rownames(myeloma_pc)) {
  hits <- grep("^cd48$", rownames(myeloma_pc), ignore.case = TRUE, value = TRUE)
  if (length(hits) > 0) cd48_gene <- hits[1] else stop("Cd48 not found.")
}
cat("\nUsing gene:", cd48_gene, "\n")

# ============================================================
# 3. FETCH DATA
# ============================================================
df_myeloma <- FetchData(myeloma_pc, vars = c("group", cd48_gene))
colnames(df_myeloma) <- c("group", "Cd48")
df_myeloma$group <- factor(df_myeloma$group, levels = myeloma_levels)

df_normal <- FetchData(normal_b, vars = c("group", cd48_gene))
colnames(df_normal) <- c("group", "Cd48")
df_normal$group <- factor(df_normal$group,
                          levels = c("Control","ND","MM","AMG"))

cat("Cd48 summary (MYC+ PCs):\n")
print(summary(df_myeloma$Cd48))
cat("% expressing:", round(100 * mean(df_myeloma$Cd48 > 0), 1), "%\n\n")

# ============================================================
# PLOT 1 — Stacked bar: Negative / Low / High per group
# ============================================================
cd48_expressing <- df_myeloma$Cd48[df_myeloma$Cd48 > 0]
threshold_high  <- median(cd48_expressing)
cat("Cd48-high threshold (median of expressors):",
    round(threshold_high, 3), "\n")

df_myeloma <- df_myeloma %>%
  mutate(
    cd48_status = case_when(
      Cd48 == 0             ~ "Negative",
      Cd48 < threshold_high ~ "Low",
      TRUE                  ~ "High"
    ),
    cd48_status = factor(cd48_status, levels = c("Negative","Low","High"))
  )

cd48_composition <- df_myeloma %>%
  group_by(group, cd48_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_stacked <- ggplot(cd48_composition,
                    aes(x = group, y = pct, fill = cd48_status)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = ifelse(pct > 5, sprintf("%.0f%%", pct), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5, color = "white", fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("Negative" = "grey85", "Low" = "#FCBBA1", "High" = "#B2182B"),
    name   = "Cd48 status"
  ) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(
    title    = "Cd48 Expression Status — MYC+ Myeloma Plasma Cells",
    subtitle = paste0("High = above median of expressors (>",
                      round(threshold_high, 2), " normalised units)"),
    x = NULL, y = "% of cells"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title          = element_text(face = "bold", hjust = 0.5),
    plot.subtitle       = element_text(color = "grey40", size = 9, hjust = 0.5),
    legend.position     = "right",
    panel.grid.major.x  = element_blank()
  )

print(p_stacked)
# ggsave("PC_Cd48_stacked.pdf", p_stacked, width = 6, height = 5)

# ============================================================
# PLOT 2 — Violin: expressing cells only (removes zero inflation)
# ============================================================
df_expressing_only <- df_myeloma %>% filter(Cd48 > 0)

cat("\nExpressing cells per group:\n")
print(table(df_expressing_only$group))

p_violin_expr <- ggplot(df_expressing_only,
                        aes(x = group, y = Cd48, fill = group)) +
  geom_violin(width = 0.85, alpha = 0.5, color = NA, trim = TRUE) +
  geom_jitter(size = 0.3, width = 0.18, alpha = 0.25, color = "black") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.4, color = "black", linewidth = 0.5) +
  stat_compare_means(
    comparisons     = comparisons_myeloma,
    method          = "wilcox.test",
    p.adjust.method = "BH",
    label           = "p.signif",
    hide.ns         = TRUE,
    tip.length      = 0.01,
    bracket.size    = 0.3,
    size            = 3.5
  ) +
  scale_fill_manual(values = myeloma_colors) +
  scale_y_log10() +
  annotation_logticks(sides = "l", size = 0.3) +
  labs(
    title    = "Cd48 Expression Level — MYC+ Myeloma PCs",
    subtitle = "Zero-expressing cells excluded | BH-corrected Wilcoxon\n*** p<0.001  ** p<0.01  * p<0.05",
    x = NULL,
    y = "Normalised expression (log10)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(color = "grey40", size = 9, hjust = 0.5),
    legend.position = "none",
    axis.text.x     = element_text(size = 12)
  )

print(p_violin_expr)
# ggsave("PC_Cd48_violin_expressing.pdf", p_violin_expr, width = 5.5, height = 6)

# ============================================================
# PLOT 3 — % expressing: MYC+ vs MYC- side by side
# ============================================================
pct_compare <- bind_rows(
  df_myeloma %>% mutate(cell_type = "MYC+ Myeloma PC"),
  df_normal  %>% mutate(cell_type = "MYC- Normal B cell")
) %>%
  group_by(cell_type, group) %>%
  summarise(
    n_cells  = n(),
    pct_expr = 100 * mean(Cd48 > 0),
    .groups  = "drop"
  ) %>%
  mutate(cell_type = factor(cell_type,
                            levels = c("MYC+ Myeloma PC",
                                       "MYC- Normal B cell")))

p_pct_compare <- ggplot(pct_compare,
                        aes(x = group, y = pct_expr, fill = group)) +
  geom_bar(stat = "identity", width = 0.65,
           color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.0f%%\nn=%d", pct_expr, n_cells)),
    vjust = -0.4, size = 3, color = "grey30"
  ) +
  facet_wrap(~ cell_type, ncol = 2) +
  scale_fill_manual(values = all_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)),
                     labels = percent_format(scale = 1)) +
  labs(
    title    = "% Cells Expressing Cd48: MYC+ Myeloma PC vs MYC- Normal B Cells",
    subtitle = "Cells with normalised Cd48 > 0",
    x = NULL, y = "% Expressing Cd48"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title          = element_text(face = "bold", hjust = 0.5),
    plot.subtitle       = element_text(color = "grey40", size = 9, hjust = 0.5),
    legend.position     = "none",
    strip.text          = element_text(face = "bold", size = 11),
    panel.grid.major.x  = element_blank(),
    axis.text.x         = element_text(size = 11)
  )

print(p_pct_compare)
# ggsave("PC_Cd48_pct_compare.pdf", p_pct_compare, width = 10, height = 5)

# ============================================================
# PLOT 4 — Combined panel
# ============================================================
panel_final <- (p_stacked | p_violin_expr) /
  p_pct_compare +
  plot_annotation(
    title    = "Cd48 in vk*MYC Myeloma Plasma Cells",
    subtitle = "Top left: expression composition | Top right: level in expressors | Bottom: prevalence",
    theme    = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9,  color = "grey40", hjust = 0.5)
    )
  )

print(panel_final)
# ggsave("PC_Cd48_final_panel.pdf", panel_final, width = 12, height = 11)

# ============================================================
# STATISTICS
# ============================================================
cat("\n--- Chi-square: Cd48 expression prevalence across myeloma groups ---\n")
cont_table <- df_myeloma %>%
  mutate(expressing = Cd48 > 0) %>%
  group_by(group, expressing) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = expressing, values_from = n, values_fill = 0) %>%
  column_to_rownames("group") %>%
  as.matrix()

print(cont_table)
print(chisq.test(cont_table))

cat("\n--- Pairwise Wilcoxon BH-corrected (expressing cells only) ---\n")
pw <- pairwise.wilcox.test(
  x               = df_expressing_only$Cd48,
  g               = df_expressing_only$group,
  p.adjust.method = "BH"
)
print(pw$p.value)

cat("\n--- Summary per group (expressing cells only) ---\n")
df_expressing_only %>%
  group_by(group) %>%
  summarise(
    n_expressing = n(),
    mean_expr    = round(mean(Cd48),   4),
    median_expr  = round(median(Cd48), 4),
    sd_expr      = round(sd(Cd48),     4),
    .groups      = "drop"
  ) %>%
  print()

cat("\n--- Analysis complete ---\n")