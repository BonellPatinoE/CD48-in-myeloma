# ================================================
# AP-1 TRANSCRIPTION FACTOR ANALYSIS IN MURINE NK CELLS
# Comparing Control / ND / MM / AMG groups
#
# Human AP-1 paralogs → Murine orthologs:
#   FOS   → Fos
#   FOSB  → Fosb
#   FOSL1 → Fosl1
#   FOSL2 → Fosl2
#   JUN   → Jun
#   JUNB  → Junb
#   JUND  → Jund
#   ATF3  → Atf3   (stress-induced AP-1 family member)
# ================================================

# ================================================
# 0. SETUP
# ================================================
set.seed(1234)

library(Seurat)
library(AUCell)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tibble)

# ================================================
# 1. LOAD YOUR ANNOTATED SEURAT OBJECT AND SUBSET NK CELLS
# ================================================
seurat_obj <- readRDS("merged_harmony_annotated_filtered.rds")

# Confirm required columns exist
stopifnot("broad_celltype" %in% colnames(seurat_obj@meta.data))
stopifnot("group"          %in% colnames(seurat_obj@meta.data))

# Subset to NK cells only
nk <- subset(seurat_obj, subset = broad_celltype == "NK cells")

# Set group factor order
nk$group <- factor(nk$group, levels = c("Control", "ND", "MM", "AMG"))

cat("NK cells per group:\n")
print(table(nk$group))

# Join layers if using Seurat v5 (required after subsetting)
nk <- JoinLayers(nk)

# ================================================
# 2. DEFINE MURINE AP-1 GENE SET AND CHECK AVAILABILITY
# ================================================
ap1_candidates <- c(
  "Fos", "Fosb", "Fosl1", "Fosl2",   # Fos family
  "Jun", "Junb", "Jund",              # Jun family
  "Atf3"                              # ATF/CREB stress-related AP-1
)

# Only keep genes actually present in this dataset
ap1_genes <- intersect(ap1_candidates, rownames(nk))

missing <- setdiff(ap1_candidates, rownames(nk))

cat("\n--- AP-1 gene availability ---\n")
cat("Found     :", paste(ap1_genes, collapse = ", "), "\n")
cat("Not found :", ifelse(length(missing) == 0, "none", paste(missing, collapse = ", ")), "\n")

if (length(ap1_genes) == 0) {
  stop("No AP-1 genes found in the dataset. Check gene name capitalisation (murine = sentence case).")
}

# ================================================
# 3. PANEL A — VIOLIN PLOTS FOR EACH INDIVIDUAL AP-1 GENE
# ================================================
# Shared pairwise comparisons
comparisons <- list(
  c("Control", "ND"),
  c("Control", "MM"),
  c("Control", "AMG"),
  c("ND",      "MM"),
  c("ND",      "AMG"),
  c("MM",      "AMG")
)

# Color palette (one colour per group)
group_colors <- c(
  "Control" = "#66C2A5",
  "ND"      = "#FC8D62",
  "MM"      = "#8DA0CB",
  "AMG"     = "#E78AC3"
)

# Function: single-gene violin with Wilcoxon pairwise stats
plot_gene_violin <- function(nk_obj, gene, colors = group_colors) {
  df <- FetchData(nk_obj, vars = c("group", gene))
  colnames(df) <- c("group", "expr")
  df$group <- factor(df$group, levels = c("Control", "ND", "MM", "AMG"))
  
  ggplot(df, aes(x = group, y = expr + 0.001, fill = group)) +
    geom_violin(width = 0.8, alpha = 0.4, color = NA, trim = TRUE) +
    geom_jitter(size = 0.3, width = 0.18, alpha = 0.3, color = "black") +
    stat_summary(
      fun    = mean,
      geom   = "crossbar",
      width  = 0.4,
      color  = "black",
      size   = 0.6,
      fatten = 2
    ) +
    stat_compare_means(
      comparisons     = comparisons,
      method          = "wilcox.test",
      p.adjust.method = "BH",
      label           = "p.signif",   # *** / ns
      hide.ns         = TRUE,
      tip.length      = 0.01,
      bracket.size    = 0.3
    ) +
    scale_fill_manual(values = colors) +
    scale_y_log10() +
    labs(
      title = gene,
      x     = NULL,
      y     = "Normalised expression (log10)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title      = element_text(face = "bold.italic", hjust = 0.5),
      axis.text.x     = element_text(size = 10, angle = 30, hjust = 1)
    )
}

# Generate one plot per AP-1 gene
violin_list <- lapply(ap1_genes, function(g) plot_gene_violin(nk, g))
names(violin_list) <- ap1_genes

# Print individually
for (g in ap1_genes) {
  print(violin_list[[g]])
}

# Combine all into a single figure (up to 8 panels in a 4×2 grid)
n_genes    <- length(ap1_genes)
n_cols     <- min(4, n_genes)
n_rows     <- ceiling(n_genes / n_cols)

panel_A <- wrap_plots(violin_list, ncol = n_cols) +
  plot_annotation(
    title    = "AP-1 Transcription Factor Expression in Murine NK Cells",
    subtitle = "BH-corrected Wilcoxon; * p<0.05  ** p<0.01  *** p<0.001",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

print(panel_A)
# ggsave("NK_AP1_violins.pdf", panel_A, width = 4 * n_cols, height = 4.5 * n_rows)

# ================================================
# 4. PANEL B — AUCell COMPOSITE AP-1 ACTIVITY SCORE
# ================================================

# 4.1 Build rankings from normalised data
expr_sparse <- as(
  GetAssayData(nk, assay = "RNA", layer = "data"),
  "dgCMatrix"
)

# Prune to HVGs + AP-1 genes to reduce memory
hvgs            <- VariableFeatures(nk)
features_ranked <- unique(c(hvgs, ap1_genes))
expr_pruned     <- expr_sparse[intersect(features_ranked, rownames(expr_sparse)), ]

cat("\nBuilding AUCell rankings...\n")
rankings <- AUCell_buildRankings(expr_pruned, plotStats = FALSE)

# 4.2 Calculate AUC for the AP-1 gene set
ap1_AUC  <- AUCell_calcAUC(ap1_genes, rankings)
nk[["AP1_score"]] <- as.numeric(getAUC(ap1_AUC))

# 4.3 Violin plot for composite AP-1 score
df_auc <- FetchData(nk, vars = c("group", "AP1_score"))
df_auc$group <- factor(df_auc$group, levels = c("Control", "ND", "MM", "AMG"))

panel_B <- ggplot(df_auc, aes(x = group, y = AP1_score + 1e-6, fill = group)) +
  geom_violin(width = 0.8, alpha = 0.4, color = NA, trim = TRUE) +
  geom_jitter(size = 0.3, width = 0.18, alpha = 0.3, color = "black") +
  stat_summary(
    fun    = mean,
    geom   = "crossbar",
    width  = 0.4,
    color  = "black",
    size   = 0.6,
    fatten = 2
  ) +
  stat_compare_means(
    comparisons     = comparisons,
    method          = "wilcox.test",
    p.adjust.method = "BH",
    label           = "p.signif",
    hide.ns         = TRUE,
    tip.length      = 0.01,
    bracket.size    = 0.3
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_log10() +
  labs(
    title    = "Composite AP-1 Activity Score (AUCell)",
    subtitle = paste("Gene set:", paste(ap1_genes, collapse = ", ")),
    x        = NULL,
    y        = "AUC score (log10)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 9, color = "grey40", hjust = 0.5),
    axis.text.x     = element_text(size = 12)
  )

print(panel_B)
# ggsave("NK_AP1_AUCell_score.pdf", panel_B, width = 5, height = 6)

# ================================================
# 5. PANEL C — DOT PLOT: % expressing + mean expression
#    (equivalent to a "bubble plot" summary across groups)
# ================================================
Idents(nk) <- nk$group

panel_C <- DotPlot(
  nk,
  features = ap1_genes,
  group.by = "group",
  cols     = c("lightgrey", "firebrick"),
  dot.scale = 8
) +
  coord_flip() +
  labs(
    title = "AP-1 TF Expression — % Expressing & Mean Level",
    x     = NULL,
    y     = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")
  )

print(panel_C)
# ggsave("NK_AP1_dotplot.pdf", panel_C, width = 6, height = 5)

# ================================================
# 6. PANEL D — FEATURE PLOTS ON UMAP (per key gene)
#    Show spatial distribution of top AP-1 members
# ================================================
# Focus on the core three highlighted in the human volcano: Fos, Jun, Junb
core_ap1 <- intersect(c("Fos", "Jun", "Junb"), ap1_genes)

if (length(core_ap1) > 0) {
  panel_D <- FeaturePlot(
    nk,
    features  = core_ap1,
    reduction = "umap",
    pt.size   = 0.4,
    ncol      = length(core_ap1),
    order     = TRUE,           # plot expressing cells on top
    cols      = c("lightgrey", "darkblue")
  ) &
    theme_minimal(base_size = 11)
  
  print(panel_D)
  # ggsave("NK_AP1_featureplots.pdf", panel_D, width = 4 * length(core_ap1), height = 4)
}

# ================================================
# 7. STATISTICAL SUMMARY TABLE
#    Mean AUC ± SD per group + pairwise Wilcoxon p-values
# ================================================
# 7.1. Per-group summary of AP-1 AUC
summary_tbl <- df_auc %>%
  group_by(group) %>%
  summarise(
    n        = n(),
    mean_AUC = mean(AP1_score),
    sd_AUC   = sd(AP1_score),
    median_AUC = median(AP1_score),
    .groups  = "drop"
  )

cat("\n--- AP-1 AUC score summary ---\n")
print(summary_tbl)

# 7.2. Pairwise Wilcoxon tests (BH corrected)
groups_vec <- as.character(df_auc$group)
scores_vec <- df_auc$AP1_score

pairwise_results <- pairwise.wilcox.test(
  x         = scores_vec,
  g         = groups_vec,
  p.adjust.method = "BH"
)

cat("\n--- Pairwise Wilcoxon (BH-corrected) p-values ---\n")
print(pairwise_results$p.value)

# 7.3. Per-gene mean expression per group
gene_summary <- FetchData(nk, vars = c("group", ap1_genes)) %>%
  group_by(group) %>%
  summarise(across(all_of(ap1_genes), mean), .groups = "drop")

cat("\n--- Mean normalised expression per group (AP-1 genes) ---\n")
print(gene_summary)

# ================================================
# 8. OPTIONAL: SAVE NK OBJECT WITH AP-1 SCORE
# ================================================
# saveRDS(nk, file = "nk_only_AP1_scored.rds")

cat("\n--- Analysis complete ---\n")
cat("Plots generated: Panel A (per-gene violins), Panel B (AUCell score),\n")
cat("                 Panel C (dot plot), Panel D (UMAP feature plots)\n")

# ================================================
# 8. OPTIONAL: SAVE NK OBJECT WITH AP-1 SCORE
# ================================================
# saveRDS(nk, file = "nk_only_AP1_scored.rds")

cat("\n--- Analysis complete ---\n")
cat("Plots generated: Panel A (per-gene violins), Panel B (AUCell score),\n")
cat("                 Panel C (dot plot), Panel D (UMAP feature plots)\n")

# ================================================
# 9. NK SUBCLUSTERING + RECEPTOR ANALYSIS
# Addresses Reviewer Comment 24
# ================================================

# Step 1: Reprocess NK cells for subclustering
nk <- NormalizeData(nk)
nk <- FindVariableFeatures(nk, nfeatures = 2000)
nk <- ScaleData(nk)
nk <- RunPCA(nk, npcs = 30)
nk <- FindNeighbors(nk, dims = 1:20)
nk <- FindClusters(nk, resolution = 0.4,
                   cluster.name = "nk_clusters")
nk <- RunUMAP(nk, dims = 1:20,
              reduction.name = "umap_nk")

cat("NK subclusters:\n")
print(table(nk$nk_clusters))

cat("\nSubcluster composition by group:\n")
comp <- table(nk$nk_clusters, nk$group)
print(comp)

# Step 2: UMAP colored by cluster and by group
p_umap_cluster <- DimPlot(
  nk,
  reduction = "umap_nk",
  group.by  = "nk_clusters",
  label     = TRUE,
  pt.size   = 0.5
) + ggtitle("Murine NK Subclusters")

p_umap_group <- DimPlot(
  nk,
  reduction = "umap_nk",
  group.by  = "group",
  pt.size   = 0.5,
  cols      = group_colors
) + ggtitle("NK Subclusters by Disease Stage")

print(p_umap_cluster + p_umap_group)

ggsave("NK_UMAP_subclusters_murine.pdf",
       p_umap_cluster + p_umap_group,
       width = 12, height = 5, dpi = 300)

# Split UMAP by group
p_split <- DimPlot(
  nk,
  reduction = "umap_nk",
  group.by  = "nk_clusters",
  split.by  = "group",
  label     = TRUE,
  pt.size   = 0.5,
  ncol      = 2
) + ggtitle("NK Subclusters Split by Disease Stage")

print(p_split)
ggsave("NK_UMAP_split_murine.pdf",
       p_split,
       width = 14, height = 10, dpi = 300)

# Step 3: AP-1 DotPlot by subcluster
Idents(nk) <- nk$nk_clusters

p_ap1_subcluster <- DotPlot(
  nk,
  features = ap1_genes,
  group.by = "nk_clusters",
  cols     = c("lightgrey", "firebrick"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "AP-1 Expression by NK Subcluster",
       x = NULL, y = "NK Subcluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, face = "italic")
  )

print(p_ap1_subcluster)
ggsave("NK_AP1_by_subcluster_murine.pdf",
       p_ap1_subcluster,
       width = 8, height = 6, dpi = 300)

# Step 4: Inhibitory and activating receptor genes by subcluster
activating_genes <- intersect(
  c("Klrk1",    # NKG2D
    "Ncr1",     # NKp46
    "Cd226",    # DNAM-1
    "Klrc2",    # NKG2C
    "Fcgr3"),   # CD16
  rownames(nk)
)

inhibitory_genes <- intersect(
  c("Klrc1",    # NKG2A
    "Klrg1",    # KLRG1
    "Lag3",     # LAG-3
    "Tigit",    # TIGIT
    "Pdcd1",    # PD-1
    "Havcr2"),  # TIM-3
  rownames(nk)
)

cat("Activating receptors found:", paste(activating_genes, collapse=", "), "\n")
cat("Inhibitory receptors found:", paste(inhibitory_genes, collapse=", "), "\n")

p_activating <- DotPlot(
  nk,
  features  = activating_genes,
  group.by  = "nk_clusters",
  cols      = c("lightgrey", "steelblue"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "Activating Receptors by NK Subcluster",
       x = NULL, y = "NK Subcluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

p_inhibitory <- DotPlot(
  nk,
  features  = inhibitory_genes,
  group.by  = "nk_clusters",
  cols      = c("lightgrey", "darkred"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "Inhibitory Receptors by NK Subcluster",
       x = NULL, y = "NK Subcluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_activating)
print(p_inhibitory)

ggsave("NK_activating_receptors_subcluster.pdf",
       p_activating, width = 8, height = 5, dpi = 300)
ggsave("NK_inhibitory_receptors_subcluster.pdf",
       p_inhibitory, width = 8, height = 5, dpi = 300)

# Step 5: Cluster composition bar chart
comp_df      <- as.data.frame.matrix(comp)
comp_df$Total <- rowSums(comp_df)
for(grp in c("Control","ND","MM","AMG")) {
  if(grp %in% colnames(comp_df)) {
    comp_df[[paste0(grp,"_pct")]] <-
      round(comp_df[[grp]] / comp_df$Total * 100, 1)
  }
}
cat("\nCluster composition (%):\n")
print(comp_df)

library(reshape2)
comp_long <- melt(
  as.matrix(comp),
  varnames   = c("Cluster","Group"),
  value.name = "Cells"
)
comp_long$Group <- factor(comp_long$Group,
                          levels = c("Control","ND","MM","AMG"))

p_comp <- ggplot(comp_long,
                 aes(x = factor(Cluster), y = Cells, fill = Group)) +
  geom_bar(stat = "identity", position = "fill",
           alpha = 0.8, color = "white") +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "NK Subcluster Composition by Disease Stage",
       x = "NK Subcluster", y = "% of Cluster") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_comp)
ggsave("NK_cluster_composition_murine.pdf",
       p_comp, width = 10, height = 5, dpi = 300)

# Save NK object with subclusters
saveRDS(nk, "nk_subclustered_AP1_scored.rds")
cat("Saved: nk_subclustered_AP1_scored.rds\n")
cat("Figures saved: UMAP, AP-1 dotplot, receptor dotplots, composition bar chart\n")

# AP-1 split by group — add this to your code
nk$cluster_group <- paste0(
  as.character(nk$nk_clusters), "_",
  as.character(nk$group)
)

# AP-1 DotPlot split by group
p_ap1_split <- DotPlot(
  nk,
  features = ap1_genes,
  group.by = "cluster_group",
  cols     = c("lightgrey", "firebrick"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "AP-1 Expression by NK Subcluster — split by Disease Stage",
       x = NULL, y = "Subcluster_Group") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )

print(p_ap1_split)
ggsave("NK_AP1_subcluster_split_group.pdf",
       p_ap1_split, width = 14, height = 6, dpi = 300)

# Activating receptors split by group
p_act_split <- DotPlot(
  nk,
  features = activating_genes,
  group.by = "cluster_group",
  cols     = c("lightgrey", "steelblue"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "Activating Receptors by NK Subcluster — split by Disease Stage",
       x = NULL, y = "Subcluster_Group") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )

print(p_act_split)
ggsave("NK_activating_subcluster_split_group.pdf",
       p_act_split, width = 14, height = 5, dpi = 300)

# Inhibitory receptors split by group
p_inh_split <- DotPlot(
  nk,
  features = inhibitory_genes,
  group.by = "cluster_group",
  cols     = c("lightgrey", "darkred"),
  dot.scale = 8
) +
  coord_flip() +
  labs(title = "Inhibitory Receptors by NK Subcluster — split by Disease Stage",
       x = NULL, y = "Subcluster_Group") +
  theme_minimal(base_size = 10) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  )

print(p_inh_split)
ggsave("NK_inhibitory_subcluster_split_group.pdf",
       p_inh_split, width = 14, height = 5, dpi = 300)


##################
##################

# ============================================================
# STACKED HORIZONTAL DOTPLOTS — AP-1, Activating, Inhibitory
# Clean version — run entirely from scratch
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── 1. COLORS ────────────────────────────────────────────────
group_colors <- c(
  "Control" = "#66C2A5",
  "ND"      = "#FC8D62",
  "MM"      = "#8DA0CB",
  "AMG"     = "#E78AC3"
)

# ── 2. GENE SETS ─────────────────────────────────────────────
ap1_genes <- intersect(
  c("Fos","Fosb","Fosl1","Fosl2","Jun","Junb","Jund","Atf3"),
  rownames(nk)
)

activating_all <- intersect(
  c("Klrc2","Cd160",                    # HLA-dependent
    "Ncr1","Klrk1","Crtam","Fcgr3"),   # HLA-independent
  rownames(nk)
)

inhibitory_all <- intersect(
  c("Klra1","Klra3","Klra6","Klra7","Klra8",
    "Klra9","Klra10","Klra12","Pirb","Lag3",  # HLA-dependent
    "Pdcd1","Siglecg","Cd300a","Cd96",         # HLA-independent
    grep("^Il1rapl", rownames(nk), value=TRUE),
    "Tigit","Havcr2"),
  rownames(nk)
)

cat("AP-1 genes:      ", paste(ap1_genes,       collapse=", "), "\n")
cat("Activating genes:", paste(activating_all,  collapse=", "), "\n")
cat("Inhibitory genes:", paste(inhibitory_all,  collapse=", "), "\n")

# ── 3. ORDERED LABELS ────────────────────────────────────────
group_levels   <- c("Control","ND","MM","AMG")
cluster_levels <- as.character(
  sort(unique(as.numeric(as.character(nk$nk_clusters))))
)

ordered_labels <- c()
for (cl in cluster_levels) {
  for (gr in group_levels) {
    ordered_labels <- c(ordered_labels, paste0(cl,"_",gr))
  }
}

# Add cluster_group to metadata
nk$cluster_group <- paste0(
  as.character(nk$nk_clusters), "_",
  as.character(nk$group)
)

# Keep only existing combinations
existing <- unique(as.character(nk$cluster_group))
ordered_labels_existing <- ordered_labels[ordered_labels %in% existing]

cat("\nOrdered labels:\n")
print(ordered_labels_existing)

# ── 4. DATA BUILDER FUNCTION ──────────────────────────────────
build_dotplot_data <- function(nk_obj, genes, ordered_labels) {
  results  <- list()
  all_expr <- FetchData(nk_obj,
                        vars = c("cluster_group", genes))
  
  for (label in ordered_labels) {
    cells_mask <- all_expr$cluster_group == label
    if (sum(cells_mask) == 0) next
    cell_data <- all_expr[cells_mask, , drop = FALSE]
    parts     <- strsplit(label, "_")[[1]]
    cluster   <- parts[1]
    group     <- parts[2]
    
    for (gene in genes) {
      if (!gene %in% colnames(cell_data)) next
      expr <- cell_data[[gene]]
      results[[length(results) + 1]] <- data.frame(
        label     = label,
        cluster   = cluster,
        group     = group,
        gene      = gene,
        pct_expr  = mean(expr > 0) * 100,
        mean_expr = mean(expr),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results)
}

# ── 5. BUILD DATA ─────────────────────────────────────────────
ap1_data     <- build_dotplot_data(nk, ap1_genes,
                                   ordered_labels_existing)
act_all_data <- build_dotplot_data(nk, activating_all,
                                   ordered_labels_existing)
inh_all_data <- build_dotplot_data(nk, inhibitory_all,
                                   ordered_labels_existing)

cat("\nAP-1 data:      ", nrow(ap1_data),     "rows\n")
cat("Activating data:", nrow(act_all_data), "rows\n")
cat("Inhibitory data:", nrow(inh_all_data), "rows\n")

# ── 6. PLOT FUNCTION ──────────────────────────────────────────
custom_dotplot_pub <- function(dot_data,
                               title,
                               group_colors,
                               ordered_labels,
                               show_x_labels = FALSE) {
  
  dot_data$label <- factor(dot_data$label,
                           levels = ordered_labels)
  dot_data$group <- factor(dot_data$group,
                           levels = names(group_colors))
  
  # Expression intensity → alpha
  max_expr           <- max(abs(dot_data$mean_expr), na.rm=TRUE)
  dot_data$alpha_val <- 0.15 + 0.85 *
    (pmax(dot_data$mean_expr, 0) / max_expr)
  
  p <- ggplot(dot_data,
              aes(x     = label,
                  y     = gene,
                  size  = pct_expr,
                  color = group,
                  alpha = alpha_val)) +
    geom_point() +
    scale_color_manual(values = group_colors,
                       name   = "Disease Stage") +
    scale_size_continuous(range  = c(0.5, 7),
                          name   = "% Expressing",
                          breaks = c(25, 50, 75, 100)) +
    scale_alpha_continuous(range = c(0.15, 1),
                           guide = "none") +
    scale_x_discrete(limits = ordered_labels) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title         = element_text(face  = "bold",
                                        hjust = 0.5,
                                        size  = 14),
      axis.text.y        = element_text(face  = "italic",
                                        size  = 13),
      axis.title         = element_blank(),
      panel.grid.major   = element_line(color = "grey90"),
      legend.position    = "right",
      legend.text        = element_text(size  = 11),
      legend.title       = element_text(size  = 12,
                                        face  = "bold")
    ) +
    labs(title = title)
  
  # X-axis labels only on bottom panel
  if (show_x_labels) {
    p <- p + theme(
      axis.text.x = element_text(angle = 90,
                                 hjust = 1,
                                 vjust = 0.5,
                                 size  = 11,
                                 face  = "bold")
    )
  } else {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  # Cluster separator vertical lines
  n_groups      <- length(group_colors)
  sep_positions <- seq(n_groups + 0.5,
                       length(ordered_labels) - 0.5,
                       by = n_groups)
  for (pos in sep_positions) {
    p <- p + geom_vline(xintercept = pos,
                        color     = "grey50",
                        linetype  = "dashed",
                        linewidth = 0.4)
  }
  
  return(p)
}

# ── 7. GENERATE THREE PANELS ──────────────────────────────────

# AP-1 — top panel
p_ap1_pub <- custom_dotplot_pub(
  ap1_data,
  title          = "AP-1 Family Genes",
  group_colors   = group_colors,
  ordered_labels = ordered_labels_existing,
  show_x_labels  = FALSE
)

# Activating — middle panel
p_act_pub <- custom_dotplot_pub(
  act_all_data,
  title          = "Activating NK Receptors",
  group_colors   = group_colors,
  ordered_labels = ordered_labels_existing,
  show_x_labels  = FALSE
)

# Inhibitory — bottom panel (shows x-axis labels)
p_inh_pub <- custom_dotplot_pub(
  inh_all_data,
  title          = "Inhibitory NK Receptors",
  group_colors   = group_colors,
  ordered_labels = ordered_labels_existing,
  show_x_labels  = TRUE
)

# ── 8. STACK WITH PATCHWORK ───────────────────────────────────
n_ap1 <- length(unique(ap1_data$gene))
n_act <- length(unique(act_all_data$gene))
n_inh <- length(unique(inh_all_data$gene))

cat("\nGenes per panel — AP-1:", n_ap1,
    "| Activating:", n_act,
    "| Inhibitory:", n_inh, "\n")

combined <- p_ap1_pub /
  p_act_pub /
  p_inh_pub +
  plot_layout(
    heights = c(n_ap1, n_act, n_inh),
    guides  = "collect"
  ) &
  theme(legend.position = "right")

# ── 9. SAVE ───────────────────────────────────────────────────
ggsave("NK_dotplots_stacked_final.pdf",
       combined,
       width  = 22,
       height = 18,
       dpi    = 300)

ggsave("NK_dotplots_stacked_final.png",
       combined,
       width  = 22,
       height = 18,
       dpi    = 300)

print(combined)
cat("Saved: NK_dotplots_stacked_final (PDF + PNG)\n")