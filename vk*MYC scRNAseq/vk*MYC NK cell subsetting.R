# ==============================================
# 0. SET UP: Seed & load libraries
# ==============================================
set.seed(1234)

# Seurat for single‐cell processing
library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================
# 1. LOAD YOUR ANNOTATED + FILTERED SEURAT OBJECT
# ==============================================
# (This should be the object you saved at the end of your previous pipeline, 
#   containing broad_celltype, SingleR_cluster, group, etc.)
merged.filtered <- readRDS("merged_harmony_annotated_filtered.rds")

# Quick check to confirm “broad_celltype” and “group” exist:
stopifnot("broad_celltype" %in% colnames(merged.filtered@meta.data))
stopifnot("group" %in% colnames(merged.filtered@meta.data))

table(merged.filtered$broad_celltype)
table(merged.filtered$group)

# ==============================================
# 2. SUBSET TO NK CELLS ONLY
# ==============================================
# We assume you collapsed “NK/ILC/NKT” into “NK cells” (broad_celltype).
# If you used a slightly different name, adjust accordingly.

# 2.1. Identify which metadata column holds your NK label:
#      In our saved object, it’s “broad_celltype” == "NK cells"
nk_cells <- WhichCells(merged.filtered, expression = broad_celltype == "NK cells")

# 2.2. Create a new Seurat object for NK cells only
nk <- subset(merged.filtered, cells = nk_cells)

# 2.3. Verify
cat("Number of NK cells:", length(nk_cells), "\n")
table(nk$group)  # how many NK cells per sample group

# ==============================================
# 3. DIFFERENTIAL EXPRESSION: ALL‐GROUP TEST
#     (Find markers across the four ‘group’ levels)
# ==============================================
# Set the active identity to your “group” factor
Idents(nk) <- nk$group

# 3.1. Run FindAllMarkers across the four sample groups for NK cells
#      This will give you, for each pairwise comparison among Control, ND, MM, AMG,
#      a table of DE genes. Note that by default Seurat does “each group vs all others.”
markers_all_groups <- FindAllMarkers(
  object          = nk,
  assay           = "RNA",         # use normalized data slot
  slot            = "data",        # log‐normalized values
  only.pos        = FALSE,         # get both up‐ and downregulated
  min.pct         = 0.10,          # genes expressed in ≥10% of cells in either group
  logfc.threshold = 0.25           # absolute log2FC ≥ 0.25
)

# 3.2. Inspect top markers for each group
top5_per_group <- markers_all_groups %>% 
  group_by(cluster) %>%                    # here, “cluster” actually holds the group name
  slice_max(n = 5, order_by = avg_log2FC) %>% 
  select(cluster, gene, avg_log2FC, p_val_adj)

print(top5_per_group)

# (Optional) Save this DE table to disk:
# saveRDS(markers_all_groups, file = "NK_all_groups_markers.rds")
# write.csv(top5_per_group, file = "NK_top5_per_group.csv", row.names = FALSE)

# ==============================================
# 4. PAIRWISE DE: CONTROL VS MM WITHIN NK CELLS
# ==============================================
# We want to compare NK cells from the “Control” group vs NK cells from the “MM” group:
# Idents(nk) is already “group”

# 4.1. Run FindMarkers for ident.1 = "Control", ident.2 = "MM"
nk_de_control_vs_mm <- FindMarkers(
  object          = nk,
  assay           = "RNA",
  slot            = "data",
  ident.1         = "Control",       # first group
  ident.2         = "MM",            # second group
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  test.use        = "wilcox"         # default: Wilcoxon rank sum test
)

# 4.2. Add gene names as a column (since FindMarkers returns a data.frame with rownames = genes)
nk_de_control_vs_mm <- nk_de_control_vs_mm %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(p_val_adj)

# 4.3. Inspect top DE genes (by adjusted p‐value)
head(nk_de_control_vs_mm, n = 10)

# (Optional) Save this pairwise DE table
# saveRDS(nk_de_control_vs_mm, file = "NK_DE_Control_vs_MM.rds")
# write.csv(nk_de_control_vs_mm, file = "NK_DE_Control_vs_MM.csv", row.names = FALSE)

# ==============================================
# 5. VOLCANO PLOT: CONTROL VS MM
# ==============================================
# We’ll plot log2FC (avg_log2FC) on the x‐axis vs –log10(adjusted p‐value) on the y‐axis,
# and color by significance (e.g. p_val_adj < 0.05 & |log2FC| > 0.25).

# 5.1. Prepare the DE table for plotting
# Add a column “significance” to classify each gene
volcano_df <- nk_de_control_vs_mm %>%
  mutate(
    log10_padj = -log10(p_val_adj + 1e-300),   # avoid -Inf by adding a tiny constant
    sig = case_when(
      p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25 ~ "padj<0.05 & |log2FC|≥0.25",
      TRUE                                       ~ "Not Significant"
    )
  )

# 5.2. Basic volcano plot with ggplot2
volcano_plot <- ggplot(volcano_df, aes(x = avg_log2FC, y = log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c(
      "padj<0.05 & |log2FC|≥0.25" = "firebrick",
      "Not Significant"           = "grey70"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "NK cells: Control vs MM",
    x     = "Average log2 Fold Change (Control vs MM)",
    y     = "-log10(adj. p‐value)",
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey")

# 5.3. Optionally, label the top hits (e.g. top 10 most significant)
top_genes_to_label <- volcano_df %>%
  filter(sig == "padj<0.05 & abs(avg_log2FC) >= 0.25") %>%
  slice_min(order_by = p_val_adj, n = 10) %>%
  pull(gene)

volcano_plot + 
  geom_text(
    data = subset(volcano_df, gene %in% top_genes_to_label),
    aes(label = gene),
    vjust = -0.5, 
    size = 3
  )

VlnPlot(
  nk,
  features = "tnf",
  group.by = "group",
  pt.size  = 0
) + labs(title = "hMYC expression across sample groups")

# 1.1. Within your NK object, plot B‐cell markers vs MYC
FeaturePlot(
  nk, 
  features  = c("Cd19","Ms4a1","MYC"), 
  reduction = "umap", 
  pt.size   = 0.5,
  combine   = TRUE
)


# If you don’t need labels, just print the base volcano_plot:
# print(volcano_plot)

# ==============================================
# 6. (OPTIONAL) SAVE PLOTS AND RESULTS
# ==============================================
# ggsave("NK_volcano_Control_vs_MM.pdf", volcano_plot, width = 6, height = 5)
# ggsave("NK_cell_level_by_group.pdf",         p_cell,   width = 8, height = 6)
# ggsave("NK_cluster_level_comparison.pdf",    p_cluster, width = 8, height = 6)

# A reminder: you can also save R objects:
# saveRDS(nk_de_control_vs_mm, file = "NK_DE_Control_vs_MM.rds")
# saveRDS(volcano_df,            file = "NK_Volcano_Data.rds")

# ==============================================
# END OF SCRIPT
# ==============================================

# 1. Start from your NK subset (assume you've already created `nk`)
#    If not yet, do:
#    nk <- subset(merged.filtered, subset = broad_celltype == "NK cells")

# 2. Pull raw UMI counts for MYC from the “counts” slot
myc_raw <- GetAssayData(nk, assay = "RNA", slot = "counts")["MYC", ]

# 3. Identify cells that have NO MYC UMI (i.e. myc_raw == 0)
cells_no_myc <- names(myc_raw)[ myc_raw == 0 ]

length(cells_no_myc) 
# e.g. number of NK cells with zero MYC—for example, 4800 out of 4929

# 4. Subset `nk` to keep only cells where MYC_raw = 0
nk_no_myc <- subset(nk, cells = cells_no_myc)

cat("Original NK count:", ncol(nk), "\n")
cat("After removing any MYC+ cells:", ncol(nk_no_myc), "\n")

# 5. Verify by plotting MYC on the “cleaned” object; nothing should show up
FeaturePlot(
  nk_no_myc,
  features  = "MYC",
  reduction = "umap",
  pt.size   = 0.5
) + labs(title = "MYC expression (all cells should be zero)")

# 6. Now rerun Control vs MM differential expression on nk_no_myc
Idents(nk_no_myc) <- nk_no_myc$group

nk_de_control_vs_mm_clean <- FindMarkers(
  object          = nk_no_myc,
  assay           = "RNA",
  slot            = "data",
  ident.1         = "Control",
  ident.2         = "MM",
  min.pct         = 0.10,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
) %>%
  rownames_to_column(var = "gene") %>%
  arrange(p_val_adj)

# 7. Inspect top DE genes to confirm MYC is gone
head(nk_de_control_vs_mm_clean, n = 10)

# 8. Optional: rebuild volcano plot with the cleaned DE results
volcano_df <- nk_de_control_vs_mm_clean %>%
  mutate(
    log10_padj = -log10(p_val_adj + 1e-300),
    sig = case_when(
      p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25 ~ "padj<0.05 & |log2FC|≥0.25",
      TRUE                                       ~ "Not Significant"
    )
  )

ggplot(volcano_df, aes(x = avg_log2FC, y = log10_padj)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 1) +
  scale_color_manual(
    values = c(
      "padj<0.05 & |log2FC|≥0.25" = "firebrick",
      "Not Significant"           = "grey70"
    )
  ) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05),    linetype = "dashed", color = "darkgrey") +
  theme_minimal(base_size = 14) +
  labs(
    title = "NK cells (filtered to MYC‐zero): Control vs MM",
    x     = "Average log2FC (Control / MM)",
    y     = "-log10(adj. p‐value)"
  )

VlnPlot(
  nk,
  features = "Tigit",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Tigit expression across sample groups")

VlnPlot(
  nk,
  features = "Lag3",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Lag3 expression across sample groups")

VlnPlot(
  nk,
  features = "Klrg1",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Klrg1 expression across sample groups")

VlnPlot(
  nk,
  features = "Klrc1",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Klrc1 expression across sample groups")

VlnPlot(
  nk,
  features = "Cish",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Cish expression across sample groups")

VlnPlot(
  nk,
  features = "Sh2d1a",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Sh2d1a expression across sample groups")

VlnPlot(
  nk,
  features = "Inpp5d",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Inpp5d expression across sample groups")

VlnPlot(
  nk,
  features = "Inpp5d",
  group.by = "group",
  pt.size  = 1
) + labs(title = "Inpp5d expression across sample groups")

VlnPlot(
  nk,
  features = "Cd244",
  group.by = "group",
  pt.size  = 1
) + labs(title = "CD244 expression across sample groups")

VlnPlot(
  nk,
  features = "Klrc2",
  group.by = "group",
  pt.size  = 1
) + labs(title = "CD244 expression across sample groups")

library(Seurat)
library(ggpubr)
library(dplyr)

# Extract expression data and group labels
expr_df <- FetchData(nk, vars = c("Klrc2", "group"))
colnames(expr_df) <- c("Expression", "Group")

# Basic violin plot with p-values (ANOVA + pairwise Wilcoxon)
ggviolin(expr_df, x = "Group", y = "Expression", fill = "Group",
         add = "jitter", palette = "jco", trim = TRUE) +
  stat_compare_means(method = "anova", label.y = max(expr_df$Expression) + 0.5) +  # ANOVA p-value
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(
                       c("AMG", "Control"),
                       c("AMG", "MM"),
                       c("AMG", "ND"),
                       c("Control", "MM"),
                       c("Control", "ND"),
                       c("MM", "ND")
                     ),
                     label = "p.signif") +  # or use "p.format" to show actual p-values
  labs(title = "Klrc2 (CD244) expression across groups")


VlnPlot(
  nk,
  features = "Klrc3",
  group.by = "group",
  pt.size  = 1
) + labs(title = "CD244 expression across sample groups")

VlnPlot(
  nk,
  features = "Lag3",
  group.by = "group",
  pt.size  = 1
) + labs(title = "CD244 expression across sample groups")


