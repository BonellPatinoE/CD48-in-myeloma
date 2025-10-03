# ====================================================================================
# 0. SETUP: seed for reproducibility, load libraries
# ====================================================================================
set.seed(1234)

library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)
library(ggplot2)
library(patchwork)

# ====================================================================================
# 1. LOAD THE PREVIOUSLY SAVED HARMONY‐INTEGRATED SEURAT OBJECT
# ====================================================================================
merged <- readRDS("merged_harmony_markers.rds")

# Quick check: how many cells per original Seurat cluster?
cat("Cells per cluster (pre‐filter):\n")
print(table(merged$seurat_clusters))

# ====================================================================================
# 2. RUN SINGLE R FOR MURINE IMMUNE REFERENCE (IMMGEN)
#    2.1. Prepare reference and “test” data
#    2.2. Cell‐level SingleR annotation
#    2.3. Cluster‐level SingleR annotation (then map back to cells)
# ====================================================================================
# 2.1. Load ImmGen reference and extract normalized (“data”) matrix from merged object
immgen.ref <- ImmGenData()
testdata   <- GetAssayData(merged, assay = "RNA", slot = "data")

# 2.2. Cell‐level SingleR
singleR.cell <- SingleR(
  test            = testdata,
  ref             = immgen.ref,
  labels          = immgen.ref$label.main,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)
merged$SingleR_celltype <- singleR.cell$labels

# 2.3. Cluster‐level SingleR
singleR.cluster <- SingleR(
  test            = testdata,
  ref             = immgen.ref,
  labels          = immgen.ref$label.main,
  clusters        = merged$seurat_clusters,
  assay.type.test = "logcounts",
  assay.type.ref  = "logcounts"
)

# Extract cluster‐level labels as a named vector (names = cluster IDs)
cl.labels <- as.character(singleR.cluster$labels)
names(cl.labels) <- rownames(singleR.cluster)

# Map cluster‐level label onto each cell, store in metadata
cell.labels <- cl.labels[ as.character(merged$seurat_clusters) ]
names(cell.labels) <- colnames(merged)

merged <- AddMetaData(
  object   = merged,
  metadata = cell.labels,
  col.name = "SingleR_cluster"
)

# Verify mapping:
cat("Mapping of clusters to SingleR labels:\n")
print(table(merged$SingleR_cluster, merged$seurat_clusters))

# ====================================================================================
# 3. VISUALIZE SINGLER ANNOTATIONS ON UMAP
#    3.1. Cell‐level labels, faceted by sample group (“group” metadata should already exist)
#    3.2. Side‐by‐side comparison: cell‐level vs. cluster‐level
# ====================================================================================
# 3.1. Cell‐level SingleR, faceted by “group” (e.g. Control/ND/MM/AMG)

# Reorder the levels of the 'group' factor
desired_order <- c("Control", "ND", "MM", "AMG")
merged$group <- factor(merged$group, levels = desired_order)

DimPlot(
  merged,
  reduction = "umap",
  group.by  = "SingleR_celltype",
  split.by  = "group",
  ncol      = 2,
  pt.size   = 0.5,
  label     = FALSE
) + labs(title = "SingleR cell‐level labels, faceted by sample group")

# 3.2. Side‐by‐side: cell‐ vs. cluster‐level
p_cell    <- DimPlot(merged, group.by = "SingleR_celltype", label = TRUE) + ggtitle("Cell‐level")
p_cluster <- DimPlot(merged, group.by = "SingleR_cluster", label = TRUE) + ggtitle("Cluster‐level")

p_cell + p_cluster

# ====================================================================================
# 4. EXPLORE HUMAN MYC TRANSCRIPT EXPRESSION ON UMAP
#    4.1. Confirm feature name
#    4.2. Plot raw UMAP (all cells) colored by “MYC”
#    4.3. Plot faceted by sample group
# ====================================================================================
# 4.1. Find exact “MYC” gene name(s) in rownames
myc_features <- grep("MYC", rownames(merged), ignore.case = TRUE, value = TRUE)
print(myc_features)
# We expect “MYC” (human transgene) and “Myc” (mouse ortholog)

# 4.2. FeaturePlot for human MYC on UMAP
FeaturePlot(
  merged,
  features  = "MYC",
  reduction = "umap",
  pt.size   = 0.5
) + labs(title = "Expression of human MYC transgene")

# 4.3. Facet “MYC” expression by sample group in a 2×2 layout
FeaturePlot(
  merged,
  features  = "MYC",
  reduction = "umap",
  split.by  = "group",
  ncol      = 2,
  pt.size   = 0.5
) + labs(title = "hMYC expression by sample group")

# (Optional alternative using patchwork)
plots_myc <- FeaturePlot(
  merged,
  features  = "MYC",
  reduction = "umap",
  split.by  = "group",
  ncol      = 1,
  combine   = FALSE,
  pt.size   = 0.5
)
wrap_plots(plots_myc, ncol = 2) +
  plot_annotation(title = "hMYC expression by sample group")

# 4.4. Violin plot of “MYC” across sample groups
VlnPlot(
  merged,
  features = "MYC",
  group.by = "group",
  pt.size  = 0
) + labs(title = "hMYC expression across sample groups")

# 4.5. Compare mouse Myc vs human MYC side by side (if desired)
FeaturePlot(
  merged,
  features  = c("Myc","MYC"),
  reduction = "umap",
  pt.size   = 0.5,
  ncol      = 2
) + labs(title = "Mouse Myc (left) vs. Human MYC (right)")

DotPlot(merged, features = "Myc", "MYC", "Jun", "Fos")

# ====================================================================================
# 5. CHECK AND COLLAPSE SINGLER LABELS INTO BROAD CATEGORIES
#    (e.g. B cells, T cells, NK cells, etc.)
# ====================================================================================
# 5.1. Inspect unique calls from SingleR
cat("Unique SingleR cell‐level labels:\n")
print(unique(merged$SingleR_celltype))

# 5.2. Build “broad_celltype” metadata, merging as follows:
#   • “B cells” + “B cells, pro” → “B cells”
#   • “Basophils” + “Mast cells” → “Basophils”
#   • “T cells” + “Tgd” + “NKT” → “T cells”
#   • “Macrophages” + “Microglia” → “Macrophages”
#   • “NK cells” + “ILC” → “NK cells”   (NKT already merged above)
#   • All others remain unchanged
merged$broad_celltype <- case_when(
  merged$SingleR_celltype %in% c("B cells",     "B cells, pro") ~ "B cells",
  merged$SingleR_celltype %in% c("Basophils",   "Mast cells") ~ "Basophils",
  merged$SingleR_celltype %in% c("T cells",     "Tgd",   "NKT") ~ "T cells",
  merged$SingleR_celltype %in% c("Macrophages", "Microglia") ~ "Macrophages",
  merged$SingleR_celltype %in% c("NK cells",    "ILC")       ~ "NK cells",
  merged$SingleR_celltype == "Stem cells"      ~ "Stem cells",
  merged$SingleR_celltype == "Fibroblasts"     ~ "Fibroblasts",
  merged$SingleR_celltype == "Neutrophils"     ~ "Neutrophils",
  merged$SingleR_celltype == "Monocytes"       ~ "Monocytes",
  merged$SingleR_celltype == "DC"              ~ "DC",
  merged$SingleR_celltype == "Epithelial cells"~ "Epithelial cells",
  merged$SingleR_celltype == "Eosinophils"     ~ "Eosinophils",
  merged$SingleR_celltype == "Endothelial cells" ~ "Endothelial cells",
  # Catch‐all—carry forward any other label
  TRUE                                          ~ merged$SingleR_celltype
)

# 5.3. Count cells per broad category
cat("Cells per broad_celltype (post‐collapse):\n")
print(table(merged$broad_celltype))

# ====================================================================================
# 6. FILTER OUT “ENDOTHELIAL CELLS”, “EPITHELIAL CELLS”, AND “FIBROBLASTS”
#    (i.e. drop those from the object)
# ====================================================================================
# 6.1. Define categories to drop
to_drop <- c("Endothelial cells", "Epithelial cells", "Fibroblasts")

# 6.2. Build a vector of cell barcodes to keep (all except the dropped categories)
keep_cells <- colnames(merged)[ ! merged$broad_celltype %in% to_drop ]

# 6.3. Subset the Seurat object by cell names
merged.filtered <- subset(merged, cells = keep_cells)

# 6.4. Verify that unwanted categories are gone
cat("Cells per broad_celltype (after filtering):\n")
print(table(merged.filtered$broad_celltype))

cat("Cell counts: Original =", ncol(merged), 
    "; Filtered =", ncol(merged.filtered), "\n")

# ====================================================================================
# 7. VISUALIZE “BROAD_CELLTYPE” ON UMAP (ONLY RETAINED CELLS)
#    7.1. UMAP colored by broad_celltype
#    7.2. Facet by sample group
# ====================================================================================
# 7.1. UMAP (merged.filtered) colored by broad_celltype
DimPlot(
  merged.filtered,
  reduction = "umap",
  group.by  = "broad_celltype",
  label     = TRUE,
  pt.size   = 0.5
) + labs(title = "Broad cell‐type annotations (filtered)")

# 7.2. UMAP faceted by sample group
DimPlot(
  merged.filtered,
  reduction = "umap",
  group.by  = "broad_celltype",
  split.by  = "group",
  ncol      = 2,
  pt.size   = 0.5,
  label     = FALSE
) + labs(title = "Broad cell‐type by sample group (filtered)")

# ====================================================================================
# 8. STRUCTURE FOR FUTURE ANALYSES (OPTIONAL)
#    If you wish, recluster or rerun UMAP on merged.filtered for finer embeddings:
#    
#    merged.filtered <- NormalizeData(merged.filtered)
#    merged.filtered <- FindVariableFeatures(merged.filtered, selection.method = "vst", nfeatures = 2000)
#    merged.filtered <- ScaleData(merged.filtered, vars.to.regress = c("nCount_RNA","percent.mt"))
#    merged.filtered <- RunPCA(merged.filtered, npcs = 30, verbose = FALSE)
#    merged.filtered <- RunHarmony(merged.filtered, group.by.vars="orig.ident", reduction.use="pca", dims.use=1:30)
#    merged.filtered <- FindNeighbors(merged.filtered, reduction="harmony", dims=1:30)
#    merged.filtered <- FindClusters(merged.filtered, resolution=0.5)
#    merged.filtered <- RunUMAP(merged.filtered, reduction="harmony", dims=1:30)
#
#    Then re‐plot as desired:
#    DimPlot(merged.filtered, reduction="umap_postHarmony", group.by="broad_celltype", label=TRUE)
#
# ====================================================================================
# 9. SAVE THE FINAL ANNOTATED + FILTERED SEURAT OBJECT
# ====================================================================================
saveRDS(
  object = merged.filtered,
  file   = "merged_harmony_annotated_filtered.rds"
)

cat("Finished: merged_harmony_annotated_filtered.rds written to disk.\n")
