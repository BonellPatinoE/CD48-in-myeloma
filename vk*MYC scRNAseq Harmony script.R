#Load libraries
library(Seurat)
library(harmony)
library(purrr)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)

# 0. For reproducibility
set.seed(1234)

# 1. Point to your directory of count matrices
data_dir <- "/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/vk*MYC scRNAseq/GSE134370_RAW"

# 2. List all .txt files
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)

# 3. Read each file into a Seurat object
seurat_list <- map(files, function(f) {
  df <- read.delim(f,
                   header      = TRUE,
                   sep         = "\t",
                   row.names   = 1,
                   check.names = FALSE)
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  CreateSeuratObject(
    counts       = mat,
    project      = tools::file_path_sans_ext(basename(f)),
    min.cells    = 3,
    min.features = 200
  )
})

# 4. Name your list elements so add.cell.ids works
names(seurat_list) <- tools::file_path_sans_ext(basename(files))
print(names(seurat_list))
#> [1] "Sample1" "Sample2" ... "Sample18"

# 5. Merge all samples
merged <- merge(
  x            = seurat_list[[1]],
  y            = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project      = "AllSamples"
)

# 6. QC: % mito, feature and count filters
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
VlnPlot(merged, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
merged <- subset(merged,
                 subset = nFeature_RNA > 200 &
                   nFeature_RNA < 2500 &
                   percent.mt    < 5)

# 7. Normalize & find variable features
merged <- NormalizeData(merged,
                        normalization.method = "LogNormalize",
                        scale.factor         = 1e4)
merged <- FindVariableFeatures(merged,
                               selection.method = "vst",
                               nfeatures        = 2000)

# 8. Scale (regress out UMI & mito)
merged <- ScaleData(merged,
                    vars.to.regress = c("nCount_RNA","percent.mt"))

# 9. PCA
merged <- RunPCA(merged,
                 npcs    = 30,
                 verbose = FALSE)
ElbowPlot(merged, ndims = 30)

# 10. Harmony batch correction (fully named args!)
merged <- RunHarmony(
  object           = merged,
  group.by.vars    = "orig.ident",
  reduction.use    = "pca",
  dims.use         = 1:30,
  
  # increase the number of Harmony rounds (default 10→20)
  max.iter.harmony   = 20,
  # increase the number of clustering iterations per Harmony round (default 20→50)
  max.iter.cluster   = 50,
  
  # optionally loosen or tighten convergence thresholds
  epsilon.cluster    = 1e-05,  # default
  epsilon.harmony    = 1e-04,  # default
  
  # visualize convergence to check it’s behaving
  plot_convergence    = TRUE,
  verbose             = TRUE
)

# 11. Neighbors, clustering & UMAP on the Harmony embeddings
merged <- FindNeighbors(merged,
                        reduction = "harmony",
                        dims      = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
merged <- RunUMAP(merged,
                  reduction = "harmony",
                  dims      = 1:30)

merged$group <- case_when(
  grepl("^Cont", merged$orig.ident) ~ "Control",
  grepl("^ND",   merged$orig.ident) ~ "ND",
  grepl("^MM",   merged$orig.ident) ~ "MM",
  grepl("^AMG",  merged$orig.ident) ~ "AMG",
  TRUE                               ~ NA_character_
)

# sanity check
table(merged$group)

# now plot by the 4 groups instead of orig.ident
DimPlot(merged, reduction = "umap", group.by = "group", 
        label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP colored by broad condition")

# (Optional) set your new group as the active identity
Idents(merged) <- merged$group
DimPlot(merged, reduction = "umap", label = TRUE) +
  ggtitle("UMAP with Idents set to group")

# 12. Visualization
DimPlot(merged, reduction = "umap", group.by = "group")
DimPlot(merged, reduction = "umap", group.by = "seurat_clusters")
DimPlot(merged, reduction = "umap", label = TRUE)

# show the first few rows of your metadata:
head(merged@meta.data )

# see exactly which values you have in the column you thought was `orig.ident`:
unique( merged$orig.ident )
unique( merged$group )

# 13. Marker detection
# Merge the layers into your default assay
merged <- JoinLayers(merged)

# Make sure you’re on the RNA assay
DefaultAssay(merged) <- "RNA"

# Switch identities to the cluster assignment
Idents(merged) <- "seurat_clusters"

# Inspect how many clusters you have
table(Idents(merged))

# Run FindAllMarkers on the joined data
markers <- FindAllMarkers(
  merged,
  assay            = "RNA",
  slot             = "data",    # use the normalized counts
  only.pos         = TRUE,
  min.pct          = 0.25,
  logfc.threshold  = 0.25
)

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

markers%>%filter(cluster ==1)

# Heatmap of top 10-15 per cluster
top15 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 15, order_by = avg_log2FC)

top15%>%filter(cluster == 5)%>%select(gene)

#Save the Seurat object itself
saveRDS(merged, file = "merged_harmony_markers.rds")

# 3. Also save the markers data.frame
saveRDS(markers, file = "cluster_markers.rds")
write.csv(markers, file = "cluster_markers.csv", row.names = FALSE)
