# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)

# Load data
#merged_seurat_filtered<-readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_merged_filtered_25percent.rds")
unique(merged_seurat_filtered@meta.data$Number)

# run Harmony -----------
seurat.harmony <- merged_seurat_filtered %>%
  RunHarmony(group.by.vars = 'Patient', plot_convergence = FALSE)

seurat.harmony@reductions

seurat.harmony.embed <- Embeddings(seurat.harmony, "harmony")
seurat.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
seurat.harmony <- seurat.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize 
p3 <- DimPlot(seurat.harmony, reduction = 'umap', group.by = 'Patient', raster=FALSE)
p4 <- DimPlot(seurat.harmony, reduction = 'umap', group.by = 'Number', raster=FALSE)
p5 <- DimPlot(seurat.harmony, reduction = 'umap', label = TRUE, raster=FALSE)

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
grid.arrange(p1, p3)
grid.arrange(p2, p4)
p3
p4
p5


saveRDS(seurat.harmony, file = "scRNAseq_Harmony_integrated_samples25percent.rds")
seurat.harmony
