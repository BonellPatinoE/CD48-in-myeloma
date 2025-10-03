library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(presto)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(ggplot2)

# Load dataset
seurat.harmony.annotated<-readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_CD48_seurat.harmony.annotated.rds")


unique(seurat.harmony.annotated@meta.data$Patient)
unique(seurat.harmony.annotated@meta.data$Number)


# Check which samples contain NK_Cells
table(seurat.harmony.annotated$Number[Idents(seurat.harmony.annotated) == "NK_Cells"])


# Subsetting the Seurat object to include only cells from cluster 0
cluster_of_interest <- "NK_Cells"
seurat_subset <- subset(seurat.harmony.annotated, idents = cluster_of_interest)

# Finding variable features (if not already done)
seurat_subset <- FindVariableFeatures(seurat_subset)

# Scaling the data
seurat_subset <- ScaleData(seurat_subset)

# Perform PCA on the subsetted data
seurat_subset <- RunPCA(seurat_subset, npcs = 30)

# Find neighbors and clusters within the subsetted data
seurat_subset <- FindNeighbors(seurat_subset, dims = 1:10)
seurat_subset <- FindClusters(seurat_subset, resolution = 0.5)

# Run UMAP for visualization
seurat_subset <- RunUMAP(seurat_subset, dims = 1:10)

# UMAP plot to visualize subtypes within cluster 0
DimPlot(seurat_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(seurat_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Patient") 

seurat_subset@meta.data

# Find markers for each subtype
markers <- FindAllMarkers(seurat_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# Ensure clusters are ordered from 0 to 11
markers$cluster <- factor(markers$cluster, levels = 0:12)

# Select top 10 genes per cluster
top_genes <- markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) %>%  # Arrange genes within each cluster
  ungroup()

# Reshape data for heatmap
heatmap_data <- top_genes %>%
  select(cluster, gene, avg_log2FC) %>%
  spread(cluster, avg_log2FC, fill = 0) 

# Reorder rows so genes appear stacked per cluster
heatmap_data <- heatmap_data %>%
  arrange(match(gene, unique(top_genes$gene))) %>%
  column_to_rownames(var = "gene")

# Order columns from 0 to 12
heatmap_data <- heatmap_data[, order(as.numeric(colnames(heatmap_data)))]

# Generate the heatmap
pheatmap(heatmap_data, 
         scale = "row",        # Standardizes rows (z-score)
         cluster_cols = FALSE, # Keep clusters in numeric order
         cluster_rows = FALSE, # Prevent clustering genes, keeping them grouped
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         color = colorRampPalette(c("purple", "black", "yellow"))(50))  # Blue-low, Red-high expression

table(seurat.harmony.annotated$Number)  # Before subsetting
table(seurat_subset$Number)             # After subsetting

markers%>%filter(cluster==4)
markers%>%filter(gene=="RNF213")
FeaturePlot(seurat_subset, features = c("KLRD1"), split.by = "Patient", min.cutoff = 0, max.cutoff = 4, raster=FALSE)

#"CD2", "GZMH", "KLRC2", "FCER1G" 

# Checking the datasets are being merged properly
unique(seurat_subset@meta.data$Patient)
unique(seurat_subset@meta.data$Number)

# Compare with clusters
DimPlot(seurat_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Patient") 

library(Seurat)
library(cowplot)

p <- FeaturePlot(seurat_subset, features = c("FCGR3A", "NCAM1"), 
                 split.by = "Patient", min.cutoff = 0, max.cutoff = 6, raster = FALSE, combine = FALSE)

# Arrange the plots while keeping the legend
plot_grid(plotlist = p, ncol = 2)

q<- FeaturePlot(seurat_subset, features = c("GNLY", "PRF1", "NKG7", "FGFBP2", "KLRB1", "KLRF1", "KLRC2", "KLRC3"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6, raster=FALSE, combine = FALSE)
plot_grid(plotlist = q, ncol = 8)

FeaturePlot(seurat_subset, features = c("FCGR3A", "NCAM1"), split.by = "Patient", min.cutoff = 0, max.cutoff = 4, raster=FALSE)
FeaturePlot(seurat_subset, features = c("GNLY", "PRF1", "NKG7", "FGFBP2"), split.by = "Patient", min.cutoff = 0, max.cutoff = 4, raster=FALSE)
FeaturePlot(seurat_subset, features = c("GNLY", "PRF1", "NKG7", "FGFBP2", "KLRB1", "KLRF1", "KLRC2", "KLRC3"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)
FeaturePlot(seurat_subset, features = c("KLRC1"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)


Idents(seurat_subset) <- "seurat_clusters"
Idents(seurat_subset) <- "Patient"

DotPlot(seurat_subset, features = c("GNLY", "PRF1", "NKG7", "FGFBP2", "KLRB1", "KLRF1", "KLRC2", "KLRC3"),
        cols = c("red", "blue"), dot.scale = 8, split.by = "Patient") + RotatedAxis()


FeaturePlot(seurat_subset, features = c(#"NCAM1", "GNLY", "PRF1", "NKG7", "FGFBP2",
                                        "KLRB1", "KLRF1", "FCER1G", "KLRC2", "KLRC3"
                                        ), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)

FeaturePlot(seurat_subset, features = c("FOS"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)
FeaturePlot(seurat_subset, features = c("FCGR3A", "KLRD1", "NKG7"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)

DotPlot(seurat_subset, features = c("FCGR3A", "KLRD1", "NKG7"),
        cols = c("red", "blue"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

VlnPlot(seurat_subset, features = c("FOS"), pt.size = 0.5)
VlnPlot(seurat_subset, features = c("CD244"), pt.size = 0.5)

DotPlot(seurat_subset, features = c("LDHA", "ALDOA"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()
VlnPlot(seurat_subset, features = c("LDHA", "ALDOA"), pt.size = 0.5)

DotPlot(seurat_subset, features = c("MT-ATP8", "MT-CO3", "MT-ND4L", "MT-CYB"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()
VlnPlot(seurat_subset, features = c("MT-ATP8", "MT-CO3", "MT-ND4L", "MT-CYB"), pt.size = 0.5)

DotPlot(seurat_subset, features = c("CDC42", "CD160", "SYK"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("PRF1", "GZMA", "GZMK", "GNLY"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("KLRG1", "TIGIT", "LAG3", "KLRC1", "CISH"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("LAG3", "CISH"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("JUN", "LGALS1", "KLF6"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("SKIL", "SMAD7"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

DotPlot(seurat_subset, features = c("LAG3", "CISH"),
        cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

FeaturePlot(seurat_subset, features = c("INPP5D", "SH2D1A"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6, raster=FALSE)
VlnPlot(seurat_subset, features = c("INPP5D", "SH2D1A"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
DotPlot(seurat_subset, features = c("INPP5D", "SH2D1A"),
        cols = c("red", "blue"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

FeaturePlot(seurat_subset, features = c("FOS", "JUN", "FOSB", "JUNB"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6, raster=FALSE)
VlnPlot(seurat_subset, features = c("FOS", "JUN", "FOSB", "JUNB"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
DotPlot(seurat_subset, features = c("FOS", "JUN", "FOSB", "JUNB"),
        cols = c("red", "blue"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

VlnPlot(seurat_subset, features = c("JUN", "BATF"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("CD48"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("CD244"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()

FeaturePlot(seurat_subset, features = c("S100A9", "CD14", "CD33", "LYZ"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6, raster=FALSE)
VlnPlot(seurat_subset, features = c("S100A9", "CD14", "CD33", "LYZ"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
DotPlot(seurat_subset, features = c("S100A9", "CD14", "CD33", "LYZ"),
        cols = c("red", "blue"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

#########################
Idents(seurat_subset) <- "seurat_clusters"
Idents(seurat_subset) <- "Patient"

FeaturePlot(seurat_subset, features = c("SIGLEC7"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)
VlnPlot(seurat_subset, features = c("NFATC2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
#Src, Fyn, Yes, Lck, Lyn, Hck, Fgr and Blk


#"CD2", "GZMH", "KLRC2", "FCER1G", "IL32" 
VlnPlot(seurat_subset, features = c("LAG3", "TIGIT", "KLRG1", "KLRC1", "CISH", "JUN", "LGALS1", "KLF6", "SKIL", "SMAD7"), pt.size = 0.5)
VlnPlot(seurat_subset, features = c("BATF"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("BATF3"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("DDIT3"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("TGFBR1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("TGFBR2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("TGFB1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("EOMES"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("HAVCR2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("KLRG1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("KLRC1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("NT5E"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("FOS"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() 
VlnPlot(seurat_subset, features = c("JUN"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() 
VlnPlot(seurat_subset, features = c("ENTPD1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #CD39
VlnPlot(seurat_subset, features = c("CISH"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() 
VlnPlot(seurat_subset, features = c("PDCD1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("LAG3"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("HAVCR2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("KLRC2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("KLRC3"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("KLRB1"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()

VlnPlot(seurat_subset, features = c("CD244"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()
VlnPlot(seurat_subset, features = c("CD226"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()

VlnPlot(seurat_subset, features = c("CD226"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.3, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("TGFB1"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.3, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("JUN"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("FOS"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("BATF"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("EOMES"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

### CMV ADAPTIVE, CD56-Dim NK cell population ###
VlnPlot(seurat_subset, features = c("FCER1G"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("KLRB1"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("KLRF1"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("KLRC2"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("KLRC3"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("KLRC1"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## SAP

VlnPlot(seurat_subset, features = c("INPP5D"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## SHIP1

VlnPlot(seurat_subset, features = c("INPP5D"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## SHIP1

VlnPlot(seurat_subset, features = c("SH2D1B"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## EAT-2

VlnPlot(seurat_subset, features = c("HAVCR2"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## TIM3

VlnPlot(seurat_subset, features = c("TIGIT"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()## TIGIT

VlnPlot(seurat_subset, features = c("LAG3"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## LAG3

VlnPlot(seurat_subset, features = c("CTLA4"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() ## CTLA4

VlnPlot(seurat_subset, features = c("CD96"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() #CD96

VlnPlot(seurat_subset, features = c("SIGLEC9"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() # SIGLEC9

VlnPlot(seurat_subset, features = c("KLRG1"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() #KLRG1

VlnPlot(seurat_subset, features = c("KLRC3"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means()  #KLRC1


VlnPlot(seurat_subset, features = c("CISH"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() #CISH

VlnPlot(seurat_subset, features = c("CD226"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() #CISH

VlnPlot(seurat_subset, features = c("CD244"), split.by = "Patient", pt.size = 0.5) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.5, color = "black") + 
  stat_compare_means() #CISH

VlnPlot(seurat_subset, features = c("SH2D1A"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means()  # SAP
VlnPlot(seurat_subset, features = c("SH2D1B"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #EAT-2
VlnPlot(seurat_subset, features = c("INPP5D"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() # SHIP1
VlnPlot(seurat_subset, features = c("PTPN11"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("PTPN6"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHP1
VlnPlot(seurat_subset, features = c("SH3BP2"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() # 3BP2

VlnPlot(seurat_subset, features = c("SH2D1A", "INPP5D"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() # 3BP2
DotPlot(seurat_subset, features = c("SH2D1A", "INPP5D"), cols = c("blue", "red"), dot.scale = 8, split.by = "Patient") + RotatedAxis()

seurat_subset$SAP_SHIP1_ratio <- seurat_subset@assays$RNA@data["INPP5D", ] / seurat_subset@assays$RNA@data["SH2D1A", ]
seurat_subset$SAP_SHIP2_ratio <- seurat_subset@assays$RNA@data["PTPN11", ] / seurat_subset@assays$RNA@data["SH2D1A", ]

VlnPlot(seurat_subset, features = c("SAP_SHIP1_ratio"), split.by = "Patient", pt.size = 0.5) + 
  stat_compare_means()

#VlnPlot(seurat_subset, features = c("SAP_SHIP2_ratio"), split.by = "Patient", pt.size = 0.5) + 
  #stat_compare_means()

FeaturePlot(seurat_subset, features = c("NCAM1", "FCGR3A", "CD226"), split.by = "Patient", min.cutoff = 0, max.cutoff = 2, raster=FALSE)
VlnPlot(seurat_subset, features = c("LAG3", "TIGIT", "KLRG1", "KLRC1", "CISH"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2
VlnPlot(seurat_subset, features = c("SOCS3"), split.by = "Patient", pt.size = 0.5)+ stat_compare_means() #SHIP2


#IFIT1, IFIT2, IFIT3, IFIH1

VlnPlot(seurat_subset, features = c("TGFBR2", "EOMES", "HAVCR2", "KLRG1", "KLRC1", "NT5E"), pt.size = 0.5)

## SCORE ANALYSIS ##

# Define a list of genes associated with exhaustion
exhaustion_genes <- c("LAG3", "TIGIT", "KLRG1", "KLRC1", "CISH")

# Calculate the module score for exhaustion
seurat_object <- AddModuleScore(seurat_subset, 
                                features = list(exhaustion_genes), 
                                name = "ExhaustionScore")

# Visualize the exhaustion score in a feature plot
VlnPlot(seurat_object, "ExhaustionScore1", split.by = "Patient", raster=FALSE) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()


# Define a list of genes associated with siglecs
siglec_genes <- c("SIGLEC7", "SIGLEC9")

# Calculate the module score for exhaustion
seurat_siglec <- AddModuleScore(seurat_subset, 
                                features = list(siglec_genes), 
                                name = "SiglecScore")

# Visualize the exhaustion score in a feature plot
VlnPlot(seurat_siglec, "SiglecScore1", split.by = "Patient", raster=FALSE) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()


# Define a list of genes associated with TGF-B pathway
TGFB_path <- c("JUN", "LGALS1", "KLF6", "SKIL", "SMAD7")
TGFB_path <- c("JUN", "LGALS1", "KLF6", "SKIL", "SMAD7")

# Calculate the module score for exhaustion
seurat_object <- AddModuleScore(seurat_subset, 
                                features = list(TGFB_path), 
                                name = "TGFBScore")

# Visualize the exhaustion score in a feature plot
FeaturePlot(seurat_object, features = "TGFBScore1", split.by = "Patient")
VlnPlot(seurat_object, "TGFBScore1", split.by = "Patient", raster=FALSE) + 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()



# Define gene sets
cytotoxicity_genes <- c("GZMA", "GZMB", "GZMH", "GZMM", "GZMK", "GNLY", "PRF1", "CTSW")
inflammatory_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CXCL10", "CXCL9", "IL1B", "IL6", "IL7", "IL15", "IL18")
stress_genes <- c("BAG3", "CALU", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HIF1A", "HSP90AA1", "HSP90AB1", 
                  "HSP90B1", "HSPA1A", "HSPA1B", "HSPA6", "HSPB1", "HSPH1", "IER2", "JUN", "JUNB", "NFKBIA", 
                  "NFKBIZ", "RGS2", "SLC2A3", "SOCS3", "UBC", "ZFAND2A", "ZFP36", "ZFP36L1")

HLA_dep_inh_rec_genes<- c("KIR2DL1", "KIR2DL3", "KIR3DL1", "KIR3DL2", "LILRB1", "LAG3") 
HLA_indep_inh_rec_genes<- c("PDCD1", "SIGLEC7", "CD300A", "CD96", "IL1RAPL1", "TIGIT", "HAVCR2") 
HLA_dep_act_rec_genes<-c("KIR2DL4", "CD160", "KLRC2") 
HLA_indep_act_rec_genes<- c("NCR3", "NCR1", "KLRK1", "CRTAM", "FCGR3A")

# Load AUCell and calculate scores
library(AUCell)

expression_matrix <- GetAssayData(seurat_subset, assay = "RNA", slot = "data")
cells_rankings <- AUCell_buildRankings(expression_matrix)  # Build gene expression rankings

# Calculate AUC for each gene set
cytotoxicity_scores <- AUCell_calcAUC(cytotoxicity_genes, cells_rankings)
inflammatory_scores <- AUCell_calcAUC(inflammatory_genes, cells_rankings)
stress_scores <- AUCell_calcAUC(stress_genes, cells_rankings)
HLA_dep_inh_rec_scores<-AUCell_calcAUC(HLA_dep_inh_rec_genes, cells_rankings)
HLA_indep_inh_rec_scores<-AUCell_calcAUC(HLA_indep_inh_rec_genes, cells_rankings)
HLA_dep_act_rec_scores<-AUCell_calcAUC(HLA_dep_act_rec_genes, cells_rankings)
HLA_indep_act_rec_scores<-AUCell_calcAUC(HLA_indep_act_rec_genes, cells_rankings)


# Add scores to metadata of Seurat object
seurat_subset[["cytotoxicity_score"]] <- as.numeric(cytotoxicity_scores@assays@data@listData$AUC)
seurat_subset[["inflammatory_score"]] <- as.numeric(inflammatory_scores@assays@data@listData$AUC)
seurat_subset[["stress_score"]] <- as.numeric(stress_scores@assays@data@listData$AUC)
seurat_subset[["HLA_dep_inh_rec_scores"]] <- as.numeric(HLA_dep_inh_rec_scores@assays@data@listData$AUC)
seurat_subset[["HLA_indep_inh_rec_score"]] <- as.numeric(HLA_indep_inh_rec_scores@assays@data@listData$AUC)
seurat_subset[["HLA_dep_act_rec_scores"]] <- as.numeric(HLA_dep_act_rec_scores@assays@data@listData$AUC)
seurat_subset[["HLA_indep_act_rec_scores"]] <- as.numeric(HLA_indep_act_rec_scores@assays@data@listData$AUC)

# Visualize the scores
VlnPlot(seurat_subset, features = c("cytotoxicity_score"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("inflammatory_score"), split.by = "Patient")+ 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("stress_score"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("HLA_dep_inh_rec_scores", "HLA_indep_inh_rec_score", "HLA_dep_act_rec_scores", "HLA_indep_act_rec_scores"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()


VlnPlot(seurat_subset, features = c("HLA_dep_inh_rec_scores"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("HLA_indep_inh_rec_score"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("HLA_dep_act_rec_scores"), split.by = "Patient")+ 
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

VlnPlot(seurat_subset, features = c("HLA_indep_act_rec_scores"), split.by = "Patient")+
  stat_summary(fun.data = function(x) data.frame(y = mean(x), ymin = mean(x), ymax = mean(x)), 
               geom = "crossbar", width = 0.6, color = "black") + 
  stat_compare_means()

# Check the structure of the calculated scores
str(HLA_indep_inh_rec_scores)  # Ensure it has the right structure

# Assuming `HLA_indep_inh_rec_scores` is an AUCell object and has been created correctly
# Check if the AUC is available
if ("AUC" %in% names(HLA_indep_inh_rec_scores@assays@data@listData)) {
  seurat_subset[["HLA_indep_inh_rec_score"]] <- 
    as.numeric(HLA_indep_inh_rec_scores@assays@data@listData$AUC)
} else {
  stop("AUC values not found in HLA_indep_inh_rec_scores.")
}

VlnPlot(seurat_subset, features = c("HLA_indep_inh_rec_score"), split.by = "Patient") + 
  stat_compare_means()


#########################
Idents(seurat_subset) <- "seurat_clusters"
Idents(seurat_subset) <- "Patient"

# Perform differential expression analysis between HD and MM
de_genes <- FindMarkers(seurat_subset, ident.1 = "MM", ident.2 = "HD", min.pct = 0.25, logfc.threshold = 0.25)
head(de_genes)

# Visualize DEGs with a volcano plot
library(EnhancedVolcano)

EnhancedVolcano(de_genes,
                lab = rownames(de_genes), #de_genes$label,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 3.0,
                xlim = c(-4, 4))  # Set x-axis limits

# Genes more expressed in MM patients
genes_MM <- rownames(de_genes[de_genes$avg_log2FC > 0, ])
genes_MM<-data.frame(genes_MM)
#write_csv(genes_MM, "genes_MM.csv")

# Genes more expressed in HD patients
genes_HD <- rownames(de_genes[de_genes$avg_log2FC < 0, ])
genes_HD<-data.frame(genes_HD)
#write_csv(genes_HD, "genes_HD.csv")

# View the lists
head(genes_MM)
head(genes_HD)

# Create a column for the genes
genes <- rownames(de_genes)

# Create a column for the patient group
patient_group <- ifelse(de_genes$avg_log2FC > 0, "MM", "HD")

# Create a column for log2FC values
log2FC <- de_genes$avg_log2FC

# Combine them into a data frame
gene_patient_table <- data.frame(Gene = genes, Patient = patient_group, Log2FC = log2FC)

# View the table
head(gene_patient_table)
#write_csv(gene_patient_table, "DEG in NK cells between HD and MM patients.csv")

# Genes upregulated in MM (negative avg_log2FC and p-value < 0.05)
mm_genes_sig <- rownames(de_genes[de_genes$avg_log2FC < 0 & de_genes$p_val_adj < 0.05, ])
head(mm_genes_sig)

# Visualize the top DEGs with a heatmap

# Set the identity classes to 'Patient' for comparison
Idents(seurat_subset) <- "seurat_clusters"
Idents(seurat_subset) <- "Patient"

top_de_genes <- rownames(de_genes[order(de_genes$p_val_adj), ])[1:203]
DoHeatmap(seurat_subset, features = top_de_genes) 

# Visualize the top DEGs with a heatmap based on clusters
DoHeatmap(seurat_subset, features = top_de_genes, size = 3) + 
  scale_fill_viridis_c()

# Average expression per cluster
avg_expression <- AverageExpression(seurat_subset, return.seurat = TRUE)

# Extract the averaged expression data
avg_data <- GetAssayData(avg_expression, slot = "data")

# Create the heatmap for the top DE genes
DoHeatmap(avg_expression, features = top_de_genes, size = 3) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))


seurat_subset@meta.data

# Perform differential expression analysis between HD and MM
de_genes_nk <- FindMarkers(seurat_subset, 
                           ident.1 = "HD", 
                           ident.2 = "MM", 
                           min.pct = 0.25, 
                           logfc.threshold = 0.25)

# Filter genes based on adjusted p-value threshold (e.g., p_val_adj < 0.05)
sig_genes_nk <- de_genes_nk[de_genes_nk$p_val_adj < 0.05, ]

# Count the total number of differentially expressed genes
num_de_genes <- nrow(sig_genes_nk)

# Display the count
print(paste("Total number of differentially expressed genes: ", num_de_genes))


## Get the transcripts numbers on the specific phenotype
## List of specific genes you are interested in
genes_of_interest <- c("NCAM1", "GNLY", "PRF1", "NKG7", "FGFBP2", "KLRB1", "KLRF1", "FCER1G", "KLRC2", "KLRC3")  # Replace with your gene names

# Ensure the genes are present in the counts matrix
valid_genes <- genes_of_interest[genes_of_interest %in% rownames(counts_matrix)]

# Filter the counts matrix for the genes of interest
filtered_counts <- counts_matrix[valid_genes, , drop = FALSE]

# Create an empty matrix to store the summarized counts for specific genes
cluster_counts_filtered <- Matrix(0, nrow = nrow(filtered_counts), ncol = length(levels(cluster_ids)),
                                  dimnames = list(rownames(filtered_counts), levels(cluster_ids)))

# Summarize the counts for each cluster and specific genes
for (cluster in levels(cluster_ids)) {
  cluster_cells <- which(cluster_ids == cluster)
  cluster_counts_filtered[, cluster] <- rowSums(filtered_counts[, cluster_cells, drop = FALSE])
}

# Convert the summarized counts matrix to a data frame
summary_counts_filtered_df <- as.data.frame(as.matrix(cluster_counts_filtered))

# Display the summarized counts for specific genes
print(summary_counts_filtered_df)
log2(summary_counts_filtered_df+1)


# Heatmap to identify clusters
library(pheatmap)

# Define the desired order of clusters
desired_order <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

# Reorder the columns of the data frame according to the desired order
ordered_df <- summary_counts_filtered_df[, as.character(desired_order)]

# Create the heatmap with reordered clusters
pheatmap(log2(ordered_df + 1), 
         cluster_rows = TRUE, 
         cluster_cols = FALSE,  # If you want to control column clustering
         scale = "none",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Heatmap of Gene Expression by Cluster")


