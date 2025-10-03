library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(presto)
library(devtools)

# Loading RDS
seurat.harmony <-readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_Harmony_integrated_samples25percent.rds")
seurat.harmony

# Checking the datasets are being merged properly
unique(seurat.harmony@meta.data$Patient)
unique(seurat.harmony@meta.data$Number)

# Run UMAP to see clusters
p5 <- DimPlot(seurat.harmony, reduction = 'umap', label = TRUE, raster=FALSE)
p5
seurat.harmony <- JoinLayers(seurat.harmony)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat.harmony, ident.1 = 2)
head(cluster2.markers, n = 15)

# find all markers distinguishing cluster 0 from clusters 2 and 6
cluster0.markers <- FindMarkers(seurat.harmony, ident.1 = 0, ident.2 = c(2, 6))
head(cluster0.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones

seurat.harmony.markers <- FindAllMarkers(seurat.harmony, only.pos = TRUE)
seurat.harmony.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)


saveRDS(seurat.harmony.markers, file = "/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_seurat_harmony_25percent.rds")

#seurat.harmony.markers<-readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_seurat_harmony_25percent.rds")
seurat.harmony.markers
seurat.harmony.markers%>%filter(cluster ==14)
p5

# B-Cells: CD79A, CD79B, MS4A1  
VlnPlot(seurat.harmony, features = c("CD79A", "CD79B", "MS4A1"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD79A", "CD79B", "MS4A1"))
DotPlot(seurat.harmony, features = c("CD79A", "CD79B", "MS4A1")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CD79B", "MS4A1", "CD79A"))

# CD8+ T-cells: CD8A, CD8B, CD7, CD3E# CD8+ T-cells: CD8A, CD8B, CD7, CD3ETRUE
VlnPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E"))
DotPlot(seurat.harmony, features = c("CD8A", "CD8B", "CD7", "CD3E")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CD8A", "CD8B", "CD7", "CD3E"))

# CD4+ T-cells: CD4, IL7R, CD7, CD3E 
VlnPlot(seurat.harmony, features = c("CD4", "IL7R", "CD7", "CD3E"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD4", "IL7R", "CD7", "CD3E"), min.cutoff = 1, max.cutoff = 4)
DotPlot(seurat.harmony, features = c("CD4", "IL7R", "CD7", "CD3E")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CD4", "IL7R", "CD7", "CD3E"))

# NK Cells: NKG7, GNLY 
VlnPlot(seurat.harmony, features = c("NKG7", "GNLY"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("NKG7", "GNLY"))
DotPlot(seurat.harmony, features = c("NKG7", "GNLY")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("NKG7", "GNLY"))

# PCs: MZB1, SDC1, IGHG1 
VlnPlot(seurat.harmony, features = c("MZB1", "SDC1", "IGHG1"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("MZB1", "SDC1", "IGHG1"))
DotPlot(seurat.harmony, features = c("MZB1", "SDC1", "IGHG1")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("MZB1", "SDC1", "IGHG1"))

# Macrophages: FCGR3A  
VlnPlot(seurat.harmony, features = c("FCGR3A"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("FCGR3A"))
DotPlot(seurat.harmony, features = c("FCGR3A")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("FCGR3A"))

# Monocytes: CD14, LYZ 
VlnPlot(seurat.harmony, features = c("CD14", "LYZ"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CD14", "LYZ"))
DotPlot(seurat.harmony, features = c("CD14", "LYZ")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CD14", "LYZ"))

# Dendritic cells: FCER1A, CLEC10A 
VlnPlot(seurat.harmony, features = c("FCER1A", "CLEC10A"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("FCER1A", "CLEC10A"))
DotPlot(seurat.harmony, features = c("FCER1A", "CLEC10A")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("FCER1A", "CLEC10A"))

# Erythrocytes: AHSP1, HBA, HBB 
VlnPlot(seurat.harmony, features = c("AHSP1", "HBA", "HBB"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("AHSP1", "HBA", "HBB"))
DotPlot(seurat.harmony, features = c("AHSP1", "HBA", "HBB")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("AHSP1", "HBA", "HBB"))

# pDC: CLEC4C, LILRA4 
VlnPlot(seurat.harmony, features = c("CLEC4C", "LILRA4"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("CLEC4C", "LILRA4"))
DotPlot(seurat.harmony, features = c("CLEC4C", "LILRA4")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("CLEC4C", "LILRA4"))

# Neutrophils: AZU1, MPO, ELANE (neutrophils).
VlnPlot(seurat.harmony, features = c("AZU1", "MPO", "ELANE"), pt.size = 0, log=F)
FeaturePlot(seurat.harmony, features = c("AZU1", "MPO", "ELANE"))
DotPlot(seurat.harmony, features = c("AZU1", "MPO", "ELANE")) + RotatedAxis()
seurat.harmony.markers%>%filter(gene %in% c("AZU1", "MPO", "ELANE"))

seurat.harmony.markers%>%filter(cluster==16)
seurat.harmony.markers%>%filter(gene %in% c("NKG7", "GNLY", "CD48", "CD244"))

dev.off()

# Compare features plot between patients
FeaturePlot(seurat.harmony, features = c("NKG7", "GNLY"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony.annotated, features = c("CD48", "CD244", "CD2", "CD58"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)

FeaturePlot(seurat.harmony, features = c("SDC1", "TNFRSF17"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony, features = c("GPRC5D", "TNFRSF17"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)

FeaturePlot(seurat.harmony, features = c("CD4", "CD25", "FOXP3"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
FeaturePlot(seurat.harmony, features = c("CD72", "CD19"), split.by = "Patient")

# Ridge Plot
features <- c("NKG7", "GNLY", "CD48", "CD244", "CD2", "CD58")
RidgePlot(seurat.harmony, features = features, ncol = 2)

# New annotations
new.cluster.ids <- c("CD4_T-cells", 
                     "CD4_T-cells", "CD8_T-cells", "Neutrophils", "Erythrocytes","NK_Cells",
                     "CD4_T-cells", "B-Cells", "B-Cells", "CD4_T-cells","Monocytes",
                     "Neutrophils", "B-Cells", "Plasma_Cells", "Erythrocytes", "Erythrocytes",
                     "Erythrocytes","Plasma_Cells", "pDCs", "Macrophages", "Plasma_Cells", 
                     "Plasma_Cells")
                     

names(new.cluster.ids) <- levels(seurat.harmony)
seurat.harmony.annotated <- RenameIdents(seurat.harmony, new.cluster.ids)
DimPlot(seurat.harmony.annotated, reduction = "umap", label = TRUE, pt.size = 0.3) + NoLegend()
DimPlot(seurat.harmony.annotated, reduction = "umap", label = FALSE, pt.size = 0.3, split.by = "Patient") 


#saveRDS(seurat.harmony.annotated, file = "/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_25percent.harmony.annotated.rds")
seurat.harmony.annotated<-readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_25percent.harmony.annotated.rds")

FeaturePlot(seurat.harmony.annotated, features = c("CD14", "LYZ", "FCGR3A", "CCR10"), split.by = "Patient")
FeaturePlot(seurat.harmony.annotated, features = c("CCR10"), split.by = "Patient", min.cutoff = 0, max.cutoff = 6)
DotPlot(seurat.harmony, features = c("CCR10"), split.by = "Patient") + RotatedAxis()

plots <- FeaturePlot(seurat.harmony.annotated, 
                     features = c("CCR10", "CD14", "LYZ"), 
                     split.by = "Patient", 
                     min.cutoff = 0, 
                     max.cutoff = 6)

CombinePlots(plots, legend = "right")  # Forces the legend to the right

FeaturePlot(seurat.harmony.annotated, 
            features = c("CCR10", "CD14", "LYZ", "TNFRSF17"), 
            split.by = "Patient", 
            min.cutoff = 0, 
            max.cutoff = 6)

# Add new cluster labels to metadata
seurat.harmony$celltype <- Idents(seurat.harmony.annotated)  

# Subset the dataset for monocytes only
monocyte_subset <- subset(seurat.harmony, subset = celltype == "Monocytes")

FeaturePlot(monocyte_subset, 
            features = "CCR10", 
            split.by = "Patient", 
            min.cutoff = 0, 
            max.cutoff = 6)


DotPlot(monocyte_subset, 
        features = "CCR10", 
        group.by = "Patient") + 
  RotatedAxis()  # Rotate axis labels for better visibility

# Subset the dataset for Plasma Cells only
PC_subset <- subset(seurat.harmony, subset = celltype == "Plasma_Cells")

# FeaturePlot for CCR10 and TNFRSF17
FeaturePlot(PC_subset, 
            features = c("CCR10", "TNFRSF17"),
            split.by = "Patient", 
            min.cutoff = 0, 
            max.cutoff = 6)

# DotPlot for CCR10 and TNFRSF17
DotPlot(PC_subset, 
        features = c("CCR10", "TNFRSF17"),
        group.by = "Patient") + 
  RotatedAxis()  # Rotate axis labels for better visibility
