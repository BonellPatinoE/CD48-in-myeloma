# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)

# get data location
dirs <- list.dirs(path = '/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/GSE223060_RAW/', recursive = F, full.names = F)
dirs


# load the data and assign them a SeuratObject
for (x in dirs) {
  name <- gsub("(?<=\\d)_", ".", paste0("MM_", x), perl = TRUE)
  name <- sub("MM_BM|MM_ND_", "HD_", name)
  name <- sub("MMRF_", "MMRF.", name)
  
  cts <- ReadMtx(
    mtx = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/GSE223060_RAW/', x, '/matrix.mtx'),
    features = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/GSE223060_RAW/', x, '/genes.tsv'),
    cells = paste0('/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/GSE223060_RAW/', x, '/barcodes.tsv'),
    feature.column = 1
  )
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}



# merge datasets
ls()

# MMRF samples with more that 25% plasma cells (MMRF_1325, MMRF1537, MMRF1640, MMRF_2038, MMRF_1267, MMRF_1720, MMRF_1505, MMRF_2251, MMRF_2259)



merged_seurat <- merge(HD_083017, y = c(HD_090617, HD_170531, HD_170607,
                                        HD_2, HD_4, HD_5, #MM_25183, MM_27522.1, MM_27522.2,
                                        #MM_27522.3, MM_27522.4, MM_27522.5,
                                        #MM_27522.6,MM_37692.2, MM_47491.1,
                                        #MM_47491.2, MM_56203.2, MM_57075.3, MM_58408.1, MM_58408.2, MM_59114.1,
                                        #MM_59114.4, MM_60359.1, MM_60359.2, MM_77570, MM_81012.1, MM_81012.2,
                                        #MM_83942, 
                                        #MM_MMRF.1267, 
                                        MM_MMRF.1325, #MM_MMRF.1413, MM_MMRF.1424, 
                                        MM_MMRF.1505,
                                        MM_MMRF.1537, MM_MMRF.1640, #MM_MMRF.1641, MM_MMRF.1686, MM_MMRF.1695, 
                                        MM_MMRF.1720,
                                        #MM_MMRF.1777, MM_MMRF.1941, 
                                        MM_MMRF.1967, MM_MMRF.2038, #MM_MMRF.2168, 
                                        MM_MMRF.2251,
                                        MM_MMRF.2259), 
                       #MM_MMY18273, MM_MMY21940, MM_MMY22933, MM_MMY34339, MM_MMY34600,
                       #MM_MMY40511, MM_MMY47218, MM_MMY67868, MM_MMY70893, MM_MMY74196, MM_MMY80649,
                       #MM_MMY83942, MM_MMY98423),
                       add.cell.ids = ls()[c(3, 4, 5, 6, 7, 8, 9, 33, 36, 37, 38, 42, 45, 46, 48, 49)],
                       project = 'scRNAseq CD48')

# MMRF samples with more that 25% plasma cells (MMRF_1325, MMRF1537, MMRF1640, MMRF_2038, MMRF_1267, MMRF_1720, MMRF_1505, MMRF_2251, MMRF_2259)

# QC & filtering -----------------------

View(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Number', 'Barcode'), 
                                    sep = '_')

# Checking the datasets are being merged properly
unique(merged_seurat@meta.data$Patient)
unique(merged_seurat@meta.data$Number)

# Calculate percentage of mitochondrial genes
merged_seurat[["mitoPercent"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")


# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   nFeature_RNA <= 50000 &
                                   nCount_RNA >= 1000 &
                                   nCount_RNA <= 10000 &
                                   mitoPercent < 10)

merged_seurat_filtered
merged_seurat


# Normalize data 
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)

# standard workflow steps
merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20, reduction = 'pca')

# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient', raster = FALSE) + ggtitle("UMAP by Patient")
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Number', raster = FALSE) + ggtitle("UMAP by Number")
grid.arrange(p1, p2, ncol = 2, nrow = 1)
p1
p2

merged_seurat_filtered

saveRDS(merged_seurat_filtered, file = "/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_merged_filtered_25percent.rds")
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
