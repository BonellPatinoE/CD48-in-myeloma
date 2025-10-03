# ========================================
# 0. Load Libraries
# ========================================
library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(hdf5r)
library(ggalluvial)

# ========================================
# 1. Load and Subset Seurat Object (MM only)
# ========================================
seurat.obj <- readRDS("/Users/bonellpatinoescobar/Desktop/UCSF/Wiita Lab/R/scRNAseq CD48/scRNAseq CD48/scRNAseq_25percent.harmony.annotated.rds")

# Subset MM cells
mm.cells <- rownames(seurat.obj@meta.data[grepl("MM", seurat.obj@meta.data$Patient), ])
seurat.mm <- subset(seurat.obj, cells = mm.cells)

# ========================================
# 2. Ensure proper cell type annotation
# ========================================
# Only run this if your cell types were renamed manually in a previous step
if (!"celltype" %in% colnames(seurat.mm@meta.data)) {
  seurat.mm$celltype <- Idents(seurat.mm)
}

# Set identity to cell type
Idents(seurat.mm) <- seurat.mm$celltype

# ========================================
# 3. Build CellChat Object
# ========================================
data.input <- GetAssayData(seurat.mm[["RNA"]], layer = "data")
meta <- seurat.mm@meta.data
meta$cell_type <- as.character(Idents(seurat.mm))  # Ensure it's character
meta$samples <- meta$Patient

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")
cellchat@DB <- CellChatDB.human

# ========================================
# 4. Run CellChat Analysis
# ========================================
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ========================================
# 5. Analyze CD48–CD244 Specifically
# ========================================
cd48_lr <- subsetCommunication(cellchat, signaling = "CD48", slot.name = "net")  # net includes more interactions
cd48_cd244 <- cd48_lr[cd48_lr$ligand == "CD48" & cd48_lr$receptor == "CD244", ]
print(cd48_cd244)

# Try bubble plot (from net slot)
if (nrow(cd48_cd244) > 0) {
  tryCatch({
    netVisual_bubble(cellchat,
                     sources.use = unique(cd48_cd244$source),
                     targets.use = unique(cd48_cd244$target),
                     pairLR.use = cd48_cd244[, c("ligand", "receptor")],
                     remove.isolate = FALSE)
  }, error = function(e) {
    message("⚠️ Bubble plot failed: ", e$message)
  })
}

# ========================================
# 6. Visualize CD48 Signaling (All)
# ========================================
pathways.show <- "CD48"

# Circle plot
par(mfrow = c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# ========================================
# 7. Expression of CD48 and CD244
# ========================================
# CellChat gene expression plot
plotGeneExpression(cellchat, c("CD48", "CD244"))

# Seurat violin plot
VlnPlot(seurat.mm, features = c("CD48", "CD244"), group.by = "celltype")


##########
# All inferred communications (strict slot)
all_comm <- subsetCommunication(cellchat, slot.name = "netP")

# Interactions where NK cells are the sender
nk_as_sender <- all_comm %>% filter(source == "NK_Cells")

# Interactions where NK cells are the receiver
nk_as_receiver <- all_comm %>% filter(target == "NK_Cells")

# View top interactions by probability
nk_as_sender %>% arrange(desc(prob)) %>% head(10)
nk_as_receiver %>% arrange(desc(prob)) %>% head(10)

if (nrow(nk_as_sender) > 0) {
  netVisual_bubble(cellchat,
                   sources.use = "NK_Cells",
                   targets.use = unique(nk_as_sender$target),
                   remove.isolate = FALSE)
}

if (nrow(nk_as_receiver) > 0) {
  netVisual_bubble(cellchat,
                   sources.use = unique(nk_as_receiver$source),
                   targets.use = "NK_Cells",
                   remove.isolate = FALSE)
}

# Circle plot – NK cells as receiver (you can change to sender too)
netVisual_aggregate(cellchat, sources.use = NULL, targets.use = "NK_Cells", layout = "circle")

# Heatmap – highlight how NK cells communicate with others
netVisual_heatmap(cellchat, sources.use = "NK_Cells")
netVisual_heatmap(cellchat, targets.use = "NK_Cells")

# Compute overall signaling roles
cellchat <- netAnalysis_computeCentrality(cellchat)

# Visualize signaling roles
ht1 <-netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")  # For senders
ht2 <-netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")  # For receivers

ht1 + ht2

# Barplot for NK cells only
netAnalysis_signalingRole_scatter(cellchat, signaling = NULL)  # Interactive

