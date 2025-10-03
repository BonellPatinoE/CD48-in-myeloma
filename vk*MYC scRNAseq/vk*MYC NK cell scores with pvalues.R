# ================================================
# 0. SET UP: Seed + Libraries
# ================================================
set.seed(1234)

# Core packages
library(Seurat)
library(AUCell)
library(Matrix)    # for coercing to sparse matrices
library(dplyr)
library(ggplot2)
library(ggpubr)    # for stat_compare_means on violins

# ================================================
# 1. LOAD YOUR MERGED, ANNOTATED SEURAT OBJECT
#    (Must already have "broad_celltype" with "NK cells" and "group" columns)
# ================================================
seurat_subset <- readRDS("merged_harmony_annotated_filtered.rds")

# Quick sanity checks:
print(seurat_subset)
table(seurat_subset$broad_celltype)  # should include “NK cells”
table(seurat_subset$group)           # e.g., “Control”, “ND”, “MM”, “AMG”

# ================================================
# 2. SUBSET TO NK CELLS ONLY
# ================================================
nk_only <- subset(seurat_subset, subset = broad_celltype == "NK cells")

# Verify:
table(nk_only$broad_celltype)  # should show only “NK cells”
table(nk_only$group)           # breakdown of groups in NK cells

# ================================================
# 3. DETECT WHICH “Gzm*” AND “Il1rapl*” GENES EXIST IN NK CELLS
# ================================================
available_gzms    <- grep("^Gzm", rownames(nk_only),  value = TRUE)
available_il1rapl <- grep("^Il1rapl", rownames(nk_only), value = TRUE)

message("Present granzyme‐family genes in NK cells:\n  ", paste(available_gzms, collapse = ", "))
message("Present Il1rapl‐family genes in NK cells:\n  ", paste(available_il1rapl, collapse = ", "))

# ================================================
# 4. DEFINE MURINE GENE SETS (ONLY SYMBOLS THAT EXIST)
# ================================================
# 4.1. Cytotoxicity (mouse granzymes + Perforin + Cathepsin W)
cytotoxicity_genes <- intersect(
  available_gzms,
  c("Gzma","Gzmb","Gzmc","Gzmm","Gzmk")
)
cytotoxicity_genes <- union(cytotoxicity_genes, c("Prf1","Ctsw"))

# 4.2. Inflammatory cytokines/chemokines (mouse)
inflammatory_genes <- intersect(
  rownames(nk_only),
  c("Ccl2","Ccl3","Ccl4","Ccl5","Cxcl9","Cxcl10","Il1b","Il6","Il7","Il15","Il18")
)

# 4.3. Stress response (mouse HSPs, immediate‐early, etc.)
stress_genes <- intersect(
  rownames(nk_only),
  c(
    "Bag3","Calu","Dnajb1","Dusp1","Egr1","Fos","Fosb","Hif1a",
    "Hsp90aa1","Hsp90ab1","Hsp90b1","Hspa1a","Hspa1b","Hspa8",
    "Hspb1","Hsph1","Ier2","Jun","Junb","Nfkbia","Nfkbiz",
    "Rgs2","Slc2a3","Socs3","Ubc","Zfand2a","Zfp36","Zfp36l1"
  )
)

# 4.4. HLA‐dependent inhibitory receptors (Ly49 family + Pirb + Lag3)
HLA_dep_inh_rec_genes <- intersect(
  rownames(nk_only),
  c("Klra1","Klra3", "Klra6", "Klra7","Klra8", "Klra9", "Klra10", "Klra12", "Pirb","Lag3")
)

# 4.5. HLA‐independent inhibitory receptors (mouse)
HLA_indep_inh_rec_genes <- intersect(
  rownames(nk_only),
  c("Pdcd1","Siglecg","Cd300a","Cd96", available_il1rapl, "Tigit","Havcr2")
)

# 4.6. HLA‐dependent activating receptors (Klrc2 + Cd160)
HLA_dep_act_rec_genes <- intersect(
  rownames(nk_only),
  c("Klrc2","Cd160")
)

# 4.7. HLA‐independent activating receptors (Ncr1, Klrk1, Crtam, Fcgr3)
HLA_indep_act_rec_genes <- intersect(
  rownames(nk_only),
  c("Ncr1","Klrk1","Crtam","Fcgr3")
)

# 4.8. Exhaustion gene set (murine markers)
exhaustion_genes <- intersect(
  rownames(nk_only),
  c("Cish","Lag3","Tigit","Klrc1","Klrg1")
)

# Warn if any set becomes empty after intersection:
all_sets <- list(
  cytotoxicity_genes        = cytotoxicity_genes,
  inflammatory_genes        = inflammatory_genes,
  stress_genes              = stress_genes,
  HLA_dep_inh_rec_genes     = HLA_dep_inh_rec_genes,
  HLA_indep_inh_rec_genes   = HLA_indep_inh_rec_genes,
  HLA_dep_act_rec_genes     = HLA_dep_act_rec_genes,
  HLA_indep_act_rec_genes   = HLA_indep_act_rec_genes,
  exhaustion_genes          = exhaustion_genes
)

for (set_name in names(all_sets)) {
  if (length(all_sets[[set_name]]) == 0) {
    warning("After intersection, no genes remain for: ", set_name)
  }
}

# Print final lists for verification:
message("Final cytotoxicity_genes:        ", paste(cytotoxicity_genes, collapse = ", "))
message("Final inflammatory_genes:        ", paste(inflammatory_genes, collapse = ", "))
message("Final stress_genes:              ", paste(stress_genes, collapse = ", "))
message("Final HLA_dep_inh_rec_genes:     ", paste(HLA_dep_inh_rec_genes, collapse = ", "))
message("Final HLA_indep_inh_rec_genes:   ", paste(HLA_indep_inh_rec_genes, collapse = ", "))
message("Final HLA_dep_act_rec_genes:     ", paste(HLA_dep_act_rec_genes, collapse = ", "))
message("Final HLA_indep_act_rec_genes:   ", paste(HLA_indep_act_rec_genes, collapse = ", "))
message("Final exhaustion_genes:          ", paste(exhaustion_genes, collapse = ", "))

# ================================================
# 5. EXTRACT EXPRESSION MATRIX AND CONVERT TO SPARSE
# ================================================
expr_matrix_nk        <- GetAssayData(nk_only, assay = "RNA", slot = "data")
expr_matrix_nk_sparse <- as(expr_matrix_nk, "dgCMatrix")

# ================================================
# 6. OPTIONAL: FILTER TO HVGs ∪ GENE SETS
#    (Reduces memory usage for AUCell rankings)
# ================================================
hvg_nk <- VariableFeatures(nk_only)

features_for_ranking_nk <- unique(c(
  hvg_nk,
  unlist(all_sets)  # include all gene sets
))

expr_matrix_pruned_nk <- expr_matrix_nk_sparse[features_for_ranking_nk, ]

# ================================================
# 7. BUILD AUCell RANKINGS ON THE PRUNED SPARSE MATRIX
# ================================================
cells_rankings_nk <- AUCell_buildRankings(
  expr_matrix_pruned_nk,
  plotStats = FALSE
)

# ================================================
# 8. CALCULATE AUC SCORES FOR EACH GENE SET (NK cells)
# ================================================
cyto_AUC_nk          <- AUCell_calcAUC(cytotoxicity_genes,       cells_rankings_nk)
inflammatory_AUC_nk  <- AUCell_calcAUC(inflammatory_genes,       cells_rankings_nk)
stress_AUC_nk        <- AUCell_calcAUC(stress_genes,             cells_rankings_nk)
HLA_dep_inh_AUC_nk   <- AUCell_calcAUC(HLA_dep_inh_rec_genes,    cells_rankings_nk)
HLA_indep_inh_AUC_nk <- AUCell_calcAUC(HLA_indep_inh_rec_genes,  cells_rankings_nk)
HLA_dep_act_AUC_nk   <- AUCell_calcAUC(HLA_dep_act_rec_genes,    cells_rankings_nk)
HLA_indep_act_AUC_nk <- AUCell_calcAUC(HLA_indep_act_rec_genes,  cells_rankings_nk)
exh_AUC_nk           <- AUCell_calcAUC(exhaustion_genes,         cells_rankings_nk)

# ================================================
# 9. ADD AUC SCORES BACK INTO THE NK‐ONLY METADATA
# ================================================
nk_only[["cytotoxicity_score"]]        <- as.numeric(cyto_AUC_nk@assays@data@listData$AUC)
nk_only[["inflammatory_score"]]         <- as.numeric(inflammatory_AUC_nk@assays@data@listData$AUC)
nk_only[["stress_score"]]               <- as.numeric(stress_AUC_nk@assays@data@listData$AUC)
nk_only[["HLA_dep_inh_rec_score"]]      <- as.numeric(HLA_dep_inh_AUC_nk@assays@data@listData$AUC)
nk_only[["HLA_indep_inh_rec_score"]]    <- as.numeric(HLA_indep_inh_AUC_nk@assays@data@listData$AUC)
nk_only[["HLA_dep_act_rec_score"]]      <- as.numeric(HLA_dep_act_AUC_nk@assays@data@listData$AUC)
nk_only[["HLA_indep_act_rec_score"]]    <- as.numeric(HLA_indep_act_AUC_nk@assays@data@listData$AUC)
nk_only[["exhaustion_score"]]           <- as.numeric(exh_AUC_nk@assays@data@listData$AUC)

# Verify that metadata columns were added:
head(nk_only@meta.data[, c(
  "group",
  "cytotoxicity_score",
  "inflammatory_score",
  "stress_score",
  "HLA_dep_inh_rec_score",
  "HLA_indep_inh_rec_score",
  "HLA_dep_act_rec_score",
  "HLA_indep_act_rec_score",
  "exhaustion_score"
)])

# ================================================
# 10. VISUALIZE AUC SCORES WITH VIOLIN PLOTS GROUPED BY “group”
#     Each dot = one NK cell; pt.size = 0.5; black crossbar = mean;
#     Wilcoxon tests (all 6 pairwise comparisons) via stat_compare_means()
# ================================================
plot_score_by_group <- function(seurat_obj, score_col, title_text) {
  df <- FetchData(seurat_obj, vars = c("group", score_col))
  df$group <- factor(df$group, levels = c("Control", "ND", "MM", "AMG"))
  
  ggplot(df, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_violin(width = 0.8, alpha = 0.3, color = NA) +
    geom_jitter(size = 0.5, width = 0.2, alpha = 0.4, color = "black") +
    stat_summary(
      fun       = mean,
      geom      = "crossbar",
      width     = 0.4,
      color     = "black",
      size      = 0.8,
      fatten    = 2
    ) +
    stat_compare_means(
      aes(label = paste0("p=", signif(..p.adj.., 2))),
      method       = "wilcox.test",
      comparisons  = list(
        c("Control","ND"),
        c("Control","MM"),
        c("Control","AMG"),
        c("ND","MM"),
        c("ND","AMG"),
        c("MM","AMG")
      ),
      p.adjust.method = "BH",
      hide.ns      = TRUE,
      tip.length   = 0.01,
      bracket.size = 0.3
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = title_text,
      x     = NULL,
      y     = title_text
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 12),
      axis.title.y    = element_text(size = 12)
    )+
    scale_y_log10()
}

# 10.1. Cytotoxicity score
p_cyto_nk <- plot_score_by_group(
  nk_only,
  score_col  = "cytotoxicity_score",
  title_text = "Cytotoxicity AUC"
)
print(p_cyto_nk)

# 10.2. Inflammatory score
p_inflam_nk <- plot_score_by_group(
  nk_only,
  score_col  = "inflammatory_score",
  title_text = "Inflammatory AUC"
)
print(p_inflam_nk)

# 10.3. Stress response score
p_stress_nk <- plot_score_by_group(
  nk_only,
  score_col  = "stress_score",
  title_text = "Stress response AUC"
)
print(p_stress_nk)

# 10.4. HLA‐dependent inhibitory receptor score
p_hla_dep_inh_nk <- plot_score_by_group(
  nk_only,
  score_col  = "HLA_dep_inh_rec_score",
  title_text = "HLA‐dep Inhibitory AUC"
)
print(p_hla_dep_inh_nk)

# 10.5. HLA‐independent inhibitory receptor score
p_hla_indep_inh_nk <- plot_score_by_group(
  nk_only,
  score_col  = "HLA_indep_inh_rec_score",
  title_text = "HLA‐indep Inhibitory AUC"
)
print(p_hla_indep_inh_nk)

# 10.6. HLA‐dependent activating receptor score
p_hla_dep_act_nk <- plot_score_by_group(
  nk_only,
  score_col  = "HLA_dep_act_rec_score",
  title_text = "HLA‐dep Activating AUC"
)
print(p_hla_dep_act_nk)

# 10.7. HLA‐independent activating receptor score
p_hla_indep_act_nk <- plot_score_by_group(
  nk_only,
  score_col  = "HLA_indep_act_rec_score",
  title_text = "HLA‐indep Activating AUC"
)
print(p_hla_indep_act_nk)

# 10.8. Exhaustion score (murine markers)
p_exh_nk <- plot_score_by_group(
  nk_only,
  score_col  = "exhaustion_score",
  title_text = "Exhaustion AUC"
)
print(p_exh_nk)

# ================================================
# 11. OPTIONAL: SAVE THE NK‐ONLY OBJECT WITH ALL SCORES
# ================================================
# Uncomment to save:
# saveRDS(nk_only, file = "nk_only_with_all_gene_set_scores.rds")
