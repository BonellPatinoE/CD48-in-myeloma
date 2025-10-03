# Load required libraries
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(readr)
library(pheatmap)
library(forcats)

# Step 1: Read count matrix
raw <- read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
counts <- as.matrix(raw[,-1])
rownames(counts) <- raw$Gene

# Step 2: Create colData (sample metadata)
ids <- data.frame(sample_id = colnames(counts))
rownames(ids) <- ids$sample_id

# Step 3: Create DESeq2 object with no experimental design
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ids, design = ~1)

# Step 4: Normalize with VST (blind=TRUE for unbiased transformation)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Step 5: Convert Ensembl IDs to HGNC gene symbols
ens_ids <- rownames(vsd_mat)
conversion <- gconvert(query = ens_ids, organism = "hsapiens", target = "HGNC",
                       mthreshold = 1, filter_na = FALSE)

# Attach gene symbols to VST matrix
vsd_df <- as.data.frame(vsd_mat)
vsd_df$gene_symbol <- conversion$target[match(rownames(vsd_df), conversion$input)]

# Step 6: Filter for NK cell ligands
nk_ligands <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", 
                "PCNA", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", 
                "NECTIN2", "PVR", "CD48")

nk_df <- vsd_df %>% 
  filter(gene_symbol %in% nk_ligands) %>%
  select(-gene_symbol) %>%
  mutate(gene = conversion$target[match(rownames(.), conversion$input)])

# Step 7: Tidy data for plotting
nk_long <- pivot_longer(nk_df, -gene, names_to = "sample_id", values_to = "vst_expr")

# Step 8: Optional metadata (e.g., patient or visit info)
nk_long <- nk_long %>%
  mutate(patient = str_extract(sample_id, "[0-9]{4}"),
         visit = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"))

nk_long <- nk_long %>%
  group_by(gene) %>%
  mutate(gene = fct_reorder(gene, vst_expr, .fun = median, .desc = TRUE)) %>%
  ungroup()

# Step 9: Plot expression

ggplot(nk_long, aes(x = fct_reorder(gene, vst_expr, .fun = median, .desc = TRUE), y = vst_expr)) +
  geom_boxplot(outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.2, alpha = 0.1, color = "red", size = 1) +
  theme_classic() +
  labs(
    title = "VST-normalized Expression of NK Cell Ligands (Myeloma Samples)",
    x = "NK Ligand Gene (ordered by median expression)",
    y = "VST Expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Step 10: Prepare matrix for heatmap
# Create NK-ligand-only matrix (rows = ligands, columns = samples)
nk_matrix <- vsd_df %>%
  filter(gene_symbol %in% nk_ligands) %>%
  select(all_of(colnames(vsd_mat)))  # only sample columns
rownames(nk_matrix) <- vsd_df$gene_symbol[vsd_df$gene_symbol %in% nk_ligands]

# Scale across rows (z-score per gene)
nk_matrix_scaled <- t(scale(t(nk_matrix)))

# Step 11: Create sample annotation (e.g., visit)
# Start from your 'ids' data frame
annotation_col <- ids %>%
  rownames_to_column("tmp") %>%  # remove existing rownames
  mutate(visit = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"),
         disease_stage = case_when(
           visit == "1" ~ "Newly Diagnosed",
           visit %in% c("2", "3", "4", "5", "6") ~ "Relapsed",
           TRUE ~ NA_character_
         )) %>%
  select(sample_id, disease_stage) %>%
  column_to_rownames("sample_id")

############
############
pheatmap(nk_matrix_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         fontsize_row = 9,
         main = "NK Ligand Expression (Newly Diagnosed vs Relapsed)")

cluster_rows = TRUE
cluster_cols = TRUE

pheatmap(nk_matrix_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,      # cluster genes (ligands)
         cluster_cols = TRUE,      # cluster samples
         cutree_cols = 4,          # optional: force 4 clusters for samples
         show_colnames = FALSE,
         fontsize_row = 9,
         main = "NK Ligand Expression (Clustered Samples)",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2")

nk_matrix_scaled_capped <- nk_matrix_scaled
nk_matrix_scaled_capped[nk_matrix_scaled > 2] <- 2
nk_matrix_scaled_capped[nk_matrix_scaled < -2] <- -2

custom_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(nk_matrix_scaled_capped,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_cols = 3,
         show_colnames = FALSE,
         fontsize_row = 9,
         main = "NK Ligand Expression (Z-score capped at ±3)",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",
         color = custom_palette)

################
###############


### Cytogenetic risk

translocations_df<-read_tsv('SeqFISH Files_MMRF_CoMMpass_IA19_genome_tumor_only_mm_igtx_pairoscope.tsv') 
translocations_df <- data.frame(translocations_df)

translocations_df2 <- data.frame(
  sample_id = translocations_df$SAMPLE,
  patient   = str_extract(translocations_df$SAMPLE, "[0-9]{4}"),
  visit     = str_extract(translocations_df$SAMPLE, "_[0-9]_") %>% str_extract("[0-9]"),
  type      = str_extract(translocations_df$SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
  NDS2  = translocations_df$NSD2_CALL,
  CCND3 = translocations_df$CCND3_CALL,
  MYC   = translocations_df$MYC_CALL,
  MAFA  = translocations_df$MAFA_CALL,
  CCND1 = translocations_df$CCND1_CALL,
  CCND2 = translocations_df$CCND2_CALL,
  MAF   = translocations_df$MAF_CALL,
  MAFB  = translocations_df$MAFB_CALL
)

translocations_df2 <- translocations_df2 %>%
  filter(type == "BM")

annotation_col <- translocations_df2 %>%
  select(sample_id, NDS2, MYC, MAF, CCND1, MAFB) %>%
  filter(sample_id %in% colnames(nk_matrix_scaled)) %>%
  column_to_rownames("sample_id")

###################
# Match sample order
nk_matrix_scaled <- nk_matrix_scaled[, rownames(annotation_col)]

nk_matrix_scaled_capped <- nk_matrix_scaled
nk_matrix_scaled_capped[nk_matrix_scaled > 3] <- 3
nk_matrix_scaled_capped[nk_matrix_scaled < -3] <- -3

custom_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

pheatmap(nk_matrix_scaled_capped,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_cols = 4,
         show_colnames = FALSE,
         fontsize_row = 9,
         color = custom_palette,
         main = "NK Ligand Expression (VST-normalized, stratified by translocations)")

########
seqFISH_raw<-read_tsv('SeqFISH Files_MMRF_CoMMpass_IA19_genome_gatk_cna_seqFISH.tsv')
seqFISH_raw <- data.frame(seqFISH_raw)

seqFISH <- seqFISH_raw %>%
  transmute(
    sample_id  = SAMPLE,
    patient    = str_extract(SAMPLE, "[0-9]{4}"),
    visit      = str_extract(SAMPLE, "_[0-9]_") %>% str_extract("[0-9]"),
    type       = str_extract(SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
    del17p13   = SeqWGS_Cp_17p13_20percent,
    TP53       = SeqWGS_Cp_TP53_20percent,
    Gain1q21   = SeqWGS_Cp_1q21_20percent,
    del1p22    = SeqWGS_Cp_1p22_20percent,
    MYC_CNA    = SeqWGS_Cp_MYC_20percent
  ) %>%
  filter(type == "BM")

# Combine the two datasets by sample_id
cyto_combined <- translocations_df2 %>%
  left_join(seqFISH, by = "sample_id", suffix = c("_trans", "_cna"))

annotation_col <- cyto_combined %>%
  select(sample_id, NDS2, MYC, MAF, CCND1, MAFB, del17p13, Gain1q21, del1p22) %>%
  column_to_rownames("sample_id")

annotation_col <- cyto_combined %>%
  select(sample_id, NDS2, MYC, MAF, CCND1, MAFB, del17p13, Gain1q21, del1p22) %>%
  column_to_rownames("sample_id")

#################

# Subset NK matrix to samples in annotation
common_samples <- intersect(colnames(nk_matrix_scaled), rownames(annotation_col))
nk_matrix_scaled <- nk_matrix_scaled[, common_samples]
annotation_col <- annotation_col[common_samples, ]

annotation_col[] <- lapply(annotation_col, as.factor)

custom_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

###
pheatmap(nk_matrix_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_cols = 4,
         show_colnames = FALSE,
         fontsize_row = 9,
         color = custom_palette,
         main = "NK Ligand Expression (VST, Stratified by Cytogenetics)")
###

annotation_colors <- list(
  NDS2      = c("0" = "white", "1" = "#e41a1c"),  # red
  MYC       = c("0" = "white", "1" = "#377eb8"),  # blue
  MAF       = c("0" = "white", "1" = "#4daf4a"),  # green
  CCND1     = c("0" = "white", "1" = "#984ea3"),  # purple
  MAFB      = c("0" = "white", "1" = "#ff7f00"),  # orange
  del17p13  = c("0" = "white", "1" = "#a65628"),  # brown
  Gain1q21  = c("0" = "white", "1" = "#f781bf"),  # pink
  del1p22   = c("0" = "white", "1" = "#999999")   # gray
)

annotation_col[] <- lapply(annotation_col, function(x) factor(x, levels = c("0", "1")))

###
pheatmap(nk_matrix_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_cols = 4,
         show_colnames = FALSE,
         fontsize_row = 9,
         color = custom_palette,
         main = "NK Ligand Expression with Cytogenetic Annotations")
#####

nk_matrix_scaled[nk_matrix_scaled > 4] <- 4
nk_matrix_scaled[nk_matrix_scaled < -4] <- -4

breaks <- seq(-6, 6, length.out = 101)

pheatmap(nk_matrix_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_cols = 4,
         show_colnames = FALSE,
         fontsize_row = 9,
         color = custom_palette,
         breaks = breaks,
         main = "NK Ligand Expression (Z-score capped at ±3)")

heat <- pheatmap(nk_matrix_scaled,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 cutree_cols = 4,
                 show_colnames = FALSE,
                 fontsize_row = 9,
                 color = custom_palette,
                 main = "NK Ligand Expression (Z-score capped at ±3)")

find("pheatmap")

sample_clusters <- cutree(heat$tree_col, k = 4)
cluster_df <- data.frame(
  sample_id = names(sample_clusters),
  cluster = LETTERS[sample_clusters]
)

nk_long <- nk_matrix_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "zscore") %>%
  left_join(cluster_df, by = "sample_id")

library(ggplot2)

ggplot(nk_long, aes(x = cluster, y = zscore, fill = cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.7) +
  facet_wrap(~gene, scales = "free_y") +
  theme_classic() +
  labs(title = "NK Ligand Expression by Heatmap Cluster",
       x = "Sample Cluster", y = "Z-scored Expression") +
  scale_fill_brewer(palette = "Dark2")

annotation_col$Cluster <- cluster_df$cluster[match(rownames(annotation_col), cluster_df$sample_id)]
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = LETTERS[1:4])


## 3 CLUSTERS
library(RColorBrewer)

sample_clusters <- cutree(heat$tree_col, k = 3)
cluster_df <- data.frame(
  sample_id = names(sample_clusters),
  cluster = LETTERS[sample_clusters]  # A, B, C
)

annotation_col$Cluster <- cluster_df$cluster[match(rownames(annotation_col), cluster_df$sample_id)]
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = LETTERS[1:3])  # A, B, C

annotation_colors$Cluster <- setNames(brewer.pal(2, "Dark2"), LETTERS[1:3])
########

# --- 1) Make sure the column is named "cluster"
colnames(cluster_df)
# Should print: "sample_id" "cluster"

# --- 2) Turn that into an ordered factor A→B→C
cluster_df$cluster <- factor(cluster_df$cluster, levels = c("A","B","C"))

# --- 3) Now get your columns in A→B→C order
ordered_samples <- cluster_df %>% 
  arrange(cluster) %>% 
  pull(sample_id)

# --- 4) How many samples per cluster?
cluster_sizes <- as.numeric(table(cluster_df$cluster))
# e.g. c( n_A, n_B, n_C )

# --- 5) Compute the “gaps” positions between clusters
gaps_col <- cumsum(cluster_sizes)[1:(length(cluster_sizes)-1)]
# That gives something like c(n_A, n_A + n_B)

# --- 6) Finally draw the heatmap with separators
pheatmap(
  nk_matrix_scaled[, ordered_samples],
  annotation_col    = annotation_col[ordered_samples, , drop = FALSE],
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,   # keep your A→B→C ordering
  gaps_col          = gaps_col, 
  show_colnames     = FALSE,
  fontsize_row      = 9,
  color             = custom_palette,
  main              = "NK Ligand Expression (clusters A→B→C)"
)


#########

nk_long <- nk_matrix_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "zscore")

nk_long <- left_join(nk_long, cluster_df, by = "sample_id")

ggplot(nk_long, aes(x = cluster, y = zscore, fill = cluster)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.7) +
  facet_wrap(~gene, scales = "free_y") +
  theme_classic() +
  labs(
    title = "NK Ligand Expression by Heatmap Cluster",
    x = "Sample Cluster", y = "Z-scored Expression"
  ) +
  scale_fill_brewer(palette = "Dark2")


nk_pca_input <- t(nk_matrix_scaled)

pca_res <- prcomp(nk_pca_input, scale. = FALSE)
pca_df <- as.data.frame(pca_res$x)  # PC1, PC2, etc.
pca_df$sample_id <- rownames(pca_df)

pca_df <- left_join(pca_df, cluster_df, by = "sample_id")

ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "PCA of NK Ligand Expression",
       subtitle = "Colored by Heatmap-Derived Cluster (A–C)",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"))

ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, type = "norm") +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    title = "PCA of NK Ligand Expression",
    subtitle = "Ellipses show 95% confidence intervals per cluster",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  )

###################################

# Step 1: Annotate relapse status from 'ids'
ids <- ids %>%
  mutate(
    visit = str_extract(sample_id, "_[0-9]_") %>% str_extract("[0-9]"),
    relapsed = case_when(
      visit == "1" ~ "Newly Diagnosed",
      visit %in% c("2", "3", "4", "5", "6") ~ "Relapsed",
      TRUE ~ NA_character_
    )
  )

# Step 2: Merge with cyto_combined to build full annotation matrix
annotation_col <- cyto_combined %>%
  left_join(ids %>% select(sample_id, relapsed), by = "sample_id") %>%
  filter(sample_id %in% colnames(nk_matrix_scaled)) %>%
  select(sample_id, relapsed, NDS2, MYC, MAF, CCND1, MAFB, del17p13, Gain1q21, del1p22) %>%
  column_to_rownames("sample_id")

# Step 3: Convert all annotations to factors
annotation_col[] <- lapply(annotation_col, as.factor)

# Step 4: Define color mapping for each annotation
annotation_colors <- list(
  relapsed   = c("Newly Diagnosed" = "white", "Relapsed" = "black"),
  NDS2       = c("0" = "white", "1" = "#e41a1c"),   # red
  MYC        = c("0" = "white", "1" = "#377eb8"),   # blue
  MAF        = c("0" = "white", "1" = "#4daf4a"),   # green
  CCND1      = c("0" = "white", "1" = "#984ea3"),   # purple
  MAFB       = c("0" = "white", "1" = "#ff7f00"),   # orange
  TP53       = c("0" = "white", "1" = "#66c2a5"),   # turquoise
  del17p13   = c("0" = "white", "1" = "#a65628"),   # brown
  Gain1q21   = c("0" = "white", "1" = "#f781bf"),   # pink
  del1p22    = c("0" = "white", "1" = "#999999")    # gray
)


# If you've already clustered with cutree_cols = 3:
heat <- pheatmap(nk_matrix_scaled,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 cutree_cols = 3,
                 show_colnames = FALSE,
                 fontsize_row = 9,
                 color = custom_palette,
                 main = "NK Ligand Expression")

# Extract cluster assignments (3 clusters)
sample_clusters <- cutree(heat$tree_col, k = 3)

# Convert to letters A, B, C
cluster_df <- data.frame(
  sample_id = names(sample_clusters),
  Cluster = LETTERS[sample_clusters]
)

annotation_col$Cluster <- cluster_df$Cluster[match(rownames(annotation_col), cluster_df$sample_id)]
annotation_col$Cluster <- factor(annotation_col$Cluster, levels = LETTERS[1:3])

library(RColorBrewer)

annotation_colors$Cluster <- setNames(brewer.pal(3, "Dark2"), LETTERS[1:3])
# Should return TRUE
all(levels(annotation_col$Cluster) %in% names(annotation_colors$Cluster))

########
# --- 1) Make sure the column is named "cluster"
colnames(cluster_df)
# Should print: "sample_id" "cluster"

# --- 2) Turn that into an ordered factor A→B→C
cluster_df$cluster <- factor(cluster_df$Cluster, levels = c("A","B","C"))

# --- 3) Now get your columns in A→B→C order
ordered_samples <- cluster_df %>% 
  arrange(cluster) %>% 
  pull(sample_id)

# --- 4) How many samples per cluster?
cluster_sizes <- as.numeric(table(cluster_df$cluster))
# e.g. c( n_A, n_B, n_C )

# --- 5) Compute the “gaps” positions between clusters
gaps_col <- cumsum(cluster_sizes)[1:(length(cluster_sizes)-1)]
# That gives something like c(n_A, n_A + n_B)

# --- 6) Finally draw the heatmap with separators
pheatmap(
  nk_matrix_scaled[, ordered_samples],
  annotation_col    = annotation_col[ordered_samples, , drop = FALSE],
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,   # keep your A→B→C ordering
  gaps_col          = gaps_col, 
  show_colnames     = FALSE,
  fontsize_row      = 9,
  color             = custom_palette,
  main              = "NK Ligand Expression"
)


#########


#######################################
library(tidyverse)

# Define binary features to test
features_to_test <- c("relapsed", "NDS2", "MYC", "MAF", "CCND1", "MAFB", "del17p13", "Gain1q21", "del1p22")
annotation_col$relapsed <- ifelse(annotation_col$relapsed == "Relapsed", "1",
                                  ifelse(annotation_col$relapsed == "Newly Diagnosed", "0", NA))

# Run Fisher's exact test + proportions + count for each feature
feature_stats_detailed <- lapply(features_to_test, function(feature) {
  data_filtered <- annotation_col[!is.na(annotation_col[[feature]]), ]
  
  tab <- table(data_filtered[[feature]], data_filtered$Cluster)
  prop_tab <- prop.table(tab, margin = 1)  # proportions per row (0 or 1)
  
  test <- fisher.test(tab)
  
  list(
    feature = feature,
    p_value = test$p.value,
    significance = case_when(
      test$p.value < 0.001 ~ "***",
      test$p.value < 0.01 ~ "**",
      test$p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    counts = tab,
    proportions = prop_tab
  )
})

# Create summary results table with p-values and stars
cluster_assoc_results <- bind_rows(lapply(feature_stats_detailed, function(x) {
  tibble(
    Feature = x$feature,
    p_value = x$p_value,
    significance = x$significance
  )
})) %>%
  arrange(p_value)

print(cluster_assoc_results)  ### p-values per cytogenetic risk

# Build binary matrix for heatmap plotting

# Arrange samples by cluster order (A–D)
ordered_samples <- cluster_df %>%
  arrange(Cluster) %>%
  pull(sample_id)

# Extract annotation for those samples and convert to binary matrix (1 = altered)
binary_matrix <- annotation_col[ordered_samples, features_to_test]
binary_matrix <- as.data.frame(lapply(binary_matrix, function(x) as.numeric(x == "1")))
binary_matrix <- t(binary_matrix)  # Now: rows = features, cols = samples


feature_stats_detailed[[1]]$counts       # for a feature
feature_stats_detailed[[1]]$proportions  # for a feature

cluster_assoc_results <- bind_rows(lapply(feature_stats_detailed, function(x) {
  # Default NA if no "1" row exists (i.e., no positive cases)
  if ("1" %in% rownames(x$proportions)) {
    prop_pos <- x$proportions["1", ]
    enriched_cluster <- names(which.max(prop_pos))
  } else {
    enriched_cluster <- NA_character_
  }
  
  tibble(
    Feature = x$feature,
    p_value = x$p_value,
    significance = x$significance,
    enriched_cluster = enriched_cluster
  )
})) %>%
  arrange(p_value)

print(cluster_assoc_results) ## p-values per cytogenetic risk associated to each cluster

library(knitr)
kable(cluster_assoc_results, caption = "Cluster Enrichment by Cytogenetic Feature")



#########################################
###########################################
###########################################

# 1) Extract CD48 VST from your VST matrix (vsd_df):
#library(dplyr)
#library(tidyr)

#cd48_df <- vsd_df %>% 
  # vsd_df has rownames = Ensembl IDs and a column gene_symbol
#  filter(gene_symbol == "CD48") %>%
#  select(-gene_symbol) %>% 
#  rownames_to_column("ens_id") %>% 
#  pivot_longer(
#    cols = -ens_id,
#    names_to  = "sample_id",
#    values_to = "CD48_expr"
#  ) %>% 
#  select(sample_id, CD48_expr)

# 2) Join it into your cox_data (which already has a sample_id column):
#cox_data2 <- cox_data %>%
#  left_join(cd48_df, by = "sample_id")

# 3) Now you can safely coerce or scale CD48_expr:
#cox_data2 <- cox_data2 %>%
#  mutate(
#    CD48_expr = as.numeric(CD48_expr),
#    Cluster   = factor(Cluster, levels = c("A","B","C")),
    # …and the rest of your factors…
#  )

# 4) Fit your full model including CD48_expr:
#cox_full <- coxph(
#  Surv(time, status) ~ CD48_expr + Cluster + relapsed + 
#    NDS2 + MYC + MAF + MAFB + del1p22 + 
#    del17p13 + Gain1q21,
#  data = cox_data2
#)

# 5) Forest‐plot it:
#library(forestmodel)
#forest_model(cox_full) +
#  ggtitle("Multivariate Cox Model\nCD48 + Clusters + Cytogenetics") +
#  theme(plot.title = element_text(hjust = 0.5))


###########################################
###########################################
###########################################


library(ggplot2)

# Make sure cluster column exists and is named "cluster"
pca_df <- pca_df %>%
  rename(cluster = Cluster) %>%
  mutate(cluster = as.factor(cluster))  # ensure it's a factor

# PCA plot with ellipses
ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster, fill = cluster)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, type = "norm") +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    title = "PCA of NK Ligand Expression",
    subtitle = "Ellipses show 95% confidence intervals per cluster",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")
  )


library(tidyverse)
library(ggpubr)
library(rstatix)

# 1) Build named cluster lookup (A–C)
cluster_letters <- setNames(LETTERS[sample_clusters], names(sample_clusters))

# 2) Prepare the base long table
nk_long_base <- nk_matrix_scaled %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "zscore")

# 0) prerequisites
library(ggplot2)
library(rstatix)
library(ggpubr)

# assume nk_long_base and cluster_letters are defined as before…

# 3) Loop and print
for(g in unique(nk_long_base$gene)) {
  df <- nk_long_base %>%
    filter(gene == g) %>%
    mutate(Cluster = factor(cluster_letters[sample_id], levels = LETTERS[1:3]))
  
  # per-gene Wilcoxon tests
  stat_test <- df %>%
    wilcox_test(zscore ~ Cluster, p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Cluster")
  
  # build the plot from df and stat_test
  p <- ggplot(df, aes(x = Cluster, y = zscore, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 2) +
    facet_wrap(~gene, scales = "free_y") +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = paste0(g, " expression by cluster"),
      x = NULL, y = "Z-score"
    ) +
    stat_pvalue_manual(
      stat_test,
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = TRUE,
      inherit.aes = FALSE
    )
  
  print(p)     # <-- this actually draws it
  Sys.sleep(0.2)  # slight pause so you can see it
}


# 1. Install + load patchwork if you haven’t yet
# install.packages("patchwork")
library(patchwork)

# 2. Build a list of plots instead of printing each one
plot_list <- lapply(unique(nk_long_base$gene), function(g) {
  df <- nk_long_base %>%
    filter(gene == g) %>%
    mutate(Cluster = factor(cluster_letters[sample_id], levels = LETTERS[1:3]))
  
  stat_test <- df %>%
    wilcox_test(zscore ~ Cluster, p.adjust.method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Cluster")
  
  ggplot(df, aes(x = Cluster, y = zscore, fill = Cluster)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
    theme_classic() +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = g, x = NULL, y = "Z-score") +
    stat_pvalue_manual(
      stat_test,
      label = "p.adj.signif",
      tip.length = 0.01,
      hide.ns = TRUE,
      inherit.aes = FALSE
    )
})

# 3. Combine them in a grid, e.g. 3 columns:
combined <- wrap_plots(plot_list, ncol = 3)

# 4. Draw it
print(combined)

