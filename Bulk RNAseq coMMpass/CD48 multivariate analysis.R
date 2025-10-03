#––– libraries
library(tidyverse)
library(DESeq2)
library(gprofiler2)
library(survival)
library(survminer)
library(forestmodel)

#––– 1) read count matrix
raw_counts <- read_tsv("Expression Estimates - Gene Based_MMRF_CoMMpass_IA19_star_geneUnstranded_counts.tsv")
counts_mat <- as.matrix(raw_counts %>% select(-Gene))
rownames(counts_mat) <- raw_counts$Gene

#––– 2) read and build cytogenetics + relapse metadata
trans_df <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_tumor_only_mm_igtx_pairoscope.tsv")
cna_df   <- read_tsv("SeqFISH Files_MMRF_CoMMpass_IA19_genome_gatk_cna_seqFISH.tsv")

trans_meta <- trans_df %>%
  transmute(
    sample_id = SAMPLE,
    patient   = str_extract(SAMPLE, "[0-9]{4}") %>% paste0("MMRF_", .),
    visit     = str_extract(SAMPLE, "_[0-9]_") %>% str_extract("[0-9]"),
    type      = str_extract(SAMPLE, "_[A-Z]{2}_") %>% str_extract("[A-Z]+"),
    relapsed  = if_else(visit == "1", 0, 1),
    NSD2      = NSD2_CALL,
    CCND1     = CCND1_CALL
  )

cna_meta <- cna_df %>%
  transmute(
    sample_id = SAMPLE,
    del17p13  = SeqWGS_Cp_17p13_20percent,
    Gain1q21  = SeqWGS_Cp_1q21_20percent
  )

cyto_combined <- trans_meta %>%
  filter(type == "BM") %>%
  left_join(cna_meta, by = "sample_id") %>%
  mutate_at(vars(relapsed, NSD2, CCND1, del17p13, Gain1q21),
            ~ factor(.x, levels = c(0,1)))

#––– 3) subset your count matrix to those same samples
common <- intersect(colnames(counts_mat), cyto_combined$sample_id)
counts_mat     <- counts_mat[, common]
cyto_combined  <- cyto_combined %>% filter(sample_id %in% common) %>%
  column_to_rownames("sample_id")

#––– 4) make a DESeq2 object + VST
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = cyto_combined,
                              design    = ~1)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

#––– 5) convert ENSG → HGNC and build a little data-frame
mapping <- gconvert(query       = rownames(vst_mat),
                    organism    = "hsapiens",
                    target      = "HGNC",
                    mthreshold  = 1,
                    filter_na   = FALSE)

vsd_df <- as.data.frame(vst_mat) %>%
  rownames_to_column("ensembl_id") %>%
  mutate(gene_symbol = mapping$target[match(ensembl_id, mapping$input)])

#––– 6) pull out CD48
if (!"CD48" %in% vsd_df$gene_symbol) {
  stop("CD48 not found in your VST matrix via the gconvert mapping.")
}
cd48_df <- vsd_df %>%
  filter(gene_symbol == "CD48") %>%
  select(-ensembl_id, -gene_symbol) %>%
  pivot_longer(cols       = everything(),
               names_to   = "sample_id",
               values_to  = "CD48_expr")

#––– 7) read survival and clean it
surv <- read_tsv("MMRF_CoMMpass_IA19_STAND_ALONE_SURVIVAL.tsv") %>%
  transmute(
    patient = PUBLIC_ID,
    status  = if_else(!is.na(deathdy) & deathdy > 1, 1, 0),
    time    = coalesce(deathdy, lstalive, lvisitdy)
  )

#––– 8) build final cox data-frame
cox_df <- cd48_df %>%
  left_join(cyto_combined %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  left_join(surv, by = "patient") %>%
  select(time, status, CD48_expr, relapsed, NSD2, CCND1, Gain1q21, del17p13) %>%
  drop_na()

#––– 9) fit and summarize Cox model
cox_mod <- coxph(
  Surv(time, status) ~ CD48_expr + relapsed + NSD2 + CCND1 + Gain1q21 + del17p13,
  data = cox_df
)
summary(cox_mod)

# --- 10) Plot the model


forest_model(cox_mod) +
  ggtitle("Multivariate Cox Model")

### Interaction between CD48 and Gain1q
library(survival)
library(survminer)

# 1) Fit your interaction model (you already have this):
cox_int <- coxph(
  Surv(time, status) ~ CD48_expr * Gain1q21 +
    relapsed + NSD2 + CCND1 + del17p13,
  data = cox_df
)

# 2) Build a newdata with four rows:
#    – CD48 “low” = median–1, “high” = median+1
#    – 1q21 gain = 0 or 1
#    – all other factors held at their reference level (“0”)
med_cd48 <- median(cox_df$CD48_expr, na.rm = TRUE)

newdata <- expand.grid(
  CD48_expr = c(med_cd48 - 1, med_cd48 + 1),
  Gain1q21  = levels(cox_df$Gain1q21),
  relapsed  = levels(cox_df$relapsed)[1],
  NSD2      = levels(cox_df$NSD2)[1],
  CCND1     = levels(cox_df$CCND1)[1],
  del17p13  = levels(cox_df$del17p13)[1]
)

# sanity check: should be 4 × 6
print(newdata)

# 3) Add a human‐friendly “group” column for the legend
newdata$group <- rep(c("CD48 low, 1q–",
                       "CD48 high,1q–",
                       "CD48 low, 1q+",
                       "CD48 high,1q+"),
                     each = 1)

# Make sure all of these columns are factors with exactly the same levels
for(v in c("Gain1q21","relapsed","NSD2","CCND1","del17p13")) {
  newdata[[v]] <- factor(newdata[[v]],
                         levels = levels(cox_df[[v]]))
}

# 4) Re-fit the survival curves at those four profiles
fit_int <- survfit(cox_int, newdata = newdata)

# 2) plot, *including* data=newdata
ggsurvplot(
  fit_int,
  data        = newdata,                    # <— this line is the fix
  palette     = c("#1b9e77","#d95f02","#7570b3","#e7298a"),
  legend.labs = newdata$group,
  conf.int    = TRUE,
  xlab        = "Days",
  ylab        = "Survival Probability",
  title       = "Interaction: CD48 × Gain1q21"
)

# look at the 4 predicted survival sets
summary(fit_int, newdata = newdata)



## Interacion CD48 and NSD2
library(survival)
library(dplyr)
library(survminer)

# — 1) Fit the interaction model
cox_int_NSD2 <- coxph(
  Surv(time, status) ~ CD48_expr * NSD2 +
    relapsed + CCND1 + Gain1q21 + del17p13,
  data = cox_df
)
summary(cox_int_NSD2)

# — 2) Compare to the main‐effects‐only model
cox_main_NSD2 <- coxph(
  Surv(time, status) ~ CD48_expr + NSD2 + relapsed + CCND1 + Gain1q21 + del17p13,
  data = cox_df
)
anova(cox_main_NSD2, cox_int_NSD2, test = "Chisq")

# — 3) Extract CD48 HRs at NSD2 = 0 vs 1
coefs <- coef(cox_int_NSD2)
vc   <- vcov(cox_int_NSD2)
se    <- sqrt(diag(vc))

# CD48 HR when NSD2 = 0
hr0    <- exp(coefs["CD48_expr"])
ci0    <- exp(coefs["CD48_expr"] + c(-1,1) * 1.96 * se["CD48_expr"])

# CD48 HR when NSD2 = 1
beta_int <- coefs["CD48_expr:NSD21"]
# var(beta0 + beta_int) = Var(beta0) + Var(beta_int) + 2 Cov(beta0,beta_int)
var_sum  <- vc["CD48_expr","CD48_expr"] +
  vc["CD48_expr:NSD21","CD48_expr:NSD21"] +
  2 * vc["CD48_expr","CD48_expr:NSD21"]
se_sum   <- sqrt(var_sum)
hr1      <- exp(coefs["CD48_expr"] + beta_int)
ci1      <- exp((coefs["CD48_expr"] + beta_int) + c(-1,1) * 1.96 * se_sum)

tibble(
  NSD2    = c("0","1"),
  HR      = c(hr0, hr1),
  CI_lo   = c(ci0[1], ci1[1]),
  CI_hi   = c(ci0[2], ci1[2])
)

# — 4) Visualize with simulated “low”/“high” CD48 in each NSD2 group
med_cd48 <- median(cox_df$CD48_expr, na.rm = TRUE)
newdata_NSD2 <- expand.grid(
  CD48_expr = c(med_cd48 - 1, med_cd48 + 1),
  NSD2       = factor(c("0","1")),
  relapsed   = "0",  # set other covariates to reference
  CCND1      = "0",
  Gain1q21   = "0",
  del17p13   = "0"
) %>%
  mutate(group = paste(
    ifelse(CD48_expr < med_cd48, "CD48 low", "CD48 high"),
    ifelse(NSD2 == "0", "NSD2–", "NSD2+"),
    sep = ", "
  ))

fit_NSD2 <- survfit(cox_int_NSD2, newdata = newdata_NSD2)

ggsurvplot(
  fit_NSD2,
  data        = newdata_NSD2,
  legend.labs = newdata_NSD2$group,
  palette     = c("#1b9e77","#d95f02","#7570b3","#e7298a"),
  conf.int    = FALSE,
  xlab        = "Days",
  ylab        = "Survival Probability",
  title       = "Interaction: CD48 × NSD2 (t(4;14))"
)

## THREE WAY INTERACTION

library(survival)
library(survminer)
library(dplyr)

# 1) Fit the 3‐way interaction Cox model (if you haven’t already)
cox_int3 <- coxph(
  Surv(time, status) ~ CD48_expr * NSD2 * Gain1q21 +
    relapsed + CCND1 + del17p13,
  data = cox_df
)

# 2) Build a “newdata” grid with +/-1 SD (or 1 unit) around median CD48 and all 8 combos of NSD2/Gain1q21
med_cd48 <- median(cox_df$CD48_expr, na.rm = TRUE)
plotdata <- expand.grid(
  CD48_expr = c(med_cd48 - 1, med_cd48 + 1),
  NSD2       = factor(c("0","1")),
  Gain1q21   = factor(c("0","1")),
  relapsed   = factor("0"),
  CCND1      = factor("0"),
  del17p13   = factor("0")
) %>%
  mutate(
    CD48_lvl   = ifelse(CD48_expr < med_cd48, "LowCD48",  "HighCD48"),
    NSD2_lbl   = ifelse(NSD2   == "1",       "NSD2+",     "NSD2-"    ),
    Gain1q_lbl = ifelse(Gain1q21 == "1",     "Gain1q+",   "Gain1q-"  ),
    group      = paste(CD48_lvl, NSD2_lbl, Gain1q_lbl, sep = " ")
  )

# 3) Generate the survival curves at those settings
fit3 <- survfit(cox_int3, newdata = plotdata)

# 4) Plot with your custom legend labels
ggsurvplot(
  fit3,
  data        = plotdata,
  legend.labs = plotdata$group,
  palette     = RColorBrewer::brewer.pal(8, "Dark2"),
  conf.int    = FALSE,
  xlab        = "Days",
  ylab        = "Survival Probability",
  title       = "Interaction: CD48 × NSD2 × Gain1q21"
)

