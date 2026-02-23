#! usr/bin/Rscript

# Install and load packages
#install.packages("Downloads/MethylCIBERSORT_Release_0.2.1/MethylCIBERSORT_0.2.1.tar.gz", repos = NULL, type = "source")
#BiocManager::install("EpiDISH")
# devtools::install_github("Moonerss/CIBERSORT")

setwd("/Users/isasiain/PhD/Projects/project_3/analysis/cell_type_deconvolution")

library(MethylCIBERSORT) # v0.2.1
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(patchwork)

#
# LOADING DATA
#

# Loading corrected and uncorrected betas and annotations
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

#
# NORMAL TISSUE
#

# Loading data
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Beta.RData")
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Annotations.RData")

normal_tissue_methylation_96 <- beta

# Imputation (median)
normal_tissue_methylation_96_imputed <- beta
na_cols <- colnames(beta)[colSums(is.na(beta)) > 0]

for (col in na_cols) {
  med <- median(beta[, col], na.rm = TRUE)
  missing_idx <- is.na(beta[, col])
  normal_tissue_methylation_96_imputed[missing_idx, col] <- med
}


#
# RUNNING METHYLCIBERSORT
#

# RUNNING METHYLCIBERSORT IN NORMAL TISSUE

# Filtering signatures
Int <- intersect(rownames(normal_tissue_methylation_96_imputed), rownames(Stromal_v2))
normal_tissue_methylation_96_imputed <- normal_tissue_methylation_96_imputed[match(Int, rownames(normal_tissue_methylation_96_imputed)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

# Redefine signatures
RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno

Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = "MyReference_NOR",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)


Prep.CancerType(Beta = normal_tissue_methylation_96, Probes = rownames(Signature$SignatureMatrix), fname = "MixtureMatrix_NOR")

### ... Running deconvolution in CIBERSORTx ...


# RUNNING METHYLCIBERSORT IN TUMOUR TISSUE

# Loading reference signatures
data("V2_Signatures")

# Getting brca signature
signature_brca <- Signatures$breast_v2_Signature.txt

# Saving brca methylartion mixture
Prep.CancerType(Beta = betaNew, 
                Probes = rownames(Signature$SignatureMatrix), 
                fname = "MixtureMatrix_TUM")

Prep.CancerType(Beta = betaAdj, 
                Probes = rownames(Signature$SignatureMatrix), 
                fname = "MixtureMatrix_TUM_adjusted")

# Saving signature
write.table(signature_brca, "Signature_TUM.txt", sep = "\t", row.names = FALSE, quote = FALSE )


### ... Running deconvolution in CIBERSORTx ...


#
# ANALYSISNG DECONVOLUTION OUTPUT
#

# Reading output
tum_cibersort <- read.csv("CIBERSORTx_tumour.csv")
tum_cibersort_adjusted <- read.csv("CIBERSORTx_tumour_adjusted.csv")

rownames(tum_cibersort) <- tum_cibersort$Mixture

tum_cibersort$SampleID <- rownames(tum_cibersort)
tum_cibersort$Mixture <- NULL
tum_cibersort$P.value <- NULL
tum_cibersort$Correlation <- NULL
tum_cibersort$RMSE <- NULL


rownames(tum_cibersort_adjusted) <- tum_cibersort_adjusted$Mixture

tum_cibersort_adjusted$SampleID <- rownames(tum_cibersort_adjusted)
tum_cibersort_adjusted$Mixture <- NULL
tum_cibersort_adjusted$P.value <- NULL
tum_cibersort_adjusted$Correlation <- NULL
tum_cibersort_adjusted$RMSE <- NULL

# Boxplots of tumour composition before and after PureBeta

tum_cibersort$SampleID <- rownames(tum_cibersort)
tum_cibersort_adjusted$SampleID <- rownames(tum_cibersort_adjusted)

# Long format for each dataframe
long_before <- tum_cibersort %>%
  pivot_longer(-SampleID, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Adjustment = "Before\nPureBeta")

long_after <- tum_cibersort_adjusted %>%
  pivot_longer(-SampleID, names_to = "CellType", values_to = "Proportion") %>%
  mutate(Adjustment = "After\nPureBeta")


# Plot 1: Original proportions. Boxplot
ggplot(long_before, aes(x = CellType, y = Proportion, fill = CellType)) +
  geom_boxplot(fill = "gray90") +
  geom_jitter(pch = 16, cex = 0.3, width = 0.2) +
  theme_bw() +
  ylab("Cell Fraction") +
  xlab(NULL) +
  theme(
    legend.position = "none",              # Removes the legend
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotates x-axis labels
  )



# Combine
long_combined <- bind_rows(long_after, long_before)

# Group non-malignnat
long_combined <- long_combined %>%
  mutate(Group = if_else(CellType == "Cancer", "Malignant", "Non-malignant"))

summarized <- long_combined %>%
  group_by(SampleID, Adjustment, Group) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop")

# Setting factor order 
summarized$Adjustment <- factor(summarized$Adjustment, levels = c("Before\nPureBeta", "After\nPureBeta"))

# Defining plotting functions
plot_boxplot_group <- function(df, group_name) {
  ggplot(df %>% filter(Group == group_name), aes(x = Adjustment, y = Proportion)) +
    geom_violin(aes(group = Adjustment), fill = "black") +
    geom_boxplot(aes(group = Adjustment), outlier.shape = NA, fill = "gray90",  col="gray60", width=0.25) +
    labs(y = "Proportion", x = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none")
}

plot_change_group <- function(df, group_name) {
  
  df_delta <- df %>% filter(Group == group_name) %>%
    pivot_wider(names_from = Adjustment, values_from = Proportion) %>%
    mutate(delta = `After\nPureBeta` - `Before\nPureBeta`)
  
  ggplot(df_delta, aes(x = "", y = delta)) +
    geom_quasirandom(size = 1, alpha = 0.5, color = "steelblue") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Δ Proportion", x = NULL) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}


# Malignant
p_malig_box <- plot_boxplot_group(summarized, "Malignant")
p_malig_change <- plot_change_group(summarized, "Malignant")

# Non-malignant
p_nonmalig_box <- plot_boxplot_group(summarized, "Non-malignant")
p_nonmalig_change <- plot_change_group(summarized, "Non-malignant")

p_malig_box | p_malig_change

p_nonmalig_box | p_nonmalig_change
