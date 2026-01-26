#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(dplyr)
library(patchwork)

#
# LOADING DATA
#

# Loading promoter cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")


# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]


# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID


# Loading data (GSE67919)
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Beta.RData")
load("/Volumes/Data/Project_3/normal_breast_methylation/GSE67919/GSE67919_Annotations.RData")

normal_tissue_methylation_96 <- beta

# Loading data (GSE225845)

epic_normal_betas <- readRDS("/Users/in2245sa/PhD/Projects/project_3/data/GSE225845_normal_epic/GSE225845_normal_normalized_betas.rds")

# Formatting data
epic_normal_betas_rownames <- rownames(epic_normal_betas)
epic_normal_betas <- as.data.frame(epic_normal_betas)
rownames(epic_normal_betas) <- epic_normal_betas_rownames


#
# PLOTTING GSE67919
#

# Prepare long-format data
selected_genes <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")

plot_data <- lapply(selected_genes, function(gene) {
  cpgs <- names(genes)[genes == gene]
  cpgs <- cpgs[cpgs %in% rownames(normal_tissue_methylation_96)]
  data <- normal_tissue_methylation_96[cpgs, , drop = FALSE]
  df <- as.data.frame(t(data))
  df$Sample <- rownames(df)
  df_long <- pivot_longer(df, -Sample, names_to = "CpG", values_to = "Beta")
  df_long$Gene <- gene
  df_long
}) %>% bind_rows()

# Removing cpgs that were not in the detected cassette.
plot_data <- plot_data[plot_data$CpG %in% names(promoter_10$colors)[promoter_10$colors == 10],]

# Plot
ggplot(plot_data, aes(x = CpG, y = Beta)) +
  geom_violin(fill="black") +
  geom_boxplot(fill = "grey90", col= "grey50", outlier.size = 0.5, width=0.2) +
  facet_grid(.~ Gene,scales = "free", space='free') +
  theme_bw(base_size = 12) +
  ylim(0,1) +
  ylab("Beta Value") +
  xlab("CpG Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#
# PLOTTING GSE225845
#
  
normal_tissue_methylation_113 <- t(epic_normal_betas)

# Prepare long-format data
selected_genes <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")

plot_data <- lapply(selected_genes, function(gene) {
  cpgs <- names(genes)[genes == gene]
  cpgs <- cpgs[cpgs %in% rownames(normal_tissue_methylation_113)]
  data <- normal_tissue_methylation_113[cpgs, , drop = FALSE]
  df <- as.data.frame(t(data))
  df$Sample <- rownames(df)
  df_long <- pivot_longer(df, -Sample, names_to = "CpG", values_to = "Beta")
  df_long$Gene <- gene
  df_long
}) %>% bind_rows()

# Removing cpgs that were not in the detected cassette.
plot_data <- plot_data[plot_data$CpG %in% names(promoter_10$colors)[promoter_10$colors == 10],]

# Plot
ggplot(plot_data, aes(x = CpG, y = Beta)) +
  geom_violin(fill="black") +
  geom_boxplot(fill = "grey90", col= "grey50", outlier.size = 0.5, width=0.2) +
  facet_grid(.~ Gene,scales = "free", space='free') +
  theme_bw(base_size = 12) +
  ylim(0,1) +
  ylab("Beta Value") +
  xlab("CpG Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
