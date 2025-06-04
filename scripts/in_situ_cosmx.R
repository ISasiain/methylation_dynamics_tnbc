#! usr/bin/Rscript

library(dplyr)

#
# LOADING DATA
#

# CosMs data
oas2_per_cell_12 <-  read.csv("/Users/isasiain/PhD/Projects/project_3/data/block12_OAS2_in_tumour.csv")

# Pdid to cores mapping
core_to_pdid <- read.table("/Volumes/Data/CosMx/mapping.txt", sep = "\t", header = T)
core_to_pdid <- unique(core_to_pdid[,c("tmaID", "PDid")])
rownames(core_to_pdid) <- core_to_pdid$tmaID

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

#
# PLOT MEAN EXPRESSION OF ALL GENES (OAS2 + CONTROLS)
#

# Genes to plot
my_genes <- c("OAS2_Count", "BRCA1_Count", "SMO_Count", "GSTP1_Count" , "AR_Count")

# Reshape to long format
df_long <- oas2_per_cell_12 %>%
  filter(Cell_Type == "Tumour") %>%
  pivot_longer(cols = all_of(my_genes), names_to = "Gene", values_to = "Expression")

# Calculate mean expression per gene
mean_expr <- df_long %>%
  group_by(Gene) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE))

# Plot
ggplot(mean_expr, aes(x = Gene, y = Mean_Expression, fill = Gene)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(title = "Mean Expression per Tumour Cell",
       y = "Mean Expression", x = "Gene") +
  theme_minimal() +
  theme(legend.position = "none")


#
# DEFINING METHYLATION STATE
#

current_gene_id = "OAS2"

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                 annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])

cluster_promoter <- kmeans(t(betaAdj[names(cpgs)[cpgs=="promoter"],]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==1]) >
                           mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


#
# FILTERING AND SUMMARIZING TUMOR CELLS
#

### BLOCK 1 and 2

par(mfrow=c(1,1))

# Plotting panCK intensity vs cell type. Group 1 corresponds to tumour
boxplot(oas2_per_cell_12$Mean.PanCK ~ oas2_per_cell_12$Cell_Type)

# Remove non-tumour cells from df
oas2_per_cell_12 <- oas2_per_cell_12[oas2_per_cell_12$Cell_Type == "Tumour", ]


# Adding pdid
oas2_per_cell_12$"PDid" <- sapply(
  oas2_per_cell_12$Tissue_ID, function(tissue_name) {
    core_to_pdid[strsplit(tissue_name, "-")[[1]][1], "PDid"]
  }
)

# Methylation state
oas2_per_cell_12$"Methylation" <- sapply(
  oas2_per_cell_12$PDid, function(pdid) {
    promoter_state[pdid]
  }
)

# SUMMARIZING

# Summarizing per Core (Tissue ID)
summary_per_tissue_12 <- oas2_per_cell_12 %>%
  group_by(Tissue_ID) %>%
  summarise(
    PDid = first(PDid),
    Methylation = first(Methylation),
    
    mean_OAS2 = mean(OAS2_Count, na.rm = TRUE),
    prop_OAS2 = mean(OAS2_Count > 0, na.rm = TRUE),
    
    mean_BRCA1 = mean(BRCA1_Count, na.rm = TRUE),
    prop_BRCA1 = mean(BRCA1_Count > 0, na.rm = TRUE),
    
    mean_SMO = mean(SMO_Count, na.rm = TRUE),
    prop_SMO = mean(SMO_Count > 0, na.rm = TRUE),
    
    mean_GSTP1 = mean(GSTP1_Count, na.rm = TRUE),
    prop_GSTP1 = mean(GSTP1_Count > 0, na.rm = TRUE),
    
    mean_AR = mean(AR_Count, na.rm = TRUE),
    prop_AR = mean(AR_Count > 0, na.rm = TRUE)
  )


# Adding PAM5O annotations
summary_per_tissue_12$PAM50 <- x[summary_per_tissue_12$PDid, "PAM50_Basal_NCN"]


### PLOTTING

## BLOCK 1 and 2

# Mean OAS2 expression

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$Methylation)
counts_samples <- c("Hypermethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypermethylated","PDid"]))),
                    "Hypomethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypomethylated","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)], 
  "\ns=", counts_samples[names(counts_cores)]
)


# Draw boxplot with custom x-axis labels
boxplot(summary_per_tissue_12$mean_OAS2 ~ summary_per_tissue_12$Methylation,
        ylab = "Mean Expression in Tumour cells",
        xlab = "Promoter Methylation",,
        ylim = c(0, 1.5),
        frame = FALSE)


# Add the annotation below the x-axis
text(x = 1:2, y = 1.35, labels = labels_with_counts, xpd = TRUE, cex = 0.8)


# Add jittered points
stripchart(summary_per_tissue_12$mean_OAS2 ~ summary_per_tissue_12$Methylation,
           method = "jitter", 
           pch = 16,
           cex = 0.6,
           col = rgb(0, 0, 0, 0.5),
           vertical = TRUE,
           add = TRUE)

# Perform Wilcoxon test
wilcox_res <- wilcox.test(summary_per_tissue_12$mean_OAS2 ~ summary_per_tissue_12$Methylation)

# Add p-value to plot
pval <- wilcox_res$p.value
text(x = 1.1, 
     y = 1, 
     labels = paste0("Mann-Whitney p = ", signif(pval, 3)),
     pos = 3, cex = 0.9)


# Proportion of cells with detected expressin of OAS2

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$Methylation)
counts_samples <- c("Hypermethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypermethylated","PDid"]))),
                    "Hypomethylated" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$Methylation == "Hypomethylated","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)], 
  "\ns=", counts_samples[names(counts_cores)]
)

# Draw boxplot with custom x-axis labels
boxplot(summary_per_tissue_12$prop_OAS2 ~ summary_per_tissue_12$Methylation,
        ylab = "Proportion of Tumour cells expressing",
        ylim=c(0, 0.4),
        xlab = "Promoter Methylation",
        frame = FALSE)

# Add the annotation below the x-axis
text(x = 1:2, y = 0.35, labels = labels_with_counts, xpd = TRUE, cex = 0.8)


# Add jittered points
stripchart(summary_per_tissue_12$prop_OAS2 ~ summary_per_tissue_12$Methylation,
           method = "jitter", 
           pch = 16,
           cex = 0.6,
           col = rgb(0, 0, 0, 0.5),
           vertical = TRUE,
           add = TRUE)

# Perform Wilcoxon test
wilcox_res <- wilcox.test(summary_per_tissue_12$prop_OAS2 ~ summary_per_tissue_12$Methylation)

# Add p-value to plot
pval <- wilcox_res$p.value
text(x = 1.1, 
     y = 0.25, 
     labels = paste0("Mann-Whitney p = ", signif(pval, 3)),
     pos = 3, cex = 0.9)



### CONTROL. Plotting AR vs Basal/NonBasal

# Mean AR expression

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$PAM50)
counts_samples <- c("Basal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "Basal","PDid"]))),
                    "nonBasal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "nonBasal","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)], 
  "\ns=", counts_samples[names(counts_cores)]
)


# Draw boxplot with custom x-axis labels
boxplot(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50,
        ylab = "Mean Expression in Tumour cells",
        xlab = "Promoter Methylation",,
        ylim = c(0, 2),
        frame = FALSE)

# Add the annotation below the x-axis
text(x = 1:2, y = 1.7, labels = labels_with_counts, xpd = TRUE, cex = 0.8)


# Add jittered points
stripchart(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50,
           method = "jitter", 
           pch = 16,
           cex = 0.6,
           col = rgb(0, 0, 0, 0.5),
           vertical = TRUE,
           add = TRUE)

# Perform Wilcoxon test
wilcox_res <- wilcox.test(summary_per_tissue_12$mean_AR ~ summary_per_tissue_12$PAM50)

# Add p-value to plot
pval <- wilcox_res$p.value
text(x = 1.1, 
     y = 1, 
     labels = paste0("Mann-Whitney p = ", signif(pval, 3)),
     pos = 3, cex = 0.9)


# Proportion of cells with detected expressin of OAS2

# Count number of data points per class
counts_cores <- table(summary_per_tissue_12$PAM50)
counts_samples <- c("Basal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "Basal","PDid"]))),
                    "nonBasal" = nrow(na.omit(unique(summary_per_tissue_12[summary_per_tissue_12$PAM50 == "nonBasal","PDid"]))))


labels_with_counts <- paste0(
  "n=", counts_cores[names(counts_cores)], 
  "\ns=", counts_samples[names(counts_cores)]
)


# Draw boxplot with custom x-axis labels
boxplot(summary_per_tissue_12$prop_AR ~ summary_per_tissue_12$PAM50,
        ylab = "Proportion of Tumour cells expressing",
        xlab = "Promoter Methylation",,
        ylim = c(0, 0.4),
        frame = FALSE)

# Add the annotation below the x-axis
text(x = 1:2, y = 0.35, labels = labels_with_counts, xpd = TRUE, cex = 0.8)


# Add jittered points
stripchart(summary_per_tissue_12$prop_AR ~ summary_per_tissue_12$PAM50,
           method = "jitter", 
           pch = 16,
           cex = 0.6,
           col = rgb(0, 0, 0, 0.5),
           vertical = TRUE,
           add = TRUE)

# Perform Wilcoxon test
wilcox_res <- wilcox.test(summary_per_tissue_12$prop_AR ~ summary_per_tissue_12$PAM50)

# Add p-value to plot
pval <- wilcox_res$p.value
text(x = 1.1, 
     y = 0.25, 
     labels = paste0("Mann-Whitney p = ", signif(pval, 3)),
     pos = 3, cex = 0.9)
