#! usr/bin/Rscript

library(ComplexHeatmap)
library(survival)
library(ggplot2)
library(patchwork)

# Load data
dist_summary_7 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_7.csv")
rownames(dist_summary_7) <- dist_summary_7$Cassette
dist_summary_7$Cassette <- NULL

dist_cassette_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_7.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]


#
# PLOTTING HEATMP OF CASSETTE 1 + PC1
# 

# Define column annotation
pc1_annotation <- HeatmapAnnotation(
  PC1 = as.numeric(dist_summary_7["1", colnames(betaAdj)]),
  col = list(
    PC1 = colorRamp2(
      c(min(as.numeric(dist_summary_7["1", colnames(betaAdj)])), 0, max(as.numeric(dist_summary_7["1", colnames(betaAdj)]))),
      c("#A6C8E4", "#FFFFFF", "#F4A6A6")  # light red → white → light blue
    )
  )
)

# Plotting heatmap
Heatmap(betaAdj[names(dist_cassette_7$colors)[dist_cassette_7$colors == 1],], 
        col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
        use_raster = F,
        show_row_dend = F, 
        show_column_names = F, 
        show_row_names = F,
        bottom_annotation = pc1_annotation,
        name="Tumor beta")


# Change due to purity adjustment

# BETA VALUE
plot(c(betaAdj[names(dist_cassette_7$colors)[dist_cassette_7$colors == 1],]),
     c(betaNew[names(dist_cassette_7$colors)[dist_cassette_7$colors == 1],]),
     pch=16,
     cex=0.1,
     xlab="Adjusted beta",
     ylab="Original beta")


# VARIANCE
var_adj <- apply(betaAdj[names(dist_cassette_7$colors)[dist_cassette_7$colors == 1],], MARGIN = 1, FUN = var)
var_orig <- apply(betaNew[names(dist_cassette_7$colors)[dist_cassette_7$colors == 1],], MARGIN = 1, FUN = var)

plot(var_adj,
     var_orig,
     pch=16,
     cex=0.5,
     xlab="Adjusted beta",
     ylab="Original beta")

#
# PLOTTING PC1 + BIOLOGICAL FEATURES
# 

# Prepare the data
plot_df <- data.frame(
  PC1_Cassette1 = as.numeric(dist_summary_7["1", colnames(betaAdj)]),
  PAM50_Basal_NCN = x$PAM50_Basal_NCN,
  PAM50_NCN = x$PAM50_NCN,
  HRD.2.status = x$HRD.2.status,
  HRD.3.status = x$HRD.3.status,
  TNBCtype4_n235_notPreCentered = x$TNBCtype4_n235_notPreCentered,
  ASCAT_TUM_FRAC = x$ASCAT_TUM_FRAC,
  TILs = x$TILs,
  IM = x$TNBCtype_IMpositive
)

# Clean the IM variable
plot_df$IM <- sapply(plot_df$IM, function(x) {
  if (is.na(x)) { NA } else if (x == 0) { "Negative" } else if (x == 1) { "Positive" }
})

# Plot 1: PC1 vs PAM50
p1 <- ggplot(plot_df, aes(x = PAM50_Basal_NCN, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "PAM50 Basal", y = "PC1 Distal Cassette 1") +
  theme(axis.text = element_text(size = 12))

pval_pam50 <-  wilcox.test(as.numeric(plot_df$PC1_Cassette1) ~ plot_df$PAM50_Basal_NCN)$p.value

# Plot 2: PC1 vs HRD.2.status
p2 <- ggplot(plot_df, aes(x = HRD.2.status, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "HRD Status", y = NULL) +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pval_hrd <- wilcox.test(plot_df$PC1_Cassette1 ~ plot_df$HRD.2.status)$p.value

# Plot 3: PC1 vs Lehmann
p3 <- ggplot(plot_lehmann <- plot_df[!is.na(plot_df$TNBCtype4_n235_notPreCentered), ], 
             aes(x = TNBCtype4_n235_notPreCentered, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "Lehmann 5", y = NULL) +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

pval_lehmann <- kruskal.test(plot_lehmann$PC1_Cassette1 ~ plot_lehmann$TNBCtype4_n235_notPreCentered)$p.value

# PC1 VS ASCAT

corr_val_ascat <- cor(plot_ASCAT$PC1_Cassette1, plot_ASCAT$ASCAT_TUM_FRAC, method="spearman")
p_val_ascat <- cor.test(plot_ASCAT$PC1_Cassette1, plot_ASCAT$ASCAT_TUM_FRAC, method="spearman", exact = FALSE)$p.value

p5 <- ggplot(plot_ASCAT <- plot_df[!is.na(plot_df$ASCAT_TUM_FRAC), ], 
             aes(x = ASCAT_TUM_FRAC, y = PC1_Cassette1)) +
  geom_point(
    color = "black",
    cex = 0.8
  ) +
  theme_bw() + 
  labs(x = "ASCAT Purity", y = "PC1 Cassette 1") +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", 
           x = 0.35 * max(plot_df$ASCAT_TUM_FRAC, na.rm = TRUE), 
           y = min(plot_df$PC1_Cassette1, na.rm = TRUE), 
           label = paste0("Cor. = ", round(corr_val_ascat, 3)), 
           hjust = 0, vjust = 0, size = 5)

# PC1 VS TILs

corr_val_tils <- cor(plot_tils$PC1_Cassette1, plot_tils$TILs, method="spearman")
p_val_tils <- cor.test(plot_tils$PC1_Cassette1, plot_tils$TILs, method="spearman", exact = FALSE)$p.value

p6 <- ggplot(plot_tils <- plot_df[!is.na(plot_df$TILs), ], aes(x = TILs, y = PC1_Cassette1)) +
  geom_point(
    color = "black",
    cex = 0.8
  ) +
  theme_bw() + 
  labs(x = "TILs", y = "") +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", 
            x = 0.3 * max(plot_df$TILs, na.rm = TRUE), 
           y = min(plot_df$PC1_Cassette1, na.rm = TRUE), 
           label = paste0("Cor. = ", round(corr_val_tils, 3)), 
           hjust = 0, vjust = 0, size = 5)



# Combine plots 
combined_plot <- p1 + p2 + p3 + p5 + p6 + plot_layout(ncol = 5, guides = "collect")

combined_plot

# Adjusting p values for multiple testing
p.adjust(c(pval_pam50, pval_hrd, pval_lehmann, p_val_ascat, p_val_tils), method = "BH")

