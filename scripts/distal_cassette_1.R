#! usr/bin/Rscript

library(ComplexHeatmap)
library(survival)
library(ggplot2)
library(patchwork)

# Load data
dist_summary_10 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(dist_summary_10) <- dist_summary_15$Cassette
dist_summary_10$Cassette <- NULL

dist_cassette_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]


#
# PLOTTING HEATMP OF CASSETTE 1 + PC1
# 

# Define column annotation
pc1_annotation <- HeatmapAnnotation(
  "PC1" = as.numeric(dist_summary_10["1", colnames(betaAdj)])
)

# Plotting heatmap
Heatmap(betaAdj[names(dist_cassette_10$colors)[dist_cassette_10$colors == 1],], 
        use_raster = F,
        show_row_dend = F, 
        show_column_names = F, 
        show_row_names = F,
        bottom_annotation = pc1_annotation,
        name="Tumor beta")

#
# PLOTTING HEATMP OF CASSETTE 1 + PC1
# 

# Prepare the data
plot_df <- data.frame(
  PC1_Cassette1 = as.numeric(dist_summary_10["1", colnames(betaAdj)]),
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

# Plot 2: PC1 vs HRD.2.status
p2 <- ggplot(plot_df, aes(x = HRD.2.status, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "HRD Status", y = NULL) +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Plot 3: PC1 vs Lehmann
p3 <- ggplot(filter(plot_df, !is.na(TNBCtype4_n235_notPreCentered)), 
             aes(x = TNBCtype4_n235_notPreCentered, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "Lehmann 5", y = NULL) +
  theme(axis.text = element_text(size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# Plot 4: PC1 vs IM
p4 <- ggplot(filter(plot_df, !is.na(IM)), aes(x = IM, y = PC1_Cassette1)) +
  geom_boxplot(color = "black", fill = "lightgrey", outlier.shape = 16, outlier.color = "red", notch = FALSE) +
  theme_bw() + 
  labs(x = "IM Status", y = NULL) +
  theme(axis.text = element_text(size = 12), 
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())

# PC1 VS ASCAT

corr_val <- cor(as.numeric(dist_summary_10["1", colnames(betaAdj)]), x$ASCAT_TUM_FRAC, method="spearman")

p5 <- ggplot(filter(plot_df, !is.na(plot_df$ASCAT_TUM_FRAC)), aes(x = ASCAT_TUM_FRAC, y = PC1_Cassette1)) +
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
           label = paste0("Cor. = ", round(corr_val, 2)), 
           hjust = 0, vjust = 0, size = 5)


# PC1 VS TILs
corr_val <- cor(as.numeric(dist_summary_10["1", colnames(betaAdj)]), x$TILs, method="spearman", use = "complete.obs")

p6 <- ggplot(filter(plot_df, !is.na(plot_df$TILs)), aes(x = TILs, y = PC1_Cassette1)) +
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
           label = paste0("Cor. = ", round(corr_val, 2)), 
           hjust = 0, vjust = 0, size = 5)


# Combine plots 
combined_plot <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 6, guides = "collect")

combined_plot

# PC1 and survival. NO EFFECT
cox_model <- coxph(Surv(x$IDFS, x$IDFSbin) ~ as.numeric(dist_summary_10["1", colnames(betaAdj)]))
summary(cox_model)

cox_model <- coxph(Surv(x$OS, x$OSbin) ~ as.numeric(dist_summary_10["1", colnames(betaAdj)]))
summary(cox_model)


# Kaplan meier plot of based on methylation state

# Create the survival object
surv_obj <- Surv(time = IDFS, event = IDFS_bin)

# Fit Kaplan-Meier survival curves stratified by promoter_state
fit <- survfit(surv_obj ~ promoter_state)

# Plot Kaplan-Meier curves
ggsurvplot(
  fit,
  data = data.frame(IDFS = IDFS, IDFS_bin = IDFS_bin, promoter_state = promoter_state),
  pval = TRUE,                  # Add p-value from log-rank test
  conf.int = TRUE,              # Add confidence interval
  risk.table = TRUE,            # Show number at risk table
  legend.title = "Promoter State",
  xlab = "Time (IDFS)",
  ylab = "Survival Probability",
  palette = "Dark2"             
)


