#! usr/bin/Rscript

library(ComplexHeatmap)
library(survival)

# Load data
dist_summary_15 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_15.csv")
rownames(dist_summary_15) <- dist_summary_15$Cassette
dist_summary_15$Cassette <- NULL

dist_cassette_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")

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
  "PC1" = as.numeric(dist_summary_15["1", colnames(betaAdj)])
)

# Plotting heatmap
Heatmap(betaAdj[names(dist_cassette_15$colors)[dist_cassette_15$colors == 1],], 
        use_raster = F,
        show_row_dend = F, 
        show_column_names = F, 
        show_row_names = F,
        bottom_annotation = pc1_annotation)

#
# PLOTTING HEATMP OF CASSETTE 1 + PC1
# 

# PC1 VS PAM50 2
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$PAM50_Basal_NCN,
        xlab = "",
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across PAM50 Subtypes",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS PAM50 6
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ ifelse(x$PAM50_NCN == "unclassified", "Uncl.", x$PAM50_NCN),
        xlab = NULL,
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 2,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across PAM50 Subtypes",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS HRD 2
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$HRD.2.status,
        xlab = NULL,
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across HRD 2",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS HRD 3
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$HRD.3.status,
        xlab = NULL,
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across HRD 3",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS LEHMANN 5
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$TNBCtype4_n235_notPreCentered,
        xlab = NULL,
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across Lehmann 5",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS IM
boxplot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ ifelse(x$TNBCtype_IMpositive == 1, "IM +", "IM -"),
        xlab = NULL,
        ylab = "PC1 Cassette 1",
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across IM",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)

# PC1 VS IM
plot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$ASCAT_TUM_FRAC,
     pch = 16,              # Solid circles for points
     xlab = "ASCAT Tumor Fraction",  # Add clear axis labels
     ylab = "PC1 Cassette 1",        
     main = "PC1 Cassette 1 vs. ASCAT Tumor Fraction", # Add a title
     cex = 1.2,             # Increase point size
     cex.axis = 1.3,        # Increase axis text size
     cex.lab = 1.4,         # Increase label text size
     frame = FALSE)         # Remove default box

# PC1 VS TILs
plot(as.numeric(dist_summary_15["1", colnames(betaAdj)]) ~ x$TILs,
     pch = 16,              # Solid circles for points
     xlab = "TILs",  # Add clear axis labels
     ylab = "PC1 Cassette 1",        
     main = "PC1 Cassette 1 vs. TILs", # Add a title
     cex = 1.2,             # Increase point size
     cex.axis = 1.3,        # Increase axis text size
     cex.lab = 1.4,         # Increase label text size
     frame = FALSE)         # Remove default box

# PC1 and survival. NO EFFECT
cox_model <- coxph(Surv(x$IDFS, x$IDFSbin) ~ as.numeric(dist_summary_15["1", colnames(betaAdj)]))
summary(cox_model)










