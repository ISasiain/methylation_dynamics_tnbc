#! usr/bin/Rscript

# Install and load packages
#install.packages("Downloads/MethylCIBERSORT_Release_0.2.1/MethylCIBERSORT_0.2.1.tar.gz", repos = NULL, type = "source")
#BiocManager::install("EpiDISH")
# devtools::install_github("Moonerss/CIBERSORT")

setwd("/Users/isasiain/PhD/Projects/project_3/analysis/cell_type_deconvolution")

library(MethylCIBERSORT)
library(EpiDISH)
library(VIM)


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


# RUNNING EPIDISH IN NORMAL TISSUE

cell_types_normal_tissue <- epidish(beta.m = normal_tissue_methylation_96,
                                    ref.m = centEpiFibFatIC.m,
                                    method = "RPC")




# Boxplot with improvements
boxplot(cell_types_normal_tissue$estF, 
        ylab = "Estimated cell fraction", 
        ylim=c(0,1),
        col = "lightblue", 
        border = "darkblue", 
        notch = FALSE, 
        outline = TRUE, 
        horizontal = FALSE, 
        pch = 16, 
        cex = 1.5)


# Barplot with improvements
# Barplot with improvements and a legend
barplot(t(cell_types_normal_tissue$estF), 
        ylab = "Value", 
        col = c("grey", "blue", "yellow", "red"), 
        border = "black", 
        cex.names = 0.8, 
        las = 2)  # Rotate the axis labels

# Add a legend
legend("right", 
       legend = c("Epithelial", "Fibroblasts", "Fatty cells", "Immune cells"), 
       fill = c("grey", "blue", "yellow", "red"),
       cex = 0.8)



# Step 1: Define X = cell type fractions (design matrix)

gbp4_cpgs <- names(genes)[genes == "GBP4"]
gbp4_cpgs <- gbp4_cpgs[gbp4_cpgs %in% rownames(normal_tissue_methylation_96)]
data_subset <- normal_tissue_methylation_96[gbp4_cpgs,]

par(mfrow = c(1, 2))

# Loop over each CpG (assuming you have 6 rows in data_subset)
for (i in 1:length(gbp4_cpgs)) {
  
  # Get CpG id
  cpg <- rownames(data_subset)[i]
  
  # Getting data for the ith CpG
  Y <- as.numeric(data_subset[i, colnames(data_subset)])  # Methylation values for CpG i
  X <- cell_types_normal_tissue$estF[colnames(data_subset), c("Epi", "Fib", "Fat", "IC")]  # Cell type fractions
  
  # Step 3: Fit linear model and estimate coefficients (no intercept)
  fit <- lm(Y ~ 0 + as.matrix(X), na.action = na.exclude)  # Fit the linear model with no intercept
  
  # Predicted vs observed plot
  plot(Y, predict(fit), 
       main = paste(cpg),
       xlab = "Observed Methylation",
       ylab = "Predicted Methylation",
       pch = 16, cex = 0.7, col = "blue")
  abline(0, 1, col = "red", lty = 2)  # Add a reference line (y = x)
  
  # Step 4: Barplot of coefficients for the current CpG
  coef_values <- coef(fit)  # Get the coefficients for the linear model
  
  barplot(coef_values, 
          main = paste(cpg),
          ylab = "Coefficient Value",
          xlab = "Cell Types",
          names.arg = c("Epi", "Fib", "Fat", "IC"),
          col = c("lightblue", "lightgreen", "khaki", "lightpink"),
          las = 1)  # Rotate labels on x-axis
  
  readline()
  
}

par(mfrow = c(1, 1))


#
# CANCER TISSUE
#

# Estimating cell types in tumour
cell_types_cancer_tissue <- epidish(beta.m = betaNew,
                                    ref.m = centEpiFibFatIC.m,
                                    method = "RPC")
par(mfrow = c(1, 1), mar=c(4,4,4,4))

# Boxplot with improvements
boxplot(cell_types_cancer_tissue$estF, 
        ylab = "Estimated cell fraction", 
        ylim=c(0,1),
        col = "lightblue", 
        border = "darkblue", 
        notch = FALSE, 
        outline = TRUE, 
        horizontal = FALSE, 
        pch = 16, 
        cex = 1.5)

barplot(t(cell_types_cancer_tissue$estF), 
        ylab = "Value", 
        col = c("grey", "blue", "yellow", "red"), 
        border = "black", 
        cex.names = 0.8, 
        las = 2)  # Rotate the axis labels

# Add a legend
legend("bottomleft", 
       legend = c("Epithelial", "Fibroblasts", "Fatty cells", "Immune cells"), 
       fill = c("grey", "blue", "yellow", "red"),
       cex = 0.8)



#
# METHYLCIBERSORT
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

# RUNNING METHYLCIBERSORT IN TUMOUR TISSUE

# Filtering signatures
Int <- intersect(rownames(betaNew), rownames(Stromal_v2))
betaNew <- betaNew[match(Int, rownames(betaNew)),]
Stromal_v2 <- Stromal_v2[match(Int, rownames(Stromal_v2)),]

# Redefine signatures
RefData <- Stromal_v2
RefPheno <- Stromal_v2.pheno

Signature <- FeatureSelect.V4(CellLines.matrix = NULL,
                              Heatmap = FALSE,
                              export = TRUE,
                              sigName = "MyReference_TUM",
                              Stroma.matrix = RefData,
                              deltaBeta = 0.2,
                              FDR = 0.01,
                              MaxDMRs = 100,
                              Phenotype.stroma = RefPheno)


Prep.CancerType(Beta = normal_tissue_methylation_96, Probes = rownames(Signature$SignatureMatrix), fname = "MixtureMatrix_TUM")
