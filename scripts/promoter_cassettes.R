#! usr/bin/Rscript

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("WGCNA")
}

library(WGCNA)


#
# Loading the data
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

#
# Preprocessing
#

# CpG context

# Getting CpGs belonging to promoter region
promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass=="promoter")]
promoter_betas <- betaNew[rownames(betaNew) %in% promoter_cpgs, ]

# Filtering based on variance

# Getting most variables CpGs
variance_prom <- sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])})

# Plotting variance
plot(density(variance_prom))
abline(v=0.05)

# Filtering data
selected_var <- sort(variance_prom, decreasing = T)[17725] # Using this to find an equivalent variance to the selected one in adjsuted betas
prom_to_analyse <- t(promoter_betas[variance_prom > selected_var,])


#
# Running WGCNA
#


# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# 
# # Call the network topology analysis function
# sft = pickSoftThreshold(
#   prom_to_analyse,             # <= Input data
#   powerVector = powers,
#   verbose = 5
# )
# 
# # Plotting
# par(mfrow = c(1,2))
# cex1 = 0.9
# 
# plot(sft$fitIndices[, 1],
#      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      xlab = "Soft Threshold (power)",
#      ylab = "Scale Free Topology Model Fit, signed R^2",
#      main = paste("Scale independence")
# )
# text(sft$fitIndices[, 1],
#      -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
#      labels = powers, cex = cex1, col = "red"
# )
# abline(h = 0.90, col = "red")
# plot(sft$fitIndices[, 1],
#      sft$fitIndices[, 5],
#      xlab = "Soft Threshold (power)",
#      ylab = "Mean Connectivity",
#      type = "n",
#      main = paste("Mean connectivity")
# )
# text(sft$fitIndices[, 1],
#      sft$fitIndices[, 5],
#      labels = powers,
#      cex = cex1, col = "red")

# Running WGCNA
betas <- c(5,8,10,15,20,25)
cor = WGCNA::cor

for (beta in betas) {

  # Calculating CpG casettes based on WGCNA
  
  netwk <- blockwiseModules(prom_to_analyse,               
                            corrType="spearman", # Using biweight midcorrelation 
                            nThreads = 10,
                            
                            # == Adjacency Function ==
                            power = beta,             
                            networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
                            
                            # == Tree and Block Options ==
                            deepSplit = 2,
                            pamRespectsDendro = F,
                            # detectCutHeight = 0.75,
                            minModuleSize = 3,
                            maxBlockSize = 6000,
                            
                            # == Module Adjustments ==
                            reassignThreshold = 0,
                            mergeCutHeight = 0.25,
                            
                            # == TOM == Archive the run results in TOM file (saves time)
                            saveTOMs = T,
                            saveTOMFileBase = "ER",
                            
                            # == Output Options
                            numericLabels = T,
                            verbose = 3)
  
  # Saving network
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_", beta, "_not_purity_adjusted.rds" )
  saveRDS(netwk, file = my_filename)
  
}
