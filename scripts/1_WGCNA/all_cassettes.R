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

# Filtering based on variance

# Getting most variables CpGs
variance_betas_unadj <- sapply(1:nrow(betaNew), FUN = function(row) {var(betaNew[row,])})
variance_betas_adj <- sapply(1:nrow(betaAdj), FUN = function(row) {var(betaAdj[row,])})

# Filtering data


# # Selecting variance thershold in unadjusted to keep the same number of CpG s analysed
#selected_var <- variance_betas_unadj[order(variance_betas_unadj, decreasing = TRUE)][sum(variance_betas_adj > 0.05)]

selected_var <- 0.05
cpgs_to_analyse <- t(betaNew[variance_betas_adj > selected_var,])

dim(cpgs_to_analyse)

# Converting to M values
eps <- 0.01
cpgs_to_analyse <- pmin(pmax(cpgs_to_analyse, eps), 1 - eps)

# Transform to M-values
cpgs_to_analyse <- log2(cpgs_to_analyse / (1 - cpgs_to_analyse))


#
# Running WGCNA
#

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  cpgs_to_analyse,
  powerVector = powers,
  verbose = 5, 
  networkType = "unsigned"
)

# Plotting
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (beta)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (beta)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")


# Running WGCNA
betas <- c(5,6,7,8,10,15)
cor = WGCNA::cor

for (beta in betas) {
  
  print(beta)
    
    # Calculating CpG cassettes based on WGCNA
    
    netwk <- blockwiseModules(cpgs_to_analyse,               
                              corType="pearson",
                              nThreads = 10,
                              
                              # == Adjacency Function ==
                              power = beta,             
                              networkType = "unsigned", 
                              
                              # == Tree and Block Options ==
                              deepSplit = 2,
                              pamRespectsDendro = T,
                              minModuleSize = 15,
                              maxBlockSize = 10000,
                              
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
    my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/all/unadjusted_cassettes_beta_", beta, ".rds" )
    saveRDS(netwk, file = my_filename)
    
  }
