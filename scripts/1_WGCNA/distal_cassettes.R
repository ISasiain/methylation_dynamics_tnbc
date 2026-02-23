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

# Getting distal CpGs

# All
distal_cpgs <- annoObj$illuminaID[which((annoObj$featureClass == "distal" | annoObj$featureClass == "distal body"))]

 # ATAC only
#distal_cpgs <- annoObj$illuminaID[which((annoObj$hasAtacOverlap == 1) & 
#                            (annoObj$featureClass == "distal" | annoObj$featureClass == "distal body"))]


distal_betas <- betaAdj[rownames(betaAdj) %in% distal_cpgs, ]

# Filtering based on variance

# Getting most variables CpGs
variance_dis <- sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])})

# Plotting variance
plot(density(variance_dis))
abline(v=0.05)

# Filtering data
#selected_var <- sort(variance_dis, decreasing = T)[7870] # Using this to find an equivalent variance to the selected one in non-adjusted betas
selected_var <- 0.05
dis_to_analyse <- t(distal_betas[variance_dis > selected_var,])

# Converting to M values
eps <- 0.01
dis_to_analyse <- pmin(pmax(dis_to_analyse, eps), 1 - eps)

# Transform to M-values
dis_to_analyse <- log2(dis_to_analyse / (1 - dis_to_analyse))

#
# Running WGCNA
#

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  dis_to_analyse,
  powerVector = powers,
  verbose = 5, 
  networkType = "unsigned")
  

# Plotting
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
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
     xlab = "Soft Threshold (power)",
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
    
    netwk <- blockwiseModules(dis_to_analyse,               
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
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/distal/atac_cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
}
