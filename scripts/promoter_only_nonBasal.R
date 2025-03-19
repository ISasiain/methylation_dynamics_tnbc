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

# Distal cassettes
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")

# Summary distal cassettes
summary_prom15 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_15.csv")
rownames(summary_prom15) <- as.character(summary_prom15$Cassette)
summary_prom15$Cassette <- NULL

# Getting promoter CpGs and betas
promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass=="promoter")]
promoter_betas <- betaAdj[rownames(betaAdj) %in% promoter_cpgs, ]

#
# Preprocessing
#

# GETTING ONLY BASAL SAMPLES

# Getting cassettes linked to basal/nonBasal split
my_cpgs_dis <-  c(
  names(promoter_15$colors)[promoter_15$colors == "1"]
)

# Geenerate data frame to store groups
groupings_df <- data.frame(matrix(nrow = length(colnames(summary_dis15)), ncol = 1))
rownames(groupings_df) <- colnames(summary_dis15)
colnames(groupings_df) <- "group_prom"


# Clustering in two groups
distance_matrix <- dist(t(betaAdj[my_cpgs_dis,]))
hc <- hclust(distance_matrix)
groupings_df$group_prom <- cutree(hc, k = 2)

# Getting basal group
promoter_betas <- promoter_betas[,rownames(groupings_df)[groupings_df==2]]


# FILTERING BASED ON VARIANCE

# Getting most variables CpGs
variance_dis <- sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])})

# Plotting variance
plot(density(variance_dis))
abline(v=0.05)

# Filtering data
dis_to_analyse <- t(promoter_betas[variance_dis > 0.05,])

#
# Running WGCNA
#

# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# 
# # Call the network topology analysis function
# sft = pickSoftThreshold(
#   dis_to_analyse,             # <= Input data
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
  
  netwk <- blockwiseModules(dis_to_analyse,               
                            corrType="bicor", # Using biweight midcorrelation 
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
                            saveTOMs = F,
                            saveTOMFileBase = "ER",
                            
                            # == Output Options
                            numericLabels = T,
                            verbose = 3)
  
  
  # Saving network
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/promoter/only_nonBasal_cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
  
}
