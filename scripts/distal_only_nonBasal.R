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
distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")

# Summary distal cassettes
summary_dis15 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_15.csv")
rownames(summary_dis15) <- as.character(summary_dis15$Cassette)
summary_dis15$Cassette <- NULL

# Getting distal CpGs and betas
distal_cpgs <- annoObj$illuminaID[which( ( (annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body") ) )]
distal_betas <- betaAdj[rownames(betaAdj) %in% distal_cpgs, ]

#
# Preprocessing
#

# GETTING ONLY BASAL SAMPLES

# Getting cassettes linked to basal/nonBasal split
my_cpgs_dis <-  c(
  names(distal_15$colors)[distal_15$colors == "2"],
  names(distal_15$colors)[distal_15$colors == "4"]
)

# Geenerate data frame to store groups
groupings_df <- data.frame(matrix(nrow = length(colnames(summary_dis15)), ncol = 1))
rownames(groupings_df) <- colnames(summary_dis15)
colnames(groupings_df) <- "group_dis"


# Clustering in two groups
distance_matrix <- dist(t(betaAdj[my_cpgs_dis,]))
hc <- hclust(distance_matrix)
groupings_df$group_dis <- cutree(hc, k = 2)

# Getting basal group
distal_betas <- distal_betas[,rownames(groupings_df)[groupings_df==1]]


# FILTERING BASED ON VARIANCE

# Getting most variables CpGs
variance_dis <- sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])})

# Plotting variance
plot(density(variance_dis))
abline(v=0.1)

# Filtering data
dis_to_analyse <- t(distal_betas[variance_dis > 0.1,])

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
bicor = WGCNA::cor

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
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/distal/only_nonBasal_cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
  
}

