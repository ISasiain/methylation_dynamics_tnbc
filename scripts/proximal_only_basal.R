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

# proximal cassettes
proximal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_15.rds")

# Summary proximal cassettes
summary_prox15 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_15.csv")
rownames(summary_prox15) <- as.character(summary_prox15$Cassette)
summary_prox15$Cassette <- NULL

# Getting proximal CpGs and betas
proximal_cpgs <- annoObj$illuminaID[which( ( (annoObj$featureClass=="proximal up") | (annoObj$featureClass=="proximal dn") ) )]
proximal_betas <- betaAdj[rownames(betaAdj) %in% proximal_cpgs, ]

#
# Preprocessing
#

# GETTING ONLY BASAL SAMPLES

# Getting cassettes linked to basal/nonBasal split
my_cpgs_prox <-  c(
  names(proximal_15$colors)[proximal_15$colors == "1"],
  names(proximal_15$colors)[proximal_15$colors == "2"]
)

# Geenerate data frame to store groups
groupings_df <- data.frame(matrix(nrow = length(colnames(summary_prox15)), ncol = 1))
rownames(groupings_df) <- colnames(summary_prox15)
colnames(groupings_df) <- "group_prox"


# Clustering in two groups
distance_matrix <- dist(t(betaAdj[my_cpgs_prox,]))
hc <- hclust(distance_matrix)
groupings_df$group_prox <- cutree(hc, k = 2)

# Getting basal group
proximal_betas <- proximal_betas[,rownames(groupings_df)[groupings_df==1]]


# FILTERING BASED ON VARIANCE

# Getting most variables CpGs
variance_prox <- sapply(1:nrow(proximal_betas), FUN = function(row) {var(proximal_betas[row,])})

# Plotting variance
plot(density(variance_prox))
abline(v=0.05)

# Filtering data
prox_to_analyse <- t(proximal_betas[variance_prox > 0.05,])

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
  
  netwk <- blockwiseModules(prox_to_analyse,               
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
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/proximal/only_basal_cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
  
}


#
# PLOTTING HEATMAPS
#

# List all files
proximal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/proximal/", full.names = TRUE, pattern = "*only_basal*")

# Initialize an empty data frame
summary_df <- data.frame(beta = numeric(), num_cassettes = numeric(), mean_cassette_length = numeric())

for (file in proximal_files) {
  
  # Getting beta
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+)\\.rds", "\\1", file))
  
  # Analysing cassettes
  my_data <- readRDS(file)$colors
  
  # Exclude cassette 0
  my_data <- my_data[my_data != 0]
  
  num_cassettes <- length(unique(unname(my_data)))
  mean_cassette_length <- mean(table(my_data))
  
  # Appending to df
  summary_df <- rbind(summary_df, data.frame(beta = beta, num_cassettes = num_cassettes, mean_cassette_length = mean_cassette_length))
}

# Convert beta to a factor with levels in the desired order
summary_df$beta <- factor(summary_df$beta, levels = c(5, 8, 10, 15, 20, 25))

# Create the plot
ggplot(summary_df, aes(x = beta)) +
  geom_bar(aes(y = num_cassettes), stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(aes(y = mean_cassette_length * 25), color = "red", size = 3) +  # Adjust scaling factor as needed
  geom_line(aes(y = mean_cassette_length * 25, group = 1), color = "red", size = 1) +  # Adjust scaling factor as needed
  labs(title = "proximal cassettes (Var > 0.1)",
       x = "Beta",
       y = "Number of Cassettes") +
  scale_y_continuous(limits = c(0, 1050), sec.axis = sec_axis(~ . / 25, name = "Mean Cassette Length")) +  # Adjust scaling factor as needed
  theme_classic()


# Plotting first cassettes with annotatios. 

# Loop through each file
for (file in proximal_files) {
  # Extract the beta value from the filename
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+)\\.rds", "\\1", file))
  
  # Load the corresponding data
  proximal <- readRDS(file)
  
  # Define the cassettes to include
  selected_cassettes <- 1:13
  
  # Extract CpGs belonging to each cassette
  cpg_list <- lapply(selected_cassettes, function(cassette) {
    names(proximal$colors[proximal$colors == cassette])
  })
  names(cpg_list) <- selected_cassettes  # Assign cassette names
  
  # Flatten the list to get selected CpGs
  selected_CpGs <- unlist(cpg_list)
  
  # Subset the adjusted beta values
  beta_subset <- proximal_betas[selected_CpGs, ]
  
  # Create a cassette grouping factor
  cassette_factor <- factor(proximal$colors[selected_CpGs], levels = selected_cassettes)
  
  # Compute CpG counts per cassette
  cassette_counts <- sapply(cpg_list, length)
  row_labels <- paste0(selected_cassettes)
  
  
  #Define ATAC annotation
  atac_annotation <- rowAnnotation(df=data.frame("ATAC"=as.factor(annoObj$hasAtacOverlap[annoObj$illuminaID %in% unname(selected_CpGs)])),
                                   col = list(ATAC = c("0" = "white", "1" = "black")))
  
  
  # Epitype annotations
  pam50_annotations <- my_annotations[colnames(beta_subset), "PAM50"]
  tnbc_annotation <- my_annotations[colnames(beta_subset), "TNBC"]
  HRD_annotation <- my_annotations[colnames(beta_subset), "HRD"]
  epi_annotation <- my_annotations[colnames(beta_subset), "NMF_atacDistal"]
  im_annotation <- my_annotations[colnames(beta_subset), "IM"]
  tils_annotation <- as.numeric(x[colnames(beta_subset), "TILs"])
  
  # Create column annotation object
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                         TNBC = tnbc_annotation,
                                         HRD = HRD_annotation,
                                         Epitype = epi_annotation,
                                         IM = im_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                           "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                           "IM"=c("Negative"="grey", "Positive"="black"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate heatmap with row annotation
  heatmap <- Heatmap(beta_subset, 
                     cluster_rows = FALSE, 
                     cluster_columns = TRUE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     row_split = cassette_factor, 
                     row_title = row_labels,  
                     column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
                     top_annotation = column_annotation,
                     right_annotation = atac_annotation,
                     use_raster = FALSE)
  
  # Save the heatmap to a file with double size
  pdf(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/diff_betas/only_basal_heatmap_beta_", beta, ".pdf"), width = 14, height = 10)  # Adjust width and height as needed
  draw(heatmap)
  dev.off()
}

