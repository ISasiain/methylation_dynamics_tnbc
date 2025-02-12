#! usr/bin/Rscript

library(WGCNA)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)


#
# Loading the data
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Loading gene expression
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# CpG context

# Getting distal CpGs
distal_cpgs <- annoObj$illuminaID[which( ( (annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body") ) )]

distal_betas <- betaAdj[rownames(betaAdj) %in% distal_cpgs, ]


plot(density(sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])})))
abline(v=0.05)
sum(sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])}) > 0.1)



# Testing WGCNA
par(mfrow = c(1,1))

# 1. DISTAL CPGS

# Getting most variables CpGs
variance_dis <- sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])})

# Plotting variance
plot(density(variance_dis))
abline(v=0.05)

# Filtering data
dis_to_analyse <- t(distal_betas[variance_dis > 0.1,])

# Running WGCNA
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  dis_to_analyse,             # <= Input data
  powerVector = powers,
  verbose = 5
)

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
abline(h = 0.90, col = "red")
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
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
  
  
}

# Calculating CpG casettes based on WGCNA





# Heatmap of PC1


# Determining values for meta-CpGs using PCA

pca_var_distal <- list()

for (module in unique(module_df_dis$labs)) {
  
  # Subsetting
  betas <- betaAdj[module_df_dis[module_df_dis$labs== module,"gene_id"],]
  
  my_pc1 <- prcomp(t(betas))$x[,1]
  
  # Store results
  pca_var_distal[[as.character(module)]] <- list("PC1"= my_pc1,
                                                   "Number_of_CpGs"=nrow(betas)) # Proportion of variance by PC1
  
}

# Create a data frame to store PC1 values
modules <- names(pca_var_distal)
pc1_matrix <- do.call(cbind, lapply(pca_var_distal, function(x) x$PC1))

# Transpose and assign module names as row names
pc1_df <- as.data.frame(t(pc1_matrix))
rownames(pc1_df) <- modules

# View the resulting data frame
head(pc1_df)


# Create annotation object
column_annotation <- HeatmapAnnotation(df = data.frame(PAM50 = pam50_annotations,
                                                       TNBC = tnbc_annotation),
                                       HRD = HRD_annotation,
                                       epitype = epi_annotation,
                                       IM = im_annotation)

# Plot the heatmap with column annotations
Heatmap(as.matrix(pc1_df[-c(1,2),]),  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations

# Plot the heatmap with column annotations
Heatmap(betaAdj[module_df_dis[module_df_dis$labs==3 | module_df_dis$labs==4, "gene_id"],],  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations

#
# CASSETTE ANALYSIS ONLY IN BASAL SAMPLES
#

basal_to_analyse <- dis_to_analyse[pam50_annotations == "Basal",]

# Running WGCNA
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  basal_to_analyse,             # <= Input data
  powerVector = powers,
  verbose = 5
)

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
abline(h = 0.90, col = "red")
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
picked_power = 10
cor <- WGCNA::cor 

# Calculating CpG casettes based on WGCNA

netwk <- blockwiseModules(basal_to_analyse,               
                          corrType="bicor", # Using biweight midcorrelation 
                          nThreads = 10,
                          
                          # == Adjacency Function ==
                          power = picked_power,             
                          networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 8,
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

# Sumarizing results
module_df_dis <- data.frame(
  gene_id = names(netwk$colors),
  labs = netwk$colors
)

# Determining values for meta-CpGs using PCA
pca_var_distal <- list()

for (module in unique(module_df_dis$labs)) {
  
  # Subsetting
  betas <- betaAdj[module_df_dis[module_df_dis$labs== module,"gene_id"], pam50_annotations=="Basal"]
  
  my_pc1 <- prcomp(t(betas[]))$x[,1]
  
  # Store results
  pca_var_distal[[as.character(module)]] <- list("PC1"= my_pc1,
                                                 "Number_of_CpGs"=nrow(betas)) # Proportion of variance by PC1
  
}

# Create a data frame to store PC1 values
modules <- names(pca_var_distal)
pc1_matrix <- do.call(cbind, lapply(pca_var_distal, function(x) x$PC1))

# Transpose and assign module names as row names
pc1_df <- as.data.frame(t(pc1_matrix))
rownames(pc1_df) <- modules

# View the resulting data frame
head(pc1_df)


# Create annotation object
column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations[pam50_annotations=="Basal"],
                                       TNBC = tnbc_annotation[pam50_annotations=="Basal"],
                                       HRD = HRD_annotation[pam50_annotations=="Basal"],
                                       epitype = epi_annotation[pam50_annotations=="Basal"],
                                       IM = im_annotation[pam50_annotations=="Basal"])

# Plot the heatmap with column annotations
Heatmap(as.matrix(pc1_df),  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations

Heatmap(betaAdj[module_df_dis[module_df_dis$labs == "1", "gene_id"], pam50_annotations=="Basal"],
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        top_annotation = column_annotation,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2")


plot(scale(t(fpkm_data[c("IL32", "TNIP1"),])))
plot(density(as.numeric(fpkm_data["IL32",])))

