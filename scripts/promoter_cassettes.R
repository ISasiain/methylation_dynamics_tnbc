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

# Getting CpGs belonging to each context

promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass=="promoter")]
promoter_betas <- betaAdj[rownames(betaAdj) %in% promoter_cpgs, ]


plot(density(sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])})))
abline(v=0.05)
sum(sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])}) > 0.05)


# Testing WGCNA


par(mfrow = c(1,1))
# 1. PROMOTER CPGS

# Getting most variables CpGs
variance_prom <- sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])})

# Plotting variance
plot(density(variance_prom))
abline(v=0.05)

# Filtering data
prom_to_analyse <- t(promoter_betas[variance_prom > 0.05,])

# Running WGCNA
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  prom_to_analyse,             # <= Input data
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
  my_filename <- paste0("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_", beta, ".rds" )
  saveRDS(netwk, file = my_filename)
  
}


#
# Determining values for meta-CpGs using PCA
#

pca_var_promoter <- list()

for (module in unique(module_df_prom$labs)) {
  
  # Subsetting
  betas <- betaAdj[module_df_prom[module_df_prom$labs== module,"gene_id"],]
  
  my_pca <- prcomp(t(betas))
  
  # Store results
  pca_var_promoter[[as.character(module)]] <- list("PC1"= my_pca$x[,1],
                                                   "Number_of_CpGs"=nrow(betas),
                                                   "Proportion_of_variance"=my_pca$importance[1,2]) 
  
}

# Create a data frame to store PC1 values
modules <- names(pca_var_promoter)
pc1_matrix <- do.call(cbind, lapply(pca_var_promoter, function(x) x$PC1))

# Transpose and assign module names as row names
pc1_df <- as.data.frame(t(pc1_matrix))
rownames(pc1_df) <- modules

# View the resulting data frame
head(pc1_df)


# Create annotation object
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

# Plot the heatmap with column annotations
Heatmap(as.matrix(pc1_df),  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations



# Heatmap  of only PAM50 Basal
pc1_basal <- pc1_df[pam50_annotations == "Basal"]

# Create annotation object
column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations[pam50_annotations == "Basal"],
                                       TNBC = tnbc_annotation[pam50_annotations == "Basal"],
                                       HRD = HRD_annotation[pam50_annotations == "Basal"],
                                       Epitype = epi_annotation[pam50_annotations == "Basal"],
                                       IM = im_annotation[pam50_annotations == "Basal"],
                                       col = list(
                                         "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                         "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                         "IM"=c("Negative"="grey", "Positive"="black"),
                                         "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                         "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                     "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
                                       )

Heatmap(as.matrix(pc1_basal),  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations


# Heatmap  of only PAM50 NonBasal
pc1_nonbasal <- pc1_df[pam50_annotations != "Basal"]

# Create annotation object
column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations[pam50_annotations != "Basal"],
                                       TNBC = tnbc_annotation[pam50_annotations != "Basal"],
                                       HRD = HRD_annotation[pam50_annotations != "Basal"],
                                       Epitype = epi_annotation[pam50_annotations != "Basal"],
                                       IM = im_annotation[pam50_annotations != "Basal"],
                                       col = list(
                                         "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                         "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                         "IM"=c("Negative"="grey", "Positive"="black"),
                                         "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                         "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                     "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
)

Heatmap(as.matrix(pc1_nonbasal),  # Exclude 'Number_of_CpGs' column for heatmap
        name = "PC1 of CpG casettes", 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = column_annotation)  # Add column annotations


#
# DETECTION OF CASSETTES LINKED TO DIFFERENTIAL EXPRESSION
#

set.seed(123)

# Getting cassettes
cassettes <- sort(unique(module_df_prom$labs))

#Getting genes and linked CpGs
#genes <- annoObj$nameUCSCknownGeneOverlap


genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Initialize an empty dataframe to store results
results <- data.frame(Cassette = integer(), Gene = character(), P_Value = numeric(), stringsAsFactors = FALSE)

# Iterate over cassettes starting from cassette 5
for (cassette in cassettes[-c(1)]) {
  # Extract CpGs and their linked genes for the current cassette
  cassette_genes <- unique(genes[module_df_prom$gene_id[module_df_prom$labs == cassette]])
  cassette_genes <- cassette_genes[cassette_genes != ""]
  
  # Split gene names with '-' or '_' and check individual parts
  valid_genes <- unique(unlist(lapply(cassette_genes, function(gene) {
    parts <- unlist(strsplit(gene, "[-_]"))
    parts[parts %in% rownames(fpkm_data)] # Retain only valid parts
  })))
  
  # Skip if no valid genes or parts are found
  if (length(valid_genes) == 0) next
  
  # Adjusted beta values for CpGs in the current cassette
  cassette_beta <- betaAdj[module_df_prom$gene_id[module_df_prom$labs == cassette], ]
  
  # Perform k-means clustering
  clusters <- kmeans(t(cassette_beta), centers = 2)$cluster
  
  # Calculate p-values for each valid gene
  cassette_results <- data.frame(
    Cassette = cassette,
    Gene = valid_genes,
    P_Value = sapply(valid_genes, function(gene) {
      cluster1 <- fpkm_data[gene, names(clusters)[clusters == 1]]
      cluster2 <- fpkm_data[gene, names(clusters)[clusters == 2]]
      if (length(cluster1) > 1 && length(cluster2) > 1) {
        wilcox.test(as.numeric(cluster1), as.numeric(cluster2))$p.value
      } else {
        NA
      }
    }),
    stringsAsFactors = FALSE
  )
  
  # Append results for the current cassette
  results <- rbind(results, cassette_results)
}

# Adjust p-values for multiple testing
cassette_results$Adjusted_P_Value <- p.adjust(cassette_results$P_Value, method = "bonferroni")

# Sort the results by P_Value
results <- results[order(results$P_Value), ]

# Get matches to include p < 1e-6
number_to_include <- sum(results$P_Value < 1e-5)

# Select the top ignificant matches
top_matches<- head(results, number_to_include)
n <- 0

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

pam50_annotations <- my_annotations[colnames(pc1_df), "PAM50"]
tnbc_annotation <- my_annotations[colnames(pc1_df), "TNBC"]
HRD_annotation <- my_annotations[colnames(pc1_df), "HRD"]
epi_annotation <- my_annotations[colnames(pc1_df), "NMF_atacDistal"]
im_annotation <- my_annotations[colnames(pc1_df), "IM"]
tils_annotation <- as.numeric(x[colnames(pc1_df), "TILs"])


# Loop through the top 100 patches
for (i in 1:nrow(top_matches)) {
  # Extract gene_id and cassette for the current patch
  current_gene_id <- top_matches$Gene[i]
  current_cassette <- top_matches$Cassette[i]
  current_p_value <- top_matches$P_Value[i]  # Extract p-value
  
  cpgs_in_cassette <- nrow(betaAdj[module_df_prom[module_df_prom$labs == current_cassette, "gene_id"], ])
  
  # Extract CpGs and their linked genes for the current cassette
  cassette_genes <- genes[module_df_prom$gene_id[module_df_prom$labs == current_cassette]]
  
  # Split gene names with '-' or '_' and check individual parts
  valid_genes <- unlist(lapply(cassette_genes, function(gene) {
    parts <- unique(unlist(strsplit(gene, "[-_]")))
    parts[parts %in% rownames(fpkm_data)] # Retain only valid parts
  }))
  
  # Get the number of CpGs linked to gene
  cpgs_linked_to_gene <- sum(valid_genes == current_gene_id)
  
  n <- n + 1
  
  # Create a filename for the output
  output_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/plots/promoter_cassettes2/", n, "_", current_cassette, "_", current_gene_id, ".pdf")

    # Create top anotation
  top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                      TNBC = tnbc_annotation,
                                      HRD = HRD_annotation,
                                      IM = im_annotation,
                                      Epitype = epi_annotation,
                                      TILs = anno_points(tils_annotation,
                                                         ylim=c(0,100),
                                                         size=unit(0.75, "mm"),
                                                         axis_param = list(
                                                           side="left",
                                                           at=c(0,25,50,75,100),
                                                           labels=c("0","25","50","75","100")
                                                         )),
                                      col = list(
                                        "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                        "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                        "IM"=c("Negative"="grey", "Positive"="black"),
                                        "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                        "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                    "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                        )
                                      )
  # Generate bottom annotation
  bottom_annotation <- HeatmapAnnotation(
    "FPKM" = anno_barplot(as.numeric(fpkm_data[current_gene_id, colnames(betaAdj)]))
  )
  
  # Updated left_annotation with color scale
  left_annotation <- rowAnnotation(
    "Normal beta" = rowMeans(betaNorm[module_df_prom[module_df_prom$labs == current_cassette, "gene_id"], , drop = FALSE]),
    col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")))
    )

  # Save as pdf
  pdf(output_filename, width = 6, height = 6)
  
  # Generate the heatmap
  my_plot <- Heatmap(
      betaAdj[module_df_prom[module_df_prom$labs == current_cassette, "gene_id"], ],
      column_km = 2,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      show_row_dend = FALSE, 
      top_annotation = top_annotation,
      bottom_annotation = bottom_annotation,
      right_annotation = left_annotation,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = "ward.D2",
      column_title = paste0(current_gene_id, ", N_of_CpGs = " , cpgs_linked_to_gene,"/", cpgs_in_cassette, ", p = ", sprintf("%.3e", current_p_value)),
      name = "Tumor beta"
    )
  
  
  # Draw and save the plot
  draw(my_plot)
  
  dev.off()
  
  # Print a message to indicate progress
  cat("Saved heatmap for Gene:", current_gene_id, "Cassette:", current_cassette, "to", output_filename, "\n")
}


# Analysinsg GBP4 cassette
current_gene_id <- "GBP4"
current_cassette <- "150"

groups_methylation <- kmeans(t(betaAdj[module_df_prom[module_df_prom$labs == current_cassette, "gene_id"], ]), centers = 2)$cluster
names(groups_methylation) <- colnames(betaAdj)

boxplot(tils_annotation ~ factor(groups_methylation, 
                                 levels = c("1", "2"), 
                                 labels = c("Hypermethylated", "Hypomethylated")), 
        ylab = "TILs (%)",
        xlab=NULL)

boxplot(as.numeric(fpkm_data["GBP4", names(groups_methylation)]) ~ factor(groups_methylation, 
                                 levels = c("1", "2"), 
                                 labels = c("Hypermethylated", "Hypomethylated")), 
        ylab = "GBP4 FPKM",
        xlab=NULL)



wilcox.test(tils_annotation ~ factor(groups_methylation, 
                                     levels = c("1", "2"), 
                                     labels = c("Hypermethylated", "Hypomethylated")))

plot(log(as.numeric(fpkm_data["ZBP1",x$PD_ID])), as.numeric(x$ASCAT_TUM_FRAC), xlab="log(GBP4 FPKM)", ylab="ASCAT Purity", pch=16, cex=0.8)
plot(log(as.numeric(fpkm_data["GBP4",x$PD_ID])), as.numeric(x$TILs), xlab="log(GBP4 FPKM)", ylab="TILs (%)", pch=16, cex=0.8, col=groups_methylation[x$PD_ID])


# Determining significant cassettes of Basal and non-Basal samples

# Divide samples into basal and non-basal groups
basal_samples <- colnames(pc1_df)[pam50_annotations == "Basal"]
non_basal_samples <- colnames(pc1_df)[pam50_annotations != "Basal"]

im_pos_samples <- na.omit(colnames(pc1_df)[im_annotation == "Positive"])
im_neg_samples <- na.omit(colnames(pc1_df)[im_annotation == "Negative"])

# Subset PC1 values for basal and non-basal groups
basal_pc1 <- pc1_df[, im_pos_samples]
non_basal_pc1 <- pc1_df[, im_neg_basal_samples]

# Initialize a data frame to store results
results <- data.frame(CpG_Cassette = rownames(pc1_df),
                      p_value = NA)

# Perform statistical tests for each cassette
for (i in 1:nrow(pc1_df)) {
  results$p_value[i] <- wilcox.test(as.numeric(basal_pc1[i, ])^(1/3), as.numeric(non_basal_pc1[i, ])^(1/3))$p.value
}

# Correct p-values using boferroni
results$adjusted_p_value <- p.adjust(results$p_value, method = "bonferroni")

# Add significance labels
results$significant <- results$adjusted_p_value < 0.01


# Subset the significant cassettes
top_cassettes <- results$CpG_Cassette[results$significant]

# Assuming the dataframe is named 'result'
top_10_cassettes <- results[order(results$adjusted_p_value), ][1:10, ]
print(top_10_cassettes)

# Extract PC1 values for the top cassettes
heatmap_data <- as.matrix(pc1_df[top_cassettes, ])

# Add annotations for PAM50 groups
pam50_colors <- c("Basal" = "red", "Non-Basal" = "blue")
im_colors <- c("Positive" = "green", "Negative"="grey")
pam50_annotations <- sapply(pam50_annotations, function(x) {if (x == "Basal") {"Basal"} else {"Non-Basal"}})


column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                       IM = im_annotation,
                                       col = list(PAM50 = pam50_colors,
                                                  IM=im_colors))

# Plot heatmap
Heatmap(heatmap_data,
        name = "PC1",
        top_annotation = column_annotation,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE)


results[results$Cassette == 1 & results$Cassette > 1e-2,2]


# Heatmap of most significant cassette

# Generate the bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "PC1" = as.numeric(pc1_df["1", ]),
  col = list(PC1 = colorRamp2(c(-30, 0, 10), c("yellow", "white", "purple")))
  )

Heatmap(betaAdj[module_df_prom[module_df_prom$labs == "24", "gene_id"], ],
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        top_annotation = column_annotation,
        bottom_annotation = bottom_annotation,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2"
)  


# Genes associated to cassette 1

genes <- annoObj$nameUCSCknownGeneOverlap
names(genes) <- annoObj$illuminaID

genes[module_df_prom$"gene_id"]


sort(table(genes[module_df_prom[module_df_prom$labs == "14", "gene_id"]]))


# Generate the bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM_exp" = anno_barplot(as.numeric(fpkm_data["IRX4",colnames(betaAdj)])),
  "Cluster" = kmeans(t(betaAdj[module_df_prom[module_df_prom$labs == "14", "gene_id"],]), centers = 2)$cluster
)
Heatmap(betaAdj[module_df_prom[module_df_prom$labs == "14", "gene_id"], ],
        cluster_rows = TRUE, 
        cluster_columns = TRUE, 
        show_row_names = FALSE, 
        show_column_names = FALSE, 
        top_annotation = column_annotation,
        bottom_annotation = bottom_annotation,
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2"
)  




par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(t(as.matrix(pc1_df[c("6","5"),])), xlab="Casette 6 (EYA4)", ylab="Casette 5 (BRCA1)")
plot(t(as.matrix(pc1_df[c("7","5"),])), xlab="Casette 7 (GSX1)", ylab="Casette 5 (BRCA1)")


heatmap(as.matrix(pc1_df[c("5","6", "7", "192"),]))


# Plotting casette 10 vs epitypes
boxplot(as.numeric(pc1_df["128", ]) ~ pam50_annotations,
        main = "Casette 1 (Basal/Non-basal)",
        xlab = "Epi Annotation",
        ylab = "PC1 Values",
        col = "lightblue")

# Boxplot with jittered scatter plot overlay
boxplot(as.numeric(pc1_df["100", ]) ~ epi_annotation,
        main = "Casette 2 (Immune infiltration?)",
        xlab = "Epi Annotation",
        ylab = "PC1 Values",
        col = "lightblue",
        border = "darkblue",
        outline = FALSE)  # Suppress outlier points to avoid duplication

# Add jittered points
points(jitter(as.numeric(as.factor(epi_annotation))), 
       as.numeric(pc1_df["100", ]),
       pch = 16,  # Solid points
       col = "red",
       cex = 0.7)  # Adjust point size

boxplot(as.numeric(pc1_df["2", ]) ~ epi_annotation,
        main = "Casette 10 (PPT2)",
        xlab = "Epi Annotation",
        ylab = "PC1 Values",
        col = "lightblue")

boxplot(as.numeric(pc1_df["5", ]) ~ epi_annotation,
        main = "Casette 11 (CYBA)",
        xlab = "Epi Annotation",
        ylab = "PC1 Values",
        col = "lightblue")


# Define the HRD and IM status annotations
HRD_annotation <- epi_annotations[colnames(pc1_df), "HRD"]
im_annotation <- epi_annotations[colnames(pc1_df), "IM"]

cassettes <- c("2")

for (casette in  cassettes) {
  
  print(casette)
  # Extract PC1 values for the specified cassette (row "103")
  pc1_values <- as.numeric(pc1_df[casette, ])
  
  # Create a data frame with the relevant information
  plot_data <- data.frame(
    PC1 = pc1_values,
    HRD = HRD_annotation,
    IM = im_annotation
  )
  
  # Subset data for HRD+ and HRD- samples
  hrd_positive_data <- plot_data[plot_data$HRD == "High", ]
  hrd_negative_data <- plot_data[plot_data$HRD == "Low/Inter", ]
  
  # Plot 1: HRD+ samples split by IM status
  hrd_positive_plot <- ggplot(hrd_positive_data, aes(x = IM, y = PC1, fill = IM)) +
    geom_boxplot(color = "black") +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, aes(color = IM)) +
    labs(
      title = "PC1 Values in HRD+ Samples, Split by IM Status",
      x = "Immune Status (IM)",
      y = "PC1 Values"
    ) +
    scale_fill_manual(values = c("lightblue", "lightcoral")) +  # Customize colors for IM status
    scale_color_manual(values = c("blue", "red")) +  # Customize colors for IM status
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Plot 2: HRD- samples split by IM status
  hrd_negative_plot <- ggplot(hrd_negative_data, aes(x = IM, y = PC1, fill = IM)) +
    geom_boxplot(color = "black") +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6, aes(color = IM)) +
    labs(
      title = "PC1 Values in HRD- Samples, Split by IM Status",
      x = "Immune Status (IM)",
      y = "PC1 Values"
    ) +
    scale_fill_manual(values = c("lightblue", "lightgreen")) +  # Customize colors for IM status
    scale_color_manual(values = c("blue", "green")) +  # Customize colors for IM status
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Display the plots side by side
  library(gridExtra)
  grid.arrange(hrd_positive_plot, hrd_negative_plot, ncol = 2)
  
  readline()
  
}

cassettes <- c("2", "65", "66", "4", "188", "243", "225", "98", "246", "215", "237", "103", "276", "11", "183", "260", "58")


# Cluster samples based on casette 2
# Perform k-means clustering for cassette "2" with 2 clusters
c2_clusters <- kmeans(t(as.matrix(pc1_df["2", ])), centers = 2)

# Create a data frame for plotting
plot_data <- data.frame(
  TILs = as.numeric(x[colnames(pc1_df), "TILs"]),
  PC1 = as.numeric(pc1_df["2", ]),
  Cluster = as.factor(c2_clusters$cluster)
)


ggplot(plot_data, aes(x = Cluster, y = TILs, fill = Cluster)) +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6, aes(color = Cluster)) +
  labs(
    title = "TILs by K-means Cluster for Cassette 2",
    x = "Cluster",
    y = "TILs Values"
  ) +
  scale_fill_manual(values = c("lightblue", "lightcoral")) +  # Customize colors for clusters
  scale_color_manual(values = c("blue", "red")) +  # Customize jitter point colors
  theme_minimal() +
  theme(legend.title = element_blank())


# Create the ggplot
ggplot(plot_data, aes(x = TILs, y = PC1, color = Cluster)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = "TILs vs PC1",
    x = "TILs",
    y = "PC1",
    color = "Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )






