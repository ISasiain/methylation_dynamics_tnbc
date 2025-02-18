#! usr/bin/Rscript

library(mclust)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

#
# DEFINING FUNCTIONS
#


calculate_logFC <- function(group1, group2, base = 2, method = "median", log_transform = FALSE) {
  if (!is.numeric(group1) || !is.numeric(group2)) {
    stop("Both inputs must be numeric vectors.")
  }
  
  if (log_transform) {
    group1 <- log2(group1 + 1)
    group2 <- log2(group2 + 1)
  }
  
  # Compute central tendency based on chosen method
  if (method == "mean") {
    stat1 <- mean(group1, na.rm = TRUE)
    stat2 <- mean(group2, na.rm = TRUE)
  } else if (method == "median") {
    stat1 <- median(group1, na.rm = TRUE)
    stat2 <- median(group2, na.rm = TRUE)
  } else if (method == "trimmed_mean") {
    stat1 <- mean(group1, trim = trim, na.rm = TRUE)
    stat2 <- mean(group2, trim = trim, na.rm = TRUE)
  } else {
    stop("Invalid method. Use 'mean', 'median', or 'trimmed_mean'.")
  }
  
  # Compute log fold change with small offset
  epsilon <- 1
  logFC <- log2((stat2 + epsilon) / (stat1 + epsilon)) / log2(base)
  
  return(logFC)
}




#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Example cassette file
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")

# Summary cassette file
summary_prom15 <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_promoter_beta_15.csv")
rownames(summary_prom15) <- as.character(summary_prom15$Cassette)
summary_prom15$Cassette <- NULL

# FPKM counts
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")


#
# CLUSTERING SAMPLES BASED ON BASAL/NON-BASAL SPLIT. CASSETTE 1
#

# Perform GMM clustering
gmm_model <- Mclust(summary_prom15["1",], G = 2)  

# Get the cluster assignments
cluster_assignments <- kmeans(t(summary_prom15["1",]), centers = 2)$cluster




#
# PLOTTING 
#

# Extracting the subset data
fpkm_subset_1 <- fpkm_data[, names(cluster_assignments)[cluster_assignments == 1]]
fpkm_subset_2 <- fpkm_data[, names(cluster_assignments)[cluster_assignments == 2]]

# Identifying genes in cassette 1
cpgs_in_cassette <- names(promoter_15$colors[promoter_15$colors == 1])
genes_in_cassette <- unique(genes[cpgs_in_cassette])

# Count the number of CpGs linked to each gene
cpg_counts <- table(genes[cpgs_in_cassette])

# Initialize results dataframe
results <- data.frame(
  Gene = genes_in_cassette,
  Wilcoxon_p = NA,
  Bonferroni_p = NA,
  Log2_fold_change = NA,
  CpG_Count = NA
)

# Compute statistics
for (gene in genes_in_cassette) {
  if (gene %in% rownames(fpkm_data)) {
    expr_group1 <- as.numeric(fpkm_subset_1[gene, ])
    expr_group2 <- as.numeric(fpkm_subset_2[gene, ])
    
    # Wilcoxon test
    test_result <- wilcox.test(expr_group1, expr_group2, exact = FALSE)
    
    # Compute median difference
    log_fold_change <- calculate_logFC(expr_group1, expr_group2)
    
    # Store results
    results[results$Gene == gene, "Wilcoxon_p"] <- test_result$p.value
    results[results$Gene == gene, "Log2_fold_change"] <- log_fold_change
  }
  
  # Store CpG count for the gene
  results[results$Gene == gene, "CpG_Count"] <- cpg_counts[gene]
}

# Apply Bonferroni correction
results$Bonferroni_p <- p.adjust(results$Wilcoxon_p, method = "bonferroni")

# Sort results by Bonferroni p-value
results <- results[order(results$Bonferroni_p), ]


# VOLCANO PLOT

# Compute -log10(Bonferroni p-value)
results$neg_log10_p <- -log10(results$Bonferroni_p)

# Define color based on significance
results$color <- ifelse(results$Bonferroni_p > 0.05, "grey", "blue")

# Volcano plot
ggplot(results, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) +  # Points with transparency
  geom_text(data = subset(results, Bonferroni_p <= 0.05),  # Show labels only for significant genes
            aes(label = Gene), vjust = -1, size = 4) +
  scale_color_manual(values = c("blue", "grey")) +  # Define colors
  scale_size(range = c(3, 8)) +                     # Adjust point size range
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10(Bonferroni p-value)",
    size = "CpG Count",
    color = "Significance"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )


# HEATMAPS

current_gene_id = "ZDH1"

pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
im_annotation <- my_annotations[colnames(betaAdj), "IM"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])


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
  "Normal beta" = rowMeans(betaNorm[names(genes)[genes == current_gene_id],]),
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")))
)

# Heatmap of genes
Heatmap(
  betaAdj[names(genes)[genes == current_gene_id],],
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE, 
  column_split = factor(cluster_assignments),
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  right_annotation = left_annotation,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  name = "Tumor beta"
)

