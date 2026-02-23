#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(energy)
library(mclust)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(missMethyl)


#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# PROMOTER
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

summary_prom10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(summary_prom10) <- as.character(summary_prom10$Cassette)
summary_prom10$Cassette <- NULL

# DISTAL
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")

summary_dis10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(summary_dis10) <- as.character(summary_dis10$Cassette)
summary_dis10$Cassette <- NULL

# PROXIMAL
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")

summary_prox10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_10.csv")
rownames(summary_prox10) <- as.character(summary_prox10$Cassette)
summary_prox10$Cassette <- NULL

# FPKM counts
fpkm_data <- read.table("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Epitype annotations
my_annotations <- read.table("/Users/in2245sa//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")


#
# PLOTTING DIFFERENTIALLY EXPRESSED GENES
#

# PROMOTER

# Clustering

# Getting data of particular cassettes
my_cpgs_prom <-  c(
  names(promoter_10$colors)[promoter_10$colors == "1"]
)

# Clustering
distance_matrix <- dist(t(betaAdj[my_cpgs_prom,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
clusters_to_analyse <- cutree(hc, k = 2)


# Extracting the subset data
fpkm_subset_1 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse== 1]]
fpkm_subset_2 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse == 2]]

# Identifying genes in cassette 1
cpgs_in_cassette <- c(names(promoter_10$colors[promoter_10$colors == 1]))
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
    test_result <- wilcox.test(expr_group2, expr_group1, exact = FALSE)
    
    # Compute log FC
    log_fold_change <- mean(log2(expr_group2 + 1)) -
      mean(log2(expr_group1 + 1))
    
    # Store results
    results[results$Gene == gene, "Wilcoxon_p"] <- test_result$p.value
    results[results$Gene == gene, "Log2_fold_change"] <- log_fold_change
  }
  
  # Store CpG count for the gene
  results[results$Gene == gene, "CpG_Count"] <- cpg_counts[gene]
}

# Apply Bonferroni correction
results$Bonferroni_p <- p.adjust(results$Wilcoxon_p, method = "BH")

# Sort results by Bonferroni p-value
results <- results[order(results$Bonferroni_p), ]


# VOLCANO PLOT

# Compute -log10(Bonferroni p-value)
results$neg_log10_p <- -log10(results$Bonferroni_p)

# Define color based on significance
results$color <- ifelse(results$Bonferroni_p > 0.05, "grey", "blue")

# Volcano plot
ggplot(results, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) + 
  geom_text_repel(
    data = subset(results, Bonferroni_p <= 0.00000000001),
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("#1f78b4", "grey")) +  
  scale_size(range = c(1, 6)) +                    
  theme_classic() +
  labs(
    x = "Log2 Fold Change in Expression",
    y = "-log10(Bonferroni p-value)",
    size = "CpG Count"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme_bw()

# Perform KEGG pathway enrichment on only differentially expressed genes linked to cassette
my_genes_for_enrichment <- na.omit(results[results$Bonferroni_p < 0.05, ])$"Gene"

# Getting CpGs linked to differentially expressed genes in cassette 1
cpgs_selected <- names(genes[genes %in% my_genes_for_enrichment])
cpgs_to_test <- my_cpgs_prom[my_cpgs_prom %in% cpgs_selected]

# Performing enrichent analysis
output <- gometh(
  cpgs_to_test,
  all.cpg = annoObj[annoObj$CpG_context == "Promoter", "illuminaID"],
  collection = "kegg",
  array.type = "EPIC",
  plot.bias = FALSE,
  prior.prob = TRUE,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = "ALL",
  sig.genes = FALSE
)

kegg_df <- output[output$FDR <= 0.05,]

kegg_df <- kegg_df[order(kegg_df$FDR, decreasing = TRUE),]
kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)


kegg_cas_1 <- ggplot(kegg_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 2.5) +
  scale_size_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

kegg_cas_1

# PROXIMAL

# Getting data
my_cpgs_prox <-  c(
  names(proximal_10$colors)[proximal_10$colors == "1"]
)

# Clustering
distance_matrix <- dist(t(betaAdj[my_cpgs_prox,]), method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
clusters_to_analyse <- cutree(hc, k = 2)

# Extracting the subset data
fpkm_subset_1 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse== 1]]
fpkm_subset_2 <- fpkm_data[, names(clusters_to_analyse)[clusters_to_analyse == 2]]

# Identifying genes in cassette 1
cpgs_in_cassette <- names(proximal_10$colors[proximal_10$colors %in% c(1)])
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
    
    # Compute log FC
    log_fold_change <- mean(log2(expr_group2 + 1)) -
      mean(log2(expr_group1 + 1))
    
    # Store results
    results[results$Gene == gene, "Wilcoxon_p"] <- test_result$p.value
    results[results$Gene == gene, "Log2_fold_change"] <- log_fold_change
  }
  
  # Store CpG count for the gene
  results[results$Gene == gene, "CpG_Count"] <- cpg_counts[gene]
}

# Apply Bonferroni correction
results$Bonferroni_p <- p.adjust(results$Wilcoxon_p, method = "BH")

# Sort results by Bonferroni p-value
results <- results[order(results$Bonferroni_p), ]


# VOLCANO PLOT

# Compute -log10(Bonferroni p-value)
results$neg_log10_p <- -log10(results$Bonferroni_p)

# Define color based on significance
results$color <- ifelse(results$Bonferroni_p > 0.05, "grey", "blue")

# Volcano plot
ggplot(results, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) + 
  geom_text_repel(
    data = subset(results, Bonferroni_p <= 0.0000000000008),
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("#1f78b4", "grey")) +  
  scale_size(range = c(1, 6)) +                    
  theme_classic() +
  labs(
    x = "Log2 Fold Change in Expression",
    y = "-log10(Bonferroni p-value)",
    size = "CpG Count"
  ) +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme_bw()


# Perform KEGG pathway enrichment on only differentially expressed genes linked to cassette
my_genes_for_enrichment <- na.omit(results[results$Bonferroni_p < 0.05, ])$"Gene"

# Getting CpGs linked to differentially expressed genes in cassette 1
cpgs_selected <- names(genes[genes %in% my_genes_for_enrichment])
cpgs_to_test <- my_cpgs_prox[my_cpgs_prox %in% cpgs_selected]

# Performing enrichent analysis
output <- gometh(
  cpgs_to_test,
  all.cpg = annoObj[annoObj$CpG_context == "Proximal", "illuminaID"],
  collection = "kegg",
  array.type = "EPIC",
  plot.bias = FALSE,
  prior.prob = TRUE,
  equiv.cpg = TRUE,
  fract.counts = TRUE,
  genomic.features = "ALL",
  sig.genes = FALSE
)

kegg_df <- output[output$FDR <= 0.05,]

kegg_df <- kegg_df[order(kegg_df$FDR, decreasing = TRUE),]
kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)


kegg_cas_1 <- ggplot(kegg_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 2.5) +
  scale_size_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

kegg_cas_1

