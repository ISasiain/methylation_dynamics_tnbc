#! usr/bin/Rscript


library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(survival)
library(survminer)
library(tidyr)
library(dplyr)
library(patchwork)



#
# LOAD DATA
#

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Loading gene expression
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Reading PDL1 TLS annotation. Matrix is called u.frame
load("PhD/Projects/project_3/summarized_TPS_data_sampleLevel.RData")
rownames(u.frame) <- u.frame$PD_ID

#
# PLOTTING CPGS AFFECTING 
#

current_gene_id = "KIT"

  
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
    "FPKM" = anno_barplot(log(as.numeric(fpkm_data[current_gene_id, colnames(betaAdj)])+1))
  )
  
  # Updated left_annotation with color scale
  right_annotation <- rowAnnotation(
    "Normal beta" = rowMeans(betaNorm[names(genes)[genes == current_gene_id],]),
    "ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
    col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
               "ATAC" = c("0" = "white", "1"= "black"))
  )
  
  # CpG context annotation
  left_annotation <- rowAnnotation("Context"= annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]]
  )
  
  # Atac annotation
  
  # Cluster based on methylation
  cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                   annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])
  
  cluster_gbp4_promoter <- kmeans(t(betaAdj[names(cpgs)[cpgs=="promoter"],]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  gbp4_promoter_state <- if (mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_gbp4_promoter$cluster==1]) >
                             mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_gbp4_promoter$cluster==2])) {
    
    as.factor(ifelse(cluster_gbp4_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
    
  } else {
    
    as.factor(ifelse(cluster_gbp4_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
    
  }
  
  # Heatmap of genes
  Heatmap(
    betaAdj[names(genes)[genes == current_gene_id],],
    cluster_rows = FALSE,
    row_order = order(annoObj$start[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]]),
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    column_split = gbp4_promoter_state,
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    right_annotation = right_annotation,
    left_annotation = left_annotation,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    name = "Tumor beta"
  )
  
  
  
  #
  # PLOTTING BOXPLOTS. TILs AND GENE EXPRESSION
  #
  
  # Convert data to a dataframe for ggplot2
  plot_data <- data.frame(
    FPKM = as.numeric(fpkm_data[current_gene_id, colnames(betaAdj),]),
    TILs = as.numeric(x[colnames(betaAdj), "TILs"]),
    PDL1_CPS = as.numeric(x[colnames(betaAdj), "PDL1_CPS"]),
    PDL1_TPS = sapply(u.frame[colnames(betaAdj), "PDL1_TPS"], function(x) {
      if (is.na(x)) {NA} 
      else if (x >= 10) {2} 
      else if (x >= 1) {1} 
      else {0}
    }
    ), 
    Methylation_State = gbp4_promoter_state,
    HRD = as.factor(x[colnames(betaAdj), "HRD.2.status"])
  )
  
  # EXPRESSION BOXPLOT
  
  # Perform Wilcoxon or Kruskal-Wallis test
  stat_test <- compare_means(FPKM ~ Methylation_State, data = plot_data, method = "wilcox.test")
  
  # Create the boxplot with p-value annotation
  ggplot(plot_data, aes(x = Methylation_State, y = FPKM, fill = Methylation_State)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Removes outliers
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Jittered points for visibility
    scale_fill_manual(values = c("indianred1", "cadetblue1")) +  # Custom fill colors
    theme_classic(base_size = 14) +  # Classic theme
    labs(x = "Promoter Methylation State", y = "GBP4 FPKM") +
    theme(legend.position = "none") +  # Hide legend
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       comparisons = list(c("Hypomethylated", "Hypermethylated")), 
                       label.x = c("Hypomethylated", "Hypermethylated"))  # Add p-value label and specify the comparisons
  
  # TILs BOXPLOT
  
  ggplot(plot_data, aes(x = Methylation_State, y = TILs, fill = Methylation_State)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Removes outliers for better jitter visibility
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Adds jittered points
    scale_fill_manual(values = c("indianred1", "cadetblue1")) +  # Custom colors
    theme_classic(base_size = 14) +  # Classic theme
    labs(x = "Promoter Methylation State", y = "TILs (%)") +
    theme(legend.position = "none") +
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       comparisons = list(c("Hypomethylated", "Hypermethylated")), 
                       label.x = c("Hypomethylated", "Hypermethylated"))  # Add p-value label and specify the comparisons
  
  
  # PDL1 CPS BARPLOT
  
  ggplot(subset(plot_data, !is.na(PDL1_CPS)), aes(x = as.factor(Methylation_State), fill = as.factor(PDL1_CPS))) +
    geom_bar(position = "fill", alpha = 0.8) +  # Stacked bar normalized to proportions
    theme_classic(base_size = 14) +  # Clean theme with larger text
    labs(x = "Promoter Methylation State", 
         y = "Proportion of Samples", 
         fill = "PD-L1 CPS") +  # Clearer legend title
    scale_fill_manual(values = c("0" = "#619CFF", 
                                 "1" = "#9B59B6", 
                                 "2" = "#F8766D"), 
                      labels = c("0" = "< 1%", 
                                 "1" = "1% & < 10%", 
                                 "2" = "≥ 10%")) +  # More detailed labels
    theme(legend.position = "right", 
          legend.text = element_text(size = 12))  # Improve readability of legend labels
  
  
  # PDL1 TPS BARPLOT
  ggplot(subset(plot_data, !is.na(PDL1_TPS)), aes(x = as.factor(Methylation_State), fill = as.factor(PDL1_TPS))) +
    geom_bar(position = "fill", alpha = 0.8) +  # Stacked bar normalized to proportions
    theme_classic(base_size = 14) +  # Clean theme with larger text
    labs(x = "Promoter Methylation State", 
         y = "Proportion of Samples", 
         fill = "PD-L1 TPS") +  # Clearer legend title
    scale_fill_manual(values = c("0" = "#619CFF", 
                                 "1" = "#9B59B6", 
                                 "2" = "#F8766D"), 
                      labels = c("0" = "< 1%", 
                                 "1" = "≥ 1% & < 10%", 
                                 "2" = "≥ 10%")) +  # More detailed labels
    theme(legend.position = "right", 
          legend.text = element_text(size = 12))  # Improve readability of legend labels
  
  # HRD VS METHYLATION
  
  ggplot(subset(plot_data, !is.na(HRD)), aes(x = as.factor(Methylation_State), fill = as.factor(HRD))) +
    geom_bar(position = "fill", alpha = 0.8) +  # Stacked bar normalized to proportions
    theme_classic(base_size = 14) +  # Clean theme with larger text
    labs(x = "Promoter Methylation State", 
         y = "Proportion of Samples", 
         fill = "PD-L1 TPS") +  # Clearer legend title
    scale_fill_manual(values = c("low/inter" = "#619CFF", 
                                 "high" = "#F8766D")) + 
    theme(legend.position = "right", 
          legend.text = element_text(size = 12))  # Improve readability of legend labels
  
  
  
  # EXPRESSION VS TILS
  
  # Calculate Spearman correlation
  cor_val <- cor(plot_data$FPKM, plot_data$TILs, method = "spearman", use="pairwise.complete.obs")
  
  # Adjust ggplot code
  ggplot(plot_data, aes(y = TILs, x = log(FPKM + 1))) +
    geom_point(alpha = 1, cex = 0.8) + 
    labs(x = "log GBP4 FPKM + 1", y = "TILs (%)",) +
    annotate("text", x = 2.4, y = 90,  
             label = paste("Sp. Cor. = ", round(cor_val, 2)), 
             hjust = 1, vjust = 1, size = 5, color = "blue") +
    theme_classic(base_size = 14) 
  
  # EXPRESSION VS PDL1 TPS AND CPS
  ggplot(subset(plot_data, !is.na(PDL1_CPS)), aes(x = factor(PDL1_CPS, levels = c(0, 1, 2), 
                                                             labels = c("[0,1)", "[1,10)", "[10,100]")), 
                                                  y = FPKM)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, fill="#F8766D") +  # Boxplot without outliers
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Add jitter for points
    theme_classic(base_size = 14) +  # Classic theme with larger text
    labs(x = "PDL1 CPS (%)", y = "GBP4 expression (FPKM)") +  # Labels for axes
    theme(legend.position = "none") +  # Remove legend
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       comparisons = list(c("[0,1)", "[1,10)"), 
                                          c("[0,1)", "[10,100]"),
                                          c("[1,10)", "[10,100]")), 
                       label.x = c(1, 2, 3))  # Specify comparisons for Wilcoxon test
  
  
  
  
  ggplot(subset(plot_data, !is.na(PDL1_TPS)), aes(x = factor(PDL1_TPS, levels = c(0, 1, 2), 
                                                             labels = c("[0,1)", "[1,10)", "[10,100]")), 
                                                  y = FPKM)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, fill="#F8766D") +  # Boxplot without outliers
    geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Add jitter for points
    theme_classic(base_size = 14) +  # Classic theme with larger text
    labs(x = "PDL1 TPS (%)", y = "GBP4 expression (FPKM)") +  # Labels for axes
    theme(legend.position = "none") +  # Remove legend
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       comparisons = list(c("[0,1)", "[1,10)"), 
                                          c("[0,1)", "[10,100]"),
                                          c("[1,10)", "[10,100]")), 
                       label.x = c(1, 2, 3))  # Specify comparisons for Wilcoxon test
  


#
# PLOTTING TILEPLOT OF PROMOTER HYPERMETHYLATION OF SELECTED GENES
#

gene_ids = c("GBP4", "ZBP1", "OAS2", "CARD16", "SAMD9L")

clusters_methylation <- data.frame(matrix(ncol = ncol(betaAdj), nrow = length(gene_ids)))
colnames(clusters_methylation) <- colnames(betaAdj)
rownames(clusters_methylation) <- gene_ids

for(current_gene_id in gene_ids) {
  
  # Cluster based on methylation
  cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                   annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])
  
  cluster_promoter <- kmeans(t(betaAdj[names(cpgs)[cpgs=="promoter"],]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  promoter_state <- if (mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==1]) >
                             mean(betaAdj[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==2])) {
    
    as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
    
  } else {
    
    as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
    
  }
  
  # Adding data to dataframe
  clusters_methylation[current_gene_id,] <- promoter_state
}

# Convert the character values to numeric for the entire matrix
# "hypermethylated" = 1, "hypomethylated" = 0

numeric_matrix <- apply(clusters_methylation, c(1, 2), function(x) {
  if (x == "Hypermethylated") {
    return(1)
  } else if (x == "Hypomethylated") {
    return(0)
  } else {
    return(NA)  # Handle any other values (if any)
  }
})

agreement <- colSums(numeric_matrix)

agreement_categories <- factor(agreement, levels = c(0, 1, 2, 3, 4, 5), 
                               labels = c("5/5", "4/5", "3/5", "3/5", "4/5", "5/5"))

# Create the heatmap 

top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    IM = im_annotation,
                                    Epitype = epi_annotation,
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)

Heatmap(clusters_methylation,
        col = c("red", "blue"),
        top_annotation = top_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_split = as.factor(ifelse(pam50_annotations == "Basal", "Basal", "NonBasal")))


agreement_data <- data.frame(
  Agreement = agreement_categories
)

# Plot proportions
ggplot(agreement_data, aes(x = Agreement, y = ..prop.., group = 1)) +
  geom_bar(stat = "count",fill = "darkgrey") +  # Counts are converted to proportions automatically
  scale_y_continuous(labels = scales::percent) +  # Format y-axis as percentages
  labs(x = "Agreement in methylation state",
       y = "Proportion") +
  theme_bw()

# Simulating proportions if the methylation states were randomly generated

set.seed(123)

# Parameters
n_reps <- 100000  # Number of simulations (each sim is 235 samples)
n_per_rep <- 235
probs <- c(0.7, 0.3)
elements <- c("A", "B")

# Pre-allocate storage for results
results <- data.frame(rep = 1:n_reps, p3 = NA, p4 = NA, p5 = NA)

for (i in 1:n_reps) {
  # Simulate one batch of 235 samples
  samples <- replicate(n_per_rep, sample(elements, size = 5, replace = TRUE, prob = probs))
  max_counts <- apply(samples, 2, function(x) max(table(x)))
  
  # Store proportions
  results$p3[i] <- mean(max_counts == 3)
  results$p4[i] <- mean(max_counts == 4)
  results$p5[i] <- mean(max_counts == 5)
}

# Reshape to long format for ggplot
long_results <- results %>%
  pivot_longer(cols = c("p3", "p4", "p5"), names_to = "EqualCount", values_to = "Proportion")

# Map better labels
long_results$EqualCount <- recode(long_results$EqualCount,
                                  p3 = "3 equal", p4 = "4 equal", p5 = "5 equal")

# Plot density
ggplot(long_results, aes(x = Proportion, fill = EqualCount, color = EqualCount)) +
  geom_histogram(aes(y = ..count../sum(..count..)), bins = 60, alpha = 0.4, position = "identity") +
  labs(title = "Proportions of Equal Elements (100,000 simulations)",
       x = "Proportion of samples (N=235)",
       y = "Probability",
       fill = "Equal elements",
       color = "Equal elements") +
  theme_bw(base_size = 14)



#
# SURVIVAL 
#

library(survival)

# Extract relevant data
cox_data <- data.frame(
  time = x[colnames(clusters_methylation), "OS"],
  status = x[colnames(clusters_methylation), "OSbin"],
  age = x[colnames(clusters_methylation), "Age"],
  OAS2 = t(clusters_methylation["OAS2",]),
  GBP4 = t(clusters_methylation["GBP4",]),
  SAMD9L = t(clusters_methylation["SAMD9L",]),
  CARD16 = t(clusters_methylation["CARD16",])
)


# Fit Cox model with methylation and age as covariates
fit <- coxph(Surv(time, status) ~ SAMD9L + age, data = cox_data)

# Show summary
summary(fit)

# Compare with other genes 

genes_to_test <- c("SOX10", "WIF1", "SFRP1", "FASN", "SREBF1", "MSMO1", 
  "HMGCS1", "HMGCR", "INSIG1", "FDFT1", "SCD", "SC5D")



# Plotting
gene_ref <- "OAS2"
par(mar=c(1,1,1,1))
plot_list <- list()

for (gene in genes_to_test) {
  
  my_df <- data.frame("Gene" = unname(t(fpkm_data[gene, colnames(clusters_methylation)])),
       "Methylation" = unname(t(clusters_methylation[gene_ref,])))
  
  plot_list[[gene]] <- ggplot(my_df, aes(x = Methylation, y = Gene)) +
    geom_boxplot(fill="lightgreen") + 
    ggtitle(paste("Plot for Gene:", gene)) +
    xlab("Methylation") +
    ylab("FPKM") +
    theme_bw()
  
}

wrap_plots(plot_list)



