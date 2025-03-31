# HEATMAPS OF PARTICULAR GENES

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(ggsignif)


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


#
# PLOTTING CPGS AFFECTING 
#<




current_gene_id = "GBP4"

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
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[names(genes)[genes == current_gene_id],]),
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")))
)

# CpG context annotation
left_annotation <- rowAnnotation("Context"= annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]]
)

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
  show_row_names = FALSE,
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
  Methylation_State = as.factor(ifelse(cluster_gbp4_promoter$cluster == 1, "Hypermethylated", "Hypomethylated")
  )
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


# PDL1 CPS BOXPLOT

ggplot(na.omit(plot_data), aes(x = as.factor(Methylation_State), fill = as.factor(PDL1_CPS))) +
  geom_bar(position = "fill", alpha = 0.8) +  # Stacked bar normalized to proportions
  theme_classic(base_size = 14) +  # Clean theme with larger text
  labs(x = "Promoter Methylation State", y = "Sample proportion", fill = "PDL1 CPS") +
  theme(legend.position = "right")  # Legend on the right

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
  


