#! usr/bin/Rscript

library(ggplot2)
library(patchwork)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(gprofiler2)
library(ggrepel)


#Loading gprofiler2#Loading data

# Load data
dist_cassette_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Loading enhancer data
enhancer_general_annotation <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_elements_v5.24.txt", header = T)
enhancer_general_annotation$chr <- paste0("chr", enhancer_general_annotation$chr)

enhancer_to_genes <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.24.txt", header = T)

# Loading cpg density data
cpg_density_5000 <- readRDS("/Users/isasiain/PhD/Projects/project_3/analysis/cpg_density/cpg_density_5000_bp.rds")
cpg_density_10000 <- readRDS("/Users/isasiain/PhD/Projects/project_3/analysis/cpg_density/cpg_density_10000_bp.rds")

# FPKM counts
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Generate obkects for genes linked to CpGs
genes <- sapply(annoObj$nameGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})


names(genes) <- annoObj$illuminaID


#
# OVERLAP OF CPGS WITH ENHANCER REGIONS
#

# DISTAL 1
cpgs_1 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 1]

matches_to_enhancers_1 <- sapply(cpgs_1, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_2)


# DISTAL 2
cpgs_2 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 2]

matches_to_enhancers_2 <- sapply(cpgs_2, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2)



# DISTAL 3
cpgs_3 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 3]

matches_to_enhancers_3 <- sapply(cpgs_3, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3)


# DISTAL 4
cpgs_4 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 4]

matches_to_enhancers_4 <- sapply(cpgs_4, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4)


# DISTAL 5
cpgs_5 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 5]

matches_to_enhancers_5 <- sapply(cpgs_5, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5)

# DISTAL 6
cpgs_6 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 6]

matches_to_enhancers_6 <- sapply(cpgs_6, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6)


# DISTAL 7
cpgs_7 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 7]

matches_to_enhancers_7 <- sapply(cpgs_7, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  filtered_enhancer_data <- filtered_enhancer_data[filtered_enhancer_data$regulatory_element_type %in% c("Enhancer", "Promoter/Enhancer"),]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)


# Calculate mapping percentages
percentages <- c(
  sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_1),
  sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2),
  sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3),
  sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4),
  sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5),
  sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6),
  sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)
) * 100  

# Create a data frame for plotting
df_plot <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages
)

# Plot
p_enhancers <- ggplot(df_plot, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping\nwith enhancers (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14) +  
  theme(axis.text = element_text(size = 12), 
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())


# 
#  REPETITIVE SEQUENCES PER CASSETTE
#


# Calculate mapping percentages of repeatsöä
percentages_repeats <- c(
  sum(annoObj[cpgs_1, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_1, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_2, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_2, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_3, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_3, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_4, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_4, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_5, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_5, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_6, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_6, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_7, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_7, "hasAnyRepeatOverlap"])
) * 100  

# Create a data frame for plotting
df_plot_repeats <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages_repeats
)

# Plot
p_repetitive <- ggplot(df_plot_repeats, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping with \nrepetitive sequences (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14)

p_enhancers / p_repetitive


#
# CpG DENSITY PER CASSETTE
#

# Get CpGs to analyse
cpgs_for_density <- names(dist_cassette_10$colors)[dist_cassette_10$colors %in% c(7,6,5,4,3,2,1)]

# Filter density values of interest
rownames(cpg_density_5000) <- cpg_density_5000$cpg_id
density_of_interest <-cpg_density_5000[cpgs_for_density,]

# Adding column for cassette
density_of_interest$cassette <- dist_cassette_10$colors[density_of_interest$cpg_id]
density_of_interest$cassette <- factor(density_of_interest$cassette, levels=c(7,6,5,4,3,2,1))

ggplot(density_of_interest, aes(x = cassette, y = cpg_density)) +
  geom_boxplot(outlier.shape = NA, fill = "grey", color = "black") +
  geom_jitter(width = 0.1, size = 0.001, alpha = 0.5, color = "black") +
  theme_bw(base_size = 14) +
  ylim(0,0.005) +
  labs(x = "Cassette", y = "CpG Density\n(CpGs / bp)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_flip()


#
# TF enrichment
#

# Proportions of TFs overlapping

tf_proportions <- rbind(
  "1" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 1],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 1]),
  "2" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 2],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 2]),
  "3" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 3],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 3]),
  "4" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 4],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 4]),
  "5" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 5],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 5]),
  "6" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 6],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 6]),
  "7" = colSums(tfMat[names(dist_cassette_10$colors)[dist_cassette_10$colors == 7],], na.rm = T) / length(names(dist_cassette_10$colors)[dist_cassette_10$colors == 7])

)

# Compute max absolute deviation from "All"

# Function to compute max pairwise absolute difference for each TF
max_pairwise_diff <- apply(tf_proportions, 2, function(x) {
  max(abs(outer(x, x, "-")))
})

# Select top 50 most differential TFs based on pairwise difference
top_tfs_pw <- names(sort(max_pairwise_diff, decreasing = TRUE)[1:35])

# Subset and melt the original proportions (excluding "All")
plot_data <- as.data.frame(tf_proportions[1:7, top_tfs_pw])
plot_data$Cassette <- factor(rownames(plot_data))
plot_data_long <- melt(plot_data, id.vars = "Cassette", variable.name = "TF", value.name = "Proportion")

ggplot(plot_data_long, aes(x = TF, y = Cassette)) +
  geom_point(aes(size = Proportion, color = Proportion)) +
  scale_color_gradient(low = "lightgrey", high = "darkblue") +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw(base_size = 14) +
    labs(x = "Transcription Factor", y = "Cassette") +
  coord_flip()



#
# MAPPING ENHANCERS TO GENES. ENRICHMENT
#

# Gatting hallmarks data
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- hallmarks[, c("gs_name", "entrez_gene")]

# CASSETTE 2
my_enhancers <- enhancer_to_genes[enhancer_to_genes$GHid %in% matches_to_enhancers_2,]
gene_list_cs <- unique(my_enhancers[my_enhancers$is_elite == 1, "symbol"])

# Convert first to ensembl ID (therev are 3 times more genes not being mapped if done directly)
ensembl_ids <- gconvert(gene_list_cs)$target

# Convert to ENTREZ-ID
gene_df_cs <- bitr(ensembl_ids,
     fromType = "ENSEMBL",
     toType = "ENTREZID",
     OrgDb = org.Hs.eg.db)

# kegg
output_cas_kegg <- enrichKEGG(gene_df_cs$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "fdr")

# go
output_cas_go <- enrichGO(gene_df_cs$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        pAdjustMethod = "fdr",
                        readable = T
                          )

# hallmarks
output_cas_hallmarks <- enricher(gene_df_cs$ENTREZID, 
                                  TERM2GENE = term2gene)




# Preprocessing and plotting hallmarks enrichment output
output_cas_kegg@result <- output_cas_kegg@result[output_cas_kegg@result$p.adjust <= 0.05,]

output_cas_kegg@result$GeneRatio_numeric <- sapply(output_cas_kegg@result$GeneRatio, function(val) {
  parts <- strsplit(val, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})

formatted_description <- tolower(gsub("_", " ", gsub("^HALLMARK_", "", output_cas_kegg@result$Description)))
formatted_description <- sub("^([a-z])", "\\U\\1", formatted_description, perl = TRUE)

output_cas_kegg@result$Description <- formatted_description

output_cas_kegg@result <- output_cas_kegg@result[order(output_cas_kegg@result$p.adjust, decreasing = TRUE),]
output_cas_kegg@result$Description <- factor(output_cas_kegg@result$Description, levels = output_cas_kegg@result$Description)

output_cas_kegg@result <- output_cas_kegg@result[order(output_cas_kegg@result$pvalue, decreasing = TRUE),]

plot_kegg <- ggplot(filter(output_cas_kegg@result, output_cas_kegg@result$p.adjust < 0.5), aes(x = -log10(p.adjust), y = Description)) +
  geom_point(aes(size = as.numeric(GeneRatio_numeric)), color = "black") +
  scale_size_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw() +
  xlim(1, 4) +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG",
    size = "Gene Ratio"
  )

# Preprocessing and plotting hallmarks enrichment output
output_cas_hallmarks@result <- output_cas_hallmarks@result[output_cas_hallmarks@result$p.adjust <= 0.05,]

output_cas_hallmarks@result$GeneRatio_numeric <- sapply(output_cas_hallmarks@result$GeneRatio, function(val) {
  parts <- strsplit(val, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})


formatted_description <- tolower(gsub("_", " ", gsub("^HALLMARK_", "", output_cas_hallmarks@result$Description)))
formatted_description <- sub("^([a-z])", "\\U\\1", formatted_description, perl = TRUE)

output_cas_hallmarks@result$Description <- formatted_description

output_cas_hallmarks@result <- output_cas_hallmarks@result[order(output_cas_hallmarks@result$p.adjust, decreasing = TRUE),]
output_cas_hallmarks@result$Description <- factor(output_cas_hallmarks@result$Description, levels = output_cas_hallmarks@result$Description)

output_cas_hallmarks@result <- output_cas_hallmarks@result[order(output_cas_hallmarks@result$pvalue, decreasing = TRUE),]

plot_hallmarks <- ggplot(filter(output_cas_hallmarks@result, output_cas_hallmarks@result$p.adjust < 0.5), aes(x = -log10(p.adjust), y = Description)) +
  geom_point(aes(size = as.numeric(GeneRatio_numeric)), color = "black") +
  scale_size_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw() +
  xlim(1, 4) +
  labs(
    x = "-log10(FDR p value)",
    y = "Hallmarks",
    size = "Gene Ratio"
  )


plot_kegg / plot_hallmarks +
  plot_layout(heights = c(1, 3), guides = "collect") &
  theme(legend.position = "bottom")


#
# EXPRESSION REGULATION. VIOLIN PLOT
#

# HCLUST BASED ON THE CASSETTE

# Define groups
betaAdj_cluster <- betaAdj[names(dist_cassette_10$colors)[dist_cassette_10$colors == 2], ]


# Transpose to cluster samples
d <- dist(t(betaAdj_cluster3))
hc <- hclust(d, method = "ward.D2")

# Cut into 2 groups
subgroups <- cutree(hc, k = 2)

Heatmap(betaAdj_cluster, 
        use_raster = F,
        column_split = subgroups,
        show_row_names = F,
        show_column_names = F)

# GETTING GENES WITH ENHANCERS LINKED TO THEM
my_enhancers <- enhancer_to_genes[enhancer_to_genes$GHid %in% matches_to_enhancers_2,]
gene_list_cs1 <- my_enhancers[my_enhancers$is_elite == 1, "symbol"]

genes_to_analyse <- unique(gene_list_cs1[gene_list_cs1 %in% rownames(fpkm_data)])


# ANALYSE AND PLOT

# Split samples into the two subgroups
group1_samples <- names(subgroups[subgroups == 1])
group2_samples <- names(subgroups[subgroups == 2])

pvals <- numeric(length(genes_to_analyse))
logFC <- numeric(length(genes_to_analyse))

for (i in seq_along(genes_to_analyse)) {
  g <- genes_to_analyse[i]
  g1 <- as.numeric(fpkm_data[g, group1_samples])
  g2 <- as.numeric(fpkm_data[g, group2_samples])
  
  if (length(unique(c(g1, g2))) > 1) {
    test <- wilcox.test(g1, g2)
    pvals[i] <- test$p.value
    logFC[i] <- log2(median(g1 + 1e-3) / median(g2 + 1e-3))  # median-based for Wilcoxon
  } else {
    pvals[i] <- NA
    logFC[i] <- NA
  }
}

# Getting data to plot
volcano_data <- data.frame(
  gene = genes_to_analyse,
  log2FC = logFC,
  pval = pvals,
  negLog10P = -log10(pvals)
)


# Prepare data
volcano_data$Bonferroni_p <- p.adjust(volcano_data$pval, method = "bonferroni")
volcano_data$neg_log10_p <- -log10(volcano_data$Bonferroni_p)
volcano_data$Log2_fold_change <- volcano_data$log2FC
volcano_data$Gene <- volcano_data$gene
volcano_data$CpG_Count <- sapply(volcano_data$gene, function(x) { # Number of CpGs mapped
  
  selected_enhancers <- enhancer_to_genes[enhancer_to_genes$symbol == x & enhancer_to_genes$is_elite == 1,]
  cpgs_of_interest <- names(matches_to_enhancers_2[matches_to_enhancers_2 %in% selected_enhancers$GHid])
  return(length(cpgs_of_interest))
  
})

# Set color: blue for significant, grey for non-significant
volcano_data$color <- ifelse(volcano_data$Bonferroni_p <= 0.05, "#1f78b4", "grey")

# Plot
ggplot(volcano_data, aes(x = Log2_fold_change, y = neg_log10_p, size = CpG_Count, color = color)) +
  geom_point(alpha = 0.7) + 
  geom_text_repel(
    data = subset(volcano_data, Bonferroni_p <= 1e-7),
    aes(label = Gene),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  scale_color_manual(values = c("#1f78b4", "grey")) +  
  scale_size(range = c(2, 6)) +                    
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


# Linking CpGs to genes regulated by enhancers
gene <- "FOXC1"
selected_enhancers <- enhancer_to_genes[enhancer_to_genes$symbol == gene & enhancer_to_genes$is_elite == 1,]
cpgs_of_interest <- names(matches_to_enhancers_2[matches_to_enhancers_2 %in% selected_enhancers$GHid])

# Generate bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM" = anno_barplot(as.numeric(fpkm_data[gene, colnames(betaAdj)]))
)

# Create top anotation
pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
im_annotation <- my_annotations[colnames(betaAdj), "IM"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])

top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    IM = im_annotation,
                                    Epitype = epi_annotation,
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)

# Cluster
cluster_methylation <- kmeans(t(betaAdj[cpgs_of_interest,]), centers = 2)

# Get CpGs annotated based on proximity
cpgs_from_proximity_annotation <- names(genes)[genes==gene]

# Generate CpG label
cpg_annotation_groups <- c(rep("Enhancer(s)", length(cpgs_of_interest)), 
  rep("Gene", length(cpgs_from_proximity_annotation)))

# Location of the CpG


interval_start <- 1609915
interval_end <- 1613897


distances <- sapply(annoObj[c(cpgs_of_interest, cpgs_from_proximity_annotation), "start"], 
                    function(pos, start = interval_start, end = interval_end) {
  if (pos < start) {
    return(pos - start)
  } else if (pos > end) {
    return((pos - end))
  } else {
    return(0)  # CpG is within the interval
  }
})


# Create row annotation with text
distance_annotation <- rowAnnotation(
  distance = anno_barplot(
    distances,
    gp = gpar(fill = "steelblue"),
    width = unit(2, "cm")
  )
)

# Add annotation to heatmap
Heatmap(betaAdj[c(cpgs_of_interest, cpgs_from_proximity_annotation), ],
        bottom_annotation = bottom_annotation,
        top_annotation = top_annotation,
        right_annotation = distance_annotation,
        column_split = subgroups,
        row_split = cpg_annotation_groups,
        show_column_names = FALSE,
        show_row_dend = F)

cpgs_of_interest

