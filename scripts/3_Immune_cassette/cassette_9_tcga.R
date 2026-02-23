#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(dplyr)
library(tidyverse)
library(grid)
library(UpSetR)

#
# LOADING DATA
#

set.seed(1234)

# TCGA-BRCA data
tcga <- new.env()
load("/Users/in2245sa/PhD/Projects/project_3/data/tcga_brca_withAnnotations.RData", envir = tcga)
load("/Users/in2245sa/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype.RData", envir = tcga)
tcga$annoObj[, "featureClass"] <- sapply(tcga$annoObj[, "featureClass"], function(val) {
  if (val == "distal" | val == "distal body") {"Distal"}
  else if (val == "proximal dn" | val == "proximal up") {"Proximal"}
  else {"Promoter"}
})


# SCANB TNBC discovery cohort
scanb <- new.env()
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData", envir=scanb)
scanb$annoObj[, "featureClass"] <- sapply(scanb$annoObj[, "featureClass"], function(val) {
  if (val == "distal" | val == "distal body") {"Distal"}
  else if (val == "proximal dn" | val == "proximal up") {"Proximal"}
  else {"Promoter"}
})

promoter_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_7.rds")

# Reading in silico TILs
tcga_tils <- read.csv2("/Users/in2245sa/PhD/Projects/project_3/data/tcga_tils_from_mpath.csv", skip = 1)
tcga_tils$sample_id <- sub("^([A-Z0-9-]+?-[A-Z0-9]+-[A-Z0-9]+)-.*$", "\\1", tcga_tils$Biospecimen_barcode)

# LINKING TCGA SAMPLES TO IN SILICO TILS
rownames(tcga_tils) <- tcga_tils$sample_id

tcga$sampleAnno[, "TILs"] <- tcga_tils[tcga$sampleAnno[, "bcr_patient_barcode"],"TIL_score"]

# Cassette 10
cpgs_9 <- names(promoter_7$colors)[promoter_7$colors == 9]
cpgs_9 <- cpgs_9[cpgs_9 %in% rownames(tcga$betaAdj)]

# PAM50 annotations (my.pam50.subtype)
load("/Users/in2245sa/PhD/Projects/project_3/data/jvc_PAM50_NCN_subtype_TCGA.RData")


#
# TILs VS TCGA SUBTYPE
#

# Prepare the data frame for ggplot
df <- data.frame(
  TILs = as.numeric(tcga$sampleAnno[, "TILs"]),
  Subtype = as.factor(my.pam50.subtype[rownames(tcga$sampleAnno)])
)

# Drop NAs if any
df <- na.omit(df)

# Drop unclassified samples
df <- df[!df$Subtype == "unclassified",]

# Custom colors
subtype_colors <- c(
  "Her2" = "purple",
  "LumA" = "darkblue",
  "Basal" = "indianred1",
  "LumB" = "lightblue",
  "Normal" = "green"
)

# Plot
ggplot(df, aes(x = Subtype, y = TILs, fill = Subtype)) +
  geom_boxplot(width = 0.7, outlier.shape = 21, alpha = 0.9, ) +
  scale_fill_manual(values = subtype_colors) +
  theme_bw(base_size = 11) +
  labs(
    x = "PAM50 Subtype",
    y = "TILs (%)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) + 
  stat_compare_means(
    method = "kruskal.test",
    label.x = 4.5,
    label.y = max(df$TILs, na.rm = TRUE) * 0.9,   
    label = "p.format"
  )


#
# ANALYSES IN TNBC
#


# Subset the beta matrix
beta_mat <- tcga$betaAdj[cpgs_9, ]
beta_mat <- beta_mat[,tcga$sampleAnno$TNBC]

# Cluster methylation state of cassette
clusters_cassette_9 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
cluster_means <- tapply(
  colMeans(beta_mat), 
  clusters_cassette_9,   
  mean                   
)
cluster_label <- ifelse(cluster_means[clusters_cassette_9] == min(cluster_means), "Hypo.", "Hyper.")

# COMPARE TILS BETWEEN CLUSTERS

# Prepare data for ggplot
tils_df <- data.frame(
  Sample = colnames(beta_mat),
  TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
  Cluster = cluster_label
)


ggplot(tils_df, aes(x = Cluster, y = TILs)) +
  geom_violin(, fill="black") +
  geom_boxplot(width=0.1) +
  theme_bw(base_size = 14)+
  stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Hypo.", "Hyper.")),
      label = "p.format",
      label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, 
      size = 5) +
  labs(,
    x = "Cluster",
    y = "TILs (%)"
  ) +
  ylim(0,max(tils_df$TILs, na.rm = TRUE) * 1.2) +
  theme(axis.title.y = ggtext::element_markdown())



# ANALYSE METHYLATION PATTERNS

# Generate annotations
subtypes <- factor(
  my.pam50.subtype[colnames(beta_mat)],
  levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
)

top_annotation <- HeatmapAnnotation(
  PAM50 = subtypes,
  TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
  #ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
  #PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
  #HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
  col = list(
    ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
    PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
    HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
    PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
  )
)


# Plotting heatmap
Heatmap(
  beta_mat,
  col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
  column_split = cluster_label,
  show_column_names = FALSE,
  show_row_names = FALSE,
  name = "Tumor\nBeta",
  top_annotation = top_annotation
)


# CHECK PER GENE
genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")

# Calculate and correct p values
p_vals_for_plot <- sapply(genes, function(gene) {
  
  my_cpgs <- tcga$betaAdj[tcga$annoObj[,"nameUCSCknownGeneOverlap"] == gene,]
  
  # Subset the beta matrix
  beta_mat <- my_cpgs[,tcga$sampleAnno$TNBC,drop=FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_9 <- as.factor(kmeans(t(beta_mat[rownames(beta_mat),, drop=FALSE]), 2)$cluster)
  cluster_means <- tapply(
    colMeans(beta_mat), 
    clusters_cassette_9,   
    mean                   
  )
  cluster_label <- ifelse(cluster_means[clusters_cassette_9] == min(cluster_means), "Hypo.", "Hyper.")
  
  # Analyse significance of GEX
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1,2]
  
  gex_data <- tcga$gexFpkm[grep(gene_ensembl,rownames(tcga$gexFpkm)), colnames(beta_mat)]
  
  gex_p <- wilcox.test(as.numeric(gex_data) ~ as.factor(cluster_label), na.action = na.omit)$p.value
  
  # Analyse significance of GEX
  tils_data <- tcga$sampleAnno[colnames(beta_mat), "TILs"]
  
  tils_p <- wilcox.test(as.numeric(tils_data) ~ as.factor(cluster_label), na.action = na.omit)$p.value
  
  return(c("GEX_p"=gex_p, "TILs_p"=tils_p))
  
})

# Adjust betas
p_vals_for_plot[1,] <- p.adjust(p_vals_for_plot[1,], method = "BH")
p_vals_for_plot[2,] <- p.adjust(p_vals_for_plot[2,], method = "BH")

# Plotting
list_of_tils_plots <- list()
list_of_gex_plots <- list()
list_of_heatmaps <- list()

for (gene in genes) {
  
  my_cpgs <- tcga$betaAdj[tcga$annoObj[,"nameUCSCknownGeneOverlap"] == gene,]
  
  # Subset the beta matrix
  beta_mat <- my_cpgs[,tcga$sampleAnno$TNBC,drop=FALSE]
  
  # Cluster methylation state of cassette
  clusters_cassette_9 <- as.factor(kmeans(t(beta_mat[rownames(beta_mat),, drop=FALSE]), 2)$cluster)
  cluster_means <- tapply(
    colMeans(beta_mat), 
    clusters_cassette_9,   
    mean                   
  )
  cluster_label <- ifelse(cluster_means[clusters_cassette_9] == min(cluster_means), "Hypo.", "Hyper.")
  
  # Getting gene expression (fpkm)
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1,2]
  
  gex_data <- tcga$gexFpkm[grep(gene_ensembl,rownames(tcga$gexFpkm)), colnames(beta_mat)]
  
  # Prepare data for ggplot
  tils_df <- data.frame(
    Sample = colnames(beta_mat),
    TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
    Cluster = cluster_label,
    GEX = gex_data
  )
  
  # Plotting TILS
  list_of_tils_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
    geom_violin(fill = "black") +
    geom_boxplot(width = 0.1) +
    theme_bw(base_size = 14) +
    annotate(
      "text",
      x = 1.5,
      y = max(tils_df$TILs, na.rm = TRUE) * 1.05,
      label = paste0("p = ", signif(p_vals_for_plot["TILs_p", gene], 3)),
      size = 5
    ) +
    labs(
      x = "Cluster",
      y = "TILs (%)"
    ) +
    ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.25) +
    theme(axis.title.y = ggtext::element_markdown())
  
  # Plotting Expression
  list_of_gex_plots[[gene]] <- ggplot(tils_df, aes(x = Cluster, y = GEX)) +
    geom_violin(, fill="black") +
    geom_boxplot(width=0.1) +
    theme_bw(base_size = 14)+
    annotate(
      "text",
      x = 1.5,
      y = max(tils_df$GEX, na.rm = TRUE) * 1.05,
      label = paste0("p = ", signif(p_vals_for_plot["GEX_p", gene], 3)),
      size = 5
    ) +
    labs(,
         x = "Cluster",
         y = paste0(gene, " FPKM")
    ) +
    ylim(0, max(tils_df$GEX, na.rm = TRUE) * 1.35) +
    theme(axis.title.y = ggtext::element_markdown())
  
  
  # Generate annotations
  subtypes <- factor(
    my.pam50.subtype[colnames(beta_mat)],
    levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
  )
  
  #Getting betas to plot
  cpgs_overlapping <- tcga$annoObj[tcga$annoObj$nameUCSCknownGeneOverlap == gene, "illuminaID"]
  betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_mat)]
  
  top_annotation <- HeatmapAnnotation(
    PAM50 = subtypes,
    TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
    #ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
    #PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
    #HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
    col = list(
      ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
      PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
      HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
      PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
    )
  )
  
  bottom_annotation <- HeatmapAnnotation("FPKM" = anno_barplot(tils_df[colnames(beta_mat), "GEX"]))
  

  
  right_annotation = rowAnnotation("In Cassette" = cpgs_overlapping %in% cpgs_9,
                                        col = list(
                                          "In Cassette" = c("TRUE" = "black", 
                                                            "FALSE" = "white")
                                          )
                                        )
  
  # Plotting heatmap
  list_of_heatmaps[[gene]] <- Heatmap(
    betas_to_plot,
    col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
    column_split = cluster_label,
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    left_annotation = right_annotation
  )
  
}

#Plotting tils vs methylation state
list_of_tils_plots[[1]] | list_of_tils_plots[[2]] | list_of_tils_plots[[3]] | list_of_tils_plots[[4]]

# Plotting gex vs methylation state
list_of_gex_plots[[1]] | list_of_gex_plots[[2]] | list_of_gex_plots[[3]] | list_of_gex_plots[[4]] 

#Heatmaps of genes
heatmap_grobs <- lapply(list_of_heatmaps, function(ht) {
  grid.grabExpr(draw(ht,
                     show_heatmap_legend = FALSE,
                     show_annotation_legend = FALSE))
})

wrap_plots(heatmap_grobs, ncol = 5)



#
# ANALYSES IN ALL BRCA
#

# Subset the beta matrix
beta_mat <- tcga$betaAdj[cpgs_9, ,drop=FALSE]
#beta_mat <- beta_mat[, my.pam50.subtype[colnames(beta_mat)] == "LumA",drop=FALSE]

# Cluster methylation state of cassette
clusters_cassette_9 <- as.factor(kmeans(t(beta_mat), 2)$cluster)
cluster_means <- tapply(
  colMeans(beta_mat), 
  clusters_cassette_9,   
  mean                   
)
cluster_label <- ifelse(cluster_means[clusters_cassette_9] == min(cluster_means), "Hypo.", "Hyper.")

# COMPARE TILS BETWEEN CLUSTERS

# Prepare data for ggplot
tils_df <- data.frame(
  Sample = colnames(beta_mat),
  TILs = tcga$sampleAnno[colnames(beta_mat), "TILs"],
  Cluster = cluster_label
)


ggplot(tils_df, aes(x = Cluster, y = TILs)) +
  geom_boxplot(fill="lightgrey") +
  theme_bw(base_size = 14)+
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Hypo.", "Hyper.")),
    label = "p.format",
    label.y = max(tils_df$TILs, na.rm = TRUE) * 1.05, # position above violins
    size = 5) +
  labs(,
       x = "Cluster",
       y = "TILs (%)"
  ) +
  ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.2) +
  theme(axis.title.y = ggtext::element_markdown())


# ANALYSE METHYLATION PATTERNS

# Generate annotations
subtypes <- factor(
  my.pam50.subtype[colnames(beta_mat)],
  levels = c("unclassified","Normal","LumA","LumB","Her2","Basal")
)

top_annotation <- HeatmapAnnotation(
  PAM50 = subtypes,
  TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"]), # dot size
  ER = tcga$sampleAnno[colnames(beta_mat),"ER"],
  PR = tcga$sampleAnno[colnames(beta_mat),"PR"],
  HER2 = tcga$sampleAnno[colnames(beta_mat),"HER2"],
  col = list(
    ER = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey"),
    PR = c("Negative"="white", "Positive"="black", "[Not Evaluated]"="grey", "Indeterminate"="grey"),
    HER2 = c("Negative"="white", "Positive"="black", "NA"="grey", "Indeterminate"="grey", "Equivocal"="grey"),
    PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
  )
)


# Plotting heatmap
Heatmap(
  beta_mat,
  col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
  column_split = cluster_label,
  show_column_names = FALSE,
  show_row_names = FALSE,
  name = "Tumor\nBeta",
  top_annotation = top_annotation
)


# CLUSTERING BASED ON ALL CPGS

genes <- c("GBP4", "OAS2", "ZBP1", "CARD16")
subtypes <- c("LumA", "LumB", "Her2", "Basal", "Normal")

# Analyse significance for plotting
p_vals_for_plot <- lapply(genes, function(gene) {
  
  # --- Select CpGs ---
  my_cpgs <- rownames(tcga$annoObj)[
    tcga$annoObj$nameUCSCknownGeneOverlap == gene &
      tcga$annoObj$featureClass %in% c("Promoter", "Proximal", "Distal")
  ]
  
  beta_mat <- tcga$betaAdj[my_cpgs, , drop = FALSE]
  
  # --- Cluster methylation ---
  clusters_gene <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(colMeans(beta_mat), clusters_gene, mean)
  cluster_label <- ifelse(
    cluster_means[clusters_gene] == min(cluster_means),
    "Hypo.",
    "Hyper."
  )
  names(cluster_label) <- names(clusters_gene)
  
  # --- Expression ---
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1, 2]
  
  gex_data <- tcga$gexFpkm[
    grep(gene_ensembl, rownames(tcga$gexFpkm)),
    colnames(beta_mat)
  ]
  
  tils_data <- tcga$sampleAnno[colnames(beta_mat), "TILs"]
  subtype_vec <- my.pam50.subtype[colnames(beta_mat)]
  
  # Store subtype p-values
  gex_pvals  <- numeric(length(subtypes))
  tils_pvals <- numeric(length(subtypes))
  names(gex_pvals)  <- subtypes
  names(tils_pvals) <- subtypes
  
  for (subtype in subtypes) {
    
    idx <- subtype_vec == subtype
    
    subset_gex  <- as.numeric(gex_data[idx])
    subset_tils <- as.numeric(tils_data[idx])
    subset_cluster <- as.factor(cluster_label[idx])
    
    # Only test if both groups present
    if (length(unique(subset_cluster)) == 2) {
      gex_pvals[subtype]  <- wilcox.test(subset_gex ~ subset_cluster, exact = FALSE)$p.value
      tils_pvals[subtype] <- wilcox.test(subset_tils ~ subset_cluster, exact = FALSE)$p.value
    } else {
      gex_pvals[subtype]  <- NA
      tils_pvals[subtype] <- NA
    }
  }
  
  # ---- BH correction per gene ----
  gex_pvals_adj  <- p.adjust(gex_pvals, method = "BH")
  tils_pvals_adj <- p.adjust(tils_pvals, method = "BH")
  
  return(list(
    gex_raw  = gex_pvals,
    tils_raw = tils_pvals,
    gex_adj  = gex_pvals_adj,
    tils_adj = tils_pvals_adj
  ))
})

names(p_vals_for_plot) <- genes



# To store plots
list_of_tils_plots <- list()
list_of_gex_plots <- list()
list_of_heatmaps <- list()
list_of_barplots <- list()

# To store proportion of hypomethylation
hypo_prop_summary <- data.frame()
all_data_long_df <- data.frame()

for (gene in genes) {
  
  # Select all cpgs per gene
  my_cpgs <- rownames(tcga$annoObj)[
    tcga$annoObj$nameUCSCknownGeneOverlap == gene &
      tcga$annoObj$featureClass %in% c("Promoter", "Proximal", "Distal")
  ]
  
  # Subset beta matrix
  beta_mat <- tcga$betaAdj[my_cpgs, , drop = FALSE]
  
  # Cluster methylation pattern
  clusters_gene <- as.factor(kmeans(t(beta_mat), 2)$cluster)
  cluster_means <- tapply(colMeans(beta_mat), clusters_gene, mean)
  cluster_label <- ifelse(
    cluster_means[clusters_gene] == min(cluster_means),
    "Hypo.",
    "Hyper."
  )
  names(cluster_label) <- names(clusters_gene)
  
  # Get expression data
  gene_ensembl <- bitr(
    gene,
    fromType = "SYMBOL",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db
  )[1, 2]
  
  gex_data <- tcga$gexFpkm[grep(gene_ensembl, rownames(tcga$gexFpkm)), colnames(beta_mat)]
  
  # Prepare heatmap
  cpgs_overlapping <- my_cpgs
  betas_to_plot <- tcga$betaAdj[cpgs_overlapping, colnames(beta_mat), drop = FALSE]
  
  top_annotation <- HeatmapAnnotation(
    TILs = anno_points(tcga$sampleAnno[colnames(beta_mat), "TILs"], 
                       size = unit(1, "mm")),
    PAM50 = my.pam50.subtype[colnames(beta_mat)],
    ER = tcga$sampleAnno[colnames(beta_mat), "ER"],
    PR = tcga$sampleAnno[colnames(beta_mat), "PR"],
    HER2 = tcga$sampleAnno[colnames(beta_mat), "HER2"],
    col = list(
      ER = c("Negative" = "white", "Positive" = "black", "[Not Evaluated]" = "grey"),
      PR = c("Negative" = "white", "Positive" = "black", "[Not Evaluated]" = "grey", "Indeterminate" = "grey"),
      HER2 = c("Negative" = "white", "Positive" = "black", "NA" = "grey", "Indeterminate" = "grey", "Equivocal" = "grey"),
      PAM50 = c("Her2"="purple","LumA"="darkblue","Basal"="indianred1","LumB"="lightblue","Normal"="green","unclassified"="grey")
    )
  )
  
  bottom_annotation <- HeatmapAnnotation(GEX = anno_barplot(gex_data[colnames(tcga$betaAdj)]))
  
  row_annotation <- rowAnnotation(
    "CpG in Cassette" = rownames(beta_mat) %in% cpgs_9,
    col = list(
      "CpG in Cassette" = c("TRUE" = "black", "FALSE" = "white")
    )
  )
  
  list_of_heatmaps[[gene]]<- Heatmap(
    betas_to_plot,
    col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
    column_split = cluster_label[colnames(beta_mat)],
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    name = "Tumor\nBeta",
    top_annotation = top_annotation,
    left_annotation = row_annotation,
    bottom_annotation = bottom_annotation
  )
  
  # Iterate over subtypes
  for (subtype in subtypes) {
    samples_in_subtype <- names(my.pam50.subtype[colnames(beta_mat)][my.pam50.subtype[colnames(beta_mat)] == subtype])
    
    # Prepare dataframe
    tils_df <- data.frame(
      Gene = gene,
      Subtype = subtype,
      Sample = samples_in_subtype,
      TILs = tcga$sampleAnno[samples_in_subtype, "TILs"],
      Cluster = cluster_label[samples_in_subtype],
      GEX = gex_data[samples_in_subtype]
    )
    
    all_data_long_df <- rbind(all_data_long_df, tils_df)
    
    # TILs plots
    list_of_tils_plots[[gene]][[subtype]] <- ggplot(tils_df, aes(x = Cluster, y = TILs)) +
      geom_violin(fill = "black", size = 0.1) +
      geom_boxplot(width = 0.1, size = 0.1) +
      theme_bw(base_size = 14) +
      annotate(
        "text",
        x = 1.5,
        y = max(tils_df$TILs, na.rm = TRUE) * 1.05,
        label = paste0("p = ", signif(p_vals_for_plot[[gene]]$tils_adj[subtype], 3)),
        size = 5
      ) +
      labs(x = "Cluster", y = "TILs (%)") +
      ylim(0, max(tils_df$TILs, na.rm = TRUE) * 1.2) +
      theme(axis.title.y = ggtext::element_markdown())
    
    # GEX plots
    list_of_gex_plots[[gene]][[subtype]] <- ggplot(tils_df, aes(x = Cluster, y = GEX)) +
      geom_violin(fill = "black", size = 0.1) +
      geom_boxplot(width = 0.1, size = 0.1) +
      theme_bw(base_size = 14) +
      annotate(
        "text",
        x = 1.5,
        y = max(tils_df$GEX, na.rm = TRUE) * 1.05,
        label = paste0("p = ", signif(p_vals_for_plot[[gene]]$gex_adj[subtype], 3)),
        size = 5
      ) +
      labs(x = "Cluster", y = paste0(gene, " FPKM")) +
      ylim(0, max(tils_df$GEX, na.rm = TRUE) * 1.2) +
      theme(axis.title.y = ggtext::element_markdown())
    
    
    # Proportion of hypomethylation per subtype
    n_sub <- length(samples_in_subtype)
    n_hypo <- sum(cluster_label[samples_in_subtype] == "Hypo.", na.rm = TRUE)
    prop_df <- data.frame(
      Gene = gene,
      PAM50 = subtype,
      n = n_sub,
      n_Hypo. = n_hypo,
      Hypo. = if (n_sub > 0) n_hypo / n_sub else NA,
      stringsAsFactors = FALSE
    )
    hypo_prop_summary <- rbind(hypo_prop_summary, prop_df)
  }
}

# Plotting TILs vs subtype and gene
(list_of_tils_plots$CARD16$Basal | list_of_tils_plots$CARD16$Her2 | list_of_tils_plots$CARD16$LumA  | list_of_tils_plots$CARD16$LumB | list_of_tils_plots$CARD16$Normal)/
  (list_of_tils_plots$GBP4$Basal | list_of_tils_plots$GBP4$Her2 | list_of_tils_plots$GBP4$LumA  | list_of_tils_plots$GBP4$LumB | list_of_tils_plots$GBP4$Normal)/
  (list_of_tils_plots$OAS2$Basal | list_of_tils_plots$OAS2$Her2 | list_of_tils_plots$OAS2$LumA  | list_of_tils_plots$OAS2$LumB | list_of_tils_plots$OAS2$Normal)/
  (list_of_tils_plots$ZBP1$Basal | list_of_tils_plots$ZBP1$Her2 | list_of_tils_plots$ZBP1$LumA  | list_of_tils_plots$ZBP1$LumB | list_of_tils_plots$ZBP1$Normal)


# Plotting GEX vs subtype and gene methylation
(list_of_gex_plots$CARD16$Basal | list_of_gex_plots$CARD16$Her2 | list_of_gex_plots$CARD16$LumA  | list_of_gex_plots$CARD16$LumB | list_of_gex_plots$CARD16$Normal)/
  (list_of_gex_plots$GBP4$Basal | list_of_gex_plots$GBP4$Her2 | list_of_gex_plots$GBP4$LumA  | list_of_gex_plots$GBP4$LumB | list_of_gex_plots$GBP4$Normal)/
  (list_of_gex_plots$OAS2$Basal | list_of_gex_plots$OAS2$Her2 | list_of_gex_plots$OAS2$LumA  | list_of_gex_plots$OAS2$LumB | list_of_gex_plots$OAS2$Normal)/
  (list_of_gex_plots$ZBP1$Basal | list_of_gex_plots$ZBP1$Her2 | list_of_gex_plots$ZBP1$LumA  | list_of_gex_plots$ZBP1$LumB | list_of_gex_plots$ZBP1$Normal)


# Plotting piecharts
pie_data <- hypo_prop_summary %>%
  mutate(
    Hyper. = 1 - Hypo.
  ) %>%
  tidyr::pivot_longer(cols = c("Hypo.", "Hyper."), 
                      names_to = "State", values_to = "Proportion")


ggplot(pie_data, aes(x = "", y = Proportion, fill = State)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_grid(Gene ~ PAM50) +
  scale_fill_manual(
    values = c("Hypo." = "#5F4B8B", "Hyper." = "#C9A441"),
    name = NULL
  ) +
  theme_void(base_size = 14) +
  labs(
  ) +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) 


# Analyse statistical significance
counts_df <- all_data_long_df %>%
  dplyr::count(Gene, PAM50 = Subtype, Cluster) %>%
  tidyr::pivot_wider(names_from = PAM50, values_from = n, values_fill = 0)

chi_results <- all_data_long_df %>%
  group_by(Gene) %>%
  summarise(
    p_value = tryCatch({
      tbl <- table(Cluster, Subtype)
      chisq.test(tbl)$p.value
    }, error = function(e) NA)
  )

# Adjust p values
chi_results$p_adj <- p.adjust(chi_results$p_value, method = "BH")

print(chi_results)

#Heatmaps
heatmap_grobs <- lapply(list_of_heatmaps, function(ht) {
  grid.grabExpr(draw(ht,
                show_heatmap_legend = FALSE,
                show_annotation_legend = FALSE))
})

wrap_plots(heatmap_grobs, ncol = 4)


# Plotting co-ocurrence of hypomethylation
hyper_matrix <- all_data_long_df %>%
  mutate(Hypo = ifelse(Cluster == "Hyper.", 1, 0)) %>%
  dplyr::select(Sample, Gene, Hypo) %>%
  pivot_wider(names_from = Gene, values_from = Hypo, values_fill = 0)

# Remove Sample column for UpSet
matrix_upset <- hyper_matrix %>% dplyr::select(-Sample)

c_matrix <- make_comb_mat(matrix_upset)

# Plot UpSet plot
comb_sets = lapply(comb_name(c_matrix), function(nm) extract_comb(c_matrix, nm))

cols_for_plot <- setNames(c("#C9A441",
           "#E8D89A", "#E8D89A", "#E8D89A","#E8D89A",
           "black","black","black","black","black", "black",
           "#B9B2D8", "#B9B2D8", "#B9B2D8", "#B9B2D8",
           "#5F4B8B"),
         comb_name(c_matrix))


# Observed frequencies from your UpSet combination matrix
observed <- comb_size(c_matrix)

# Expected frequencies assuming independence
gene_freq <- colMeans(matrix_upset)  # P(Hypo) for each gene

# Expected number of samples per combination (in same order as comb_name(c_matrix))
expected <- sapply(strsplit(comb_name(c_matrix), ""), function(bits) {
  probs <- ifelse(bits == "1", gene_freq, 1 - gene_freq)
  prod(probs) * nrow(matrix_upset)
})

# Avoid divisions by zero
expected[expected == 0] <- 1e-10

# Compute log2 enrichment ratio
log2_enrichment <- log2(observed / expected)

# Compute binomial p-values for each combination
pvals <- mapply(function(obs, exp) {
  binom.test(x = obs, n = nrow(matrix_upset), p = exp / nrow(matrix_upset), alternative = "two.sided")$p.value
}, observed, expected)

# Adjust for multiple testing
padj <- p.adjust(pvals, method = "BH")

# Combine results
enrichment_df <- data.frame(
  combination = comb_name(c_matrix),
  observed = observed,
  expected = round(expected, 2),
  log2_enrichment = round(log2_enrichment, 2),
  pval = signif(pvals, 3),
  padj = signif(padj, 3)
)

# Add labels for significance
enrichment_df$category <- ifelse(padj < 0.05 & log2_enrichment > 0, "enriched",
                                 ifelse(padj < 0.05 & log2_enrichment < 0, "depleted", "ns"))

# Order by enrichment
enrichment_df <- enrichment_df[order(-enrichment_df$log2_enrichment), ]
enrichment_df$stars <- cut(
  enrichment_df$padj,
  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
  labels = c("****", "***", "**", "*", "")
)

# Match stars to the combination order in c_matrix
stars_ordered <- enrichment_df$stars[match(comb_name(c_matrix), enrichment_df$combination)]

bottom_annotation = HeatmapAnnotation(
  signif = anno_text(
    stars_ordered,
    gp = gpar(col = ifelse(enrichment_df$category[match(comb_name(c_matrix), enrichment_df$combination)] == "enriched",
                           "black", "darkgrey"), fontsize = 14),
    location = 0.5
  ),
  TILs = anno_boxplot(
    lapply(seq_along(comb_sets), function(i) {
      my_samples <- hyper_matrix$Sample[comb_sets[[i]]]
      tcga$sampleAnno[my_samples, "TILs"]
    }),
    outline = FALSE
  ),
  annotation_name_side = "left"
)

UpSet(c_matrix,
      bottom_annotation=bottom_annotation,
      comb_col = cols_for_plot
)

