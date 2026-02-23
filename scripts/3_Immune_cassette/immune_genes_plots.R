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
library(corrplot)
library(rstatix)

#
# LOAD DATA
#

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Loading gene expression
fpkm_data <- read.table("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Epitype annotations
my_annotations <- read.table("/Users/in2245sa//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# Reading PDL1 TLS annotation. Matrix is called u.frame
load("/Users/in2245sa/PhD/Projects/project_3/data/summarized_TPS_data_sampleLevel.RData")
rownames(u.frame) <- u.frame$PD_ID

# Loading promoter cassettes
promoter_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_7.rds")

# Cibersort output
tum_cibersort <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/cell_type_deconvolution/CIBERSORTx_tumour.csv")
rownames(tum_cibersort) <- tum_cibersort$Mixture

# TMArQ counts
tmarq_counts <- read.csv2("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/supplData_withimages.csv")

#
# PLOTTING IMMUNE CASSETTE (Cassette 9)
#

set.seed(123)

data <- promoter_7$colors
cpgs <- names(data)[data == 9]


pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
pam50_annotations <- ifelse(pam50_annotations == "Uncl.", 
                            "Uncl.", 
                            ifelse(pam50_annotations == "Basal", 
                                   "Basal", 
                                   "Non-Basal"))
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])


# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    TILs = anno_points(tils_annotation,
                                                       ylim=c(0,100),
                                                       size=unit(0.75, "mm"),
                                                       axis_param = list(
                                                         side="left",
                                                         at=c(0,25,50,75,100),
                                                         labels=c("0","25","50","75","100")
                                                       )),
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey")
                                    )
)


# Updated left_annotation with color scale
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[cpgs,]),
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("#372B63", "#EFEFEA", "#8F7417")))
)


# Cluster based on methylation
cluster_promoter <- kmeans(t(betaAdj[cpgs,]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[cpgs,cluster_promoter$cluster==1]) >
                      mean(betaAdj[cpgs,cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Heatmap of genes
Heatmap(
  betaAdj[cpgs,],
  col=colorRamp2(c(0, 0.5, 1), c("#6A5FB5", "#FFFFFF", "#E1C04A")),
  cluster_rows = FALSE,
  row_order = order(annoObj$start[annoObj$illuminaID %in% cpgs]),
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = promoter_state,
  top_annotation = top_annotation,
  #right_annotation = right_annotation,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  name = "Tumor beta"
)

#
# GETTING GENES TO ANALYSE
#

# Getting genes to analyse
genes_in_cassette <- unname(genes[cpgs])

# Spliting gene id in two if associated to two genes
genes_in_cassette <- unlist(strsplit(genes_in_cassette, "_"))

# Getting only genes with 2 or more cpgs
genes_in_cassette <- names(table(genes_in_cassette)[table(genes_in_cassette) >= 2])

# Checking if methylation patterns are associated to gene expression changes and different TILs
results <- data.frame(
  gene = character(),
  p_gex = numeric(),
  p_tils = numeric(),
  stringsAsFactors = FALSE
)

for (gene in genes_in_cassette) {
  
  gene_symbol <- gene
  if (gene == "CYREN") gene_symbol <- "C7orf49"
  
  # Get CpGs mapping to this gene
  cpgs <- names(genes[grepl(gene, genes)])
  
  if (length(cpgs) < 2) next
  if (!gene_symbol %in% rownames(fpkm_data)) next
  
  # Cluster samples based on methylation
  km <- kmeans(t(betaAdj[cpgs, ]), centers = 2)
  clusters <- km$cluster
  
  # Match sample order
  common_samples <- intersect(names(clusters), colnames(fpkm_data))
  if (length(common_samples) < 5) next
  
  clusters <- clusters[common_samples]
  
  # --- Gene expression test ---
  gex <- as.numeric(fpkm_data[gene_symbol, common_samples])
  p_gex <- tryCatch(
    wilcox.test(gex ~ clusters)$p.value,
    error = function(e) NA
  )
  
  # --- TILs test ---
  tils <- as.numeric(x[common_samples, "TILs"])
  p_tils <- tryCatch(
    wilcox.test(tils ~ clusters)$p.value,
    error = function(e) NA
  )
  
  results <- rbind(results,
                   data.frame(gene = gene,
                              p_gex = p_gex,
                              p_tils = p_tils))
}

# Multiple testing correction
results$padj_gex  <- p.adjust(results$p_gex,  method = "bonferroni")
results$padj_tils <- p.adjust(results$p_tils, method = "bonferroni")

# Keep genes significant for BOTH
significant_genes <- results[
  results$padj_gex < 0.05 & results$padj_tils < 0.05,
]

print(significant_genes)

#
# PLOTTING CPGS AFFECTING 
#

my_genes <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")
heatmap_list <- list()

for (current_gene_id in my_genes) {
  
  # Create top anotation
  top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                      TNBC = tnbc_annotation,
                                      HRD = HRD_annotation,
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
                                        "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                        "HRD"=c("High"="black", "Low/Inter"="grey"),
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
    #"ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
    col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
               "ATAC" = c("0" = "white", "1"= "black"))
  )
  
  # CpG annotation. Context and included in cassette
  left_annotation <- rowAnnotation("CpG_in_cassette" = names(genes)[genes == current_gene_id] %in% names(promoter_7$colors)[promoter_7$colors == 9],
                                   "Context"= annoObj$CpG_context[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                                   col=list("Context"=c("Distal" = "#f8766d", 
                                                        "Promoter" = "#00ba38", 
                                                        "Proximal" = "#619cff"),
                                            "CpG_in_cassette"=c("TRUE" = "black", 
                                                                "FALSE" = "white")))
  
  # Cluster based on methylation
  cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                   annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])
  
  cluster_promoter <- kmeans(t(betaAdj[names(cpgs),]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  promoter_state <- if (mean(betaAdj[names(cpgs),cluster_promoter$cluster==1]) >
                        mean(betaAdj[names(cpgs),cluster_promoter$cluster==2])) {
    
    as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
    
  } else {
    
    as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
    
  }
  
  
  # Heatmap of genes
  heatmap_list[[current_gene_id]] <- Heatmap(
    betaAdj[names(genes)[genes == current_gene_id],],
    col=colorRamp2(c(0, 0.5, 1), c("#5F4B8B", "#F5F5F2", "#C9A441")),
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_row_dend = FALSE,
    column_split = promoter_state,
    top_annotation = top_annotation,
    bottom_annotation = bottom_annotation,
    #right_annotation = right_annotation,
    left_annotation = left_annotation,
    clustering_distance_columns =  "euclidean",
    clustering_method_columns = "ward.D2",
    clustering_distance_rows =  "euclidean",
    clustering_method_rows = "ward.D2",
    name = "Tumor beta"
  )
  
  
}

#Heatmaps
heatmap_grobs <- lapply(heatmap_list, function(ht) {
  grid.grabExpr(draw(ht,
                     show_heatmap_legend = FALSE,
                     show_annotation_legend = FALSE))
})

wrap_plots(heatmap_grobs, ncol = 5)


#
# PLOTTING BOXPLOTS. TILs AND GENE EXPRESSION
#
  
  # Convert data to a dataframe for ggplot2
  plot_data <- data.frame(
    GBP4_FPKM = as.numeric(fpkm_data["GBP4", colnames(betaAdj),]),
    ZBP1_FPKM = as.numeric(fpkm_data["ZBP1", colnames(betaAdj),]),
    OAS2_FPKM = as.numeric(fpkm_data["OAS2", colnames(betaAdj),]),
    CARD16_FPKM = as.numeric(fpkm_data["CARD16", colnames(betaAdj),]),
    SAMD9L_FPKM = as.numeric(fpkm_data["SAMD9L", colnames(betaAdj),]),
    TILs = as.numeric(x[colnames(betaAdj), "TILs"]),
    PDL1_CPS = as.numeric(x[colnames(betaAdj), "PDL1_CPS"]),
    PDL1_TPS = sapply(u.frame[colnames(betaAdj), "PDL1_TPS"], function(x) {
      if (is.na(x)) {NA}
      else if (x >= 10) {2}
      else if (x >= 1) {1}
      else {0}
    }
    ),
    GBP4_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "GBP4"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "GBP4"]])
      km <- kmeans(t(betaAdj[names(cpgs), ]), centers = 2)
      if (mean(betaAdj[names(cpgs), km$cluster == 1]) >
          mean(betaAdj[names(cpgs), km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    OAS2_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "OAS2"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "OAS2"]])
      km <- kmeans(t(betaAdj[names(cpgs), ]), centers = 2)
      if (mean(betaAdj[names(cpgs), km$cluster == 1]) >
          mean(betaAdj[names(cpgs), km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    ZBP1_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "ZBP1"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "ZBP1"]])
      km <- kmeans(t(betaAdj[names(cpgs), ]), centers = 2)
      if (mean(betaAdj[names(cpgs), km$cluster == 1]) >
          mean(betaAdj[names(cpgs), km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    CARD16_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "CARD16"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "CARD16"]])
      km <- kmeans(t(betaAdj[names(cpgs), ]), centers = 2)
      if (mean(betaAdj[names(cpgs), km$cluster == 1]) >
          mean(betaAdj[names(cpgs), km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo"))
      }
    }),
    
    SAMD9L_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "SAMD9L"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "SAMD9L"]])
      km <- kmeans(t(betaAdj[names(cpgs), ]), centers = 2)
      if (mean(betaAdj[names(cpgs), km$cluster == 1]) >
          mean(betaAdj[names(cpgs), km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }))
  
  
  # EXPRESSION BOXPLOT
  
  # Determine p values BH correction
  pvals <- sapply(my_genes, function(gene) {
    test <- wilcox.test(
      plot_data[,paste0(gene, "_FPKM")] ~ plot_data[,paste0(gene, "_Methylation_State")]
    )
    test$p.value
  })
  
  # Apply BH correction
  pvals_adj <- p.adjust(pvals, method = "BH")
  
  # Function to make plots
  make_expression_boxplot <- make_expression_boxplot <- function(gene, y_label, x_label, ymax, pval_adj) {
    ggplot(plot_data, aes_string(x = paste0(gene, "_Methylation_State"), y = paste0(gene, "_FPKM"))) +
      geom_violin(fill="black", size = 0.2) +
      geom_boxplot(outlier.shape = NA, width=0.15, fill="grey90", col="grey40", size = 0.2) +
      theme_bw(base_size = 14) +
      labs(x = x_label, y = y_label) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0)) +
      annotate("text", x = 1.5, y = ymax, label = paste0("adj.p = ", signif(pval_adj, 3)))
  }
  
  GBP4_plot <- make_expression_boxplot("GBP4", "FPKM", NULL, 54, pvals_adj["GBP4"])
  OAS2_plot <- make_expression_boxplot("OAS2", "FPKM", NULL, 98, pvals_adj["OAS2"])
  ZBP1_plot <- make_expression_boxplot("ZBP1", "FPKM", NULL, 19, pvals_adj["ZBP1"])
  CARD16_plot <- make_expression_boxplot("CARD16", "FPKM", NULL, 37, pvals_adj["CARD16"])
  SAMD9L_plot <- make_expression_boxplot("SAMD9L", "FPKM", NULL, 41, pvals_adj["SAMD9L"])

  
  # Combine plots vertically
  (GBP4_plot | OAS2_plot | ZBP1_plot | CARD16_plot | SAMD9L_plot) + plot_layout(guides = "collect")
  
  
  # TILs BOXPLOT
  
  # Compute p-values for TILs association with methylation state
  pvals_tils <- sapply(my_genes, function(gene) {
    test <- wilcox.test(
      plot_data[["TILs"]] ~ plot_data[[paste0(gene, "_Methylation_State")]]
    )
    test$p.value
  })
  
  # BH correction
  pvals_tils_adj <- p.adjust(pvals_tils, method = "BH")
  
  # Function to make TILs boxplots with adjusted p-values
  make_tils_boxplot <- function(gene, y_label, x_label, ymax, pval_adj) {
    ggplot(plot_data, aes_string(x = paste0(gene, "_Methylation_State"), y = "TILs")) +
      geom_violin(fill="black", size = 0.2) +
      geom_boxplot(outlier.shape = NA, width=0.15, fill="grey90", col="grey40", size = 0.2) +
      theme_bw(base_size = 14) +
      labs(x = x_label, y = y_label) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0)) +
      annotate("text", x = 1.5, y = ymax, label = paste0("adj.p = ", signif(pval_adj, 3)))
  }
  
  # Make plots for each gene
  GBP4_plot <- make_tils_boxplot("GBP4", "TILs (%)", NULL, 80, pvals_tils_adj["GBP4"])
  OAS2_plot <- make_tils_boxplot("OAS2", NULL, NULL, 80, pvals_tils_adj["OAS2"])
  ZBP1_plot <- make_tils_boxplot("ZBP1", NULL, NULL, 80, pvals_tils_adj["ZBP1"])
  CARD16_plot <- make_tils_boxplot("CARD16", NULL, NULL, 90, pvals_tils_adj["CARD16"])
  SAMD9L_plot <- make_tils_boxplot("SAMD9L", NULL, NULL, 80, pvals_tils_adj["SAMD9L"])
  
  # Combine horizontally
  GBP4_plot | OAS2_plot | ZBP1_plot | CARD16_plot | SAMD9L_plot

  
  
  # PDL1 CPS BARPLOT
  make_pdl1CPS_barplot <- function(gene, y_label, x_label, legend_pos, show_y_axis = TRUE) {
    gene_meth_col <- paste0(gene, "_Methylation_State")
    
    base_plot <- ggplot(subset(plot_data, !is.na(PDL1_CPS)), 
                        aes_string(x = gene_meth_col, fill = "as.factor(PDL1_CPS)")) +
      geom_bar(position = "fill", alpha = 0.8) +  # Proportional stacked bar
      theme_bw(base_size = 14) +
      labs(x = x_label, 
           y = y_label, 
           fill = "PD-L1 CPS") +
      scale_fill_manual(values = c("0" = "#619CFF", 
                                   "1" = "#9B59B6", 
                                   "2" = "#F8766D"),
                        labels = c("0" = "< 1%", 
                                   "1" = "1% & < 10%", 
                                   "2" = "≥ 10%")) +
      theme(legend.position = legend_pos, 
            legend.text = element_text(size = 12))
    
    if (!show_y_axis) {
      base_plot <- base_plot + 
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    return(base_plot)
  }
  
  GBP4_plot <- make_pdl1CPS_barplot("GBP4", "Proportion of samples", NULL, "none", show_y_axis = TRUE)
  OAS2_plot <- make_pdl1CPS_barplot("OAS2", NULL, NULL, "none", show_y_axis = FALSE)
  ZBP1_plot <- make_pdl1CPS_barplot("ZBP1", NULL, NULL, "none",  show_y_axis = FALSE)
  CARD16_plot <- make_pdl1CPS_barplot("CARD16", NULL, NULL, "none", show_y_axis = TRUE)
  SAMD9L_plot <- make_pdl1CPS_barplot("SAMD9L", NULL, NULL, "right", show_y_axis = FALSE)

  # Combine plots horizontally
  GBP4_plot | OAS2_plot | ZBP1_plot 
  CARD16_plot | SAMD9L_plot

#
# CORRELATION PLOTS OF GENE EXPRESSION
#  

# Getting correlation of lg1p transformed gex data
corr_matrix <- cor(log1p(t(fpkm_data[c("GBP4", "ZBP1", "OAS2", "CARD16", "SAMD9L"),])), method = "spearman", use = "pairwise.complete.obs")

# Plot
corrplot(corr_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", number.cex = 0.7)

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
  
  cluster_promoter <- kmeans(t(betaAdj[names(cpgs),]), centers = 2)
  
  # Determine hypo and hypermethylated cluster
  promoter_state <- if (mean(betaAdj[names(cpgs),cluster_promoter$cluster==1]) >
                             mean(betaAdj[names(cpgs),cluster_promoter$cluster==2])) {
    
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
    return(NA)
  }
})

# Getting number of hypermethylated samples
agreement_hyper <- colSums(numeric_matrix == 1)


# Create the heatmap 
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    Epitype = epi_annotation,
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)

agreement_hyper_annotation <- HeatmapAnnotation("Consenus" = anno_barplot(unname(agreement_hyper)))

# Function to convert pvals to stars
pval_to_stars <- function(pval) {
  if (pval < 0.0001) return("****")
  else if (pval < 0.001) return("***")
  else if (pval < 0.01) return("**")
  else if (pval < 0.05) return("*")
  else return("ns") 
}

# Compute raw p-values
pvals <- sapply(gene_ids, function(gene) {
  chisq.test(clusters_methylation[gene, ], pam50_annotations == "Basal")$p.value
})

# Apply BH correction
pvals_adj <- p.adjust(pvals, method = "BH")

# Convert adjusted p-values to stars
pval_stars <- sapply(pvals_adj, pval_to_stars)

chi_p_annotation <- rowAnnotation(
  pvalue = anno_text(pval_stars, gp = gpar(fontsize = 10))
)

Heatmap(clusters_methylation,
        col = c("red", "blue"),
        top_annotation = top_annotation,
        bottom_annotation = agreement_hyper_annotation,
        left_annotation = chi_p_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = as.factor(ifelse(pam50_annotations == "Basal", "Basal", "NonBasal")))



tnbc_annotation_factor <- factor(tnbc_annotation,
                          levels = c("BL1", "BL2", "M", "LAR"))

keep <- !is.na(tnbc_annotation_factor)

# Compute raw p-values
pvals <- sapply(gene_ids, function(gene) {
  chisq.test(clusters_methylation[gene, keep], tnbc_annotation[keep])$p.value
})

# Apply BH correction
pvals_adj <- p.adjust(pvals, method = "BH")

# Convert adjusted p-values to stars
pval_stars <- sapply(pvals_adj, pval_to_stars)

chi_p_annotation <- rowAnnotation(
  pvalue = anno_text(pval_stars, gp = gpar(fontsize = 10))
)

Heatmap(clusters_methylation[, keep],
        col = c("red", "blue"),
        top_annotation = top_annotation[keep],
        bottom_annotation = agreement_hyper_annotation[keep],
        left_annotation = chi_p_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = tnbc_annotation_factor[keep])

#
# PLOTTNG TMARQ CELL COUNTS VS METHYLATION STATUS
#

# Getting TMArQ cell counts from pdids of interest
cols_of_interest <- c(
  "PDid",
  "CD3_mean",
  "CD4_mean",
  "CD8_mean",
  "CD20_mean",
  "CD68_mean",
  "FOXP3_mean",
  "PDL1_SP142_mean"
)

# Getting mean counts (two cores) from PDIDs
filtered_tmarq <- unique(tmarq_counts[tmarq_counts$PDid %in% colnames(clusters_methylation), cols_of_interest])

# Using PDiDs as rownames
rownames(filtered_tmarq) <- filtered_tmarq$PDid
filtered_tmarq$PDid <- NULL

# Convert tp numeric and sqrt
filtered_tmarq[] <- lapply(filtered_tmarq, function(x) sqrt(as.numeric(as.character(x))))


# Genes of interest
genes_of_interest <- my_genes

# Create tidy data frame
plot_df <- lapply(genes_of_interest, function(gene) {
  data.frame(
    sample = colnames(clusters_methylation),
    gene   = gene,
    methylation_state = unlist(clusters_methylation[gene, colnames(clusters_methylation)]),
    filtered_tmarq[colnames(clusters_methylation), ]
  )
}) %>% bind_rows()

# Reshape counts into long format
plot_df_long <- plot_df %>%
  pivot_longer(
    cols = -c(sample, gene, methylation_state),
    names_to = "marker",
    values_to = "value"
  )

# Formatting marker name
plot_df_long$marker <- sub("_mean$", "", plot_df_long$marker)

# Compute BH-corrected p-values per gene
pvals_df <- plot_df_long %>%
  group_by(gene, marker) %>%
  summarise(
    wilcox_p = wilcox.test(value ~ methylation_state)$p.value,
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    p_adj = p.adjust(wilcox_p, method = "BH"),
    p_star = case_when(
      p_adj <= 0.0001 ~ "****",
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Merge adjusted p-values back to the plot data
plot_df_long <- plot_df_long %>%
  left_join(pvals_df, by = c("gene", "marker"))

# Now plot with manual annotation
ggplot(plot_df_long, aes(x = methylation_state, y = value)) +
  geom_violin(fill = "black") +
  geom_boxplot(width = 0.2, outlier.colour = "grey", size = 0.25) +
  facet_grid(gene ~ marker, scales = "fixed") +
  coord_cartesian(ylim = c(0, 100)) +
  geom_text(
    data = distinct(plot_df_long, gene, marker, p_star),
    aes(x = 1.5, y = 95, label = p_star),
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 9)
  ) +
  labs(x = "Methylation state", y = "Sqrt TMArQ mean counts")


#
# PLOTTNG METHYLCIBERSORT IMMUNE FRACTIONS VS METHYLATION STATUS
#

# Get data of interest from methylcibersort output
tum_cibersort_to_plot <- tum_cibersort[colnames(clusters_methylation), c("CD14","CD19","CD4_Eff","CD56","CD8","Eos","Neu","Treg")]

# Create tidy data frame
plot_df <- lapply(genes_of_interest, function(gene) {
  data.frame(
    sample = colnames(clusters_methylation),
    gene   = gene,
    methylation_state = unlist(clusters_methylation[gene, colnames(clusters_methylation)]),
    tum_cibersort_to_plot
  )
}) %>% bind_rows()


# Reshape counts into long format
plot_df_long <- plot_df %>%
  pivot_longer(
    cols = -c(sample, gene, methylation_state),
    names_to = "cell",
    values_to = "value"
  )

# Compute BH-corrected p-values per gene across cells
pvals_df <- plot_df_long %>%
  group_by(gene, cell) %>%
  summarise(
    wilcox_p = wilcox.test(value ~ methylation_state)$p.value,
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    p_adj = p.adjust(wilcox_p, method = "BH"),
    p_star = case_when(
      p_adj <= 0.0001 ~ "****",
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

# Merge adjusted p-values back to the plot data
plot_df_long <- plot_df_long %>%
  left_join(pvals_df, by = c("gene", "cell"))

# Plot with BH-adjusted significance stars
ggplot(plot_df_long, aes(x = methylation_state, y = value)) +
  geom_violin(fill = "black") +
  geom_boxplot(width = 0.2, outlier.colour = "grey", size = 0.25) +
  facet_grid(gene ~ cell, scales = "fixed") +
  coord_cartesian(ylim = c(0, 0.5)) +
  geom_text(
    data = distinct(plot_df_long, gene, cell, p_star),
    aes(x = 1.5, y = 0.48, label = p_star),
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 9)
  ) +
  labs(x = "Methylation state", y = "Methylcibersort cell fraction")

#
# PLOTTNG FREDLUND STROMA, IMMUNE AND TILS VS HYPERMETHYLATION AGREEMENT
#


# Fredlund stroma
p1 <- ggplot(data.frame(Hyper = as.factor(agreement_hyper),
                  Value = x[names(agreement_hyper), "Fredlund.rank.signatures"][, "Stroma"]),
       aes(x = Hyper, y = Value)) +
  geom_boxplot(fill = "goldenrod", alpha = 0.7) +
  labs(x = "Hypermethylated genes", y = "Fredlund stroma") +
  theme_bw(base_size = 14) 

# MethylCIBERSORT Fibroblast fraction
p2 <- ggplot(data.frame(Hyper = as.factor(agreement_hyper),
                  Value = tum_cibersort[names(agreement_hyper), "Fibroblast"]),
       aes(x = Hyper, y = Value)) +
  geom_boxplot(fill = "goldenrod", alpha = 0.7) +
  labs(x = "Hypermethylated genes", y = "MethylCIBERSORT\nFibroblast fraction") +
  theme_bw(base_size = 14)

# Combine and add shared x label
(p1 / p2) +
  plot_layout(guides = "collect")


# Fredlund immune response
p1 <- ggplot(data.frame(Hyper = as.factor(agreement_hyper),
                  Value = x[names(agreement_hyper), "Fredlund.rank.signatures"][, "ImmuneResponse"]),
       aes(x = Hyper, y = Value)) +
  geom_boxplot(fill = "goldenrod", alpha = 0.7) +
  labs(x = "Hypermethylated genes", y = "Fredlund\nimmune response") +
  theme_bw(base_size = 14)

# TILs
p2 <- ggplot(data.frame(Hyper = as.factor(agreement_hyper),
                  Value = x[names(agreement_hyper), "TILs"]),
       aes(x = Hyper, y = Value)) +
  geom_boxplot(fill = "goldenrod", alpha = 0.7) +
  labs(x = "Hypermethylated genes", y = "TILs") +
  theme_bw(base_size = 14)

# Combine and add shared x label
(p1 / p2) +
  plot_layout(guides = "collect")

