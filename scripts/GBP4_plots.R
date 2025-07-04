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

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# Reading PDL1 TLS annotation. Matrix is called u.frame
load("PhD/Projects/project_3/data/summarized_TPS_data_sampleLevel.RData")
rownames(u.frame) <- u.frame$PD_ID

# Loading promoter cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

# Cibersort output
tum_cibersort <- read.csv("PhD/Projects/project_3/analysis/cell_type_deconvolution/CIBERSORTx_tumour.csv")
rownames(tum_cibersort) <- tum_cibersort$Mixture

#
# PLOTTING CPGS AFFECTING 
#

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
    #"ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
    col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
               "ATAC" = c("0" = "white", "1"= "black"))
  )
  
  # CpG annotation. Context and included in cassette
  left_annotation <- rowAnnotation("CpG_in_cassette" = names(genes)[genes == current_gene_id] %in% names(proximal_10$colors)[proximal_10$colors == 1],
                                   "Context"= annoObj$CpG_context[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                                   col=list("Context"=c("Distal" = "#f8766d", 
                                                        "Promoter" = "#00ba38", 
                                                        "Proximal" = "#619cff"),
                                            "CpG_in_cassette"=c("TRUE" = "black", 
                                                                "FALSE" = "white")))
  
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
      km <- kmeans(t(betaAdj[names(cpgs)[cpgs == "promoter"], ]), centers = 2)
      if (mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 1]) >
          mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    OAS2_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "OAS2"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "OAS2"]])
      km <- kmeans(t(betaAdj[names(cpgs)[cpgs == "promoter"], ]), centers = 2)
      if (mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 1]) >
          mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    ZBP1_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "ZBP1"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "ZBP1"]])
      km <- kmeans(t(betaAdj[names(cpgs)[cpgs == "promoter"], ]), centers = 2)
      if (mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 1]) >
          mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    CARD16_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "CARD16"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "CARD16"]])
      km <- kmeans(t(betaAdj[names(cpgs)[cpgs == "promoter"], ]), centers = 2)
      if (mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 1]) >
          mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo"))
      }
    }),
    
    SAMD9L_Methylation_State = with(annoObj, {
      cpgs <- setNames(featureClass[illuminaID %in% names(genes)[genes == "SAMD9L"]],
                       illuminaID[illuminaID %in% names(genes)[genes == "SAMD9L"]])
      km <- kmeans(t(betaAdj[names(cpgs)[cpgs == "promoter"], ]), centers = 2)
      if (mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 1]) >
          mean(betaAdj[names(cpgs)[cpgs == "promoter"], km$cluster == 2])) {
        as.factor(ifelse(km$cluster == 1, "Hyper.", "Hypo."))
      } else {
        as.factor(ifelse(km$cluster == 2, "Hyper.", "Hypo."))
      }
    }),
    HRD = as.factor(x[colnames(betaAdj), "HRD.2.status"])
  )
  
  
  # EXPRESSION BOXPLOT
  
  # Function to make plots
  make_expression_boxplot <- function(gene, y_label, x_label, ymax) {
    ggplot(plot_data, aes_string(x = paste0(gene, "_Methylation_State"), y = paste0(gene, "_FPKM"))) +
      geom_violin(fill="black") +
      geom_boxplot(outlier.shape = NA, width=0.15, fill="grey90", col="grey40",) +
      theme_bw(base_size = 14) +
      labs(x = x_label, y = y_label) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0)) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        comparisons = list(c("Hypo.", "Hyper.")),
        label.x = c("Hypo.", "Hyper."),
        label.y = ymax
      )
  }
  
  GBP4_plot <- make_expression_boxplot("GBP4", "FPKM", NULL, 54)
  OAS2_plot <- make_expression_boxplot("OAS2", "FPKM", NULL, 98)
  ZBP1_plot <- make_expression_boxplot("ZBP1", "FPKM", NULL, 19)
  CARD16_plot <- make_expression_boxplot("CARD16", "FPKM", NULL, 37)
  SAMD9L_plot <- make_expression_boxplot("SAMD9L", "FPKM", NULL, 41)
  
  # Combine plots vertically
  (GBP4_plot | OAS2_plot | ZBP1_plot | CARD16_plot | SAMD9L_plot) + plot_layout(guides = "collect")
  
  
  # TILs BOXPLOT
  
  make_tils_boxplot <- function(gene, y_label, x_label, ymax) {
    
    ggplot(plot_data, aes_string(x = paste0(gene, "_Methylation_State"), y = "TILs", fill = paste0(gene, "_Methylation_State"))) +
      geom_violin(fill="black") +
      geom_boxplot(outlier.shape = NA, width=0.15, fill="grey90", col="grey40",) +
      theme_bw(base_size = 14) +  # Classic theme
      labs(x = x_label, y = y_label) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 0)) +
      stat_compare_means(method = "wilcox.test", label = "p.format", 
                       comparisons = list(c("Hypo.", "Hyper.")), 
                       label.x = c("Hypo.", "Hyper."),
                       label.y = ymax)
  }
  
  GBP4_plot <- make_tils_boxplot("GBP4", "TILs (%)", NULL, 80)
  OAS2_plot <- make_tils_boxplot("OAS2", NULL, NULL, 80)
  ZBP1_plot <- make_tils_boxplot("ZBP1", NULL, NULL, 80)
  CARD16_plot <- make_tils_boxplot("CARD16", NULL, NULL, 80)
  SAMD9L_plot <- make_tils_boxplot("SAMD9L", NULL, NULL, 80)
  
  # Combine plots horizontally
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
    return(NA)
  }
})

# Getting number of hypermethylated samples
agreement_hyper <- colSums(numeric_matrix == 1)


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

agreement_hyper_annotation <- HeatmapAnnotation("Consenus" = anno_barplot(unname(agreement_hyper)))

Heatmap(clusters_methylation,
        col = c("red", "blue"),
        top_annotation = top_annotation,
        bottom_annotation = agreement_hyper_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = as.factor(ifelse(pam50_annotations == "Basal", "Basal", "NonBasal")))


Heatmap(clusters_methylation,
        col = c("red", "blue"),
        top_annotation = top_annotation,
        bottom_annotation = agreement_hyper_annotation,
        show_column_names = FALSE,
        show_row_names = TRUE,
        name = "Methylation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_split = as.factor(tnbc_annotation))



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

