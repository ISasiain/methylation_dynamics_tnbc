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

# Loading promoter cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

#
# PLOTTING CPGS AFFECTING 
#

current_gene_id = "BRCA1"

pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
pam50_annotations <- ifelse(pam50_annotations == "Uncl.", 
                            "Uncl.", 
                            ifelse(pam50_annotations == "Basal", 
                                   "Basal", 
                                   "Non-Basal"))
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
pyro_annotation <- x[colnames(betaAdj), "BRCA1_Hypermethylated"]


# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    BRCA1_pyrosequencing = pyro_annotation,
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Non-Basal"="darkblue","Uncl."="grey"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "BRCA1_pyrosequencing"=c("0"="grey", "1"="black")
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
left_annotation <- rowAnnotation("CpG_in_cassette" = names(genes)[genes == current_gene_id] %in% names(promoter_10$colors)[promoter_10$colors == 11],
                                 "Context"= annoObj$CpG_context[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                                 col=list("Context"=c("Distal" = "#f8766d", 
                                                      "Promoter" = "#00ba38", 
                                                      "Proximal" = "#619cff"),
                                          "CpG_in_cassette"=c("TRUE" = "black", 
                                                              "FALSE" = "white")))

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


# Heatmap of genes
Heatmap(
  betaAdj[names(genes)[genes == current_gene_id],],
  cluster_rows = FALSE,
  row_order = order(annoObj$start[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]]),
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = promoter_state,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  right_annotation = right_annotation,
  left_annotation = left_annotation,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  name = "Tumor beta"
)

