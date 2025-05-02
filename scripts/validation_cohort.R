#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(gprofiler2)

#
# LOADING DATA
#

load("/Volumes/Data/Project_3/validation_cohort/Combined_annotations_rel4SCANB_deNovoMainNMF_distalAtac5000_supervisedSubNMF.RData")
load("/Volumes/Data/Project_3/validation_cohort/FPKM_rel4TNBC_validationCohort_n136.RData")
load("/Volumes/Data/Project_3/validation_cohort/PurBeta_adjustedTumor_betaMatrix_V1_V2_reduced_717459commonCpGs_TNBCs_n136.RData")

# Defining gene-cpg dictionary
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")

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
# CONVERT ENSEMBL IDS INTO GENE IDS
#

# Extract Ensembl IDs without versioning (i.e., without ".1", ".2", etc.)
ensembl_ids <- sapply(rownames(gex.data), function(id) {strsplit(id, "\\.")[[1]][1]})

# Convert Ensembl IDs to gene names (e.g., ENTREZGENE)
new_names <- gconvert(ensembl_ids, organism="hsapiens", target="ENTREZGENE", filter_na = F)


combined_names <- sapply(sort(unique(new_names$input_number)), function(num) {

  paste(new_names$name[grep(paste0("^", num, "\\."), new_names$target_number)], collapse = "_")
  
})

# Renaming gex 
rownames(gex.data) <- combined_names

#
# ANALYSING GENES OF INTEREST
#

# Defining genes of interest
current_gene_id <- "CARD16"

# Cluster based on methylation
cpgs <- setNames(annoObj$featureClass[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]],
                 annoObj$illuminaID[annoObj$illuminaID %in% names(genes)[genes == current_gene_id]])

cluster_promoter <- kmeans(t(beta.adjusted[names(cpgs)[cpgs=="promoter"],]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(beta.adjusted[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==1]) >
                           mean(beta.adjusted[names(cpgs)[cpgs=="promoter"],cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}

# ANNOTATION


# Generating annotation for heatmap
tnbc_annotation <- annotations[,"TNBCtype4"]
im_annotation <- annotations[,"TNBCtype_IM"]
pam50_annotation <- annotations[,"majorityClass"]
epi_annotation <- annotations[,"NMF_ATAC_finalSubClusters"]


# Generate bottom annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM" = anno_barplot(as.numeric(gex.data[current_gene_id,]))
)

column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotation,
                                       TNBC = tnbc_annotation,
                                       Epitype = epi_annotation,
                                       IM = im_annotation,
                                       col = list(
                                         "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                         "IM"=c("0"="grey", "1"="black"),
                                         "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey", "UNS"="black"),
                                         "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                     "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
)


Heatmap(beta.adjusted[names(cpgs)[names(cpgs) %in% rownames(beta.adjusted)],],
        column_split = promoter_state,
        show_column_names = FALSE,
        show_row_dend =  FALSE,
        bottom_annotation = bottom_annotation, 
        top_annotation = column_annotation)

boxplot(
    as.numeric(gex.data[current_gene_id,]) ~ promoter_state,
    ylab="FPKM"
)


