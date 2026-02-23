#! usr/bin/Rscript

#
# LOAD DATA
#

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

#
# ANALYSING ATAC PER CONTEXT
#

# Analysing atac per context
table(annoObj$hasAtacOverlap, annoObj$CpG_context)

# All ATAC +
sum(table(annoObj$hasAtacOverlap, annoObj$CpG_context)[2,])

# All ATAC -
sum(table(annoObj$hasAtacOverlap, annoObj$CpG_context)[1,])
