#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize) 

#
# LOADING DATA
#

sample_annotations <- read.table("/Volumes/Data/Project_3/expression_in_immune_cells/OneDrive_1_2025-04-03/GSE35069_Annotations_Step1.txt", 
                                 sep="\t", 
                                 header=TRUE, 
                                 stringsAsFactors=FALSE)

# The beta matrix is called beta.matrix
load("/Volumes/Data/Project_3/expression_in_immune_cells/OneDrive_1_2025-04-03/GSE35069_Beta.RData")


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
# Methylation state per cel type
#

# Reordering matrix
beta.matrix <- beta.matrix[,sample_annotations$GSMid]

# Getting cpgs of interest
gbp4_cpgs <- names(genes)[genes == "SAMD9L"]
gbp4_cpgs <- gbp4_cpgs[gbp4_cpgs %in% rownames(beta.matrix)]
data_subset <- beta.matrix[gbp4_cpgs,]

# Plotting heatmap


Heatmap(data_subset,
        column_split = as.factor(sample_annotations$Sample_source),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_column_names = FALSE,
        column_title_rot = 90,
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
)


