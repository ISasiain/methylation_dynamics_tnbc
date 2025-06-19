#! usr/bin/Rscript

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(biomaRt)
library(pathwork)

#
# READING SINGLE CELL FILES
#

#The data comes from: 10.1038/s41588-024-01688-9 

# Replace with your actual path
hbca_immune <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_immune.rds")
hbca_epithelial <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_epithelial.rds")
hbca_stroma <- readRDS("/Volumes/Data/Project_3/normal_breast_single_cell/hbca_stroma.rds")


#
# GETTING ENSEMBL IDS FOR GENES OF INTEREST
#

# Defining genes of interest
genes_of_interest <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")


# Getting ensembl ids
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_of_interest_ensembl <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",  
  values = genes_of_interest, 
  mart = ensembl  
)

# Create dictionary of gene ids and ensembl ids
genes_dict <- setNames(genes_of_interest_ensembl$hgnc_symbol, 
                       genes_of_interest_ensembl$ensembl_gene_id)


#
# ANALYSING SINGLE CELL DATA
#

# EPITHELIAL CELLS

# Getting Ensembl IDS from gene names
epi_normal <- DotPlot(hbca_epithelial, 
                      features = names(genes_dict), 
                      assay = "RNA",
                      group.by = "cell_type") + 
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(labels = NULL) +
  scale_y_discrete(labels = function(y) {
    sapply(y, function(label) {
      paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    })
  }) +
  scale_color_gradient2(
    limits = c(-2, 3),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  ) +
  scale_size(
    limits = c(0, 40),
    range = c(0, 8)   # Adjust dot size appearance as desired
  )
 
# IMMUNE CELLS

# Getting Ensembl IDS from gene names
immune_normal <- DotPlot(hbca_immune, 
        features = names(genes_dict), 
        assay = "RNA",
        group.by = "cell_type") + 
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = genes_dict) +
  scale_y_discrete(labels = function(y) {
    sapply(y, function(label) {
      paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
    })
  }) +
  scale_color_gradient2(
    limits = c(-2, 3),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  ) +
  scale_size(
    limits = c(0, 40),
    range = c(0, 8)   # Adjust dot size appearance as desired
  )


# Epi:Immune = 2:1 in vertical stacking
epi_normal / immune_normal + plot_layout(heights = c(1, 2.75)) + plot_layout(heights = c(1, 2.75), guides = "collect")
