#! usr/bin/Rscript

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(biomaRt)

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
genes_dict <- setNames(genes_of_interest_ensembl$ensembl_gene_id, 
         genes_of_interest_ensembl$hgnc_symbol)


#
# ANALYSING SINGLE CELL DATA
#

# EPITHELIAL CELLS

DimPlot(hbca_epithelial,
        group.by = "cell_type")

# Getting Ensembl IDS from gene names
DotPlot(hbca_epithelial, 
        features = unname(genes_dict), 
        assay = "RNA",
        group.by = "cell_type")


# IMMUNE CELLS

DimPlot(hbca_immune,
        group.by = "cell_type")

# Getting Ensembl IDS from gene names
DotPlot(hbca_immune, 
        features = unname(genes_dict), 
        assay = "RNA",
        group.by = "cell_type")


