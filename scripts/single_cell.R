#! usr/bin/Rscript

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

#
# Loading data
#

only_tum <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCTum.rds")
only_tum <- UpdateSeuratObject(only_tum)

only_str <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCSub.rds")
only_str <- UpdateSeuratObject(only_str)

only_tc <- readRDS("/Volumes/Data/Project_3/single_cell_brca/gse161529/GSE161529/seurat_objects/SeuratObject_TNBCTC.rds")
only_tc <- UpdateSeuratObject(only_tc)


# GENES OF INTEREST IN CANCER CELLS

# Plotting expression of genes of interest in tumour cells
DotPlot(only_tum, features = c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L", "IL18R1", "BATF2", "CD69"), group.by = "group", assay = "RNA") + RotatedAxis()


# STROMAL AND IMMUNE CELLS PER SAMPLE

cell_types <- c("0"="T cells", "1"="Macrophagues", "2"="Plasma cells", "3"="Fibroblasts", "4"="T cells","5"="B cells", "6"="Dendritic cells", "7"="Endothelial cells", "8"="Pericytes", "9"="Myeloid cells")

# Get log count data
cell_type_counts <- table(only_str$group, unname(cell_types[only_str$seurat_clusters]))
#cell_type_counts_log <- log(cell_type_counts + 1)

# Z Scale per cell type 
scaled_per_cell_type <- t(scale(cell_type_counts_log))

# Convert to long format for plotting
cell_type_counts_long <- melt(scaled_per_cell_type)

# Add original count values to the long format
cell_type_counts_long$original_counts <- mapply(function(x, y) cell_type_counts[x, y],
                                                cell_type_counts_long$Var2, 
                                                cell_type_counts_long$Var1)

# Plotting with original counts for tile annotation
ggplot(cell_type_counts_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Midpoint set to 0
  labs(x = "Group", y = "Cell Type", fill = "Scaled cell counts") +
  theme_minimal() +
  geom_text(aes(label = original_counts),  color = "black", size = 3)  +  # Use raw counts for annotation
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis

# T CELL SUBTYPES


cell_types <- c("0"="Effector T cells", "1"="Naive/resting T cells", "2"="Regulatory T cells", "3"="Cycling T cells", "4"="TR cells (memory)","5"="NK cells", "6"="Plasma")


cell_type_counts <- t(prop.table(table(only_tc$group, cell_types[only_tc$seurat_clusters]), margin=1))

# Convert to long format for plotting
cell_type_counts_long <- melt(cell_type_counts)

# Plotting
ggplot(cell_type_counts_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Midpoint set to 0
  labs(x = "Group", y = "Cell Type", fill = "Scaled cell counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels if needed
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)  # Add numbers inside the tiles


DotPlot(only_tc, features = c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "MKI67"), group.by = "group", assay = "RNA")
