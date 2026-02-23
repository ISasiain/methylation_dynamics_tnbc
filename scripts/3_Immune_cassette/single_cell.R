#! usr/bin/Rscript

library(Seurat)
library(presto)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(patchwork)

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

# Count number of cells per tumor
table(only_tum$group)
table(only_str$group)
table(only_tc$group)


# Count tumor cells per group
cell_counts <- only_tum@meta.data %>%
  group_by(group) %>%
  summarise(n_cells = n())

# Base DotPlot
gene_expression <- DotPlot(
  only_tum,
  features = c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L"),
  group.by = "group",
  assay = "RNA"
) +
  RotatedAxis() +
  xlab("Genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  scale_color_gradient2(
    limits = c(-1, 4),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  )

# STROMAL AND IMMUNE CELLS PER SAMPLE
cell_types <- c("0"="T cells", "1"="Macrophagues", "2"="Plasma cells", "3"="Fibroblasts", "4"="T cells","5"="B cells", "6"="Dendritic cells", "7"="Endothelial cells", "8"="Pericytes", "9"="Myeloid cells")

# Get normalized cell counts
cell_type_counts <- table(only_str$group, unname(cell_types[only_str$seurat_clusters]))
tumor_counts <- as.numeric(table(only_tum$group)[rownames(cell_type_counts)])
cell_type_counts <- cbind(cell_type_counts, Tumor = tumor_counts)

total_cells_per_sample <- rowSums(cell_type_counts)
cell_type_fractions <- sweep(cell_type_counts, 1, total_cells_per_sample, FUN = "/")

# Convert to long format for plotting
cell_type_counts_long <- melt(t(cell_type_fractions))

# Add original count values to the long format
cell_type_counts_long$original_counts <- mapply(function(x, y) cell_type_counts[x, y],
                                                cell_type_counts_long$Var2, 
                                                cell_type_counts_long$Var1)

# Plotting with original counts for tile annotation
composition_of_microenvironment <- ggplot(cell_type_counts_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Midpoint set to 0
  labs(x = "Cell types", y = "Sample", fill = "Proportion of cell type") +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_text(aes(label = original_counts),  color = "black", size = 4)  +  # Use raw counts for annotation
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis


# Combine the plots
(gene_expression | composition_of_microenvironment) + 
  plot_layout(widths = c(1, 1.75), guides = "collect") &
  theme(legend.position = "right")


# Base DotPlot
gene_expression <- DotPlot(
  only_tum,
  features = c("CD"),
  group.by = "group",
  assay = "RNA"
) +
  RotatedAxis() +
  xlab("Genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  scale_color_gradient2(
    limits = c(-1, 4),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  )


# Define the cell type labels
cell_types <- c(
  "0" = "T cells", 
  "1" = "Macrophages", 
  "2" = "Plasma cells", 
  "3" = "Fibroblasts", 
  "4" = "T cells",
  "5" = "B cells", 
  "6" = "Dendritic cells", 
  "7" = "Endothelial cells", 
  "8" = "Pericytes", 
  "9" = "Myeloid cells"
)

# Optionally, map the identities in Seurat object
Idents(only_str) <- factor(Idents(only_str), levels = names(cell_types), labels = cell_types)

# Create the dot plot
DotPlot(
  only_str,
  features = c("CD3E", "ITGAX", "PLEK",
    "CD4", "CD8A", "MS4A1", "FOXP3", "CD68",
    "CD163", "PECAM1", "SPP1", "NTAN1", "VSIR",
    "HAVCR2", "LAG3", "PDCD1", "CD274", "ACTA2",
    "FAP", "PDGFRA", "MKI67", "CD80", "GZMB"
  ),
  assay = "RNA"
) +
  RotatedAxis() +
  xlab("Genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_gradient2(
    limits = c(-1, 4),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  )



# Create the dot plot
DotPlot(
  only_tc,
  features = c("CD3E", "PLEK",
    "CD4", "CD8A", "MS4A1", "FOXP3", "CD68",
    "CD163", "PECAM1", "SPP1", "NTAN1", "VSIR",
    "HAVCR2", "LAG3", "PDCD1", "CD274", "ACTA2",
    "FAP", "PDGFRA", "MKI67", "CD80", "GZMB"
  ),
  assay = "RNA"
) +
  RotatedAxis() +
  xlab("Genes") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_gradient2(
    limits = c(-1, 4),
    low = "gold1",
    mid = "lightgrey",
    high = "purple",
    midpoint = 0
  )
