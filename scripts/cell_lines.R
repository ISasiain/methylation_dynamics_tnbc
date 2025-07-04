#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(corrplot)


#
# Loading data
#

# Only GEX
gex_data_25_lines <- read.csv("PhD/Projects/project_3/data/GSE202770_TNBC_Cell_line_RNAseq_Gene_FPKM.csv", header = T, row.names = 1)

# GEX + Methylation
load("PhD/Projects/project_3/data/TNBC_8cellines_GEX_Beta.RData") #celline.datd.list

promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")


#
# Plotting promoters and expression
#

cpgs_in_cassette <- names(promoter_10$colors)[promoter_10$colors == 10]


gene_of_interest <- "GBP4"
gene_of_interest_ensembl <- celline.data.list$GEXanno[celline.data.list$GEXanno$Gene.Name == gene_of_interest, "Gene.ID"]

# Generating heatmap annotations
types <- as.character(celline.data.list$cellTypes[, "Type"])
cells <- as.character(celline.data.list$cellTypes[, "Cell"])
colors <- as.character(celline.data.list$cellTypes[, "Color"])

# Define annotation with proper color mapping
top_annotation <- HeatmapAnnotation(
  TNBC = types,
  Cell_lines = cells,
  col = list(
    TNBC = c("BL1" = "red", "BL2" = "blue", "LAR" = "green", "M" = "grey", "MSL" = "pink"),
    Cell_lines = setNames(colors, cells)
  ),
  show_legend = TRUE
)


# Create barplot annotation
gex_annotation <- HeatmapAnnotation(
  log1p_FPKM = anno_barplot(log1p(celline.data.list$GEX[gene_of_interest_ensembl,]), gp = gpar(fill = "black")),
  annotation_name_side = "right",
  show_legend = FALSE
)

# Cassette annotation
cpg_in_cassette <- as.character(rownames(celline.data.list$Beta[names(genes)[genes == gene_of_interest], ]) %in% cpgs_in_cassette)

row_annot <- rowAnnotation(
  CpG_in_cassette = cpg_in_cassette,
  col = list(CpG_in_cassette = c("TRUE" = "black", "FALSE" = "white")),
  show_annotation_name = FALSE,  # This hides the label on the plot but keeps it in the legend
  show_legend = TRUE
)


Heatmap(
  celline.data.list$Beta[names(genes)[genes == gene_of_interest], ],
  column_split = types,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  top_annotation = top_annotation,
  cluster_columns = FALSE,
  bottom_annotation = gex_annotation,
  left_annotation = row_annot,
  show_heatmap_legend = TRUE)



#
# Plotting expression in cell lines and correlation
#


# Subset and transform the expression matrix
genes_of_interest <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")
expr_mat <- log1p(gex_data_25_lines[genes_of_interest, ])

# Reformatting names
colnames(expr_mat) <- gsub("_RNA", "", colnames(expr_mat))


# Convert to long format for ggplot
df_long <- melt(as.matrix(expr_mat))
colnames(df_long) <- c("Gene", "CellLine", "log1p_expression")


# Sort CellLine by mean expression
df_long$CellLine <- factor(df_long$CellLine, levels = names(sort(mean_expr, decreasing = TRUE)))

# Plot
ggplot(df_long, aes(x = CellLine, y = Gene, color = log1p_expression, size = log1p_expression)) +
  geom_point() +
  scale_color_gradient(low = "lightgray", high = "darkred") +
  scale_size(range = c(0, 5)) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = NULL, y = NULL,
    color = "log1p_FPKM", size = "log1p(expr)"
  )

# Transpose and compute correlation
cor_matrix <- cor(t(expr_mat), method = "spearman")

# Plot the correlation matrix
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45)
