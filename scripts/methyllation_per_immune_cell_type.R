#! usr/bin/Rscript

library(ComplexHeatmap)
library(circlize) 
library(dplyr)
library(tidyr)
library(ggplot2)

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

# Prepare plot_data as before
selected_genes <- c("GBP4", "OAS2", "ZBP1", "CARD16", "SAMD9L")

plot_data <- lapply(selected_genes, function(gene) {
  cpgs <- names(genes)[genes == gene]
  cpgs <- cpgs[cpgs %in% rownames(beta.matrix)]
  data <- beta.matrix[cpgs, , drop = FALSE]
  df <- as.data.frame(t(data))
  df$Sample <- rownames(df)
  df_long <- pivot_longer(df, -Sample, names_to = "CpG", values_to = "Beta")
  df_long$Gene <- gene
  df_long
}) %>% bind_rows()

# Filter CpGs
plot_data <- plot_data[plot_data$CpG %in% names(promoter_10$colors)[promoter_10$colors == 10],]

# Merge with sample annotations to get CellType
plot_data <- plot_data %>%
  left_join(sample_annotations[, c("GSMid", "Sample_source")],
            by = c("Sample" = "GSMid")) 

# Abbreviations
short_names <- gsub("Whole blood", "WB", plot_data$Sample_source)
short_names <- gsub("CD4\\+ T cells", "Th", short_names)
short_names <- gsub("CD8\\+ T cells", "Tc", short_names)
short_names <- gsub("CD14\\+ Monocytes", "Mono", short_names)
short_names <- gsub("CD19\\+ B cells", "B", short_names)
short_names <- gsub("CD56\\+ NK cells", "NK", short_names)
short_names <- gsub("Neutrophils", "Neu.", short_names)
short_names <- gsub("Eosinophils", "Eos.", short_names)
short_names <- gsub("Granulocytes", "Gran.", short_names)

plot_data$Sample_source <- short_names

# Plot with jittered points and median lines
ggplot(plot_data, aes(x = CpG, y = Beta)) +
  geom_jitter(aes(fill = "grey20"), width = 0.2, size = 0.8, alpha = 1, show.legend = FALSE) +
  facet_grid(Sample_source ~ Gene, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 12) +
  ylim(0, 1) +
  ylab("Beta Value") +
  xlab("CpG Site") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

