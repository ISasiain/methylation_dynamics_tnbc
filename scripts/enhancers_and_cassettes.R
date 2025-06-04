#! usr/bin/Rscript

library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(msigdb)

#Loadinggprofiler2#Loading data

# Load data
dist_cassette_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

x <- x[colnames(betaAdj), ]

# Loading enhancer data
enhancer_general_annotation <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_elements_v5.24.txt", header = T)
enhancer_general_annotation$chr <- paste0("chr", enhancer_general_annotation$chr)

enhancer_to_genes <- read.table("/Volumes/Data/GeneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.24.txt", header = T)

# Loading cpg density data
cpg_density_5000 <- readRDS("/Users/isasiain/PhD/Projects/project_3/analysis/cpg_density/cpg_density_5000_bp.rds")

#
# EXPLORING DISTAL CASSETTES
#

# DISTAL 1
cpgs_1 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 1]

matches_to_enhancers_1 <- sapply(cpgs_1, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_1)


# DISTAL 2
cpgs_2 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 2]

matches_to_enhancers_2 <- sapply(cpgs_2, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2)



# DISTAL 3
cpgs_3 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 3]

matches_to_enhancers_3 <- sapply(cpgs_3, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3)


# DISTAL 4
cpgs_4 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 4]

matches_to_enhancers_4 <- sapply(cpgs_4, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4)


# DISTAL 5
cpgs_5 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 5]

matches_to_enhancers_5 <- sapply(cpgs_5, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5)

# DISTAL 6
cpgs_6 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 6]

matches_to_enhancers_6 <- sapply(cpgs_6, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6)


# DISTAL 7
cpgs_7 <- names(dist_cassette_10$colors)[dist_cassette_10$colors == 7]

matches_to_enhancers_7 <- sapply(cpgs_7, function(cpg) {
  
  chr <- annoObj[cpg, "chr"]
  start <- annoObj[cpg, "start"]
  end <- annoObj[cpg, "end"]
  
  # Get chr of interest
  filtered_enhancer_data <- enhancer_general_annotation[enhancer_general_annotation$chr == chr,]
  
  # Find overlapping enhancers (if any)
  is_inside <- start >= filtered_enhancer_data$element_start & end <= filtered_enhancer_data$element_end
  
  # Return enhancer/s
  paste0(filtered_enhancer_data$GHid[is_inside], collapse = "_")
})

sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)


# Calculate mapping percentages
percentages <- c(
  sum(matches_to_enhancers_1 != "") / length(matches_to_enhancers_1),
  sum(matches_to_enhancers_2 != "") / length(matches_to_enhancers_2),
  sum(matches_to_enhancers_3 != "") / length(matches_to_enhancers_3),
  sum(matches_to_enhancers_4 != "") / length(matches_to_enhancers_4),
  sum(matches_to_enhancers_5 != "") / length(matches_to_enhancers_5),
  sum(matches_to_enhancers_6 != "") / length(matches_to_enhancers_6),
  sum(matches_to_enhancers_7 != "") / length(matches_to_enhancers_7)
) * 100  

# Create a data frame for plotting
df_plot <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages
)

# Plot
p_enhancers <- ggplot(df_plot, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping\nwith enhancers (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14) +  
  theme(axis.text = element_text(size = 12), 
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())


# 
#  REPETITIVE SEQUENCES PER CASSETTE
#



# Calculate mapping percentages of repeatsöä
percentages_repeats <- c(
  sum(annoObj[cpgs_1, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_1, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_2, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_2, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_3, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_3, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_4, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_4, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_5, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_5, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_6, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_6, "hasAnyRepeatOverlap"]),
  sum(annoObj[cpgs_7, "hasAnyRepeatOverlap"]) / length(annoObj[cpgs_7, "hasAnyRepeatOverlap"])
) * 100  

# Create a data frame for plotting
df_plot_repeats <- data.frame(
  Cassette = paste0(1:7),
  PercentMapped = percentages_repeats
)

# Plot
p_repetitive <- ggplot(df_plot_repeats, aes(x = Cassette, y = PercentMapped)) +
  geom_bar(stat = "identity", fill = "black") +
  ylab("CpGs overlapping with \nrepetitive sequences (%)") +
  xlab("Distal Cassettes") +
  theme_bw(base_size = 14)

p_enhancers / p_repetitive


#
# CpG DENSITY PER CASSETTE
#

# Get CpGs to analyse
cpgs_for_density <- names(dist_cassette_10$colors)[dist_cassette_10$colors %in% c(7,6,5,4,3,2,1)]

# Filter density values of interest
rownames(cpg_density_5000) <- cpg_density_5000$cpg_id
density_of_interest <-cpg_density_5000[cpgs_for_density,]

# Adding column for cassette
density_of_interest$cassette <- dist_cassette_10$colors[density_of_interest$cpg_id]
density_of_interest$cassette <- factor(density_of_interest$cassette, levels=c(7,6,5,4,3,2,1))

ggplot(density_of_interest, aes(x = cassette, y = cpg_density)) +
  geom_boxplot(outlier.shape = NA, fill = "#67A710", color = "black") +
  geom_jitter(width = 0.1, size = 0.001, alpha = 0.5, color = "black") +
  theme_bw(base_size = 14) +
  ylim(0,0.005) +
  labs(x = "Cassette", y = "CpG Density\n(CpGs / bp)") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) + coord_flip()


#
# MAPPING ENHANCERS TO GENES. GO ENRICHMENT
#

my_enhancers <- enhancer_to_genes[enhancer_to_genes$GHid %in% matches_to_enhancers_2,]
gene_list_cs1 <- unique(my_enhancers[my_enhancers$is_elite == 1, "symbol"])

# Convert first to ensembl ID (therev are 3 times more genes not being mapped if done directly)
ensembl_ids <- gconvert(gene_list_cs1)$target

# Convert to ENTREZ-ID
gene_df_cs1 <- bitr(ensembl_ids,
     fromType = "ENSEMBL",
     toType = "ENTREZID",
     OrgDb = org.Hs.eg.db)

output_cas1 <- enrichKEGG(gene_df_cs1$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "fdr")


output_cas1 <- enrichGO(gene_df_cs1$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        pAdjustMethod = "fdr",
                        readable = T
                          )


# 
hallmark <- msigdb(species = "Homo sapiens", category = "H")

# Enrichment using ENTREZ IDs and Hallmark gene sets
enrich_result <- enricher(gene = gene_df_cs1$ENTREZID,
                          TERM2GENE = hallmark)

# View results
head(enrich_result)


summary(output_cas1)

