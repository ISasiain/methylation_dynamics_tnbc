#! usr/bin/Rscript

library(ggplot2)
library(gridExtra)
library(dplyr)

#
# LOADING DATA
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")


# CpG density data
density_500 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_density_500_bp.rds")
density_1000 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_density_1000_bp.rds")
density_2000 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_density_2000_bp.rds")
density_5000 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_density_5000_bp.rds")
density_10000 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_density_10000_bp.rds")


# Cassettes
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")
distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
proximal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_15.rds")


#
# ANALYSE DENSITY PER CASSETTE
#

# PROMOTER

promoter_cassettes <- sort(unique(promoter_15$colors))

# Get the first 10 cassettes
cassettes_to_plot <- head(promoter_cassettes, 6)

# Density to plot
density_to_plot <- density_1000[density_1000$cpg_id %in% names(promoter_15$colors),]
density_to_plot$"Cassette" <- as.factor(promoter_15$colors[density_to_plot$cpg_id])

# Filter to include cassettes of interest
density_to_plot <- density_to_plot[density_to_plot$Cassette %in% cassettes_to_plot,]

# Generate the violin plot
ggplot(density_to_plot, aes(x = Cassette, y = cpg_density)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 1, aes(fill = Cassette)) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) + # Add boxplot inside violins
  labs(x = "Cassette", y = "CpG Density (CpGs / bp)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()


# Create Cassette-based density plot
cassette_density_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_density(alpha = 0.6, adjust = 2.5) +  # Separate densities with adjust
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none")

# Create Cassette-based histogram plot
cassette_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_histogram(aes(y = ..density..), alpha = 0.6, bins = 30, position = "dodge") +  # Ensure histogram scale is density
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none")

# Create all CpGs density plot (without faceting)
all_cpg_density_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_density(adjust = 2.5, fill = "grey") +  # Separate densities with adjust
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic()

# Create all CpGs histogram plot (without faceting)
all_cpg_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_histogram(aes(y = ..density..), alpha = 0.6, bins = 30) +  # Ensure histogram scale is density
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic()

# Arrange the plots side by side
grid.arrange(
  cassette_density_plot, 
  all_cpg_density_plot, 
  ncol = 2,
  widths = c(2, 1)
)

# Arrange the plots side by side
grid.arrange(
  cassette_histogram_plot, 
  all_cpg_histogram_plot, 
  ncol = 2,
  widths = c(2, 1)
)



# DISTAL

distal_cassettes <- sort(unique(distal_15$colors))

# Get the first 10 cassettes
cassettes_to_plot <- head(distal_cassettes, 6)

# Density to plot
density_to_plot <- density_1000[density_10000$cpg_id %in% names(distal_15$colors),]
density_to_plot$"Cassette" <- as.factor(distal_15$colors[density_to_plot$cpg_id])

# Filter to include cassettes of interest
density_to_plot <- density_to_plot[density_to_plot$Cassette %in% cassettes_to_plot,]

# Adding annotations about repetitive sequences
repeats <- as.factor(annoObj$hasAnyRepeatOverlap)
names(repeats) <- annoObj$illuminaID

density_to_plot$hasAnyRepeatOverlap <- repeats[density_to_plot$cpg_id]
levels(density_to_plot$hasAnyRepeatOverlap) <- c("False", "True")

# Create bar plot of repetitive sequences in cassette
ggplot(density_to_plot, aes(x = Cassette, fill = as.factor(hasAnyRepeatOverlap))) +
  geom_bar(position = "fill") +
  labs(x = "CpG Cassette", y = "Count", fill = "Has Repeat Overlap") +
  theme_minimal()

# Generate the violin plot
ggplot(density_to_plot, aes(x = Cassette, y = cpg_density, fill = Cassette)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) + # Add boxplot inside violins
  labs(x = "Cassette", y = "CpG Density (CpGs / bp)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()

# Create Cassette-based density plot
cassette_density_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_density(alpha = 0.6, adjust = 2.5) +  # Separate densities with adjust
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none")

# Create Cassette-based histogram plot
cassette_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "dodge") +  # Ensure histogram scale is density
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  xlim(0,0.02) +
  theme(legend.position = "none")

# Create all CpGs density plot (without faceting)
all_cpg_density_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_density(adjust = 2.5, fill = "grey") +  # Separate densities with adjust
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic()

# Create all CpGs histogram plot (without faceting)
all_cpg_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_histogram(alpha = 0.6, bins = 30) +  # Ensure histogram scale is density
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  xlim(0,0.02) +
  theme_classic()

# Arrange the plots side by side
grid.arrange(
  cassette_density_plot, 
  all_cpg_density_plot, 
  ncol = 2,
  widths = c(2, 1)
)

# Arrange the plots side by side
grid.arrange(
  cassette_histogram_plot, 
  all_cpg_histogram_plot, 
  ncol = 2,
  widths = c(2, 1)
)

# PROXIMAL

# Get the first 10 cassettes
cassettes_to_plot <- head(proximal_cassettes, 6)

# Density to plot
density_to_plot <- density_1000[density_1000$cpg_id %in% names(proximal_15$colors),]
density_to_plot$"Cassette" <- as.factor(proximal_15$colors[density_to_plot$cpg_id])

# Filter to include cassettes of interest
density_to_plot <- density_to_plot[density_to_plot$Cassette %in% cassettes_to_plot,]

# Generate the violin plot
ggplot(density_to_plot, aes(x = Cassette, y = cpg_density, fill = Cassette)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) + # Add boxplot inside violins
  labs(x = "Cassette", y = "CpG Density (CpGs / bp)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic()

# Create Cassette-based density plot
cassette_density_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_density(alpha = 0.6, adjust = 2.5) +  # Separate densities with adjust
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none")

# Create Cassette-based histogram plot
cassette_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density, fill = Cassette, color = Cassette)) +
  geom_histogram(aes(y = ..density..), alpha = 0.6, bins = 30, position = "dodge") +  # Ensure histogram scale is density
  facet_wrap(~Cassette, scales = "free_y") +  # Facet by Cassette
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic() +
  theme(legend.position = "none")

# Create all CpGs density plot (without faceting)
all_cpg_density_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_density(adjust = 2.5, fill = "grey") +  # Separate densities with adjust
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic()

# Create all CpGs histogram plot (without faceting)
all_cpg_histogram_plot <- ggplot(density_to_plot, aes(x = cpg_density)) +
  geom_histogram(aes(y = ..density..), alpha = 0.6, bins = 30) +  # Ensure histogram scale is density
  labs(x = "CpG Density (CpGs / bp)", y = "Density") +
  theme_classic()

# Arrange the plots side by side
grid.arrange(
  cassette_density_plot, 
  all_cpg_density_plot, 
  ncol = 2,
  widths = c(2, 1)
)

# Arrange the plots side by side
grid.arrange(
  cassette_histogram_plot, 
  all_cpg_histogram_plot, 
  ncol = 2,
  widths = c(2, 1)
)