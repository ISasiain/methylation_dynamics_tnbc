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
density_1000 <- readRDS("PhD/Projects/project_3/analysis/cpg_density/cpg_desnity_500_bp.rds")

# Cassettes
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")
distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
proximal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_15.rds")

#
# DENSITY VS VARIANCE
#

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

density_1000$Beta_variance <- variance_of_betas[density_1000$cpg_id]

# Create 50 groups based on CpG density
density_1000 <- density_1000 %>%
  mutate(density_group = cut(cpg_density, breaks = 50, labels = FALSE))  # Creates 50 groups

# Create 50 groups based on variance
density_1000 <- density_1000 %>%
  mutate(variance_group = cut(Beta_variance, breaks = 50, labels = FALSE))  # Creates 50 groups

# Compute mean and SD of variance per CpG density group
summary_stats <- density_1000 %>%
  group_by(density_group) %>%
  summarise(
    mean_variance = mean(Beta_variance, na.rm = TRUE),
    sd_variance = sd(Beta_variance, na.rm = TRUE),
    mean_density = mean(cpg_density, na.rm = TRUE)  # Get mean density for x-axis
  )

# Plot mean variance ± SD across CpG density groups
ggplot(summary_stats, aes(x = mean_density, y = mean_variance)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = mean_variance - sd_variance, ymax = mean_variance + sd_variance),
              alpha = 0.2, fill = "blue") +
  theme_minimal() +
  labs(x = "Mean CpG Density per Group", y = "Mean Beta Variance",
       title = "Beta Variance Change Across CpG Density Groups")

# Compute mean and SD of variance per CpG density group
summary_stats <- density_1000 %>%
  group_by(variance_group) %>%
  summarise(
    mean_density = mean(cpg_density, na.rm = TRUE),
    sd_density = sd(cpg_density, na.rm = TRUE),
    mean_variance = mean(Beta_variance, na.rm = TRUE)  # Get mean density for x-axis
  )

# Plot mean variance ± SD across CpG density groups
ggplot(summary_stats, aes(x = mean_variance, y = mean_density)) +
  geom_line(color = "blue", size = 1) +
  geom_ribbon(aes(ymin = mean_density - sd_density, ymax = mean_density + sd_density),
              alpha = 0.2, fill = "blue") +
  theme_minimal() +
  labs(x = "Mean Beta Variance Group", y = "Mean CpG Density",
       title = "CpG Density Change Across Beta Variance Groups")


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
density_to_plot <- density_1000[density_1000$cpg_id %in% names(distal_15$colors),]
density_to_plot$"Cassette" <- as.factor(distal_15$colors[density_to_plot$cpg_id])

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

