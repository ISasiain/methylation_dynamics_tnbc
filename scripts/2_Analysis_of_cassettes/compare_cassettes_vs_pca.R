#! usr/bin/Rscript

library(ggplot2)
library(patchwork)
library(energy)
library(dplyr)
library(tidyr)
library(purrr)
library(rlang)

# USER DEFINED FUNCTIONS
get_pc_variance <- function(pca_obj) {
  var_explained <- (pca_obj$sdev^2) / sum(pca_obj$sdev^2)
  
  data.frame(
    PC = paste0("PC", seq_along(var_explained)),
    Variance = var_explained,
    Percent = 100 * var_explained,
    CumulativePercent = 100 * cumsum(var_explained)
  )
}


#
# Loading the data
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

summary_distal_cassette_10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_7.csv")
summary_promoter_cassette_10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_7.csv")
summary_proximal_cassette_10 <- read.csv("/Users/in2245sa/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_7.csv")

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# 0. PREPROCESSING

# CpG context

# Getting distal CpGs
distal_cpgs <- annoObj$illuminaID[which((annoObj$featureClass == "distal" | annoObj$featureClass == "distal body"))]
distal_betas <- betaAdj[rownames(betaAdj) %in% distal_cpgs, ]

# Getting promoter CpGs
promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass == "promoter")]
promoter_betas <- betaAdj[rownames(betaAdj) %in% promoter_cpgs, ]

# Getting promoter CpGs
proximal_cpgs <- annoObj$illuminaID[which((annoObj$featureClass == "proximal up" | annoObj$featureClass == "proximal dn"))]
proximal_betas <- betaAdj[rownames(betaAdj) %in% proximal_cpgs, ]

# 1. VARIANCE BASED FILTERING

# Getting most variables CpGs
variance_dis <- sapply(1:nrow(distal_betas), FUN = function(row) {var(distal_betas[row,])})
variance_prom <- sapply(1:nrow(promoter_betas), FUN = function(row) {var(promoter_betas[row,])})
variance_prox <- sapply(1:nrow(proximal_betas), FUN = function(row) {var(proximal_betas[row,])})

# Filtering data
selected_var <- 0.05

dis_to_analyse <- t(distal_betas[variance_dis > selected_var,])
prom_to_analyse <- t(promoter_betas[variance_prom > selected_var,])
prox_to_analyse <- t(proximal_betas[variance_prox > selected_var,])

# 2. PCs OF CPGS

# Converting to M values

#PROMOTER
eps <- 0.01
prom_to_analyse <- pmin(pmax(prom_to_analyse, eps), 1 - eps)
prom_to_analyse <- log2(prom_to_analyse / (1 - prom_to_analyse))

#PROXIMAL
eps <- 0.01
prox_to_analyse <- pmin(pmax(prox_to_analyse, eps), 1 - eps)
prox_to_analyse <- log2(prox_to_analyse / (1 - prox_to_analyse))

#DISTAL
eps <- 0.01
dis_to_analyse <- pmin(pmax(dis_to_analyse, eps), 1 - eps)
dis_to_analyse <- log2(dis_to_analyse / (1 - dis_to_analyse))

# DISTAL
pc_distal <- prcomp(dis_to_analyse)
variance_per_pc_distal <- get_pc_variance(pc_distal)

# PROMOTER
pc_promoter <- prcomp(prom_to_analyse)
variance_per_pc_promoter <- get_pc_variance(pc_promoter)

# PROXIMAL
pc_proximal <- prcomp(prox_to_analyse)
variance_per_pc_proximal <- get_pc_variance(pc_proximal)

# 3. PLOTTING CASSETTES VS PCS

# DISTAL

# Prepare a tidy data frame
plot_data <- expand.grid(
  PC = 1:4,
  Cassette = 1:4
) %>%
  mutate(
    data = map2(PC, Cassette, ~{
      y_plot <- as.numeric(pc_distal$x[, .x])
      x_plot <- as.numeric(summary_distal_cassette_10[.y + 1, rownames(pc_distal$x)])
      rho <- dcor(x_plot, y_plot)
      data.frame(
        y = scale(y_plot),
        x = scale(x_plot),
        rho = round(rho, 2)
      )
    })
  ) %>%
  unnest(cols = data)

plot_data$PC <- paste0("PC ", plot_data$PC)
plot_data$Cassette <- paste0("PC1 of Cassette ", plot_data$Cassette)

# Convert PC and Cassette to factor for ordering in facets
plot_data$PC <- factor(plot_data$PC)
plot_data$Cassette <- factor(plot_data$Cassette)

# Get data data for labels

label_data <- plot_data %>%
  group_by(PC, Cassette) %>%
  summarize(
    x_label = min(x) + 0.4 * (max(x) - min(x)),  # slightly right of left edge
    y_label = max(y) - 0.9 * (max(y) - min(y)),  # slightly below top edge
    rho = unique(rho)
  )

# Plot
plot_1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(cex = 0.3) +
  facet_grid(PC ~ Cassette, scales = "free") +
  geom_label(data = label_data,
             aes(x = x_label, y = y_label, label = paste0("dCor=", rho)),
             size = 4, fontface = "bold",
             fill = NA, label.size = 0) +   # no box
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# 4. PLOTTING CASSETTES 2 AND 3 VS PC 2

# Generating dataframe
plot_basal_non_basal <- data.frame(
  "PC1_of_Cassette_2" = scale(as.numeric(summary_distal_cassette_10[3, rownames(pc_distal$x)])),
  "PC1_of_Cassette_3" = scale(as.numeric(summary_distal_cassette_10[4, rownames(pc_distal$x)])),
  "PC2" = scale(as.numeric(pc_distal$x[,2])),
  "PAM50" = x[rownames(pc_distal$x), "PAM50_Basal_NCN"])

plot_2 <- ggplot(plot_basal_non_basal, 
       aes(x = PC1_of_Cassette_2, y = PC1_of_Cassette_3)) +
  geom_point(aes(color = PC2, shape = PAM50), size = 2) +  # points
  scale_color_gradient2(
    low = "blue",      # negative values
    mid = "lightgrey",     # zero
    high = "red",      # positive values
    midpoint = 0       # center the gradient at 0
  ) +
  theme_bw() +
  labs(
    x = "PC1 of Cassette 2",
    y = "PC1 of Cassette 3",
    color = "PC2"
  )

# 5. PLOTTING VARIANCE OF PCs
subset_variance_per_pc_distal <- variance_per_pc_distal[variance_per_pc_distal$CumulativePercent < 40,]
subset_variance_per_pc_distal$PC <- factor(subset_variance_per_pc_distal$PC, levels = c("PC1", "PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17"))

plot_3 <- ggplot(data = subset_variance_per_pc_distal, aes(x=PC, y=Percent)) + 
  geom_col() +
  theme_bw() + 
  xlab("Principal Components") +
  ylab("% of variance explained") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



(plot_1 | (plot_3 / plot_2)) +
  plot_layout(
    widths = c(2, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "" 
  )

# PROMOTER

# Prepare a tidy data frame
plot_data <- expand.grid(
  PC = 1:4,
  Cassette = 1:4
) %>%
  mutate(
    data = map2(PC, Cassette, ~{
      y_plot <- as.numeric(pc_promoter$x[, .x])
      x_plot <- as.numeric(summary_promoter_cassette_10[.y + 1, rownames(pc_promoter$x)])
      rho <- dcor(x_plot, y_plot)
      data.frame(
        y = scale(y_plot),
        x = scale(x_plot),
        rho = round(rho, 2)
      )
    })
  ) %>%
  unnest(cols = data)

plot_data$PC <- paste0("PC ", plot_data$PC)
plot_data$Cassette <- paste0("PC1 of Cassette ", plot_data$Cassette)

# Convert PC and Cassette to factor for ordering in facets
plot_data$PC <- factor(plot_data$PC)
plot_data$Cassette <- factor(plot_data$Cassette)

# Get data data for labels

label_data <- plot_data %>%
  group_by(PC, Cassette) %>%
  summarize(
    x_label = min(x) + 0.4 * (max(x) - min(x)),  # slightly right of left edge
    y_label = max(y) - 0.9 * (max(y) - min(y)),  # slightly below top edge
    rho = unique(rho)
  )

# Plot
plot_1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(cex = 0.3) +
  facet_grid(PC ~ Cassette, scales = "free") +
  geom_label(data = label_data,
             aes(x = x_label, y = y_label, label = paste0("dCor=", rho)),
             size = 4, fontface = "bold",
             fill = NA, label.size = 0) +   # no box
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# 4. PLOTTING CASSETTES 2 AND 3 VS PC 2

# Generating dataframe
plot_basal_non_basal <- data.frame(
  "PC1_of_Cassette_2" = scale(as.numeric(summary_promoter_cassette_10[3, rownames(pc_promoter$x)])),
  "PC1_of_Cassette_3" = scale(as.numeric(summary_promoter_cassette_10[4, rownames(pc_promoter$x)])),
  "PC2" = scale(as.numeric(pc_promoter$x[,2])),
  "PAM50" = x[rownames(pc_promoter$x), "PAM50_Basal_NCN"])

plot_2 <- ggplot(plot_basal_non_basal, 
                 aes(x = PC1_of_Cassette_2, y = PC1_of_Cassette_3)) +
  geom_point(aes(color = PC2, shape = PAM50), size = 2) +  # points
  scale_color_gradient2(
    low = "blue",      # negative values
    mid = "lightgrey",     # zero
    high = "red",      # positive values
    midpoint = 0       # center the gradient at 0
  ) +
  theme_bw() +
  labs(
    x = "PC1 of Cassette 2",
    y = "PC1 of Cassette 3",
    color = "PC2"
  )

# 5. PLOTTING VARIANCE OF PCs
subset_variance_per_pc_promoter <- variance_per_pc_promoter[variance_per_pc_promoter$CumulativePercent < 40,]
subset_variance_per_pc_promoter$PC <- factor(subset_variance_per_pc_promoter$PC, levels = c("PC1", "PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17"))

plot_3 <- ggplot(data = subset_variance_per_pc_promoter, aes(x=PC, y=Percent)) + 
  geom_col() +
  theme_bw() + 
  xlab("Principal Components") +
  ylab("% of variance explained") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



(plot_1 | (plot_3 / plot_2)) +
  plot_layout(
    widths = c(2, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "" 
  )

# PROXIMAL

# Prepare a tidy data frame
plot_data <- expand.grid(
  PC = 1:4,
  Cassette = 1:4
) %>%
  mutate(
    data = map2(PC, Cassette, ~{
      y_plot <- as.numeric(pc_proximal$x[, .x])
      x_plot <- as.numeric(summary_proximal_cassette_10[.y + 1, rownames(pc_proximal$x)])
      rho <- dcor(x_plot, y_plot)
      data.frame(
        y = scale(y_plot),
        x = scale(x_plot),
        rho = round(rho, 2)
      )
    })
  ) %>%
  unnest(cols = data)

plot_data$PC <- paste0("PC ", plot_data$PC)
plot_data$Cassette <- paste0("PC1 of Cassette ", plot_data$Cassette)

# Convert PC and Cassette to factor for ordering in facets
plot_data$PC <- factor(plot_data$PC)
plot_data$Cassette <- factor(plot_data$Cassette)

# Get data data for labels

label_data <- plot_data %>%
  group_by(PC, Cassette) %>%
  summarize(
    x_label = min(x) + 0.4 * (max(x) - min(x)),  # slightly right of left edge
    y_label = max(y) - 0.9 * (max(y) - min(y)),  # slightly below top edge
    rho = unique(rho)
  )

# Plot
plot_1 <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(cex = 0.3) +
  facet_grid(PC ~ Cassette, scales = "free") +
  geom_label(data = label_data,
             aes(x = x_label, y = y_label, label = paste0("dCor=", rho)),
             size = 4, fontface = "bold",
             fill = NA, label.size = 0) +   # no box
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

# 4. PLOTTING CASSETTES 2 AND 3 VS PC 2

# Generating dataframe
plot_basal_non_basal <- data.frame(
  "PC1_of_Cassette_2" = scale(as.numeric(summary_proximal_cassette_10[3, rownames(pc_proximal$x)])),
  "PC1_of_Cassette_3" = scale(as.numeric(summary_proximal_cassette_10[4, rownames(pc_proximal$x)])),
  "PC2" = scale(as.numeric(pc_proximal$x[,2])),
  "PAM50" = x[rownames(pc_proximal$x), "PAM50_Basal_NCN"])

plot_2 <- ggplot(plot_basal_non_basal, 
                 aes(x = PC1_of_Cassette_2, y = PC1_of_Cassette_3)) +
  geom_point(aes(color = PC2, shape = PAM50), size = 2) +  # points
  scale_color_gradient2(
    low = "blue",      # negative values
    mid = "lightgrey",     # zero
    high = "red",      # positive values
    midpoint = 0       # center the gradient at 0
  ) +
  theme_bw() +
  labs(
    x = "PC1 of Cassette 2",
    y = "PC1 of Cassette 3",
    color = "PC2"
  )

# 5. PLOTTING VARIANCE OF PCs
subset_variance_per_pc_proximal <- variance_per_pc_proximal[variance_per_pc_proximal$CumulativePercent < 40,]
subset_variance_per_pc_proximal$PC <- factor(subset_variance_per_pc_proximal$PC, levels = c("PC1", "PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17"))

plot_3 <- ggplot(data = subset_variance_per_pc_proximal, aes(x=PC, y=Percent)) + 
  geom_col() +
  theme_bw() + 
  xlab("Principal Components") +
  ylab("% of variance explained") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



(plot_1 | (plot_3 / plot_2)) +
  plot_layout(
    widths = c(2, 1)
  ) +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = "" 
  )
