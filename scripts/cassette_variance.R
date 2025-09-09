#! usr/bin/Rscript

library(ggplot2)
library(patchwork)
library(dplyr)


#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Cassettes labels
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")$colors
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")$colors
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")$colors

#
# ANALYSING ALL CPGS REGARDLESS OF CONTEXT
#

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Merging labels
cpg_labels <- c(sapply(proximal_10, function(x) {paste0("proximal_", x)}),
                sapply(distal_10, function(x) {paste0("distal_", x)}),
                sapply(promoter_10, function(x) {paste0("promoter_", x)}))

# Initialize a dataframe to store results
cassette_stats <- data.frame(
  cassette = character(),
  mean_var = numeric(),
  sd_var = numeric(),
  stringsAsFactors = FALSE
)

for (cassette in unique(cpg_labels)) {
  # Subset rows belonging to the current cassette
  cassette_rows <- names(which(cpg_labels == cassette))
  
  if (length(cassette_rows) > 0) {
    # Compute variance of betas for the cassette
    variance_of_cassette <- apply(betaAdj[cassette_rows, , drop = FALSE], MARGIN = 1, FUN = var)
    
    # Compute mean and standard deviation of variance
    mean_var <- mean(variance_of_cassette, na.rm = TRUE)
    sd_var <- sd(variance_of_cassette, na.rm = TRUE)
    
    # Append to the dataframe
    cassette_stats <- rbind(cassette_stats, data.frame(
      cassette = cassette,
      mean_var = mean_var,
      sd_var = sd_var
    ))
  }
}

# Getting most variant CpGs
top_1000 <- names(sort(variance_of_betas, decreasing = T)[1:1000])
top_5000 <- names(sort(variance_of_betas, decreasing = T)[1:5000])
top_10000 <- names(sort(variance_of_betas, decreasing = T)[1:10000])
top_20000 <- names(sort(variance_of_betas, decreasing = T)[1:20000])
top_30000 <- names(sort(variance_of_betas, decreasing = T)[1:30000])


# Create a combined data frame for top 1000, 5000, 10000, and 17725
df_combined <- data.frame()

# Loop through different sets of CpGs (top 1000, 5000, etc.)
for (top_n in list(top_1000, top_5000, top_10000, top_20000, top_30000)) {
  # Count frequencies for each level based on the current top CpG set
  level_counts <- table(as.factor(cpg_labels[top_n]))
  
  # Identify levels with frequency less than 0.5 %
  infrequent_levels <- names(level_counts[level_counts < length(top_n) * 0.01])
  others <- sum(level_counts[infrequent_levels])
  
  # Combine infrequent levels into 'Others'
  level_counts <- level_counts[level_counts >= length(top_n) * 0.01]
  level_counts["Others"] <- others
  
  # Add data for current top_n to the combined data frame
  df_temp <- data.frame(
    Cassette = names(level_counts),
    Frequency = as.vector(level_counts),
    Top_CPG_Set = paste("Top", length(top_n))
  )
  
  df_combined <- rbind(df_combined, df_temp)
}

# Data for variance violin plots
top_sets <- list(
  `Top 1000` = top_1000,
  `Top 5000` = top_5000,
  `Top 10000` = top_10000,
  `Top 20000` = top_20000,
  `Top 30000` = top_30000
)

violin_df <- do.call(rbind, lapply(names(top_sets), function(set_name) {
  cpgs <- top_sets[[set_name]]
  data.frame(Top_CPG_Set = set_name, Variance = variance_of_betas[cpgs])
}))


# Define the order for the top CpG sets and the cassettes
top_cpg_order <- c("Top 1000", "Top 5000", "Top 10000", "Top 20000","Top 30000")

cassette_colors <- c(
  # Distal (reds/oranges)
  "distal_0" = "#a63603",
  "distal_1" = "#fee090",
  "distal_2" = "#d73027",
  "distal_3" = "#fdae61",
  "distal_4" = "#e6550d",
  "distal_5" = "#f46d43",
  "distal_6" = "#fc8d59",
  
  # Proximal (blues)
  "proximal_0" = "#084599",
  "proximal_1" = "#c6dbef",
  "proximal_2" = "#9ecae1",
  "proximal_3" = "#6baed6",
  "proximal_4" = "#4292c6",
  "proximal_5" = "#2171b5",
  "proximal_6" = "#084591",
  
  # Promoter (greens)
  "promoter_1" = "#005a39",
  "promoter_1" = "#a1d99b",
  "promoter_2" = "#74c476",
  "promoter_3" = "#41ab5d",
  "promoter_4" = "#238b45",
  "promoter_5" = "#005a31",
  
  # Other
  "Others" = "grey"
)


# First plot: Stacked bar (relative proportions)
bar_plot <- ggplot(df_combined, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), 
                                    y = Frequency, fill = factor(Cassette))) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = cassette_colors) +
  labs(x = NULL, y = "Frequency", fill = "Cassettes") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Second plot: Violin + boxplot (variance)
violin_plot <- ggplot(violin_df, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), y = Variance)) +
  geom_violin(fill = "#999999", color = "black") +
  geom_boxplot(width = 0.1) +
  labs(x = NULL, y = "CpG Variance") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


# Combine all three plots
combined_1 <- violin_plot / bar_plot + plot_layout(
  heights = unit(c(1,11),c("null","cm")))
  

# Print
combined_1


