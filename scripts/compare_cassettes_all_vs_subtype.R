#! usr/bin/Rscript

library(ggplot2)
library(ggsankey)
library(patchwork)
library(dplyr)


#
# LOADING DATA
#

all_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5.rds")
basal_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_basal_cassettes_beta_5.rds")
nonBasal_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_nonBasal_cassettes_beta_5.rds")

atac_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5_only_atac.rds")

#
# COMPARING CASSETTE ASSIGNMENTS
#


# FROM ALL TO NON BASAL

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(all_distal_5$colors)[all_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- nonBasal_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "All_Samples" = paste0("Cassette_", plot),
    "NonBasal_Samples" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("All_Samples", "NonBasal_Samples")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node, 
                      label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot


# FROM ALL TO BASAL

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(all_distal_5$colors)[all_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- basal_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "All_Samples" = paste0("Cassette_", plot),
    "NonBasal_Samples" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("All_Samples", "NonBasal_Samples")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                                           next_x = next_x, 
                                           node = node, 
                                           next_node = next_node, 
                                           label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot


# FROM NON-BASAL TO ALL

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(nonBasal_distal_5$colors)[nonBasal_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- all_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "NonBasal_Samples" = paste0("Cassette_", plot),
    "All_Samples" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("NonBasal_Samples", "All_Samples")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                                           next_x = next_x, 
                                           node = node, 
                                           next_node = next_node, 
                                           label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot


# FROM BASAL TO ALL

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(basal_distal_5$colors)[basal_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- all_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "NonBasal_Samples" = paste0("Cassette_", plot),
    "All_Samples" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("NonBasal_Samples", "All_Samples")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                                           next_x = next_x, 
                                           node = node, 
                                           next_node = next_node, 
                                           label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot



# FROM ATAC-ONLY TO ALL

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(atac_distal_5$colors)[atac_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- all_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "Atac_CpGs" = paste0("Cassette_", plot),
    "All_CpGs" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("Atac_CpGs", "All_CpGs")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                                           next_x = next_x, 
                                           node = node, 
                                           next_node = next_node, 
                                           label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot


# FROM ALL TO ATAC-ONLY

# Plotting first 10 cassettes
cassettes <- 1:10
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(all_distal_5$colors)[all_distal_5$colors == plot]
  
  
  # Prepare the data
  colors <- atac_distal_5$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.02 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Build dataframe
  df <- data.frame(
    "Atac_CpGs" = paste0("Cassette_", plot),
    "All_CpGs" = unname(colors_grouped)
  )
  
  # Convert to long format
  df_long <- df %>%
    make_long("Atac_CpGs", "All_CpGs")
  
  # Sorting nodes for plotting
  nodes <- unique(df_long$node)
  
  # Separate special values
  special_na <- nodes[nodes == "NA"]
  special_other <- nodes[nodes == "Other"]
  
  # Get everything else
  normal_nodes <- setdiff(nodes, c("NA", "Other"))
  
  # Try to sort numerics correctly (as character but numerically)
  suppressWarnings({
    numeric_order <- order(as.numeric(normal_nodes), na.last = NA)
  })
  
  # Handle non-numeric values that can't be coerced (e.g., "CD8", etc.)
  non_numeric_nodes <- normal_nodes[is.na(suppressWarnings(as.numeric(normal_nodes)))]
  numeric_nodes <- normal_nodes[!is.na(suppressWarnings(as.numeric(normal_nodes)))]
  
  # Sort each group
  sorted_numeric <- numeric_nodes[order(as.numeric(numeric_nodes), decreasing = TRUE)]
  sorted_non_numeric <- sort(non_numeric_nodes)
  
  # Combine everything
  nodes <- c(special_other, sorted_numeric, sorted_non_numeric, special_na)
  
  # Sort
  df_long$node <- factor(df_long$node,levels = nodes)
  df_long$next_node <- factor(df_long$next_node,levels = nodes)
  
  # Plot
  plot_list[[plot]] <- ggplot(df_long, aes(x = x, 
                                           next_x = next_x, 
                                           node = node, 
                                           next_node = next_node, 
                                           label = node)) +
    geom_sankey(flow.alpha = 0.6, node.color = "gray") +
    geom_sankey_label(size = 4, color = "black", fill = "white") +
    theme_sankey(base_size = 16) +
    theme(axis.title.x = element_blank())
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 5)
combined_plot


