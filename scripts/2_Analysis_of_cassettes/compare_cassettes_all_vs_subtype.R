#! usr/bin/Rscript

library(ggplot2)
library(ggsankey)
library(patchwork)
library(dplyr)


#
# LOADING DATA
#

all_distal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_7.rds")
atac_distal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/atac_cassettes_beta_7.rds")

#
# COMPARING CASSETTE ASSIGNMENTS
#


# FROM ATAC-ONLY TO ALL

# Plotting first 10 cassettes
cassettes <- 1:4
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(atac_distal_7$colors)[atac_distal_7$colors == plot]
  
  
  # Prepare the data
  colors <- all_distal_7$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.04 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Replace cassette 0 with uncassified
  colors_grouped <- ifelse(colors_grouped == "0", "Uncl.", colors_grouped)
  
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
  
  # Handle non-numeric values that can't be coerced 
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
    geom_sankey(flow.alpha = 0.7, node.color = "gray40", fill = "lightgray") +
    geom_sankey_label(size = 3.8, color = "black", fill = "white", fontface = "bold") +
    theme_sankey(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 4)
combined_plot


# FROM ALL TO ATAC-ONLY

# Plotting first 10 cassettes
cassettes <- 1:4
plot_list <- list()

# Generating plots
for (plot in cassettes) {
  
  # Getting CpGs
  cpgs_cas <- names(all_distal_7$colors)[all_distal_7$colors == plot]
  
  
  # Prepare the data
  colors <- atac_distal_7$colors[cpgs_cas]
  colors[is.na(colors)] <- "NA"
  
  # Count frequencies
  color_freq <- table(colors)
  threshold <- 0.04 * length(colors)
  
  # Replace infrequent categories with "Other"
  colors_grouped <- ifelse(color_freq[colors] < threshold, "Other", colors)
  
  # Replace cassette 0 with uncassified
  colors_grouped <- ifelse(colors_grouped == "0", "Uncl.", colors_grouped)
  
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
    geom_sankey(flow.alpha = 0.7, node.color = "gray40", fill = "lightgray") +
    geom_sankey_label(size = 3.8, color = "black", fill = "white", fontface = "bold") +
    theme_sankey(base_size = 15) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
}


# Plotting
combined_plot <- wrap_plots(plot_list, nrow = 1, ncol = 4)
combined_plot


