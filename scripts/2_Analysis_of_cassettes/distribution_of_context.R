#! usr/bin/Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

#
# LOAD DATA
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Cassettes labels
all_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/all/cassettes_beta_7.rds")$colors

# Grouping CpG contexts
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

#
# CONTEXT PER CASSETTES (using all cpgs)
#

# Defining N most variable CpGs to check
N_most_variable <- c(1000, 2500, 5000, 10000, 20000, 50000)

# Getting beta varaince
var_adj <- apply(betaAdj, MARGIN = 1, FUN=var)
var_unadj <- apply(betaNew, MARGIN = 1, FUN=var)

# Plotting adjusted

# Create a results list for different N thresholds
res_list <- lapply(N_most_variable, function(N) {
  top_cpgs <- names(sort(var_adj, decreasing = TRUE))[1:N]
  annoObj %>%
    filter(rownames(annoObj) %in% top_cpgs) %>%
    dplyr::count(CpG_context) %>%
    mutate(N = N,
           freq = n / sum(n))
})

# Combine into one data.frame
res_df <- bind_rows(res_list)
res_df$N <- factor(res_df$N, levels = rev(N_most_variable))

# Make stacked barplot (horizontal)
ggplot(res_df, aes(x = factor(N), y = freq, fill = CpG_context)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  labs(x = "Top N most variable CpGs",
       y = "Proportion",
       fill = "CpG context") +
  theme_bw(base_size = 14)

# Plotting unadjusted
res_list <- lapply(N_most_variable, function(N) {
  top_cpgs <- names(sort(var_unadj, decreasing = TRUE))[1:N]
  annoObj %>%
    filter(rownames(annoObj) %in% top_cpgs) %>%
    dplyr::count(CpG_context) %>%
    mutate(N = N,
           freq = n / sum(n))
})

# Combine into one data.frame
res_df <- bind_rows(res_list)
res_df$N <- factor(res_df$N, levels = rev(N_most_variable))


# Make stacked barplot (horizontal)
ggplot(res_df, aes(x = factor(N), y = freq, fill = CpG_context)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip(xlim = ) +
  labs(x = "Top N most variable CpGs",
       y = "Proportion",
       fill = "CpG context") +
  theme_bw(base_size = 14)


#
# CHANGE DISTRBUTION OF CPG CONTEXT IN MOST VARIABLE CPGS
#

# ALL CpGs after variance filtering

# Subset CpG contexts
cpg_subset <- annoObj[names(all_7), "CpG_context"]

# Turn into a data frame
df <- data.frame(CpG_context = cpg_subset)

# Count + proportions
df_counts <- df %>%
  dplyr::count(CpG_context) %>%
  mutate(proportion = n / sum(n))
df_counts$CpG_context<- factor(df_counts$CpG_context, levels = rev(df_counts$CpG_context))

# Plot barplot of proportions
plot1 <- ggplot(df_counts, aes(x = CpG_context, y = proportion * 100, fill = as.character(CpG_context))) +
  geom_col() +
  coord_flip(xlim = ) +
  labs(x = "", y = "Proportion (%)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

# Change per cassette

# Defining list of cpgs per cassette
cassette_list <- list("0"=names(all_7)[all_7 == 0], 
                      "1"=names(all_7)[all_7 == 1], 
                      "2"=names(all_7)[all_7 == 2],
                      "3"=names(all_7)[all_7 == 3], 
                      "4"=names(all_7)[all_7 == 4]
                      )

# Function to compute proportions per cassette
get_context_proportions <- function(cpg_ids, annoObj) {
  df <- data.frame(CpG_context = annoObj[cpg_ids, "CpG_context"])
  df %>%
   dplyr::count(CpG_context) %>%
    mutate(proportion = n / sum(n))
}

# Background proportions (all CpGs used as reference)
background <- data.frame(CpG_context = annoObj[, "CpG_context"]) %>%
  dplyr::count(CpG_context) %>%
  mutate(bg_proportion = n / sum(n))

# Calculate proportions for each cassette
cassette_props <- lapply(names(cassette_list), function(cassette) {
  get_context_proportions(cassette_list[[cassette]], annoObj) %>%
    mutate(cassette = cassette)
}) %>%
  bind_rows()

# Merge with background and calculate difference
cassette_change <- cassette_props %>%
  left_join(background, by = "CpG_context") %>%
  mutate(delta = proportion - bg_proportion)

# Plot
plot2 <- ggplot(cassette_change, aes(x = cassette, y = delta * 100, fill = CpG_context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Cassette", y = expression(Delta~"CpG Context (% points)")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank()
  )

wrap_plots(plot1, plot2, ncol = 2, widths = c(1, 3))
