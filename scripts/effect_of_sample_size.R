#! usr/bin/Rscript

library(WGCNA)
library(ComplexHeatmap)
library(pheatmap)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


#
# LOADING DATA
#

# Loading corrected betas and annotations
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation files
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

orig_epi <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Loading Distal cassettes. Beta = 10
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")


#
# CASSETTE DETECTION WITH RESAMPLINGS. All Distal
#

# Variance filtering after resampling

# cor <- WGCNA::cor
# 
# # Define number for resampling
# samples_n <- c(60, 100, 120, 140)
# 
# for (sample_n in samples_n) {
#   
#   # Define runs
#   runs <- 1:20
#   
#   for (run in runs) {
#     
#     # Defining proportions for sampling
#     proportions_of_epitypes <- prop.table(table(orig_epi$NMF_atacDistal))
#     
#     # Sample indices with probabilities based on the corresponding epitype
#     sampled_indices <- sample(
#       seq_len(nrow(orig_epi)),
#       size = sample_n,
#       replace = FALSE,
#       prob = proportions_of_epitypes[as.character(orig_epi$NMF_atacDistal)]
#     )
#     
#     # Getting beta values of sampled indices
#     sampled_betas <- betaAdj[,sampled_indices]
#     
#     # Using only distal cpgs
#     distal_cpgs <- annoObj$illuminaID[which(annoObj$featureClass == "distal" | annoObj$featureClass == "distal body")]
#     
#     sampled_betas <- sampled_betas[rownames(sampled_betas) %in% distal_cpgs, ]
#     
#     # Calculating variance of CpGs
#     variance_dis <- sapply(1:nrow(sampled_betas), FUN = function(row) {var(sampled_betas[row,])})
# 
#     # Filtering CpGs with low variance
#     sampled_betas <-  t(sampled_betas[variance_dis > 0.1,])
#     
#     # Defining betas to use 
#     betas <- c(5, 10, 15)
#     
#     # Running WGCNA
#     for (beta in betas) {
#       
#       message("Running: sample_n = ", sample_n, ", run = ", run, ", beta = ", beta)
#       
#       netwk <- blockwiseModules(sampled_betas,               
#                                 corrType="bicor", # Using biweight midcorrelation 
#                                 nThreads = 10,
#                                 
#                                 # == Adjacency Function ==
#                                 power = beta,             
#                                 networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
#                                 
#                                 # == Tree and Block Options ==
#                                 deepSplit = 2,
#                                 pamRespectsDendro = F,
#                                 minModuleSize = 3,
#                                 maxBlockSize = 6000,
#                                 
#                                 # == Module Adjustments ==
#                                 reassignThreshold = 0,
#                                 mergeCutHeight = 0.25,
#                                 
#                                 # == Output Options
#                                 numericLabels = T,
#                                 verbose = 0)
#       
#       # Saving network
#       my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/", "n_", sample_n, "_beta_", beta, "_resampling_", run, ".rds" )
#       saveRDS(netwk, file = my_filename)
#       
#     }
#   }
# }

# Variance filtering before resampling

# cor <- WGCNA::cor
# 
# # Define number for resampling
# samples_n <- c(60, 100, 140)
# 
# # Run variance based filtering before running WGCNA
# variance_dis <- sapply(1:nrow(betaAdj), FUN = function(row) {var(betaAdj[row,])})
# 
# for (sample_n in samples_n) {
# 
#   # Define runs
#   runs <- 1:20
# 
#   for (run in runs) {
#     
#     # Filtering CpGs with low variance
#     sampled_betas <-  betaAdj[variance_dis > 0.1,]
#     
#     # Defining proportions for sampling
#     proportions_of_epitypes <- prop.table(table(orig_epi$NMF_atacDistal))
# 
#     # Sample indices with probabilities based on the corresponding epitype
#     sampled_indices <- sample(
#       seq_len(nrow(orig_epi)),
#       size = sample_n,
#       replace = FALSE,
#       prob = proportions_of_epitypes[as.character(orig_epi$NMF_atacDistal)]
#     )
# 
#     # Getting beta values of sampled indices
#     sampled_betas <- sampled_betas[,sampled_indices]
#     
#     # Using only distal cpgs
#     distal_cpgs <- annoObj$illuminaID[which(annoObj$featureClass == "distal" | annoObj$featureClass == "distal body")]
# 
#     sampled_betas <- t(sampled_betas[rownames(sampled_betas) %in% distal_cpgs, ])
# 
#     # Defining betas to use
#     betas <- c(10)
# 
#     # Running WGCNA
#     for (beta in betas) {
# 
#       message("Running: sample_n = ", sample_n, ", run = ", run, ", beta = ", beta)
# 
#       netwk <- blockwiseModules(sampled_betas,
#                                 corrType="bicor", # Using biweight midcorrelation
#                                 nThreads = 10,
# 
#                                 # == Adjacency Function ==
#                                 power = beta,
#                                 networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
# 
#                                 # == Tree and Block Options ==
#                                 deepSplit = 2,
#                                 pamRespectsDendro = F,
#                                 minModuleSize = 3,
#                                 maxBlockSize = 6000,
# 
#                                 # == Module Adjustments ==
#                                 reassignThreshold = 0,
#                                 mergeCutHeight = 0.25,
# 
#                                 # == Output Options
#                                 numericLabels = T,
#                                 verbose = 0)
# 
#       # Saving network
#       my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/", "filtered_before_resampling_n_", sample_n, "_beta_", beta, "_resampling_", run, ".rds" )
#       saveRDS(netwk, file = my_filename)
# 
#     }
#   }
# }

#
# ANALYSE OUTPUT
#


### USED CPGS PER RESAMPLING

# List files of interest
my_files <- list.files("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/")
my_files_60 <- my_files[grepl("^n_60_beta_10_*", my_files)]
my_files_100 <- my_files[grepl("^n_100_beta_10_*", my_files)]
my_files_120 <- my_files[grepl("^n_120_beta_10_*", my_files)]
my_files_140 <- my_files[grepl("^n_140_beta_10_*", my_files)]

my_files_60_filbefore <- my_files[grepl("^filtered_before_resampling_n_60_beta_10_*", my_files)]
my_files_100_filbefore <- my_files[grepl("^filtered_before_resampling_n_100_beta_10_*", my_files)]
my_files_140_filbefore <- my_files[grepl("^filtered_before_resampling_n_140_beta_10_*", my_files)]



# Lists of cassettes and samples used
list_of_cassettes <- list()
list_of_samples <- list()

for (file in my_files_100) {

  # Appending cassettes
  list_of_cassettes[[file]] <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
                                              file))$colors
  # Appending samples
  list_of_samples[[file]] <- rownames(readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
                                                       file))$MEs)
}

# Comparing CpGs included after variance filtering
all_cpgs <- unique(unlist(lapply(list_of_cassettes, names)))

df_included_cpgs <- as.data.frame(matrix(ncol=length(my_files),
                                         nrow=length(all_cpgs)))

colnames(df_included_cpgs) <- my_files
rownames(df_included_cpgs) <- all_cpgs

# Counting the included CpGs
for (file_name in colnames(df_included_cpgs)) {

  df_included_cpgs[,file_name] <- all_cpgs %in% names(list_of_cassettes[[file_name]])

}

# Calculating number of kept CpGs in each resampling
used_cpgs <- rowSums(df_included_cpgs)

# Plotting hist + density
hist(rowSums(df_included_cpgs),
     freq = FALSE,
     breaks = 20,
     col = "lightblue",
     border = "black",
     xaxt = "n",  # suppress x-axis
     main = NULL,
     xlab = "Resamplings",
     ylab="Fraction of CpGs",
     ylim=c(0,0.6))

# Add your custom x-axis
axis(side = 1, at = c(1, 5, 10, 15, 20))

# Add the density line
lines(density(rowSums(df_included_cpgs)), col = "red", lwd = 2)


# ### HEATMAPS OF CASSETTES
# 
# # Plotting heatmaps
# 
# # Define the cassettes to include
# selected_cassettes <- 1:7
# file_names <- "n_60_beta_10_resampling_10.rds"
# 
# # Extract CpGs belonging to each cassette
# cpg_list <- lapply(selected_cassettes, function(cassette) {
#   names(list_of_cassettes[[file_names]][list_of_cassettes[[file_names]] == cassette])
# })
# 
# names(cpg_list) <- selected_cassettes  # Assign cassette names
# 
# # Flatten the list to get selected CpGs
# selected_CpGs <- as.character(unlist(cpg_list))
# 
# # Subset the adjusted beta values
# beta_subset <- betaAdj[selected_CpGs, list_of_samples[[file_names]]]
# 
# # Create a cassette grouping factor
# cassette_factor <- factor(list_of_cassettes[[file_names]][selected_CpGs], levels = selected_cassettes)
# 
# # Compute CpG counts per cassette
# cassette_counts <- sapply(cpg_list, length)
# row_labels <- paste0(selected_cassettes)
# 
# # Determine column annotation
# col_annotation <- HeatmapAnnotation(df=data.frame("Epitype" = orig_epi[list_of_samples[[file_names]], "NMF_atacDistal"]),
#                                     col=list("Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2",
#                                                 "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")))
# 
# Heatmap(beta_subset,
#         cluster_rows = FALSE,
#         cluster_columns = TRUE,
#         show_row_names = FALSE,
#         show_column_names = FALSE,
#         row_split = cassette_factor,
#         row_title = row_labels,
#         column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
#         top_annotation = col_annotation,
#         use_raster = FALSE)


### CONSENSUS ANALYSIS ACROSS RESAMPLINGS 

list_of_cassettes <- list()

# Getting cassetttes for n_resapmlings of interest
for (file in my_files_100) {
  
  # Appending cassettes
  list_of_cassettes[[file]] <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
                                              file))$colors
}

# Getting all CpGs used in every resampling
all_cpgs <- unique(unlist(lapply(list_of_cassettes, names)))

# Generating df to store the data
df_of_cassette_assignments <- as.data.frame(matrix(NA, nrow = length(all_cpgs), ncol = length(list_of_cassettes)))
rownames(df_of_cassette_assignments) <- all_cpgs
colnames(df_of_cassette_assignments) <- names(list_of_cassettes)

# Populating dataframe
for (n_col in names(list_of_cassettes)) {
  cassette <- list_of_cassettes[[n_col]]
  common_cpgs <- intersect(all_cpgs, names(cassette))
  df_of_cassette_assignments[common_cpgs, n_col] <- cassette[common_cpgs]
}


# Defining dataframe to save cassette assignment for the whole cohort
original_assigments_df <- data.frame("CpG" = names(distal_10$colors),
                           "Cassettes" = distal_10$colors)

# Define cassette to analyse
cassette_to_analyse <- 3

# Initialize empty list
cassette_prop_df <- list()  

for (res_num in 1:20) {
  # Calculate proportions
  proportions <- prop.table(table(
    df_of_cassette_assignments[
      rownames(original_assigments_df[original_assigments_df$Cassettes == cassette_to_analyse, ]), res_num
    ],
    useNA = "ifany"
  ))
  
  # Get top cluster (excluding "0" and NA)
  valid_names <- names(proportions)[!names(proportions) %in% c("0", NA)]
  top_cassette <- names(which.max(proportions[valid_names]))
  
  # Identify low frequency clusters to group as "Other"
  low_names <- setdiff(valid_names, top_cassette)
  
  # Grouping
  grouped <- proportions
  grouped["Others"] <- sum(grouped[low_names])
  grouped <- grouped[!(names(grouped) %in% low_names)]
  
  # Convert to data.frame
  df <- data.frame(
    Cluster = names(grouped),
    Proportion = as.numeric(grouped),
    Resamplings = colnames(df_of_cassette_assignments)[res_num]
  )
  
  # Rename clusters
  df$Cluster[df$Cluster == top_cassette] <- "Most_common_cassette"
  df$Cluster[df$Cluster == "0"] <- "Unclassified"
  df$Cluster[is.na(df$Cluster)] <- "Filtered"
  
  cassette_prop_df[[res_num]] <- df
}


# Combine all into one data frame
cassette_prop_df <- bind_rows(cassette_prop_df)

# Sort by proportion of most_common_cluster
resampling_order <- cassette_prop_df %>%
  filter(Cluster == "Most_common_cassette") %>%
  arrange(desc(Proportion))

# Plot scatterplot with lines connecting points for each Cluster
plot1 <- ggplot(cassette_prop_df, aes(x = factor(cassette_prop_df$Resampling, resampling_order$"Resamplings"), y = Proportion, color = Cluster, group = Cluster)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(y = "Cassette assignment \nper resampling", x = "") +
  theme(axis.text.x = element_blank()) +
  scale_color_brewer(palette = "Set1")

# Plotting 

list_of_most_common_cassette_composition <- list()

# Getting most common non NA, and 0 cassette in each resampling
for (resampling_file in unique(resampling_order$Resamplings)) {
  
  # Read file
  my_cassettes <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
      resampling_file))$colors
  
  # Calculating proportions
  proportions <- prop.table(table(
      my_cassettes[original_assigments_df$CpG[original_assigments_df$Cassettes == cassette_to_analyse]],
    useNA = "ifany"
  ))
  
  valid_names <- names(proportions)[!names(proportions) %in% c("0", NA)]
  top_cassette <- names(which.max(proportions[valid_names]))
  
  # Getting CpGs from the most common cassette
  cpgs_to_plot <- names(my_cassettes)[my_cassettes == top_cassette]
  
  # Saving in the list
  list_of_most_common_cassette_composition[[resampling_file]] <- original_assigments_df[cpgs_to_plot, ]
  
}

# Transforming data to list of cassette proportions
list_of_most_common_cassette_composition <- lapply(list_of_most_common_cassette_composition, function(res) {prop.table(table(unname(res$Cassettes), useNA = "ifany"))})

# Convert list of named vectors into a tidy data frame
to_plot_df <- lapply(names(list_of_most_common_cassette_composition), function(resamp) {
  props <- list_of_most_common_cassette_composition[[resamp]]
  data.frame(
    Resampling = resamp,
    Cassette = names(props),
    Proportion = as.numeric(props),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# Grouping
to_plot_df$Cassette_grouped <- ifelse(
  to_plot_df$Cassette %in% c("NA", "0"), to_plot_df$Cassette,
  ifelse(to_plot_df$Proportion < 0.05, "Others", to_plot_df$Cassette)
)

to_plot_df <- to_plot_df %>%
  group_by(Resampling, Cassette_grouped) %>%
  summarise(Proportion = sum(Proportion), .groups = "drop")

plot2 <- ggplot(to_plot_df, aes(x = factor(to_plot_df$Resampling, resampling_order$"Resamplings"), y = Proportion, fill = Cassette_grouped)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Resampling", y = "Composition of most \ncommon cassette", fill = "Cassette") +
  theme(axis.text.x = element_blank()) +
  scale_fill_brewer(palette = "Set3")


# Plotting plots
plot1 / plot2

# Why arent we detecting cassettes in some resamplings? Sample size and filtering

# Define list of variance filtered plots
list_of_variance_filtered_plots <- list()

for (cassette in 2:5) {
  
  list_of_variance_filtered_plots[[cassette]] <- list()
  
  for (file in my_files_140) {
    
    a <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/", file))
    
    subset_of_interest <- betaAdj[rownames(cassettes_df[original_assigments_df$Cassettes == cassette, ]),rownames(a$MEs)]
    
    var_of_subset <- apply(subset_of_interest, MARGIN = 1, FUN = var)
    
    # Row annotation: variance > 0.1
    high_variance_flag <- var_of_subset > 0.1
    high_variance_annotation <- rowAnnotation(
      HighVariance = high_variance_flag,
      col = list(HighVariance = c("TRUE" = "red", "FALSE" = "lightgrey")),
      width = unit(0.5, "cm")
    )
    
    # Row annotation: variance barplot
    variance_barplot_annotation <- rowAnnotation(
      Variance = anno_barplot(
        var_of_subset,
        border = FALSE,
        gp = gpar(fill = "grey"),
        axis_param = list(
          at = c(0, 0.1, max(var_of_subset)),
          labels = c("0", "0.1", round(max(var_of_subset), 2))
        )
      ),
      annotation_name_side = "top"
    )
    
    # Final heatmap with both annotations
    list_of_variance_filtered_plots[[cassette]][[file]] <- Heatmap(
      subset_of_interest,
      use_raster = FALSE,
      show_row_dend = FALSE,
      show_column_dend = FALSE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      left_annotation = high_variance_annotation,
      right_annotation = variance_barplot_annotation
    )
  }
}

