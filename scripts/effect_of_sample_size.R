#! usr/bin/Rscript

library(WGCNA)
library(ComplexHeatmap)
library(pheatmap)
library(Matrix)


#
# LOADING DATA
#

# Loading corrected betas and annotations
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation files
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

orig_epi <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")



#
# CASSETTE DETECTION WITH RESAMPLINGS. All Distal
#

cor <- WGCNA::cor

# Define number for resampling
samples_n <- c(60, 100, 120)

for (sample_n in samples_n) {
  
  # Define runs
  runs <- 1:20
  
  for (run in runs) {
    
    # Defining proportions for sampling
    proportions_of_epitypes <- prop.table(table(orig_epi$NMF_atacDistal))
    
    # Sample indices with probabilities based on the corresponding epitype
    sampled_indices <- sample(
      seq_len(nrow(orig_epi)),
      size = sample_n,
      replace = FALSE,
      prob = proportions_of_epitypes[as.character(orig_epi$NMF_atacDistal)]
    )
    
    # Getting beta values of sampled indices
    sampled_betas <- betaAdj[,sampled_indices]
    
    # Using only distal cpgs
    distal_cpgs <- annoObj$illuminaID[which(annoObj$featureClass == "distal" | annoObj$featureClass == "distal body")]
    
    sampled_betas <- sampled_betas[rownames(sampled_betas) %in% distal_cpgs, ]
    
    # Calculating variance of CpGs
    variance_dis <- sapply(1:nrow(sampled_betas), FUN = function(row) {var(sampled_betas[row,])})
    
    # Filtering CpGs with low variance
    sampled_betas <-  t(sampled_betas[variance_dis > 0.1,])
    
    # Defining betas to use 
    betas <- c(5, 10, 15)
    
    # Running WGCNA
    for (beta in betas) {
      
      message("Running: sample_n = ", sample_n, ", run = ", run, ", beta = ", beta)
      
      netwk <- blockwiseModules(sampled_betas,               
                                corrType="bicor", # Using biweight midcorrelation 
                                nThreads = 10,
                                
                                # == Adjacency Function ==
                                power = beta,             
                                networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
                                
                                # == Tree and Block Options ==
                                deepSplit = 2,
                                pamRespectsDendro = F,
                                minModuleSize = 3,
                                maxBlockSize = 6000,
                                
                                # == Module Adjustments ==
                                reassignThreshold = 0,
                                mergeCutHeight = 0.25,
                                
                                # == Output Options
                                numericLabels = T,
                                verbose = 0)
      
      # Saving network
      my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/", "n_", sample_n, "_beta_", beta, "_resampling_", run, ".rds" )
      saveRDS(netwk, file = my_filename)
      
    }
  }
}



# #
# # ANALYSE OUTPUT
# #
# 
# # List files of interest
# my_files <- list.files("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/")
# my_files_60 <- my_files[grepl("^n_100_beta_10_*", my_files)]
# my_files_100 <- my_files[grepl("^n_100_beta_10_*", my_files)]
# my_files_120 <- my_files[grepl("^n_100_beta_10_*", my_files)]
# 
# # Lists of cassettes and samples used
# list_of_cassettes <- list()
# list_of_samples <- list()
# 
# for (file in my_files_100) {
#   
#   # Appending cassettes 
#   list_of_cassettes[[file]] <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
#                                               file))$colors
#   # Appending samples
#   list_of_samples[[file]] <- rownames(readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
#                                                        file))$MEs)
# }
# 
# # Comparing CpGs included after variance filtering
# all_cpgs <- unique(unlist(lapply(list_of_cassettes, names)))
# 
# df_included_cpgs <- as.data.frame(matrix(ncol=length(my_files),
#                                          nrow=length(all_cpgs)))
# 
# colnames(df_included_cpgs) <- my_files
# rownames(df_included_cpgs) <- all_cpgs
# 
# # Counting the included CpGs
# for (file_name in colnames(df_included_cpgs)) {
#   
#   df_included_cpgs[,file_name] <- all_cpgs %in% names(list_of_cassettes[[file_name]]) 
#   
# }
# 
# # Calculating number of kept CpGs in each resampling
# used_cpgs <- rowSums(df_included_cpgs)
# 
# # Plotting hist + density
# hist(rowSums(df_included_cpgs),
#      freq = FALSE,
#      breaks = 20,
#      col = "lightblue",
#      border = "black",
#      xaxt = "n",  # suppress x-axis
#      main = NULL,
#      xlab = "Resamplings",
#      ylab="Fraction of CpGs")
# 
# # Add your custom x-axis
# axis(side = 1, at = c(1, 5, 10, 15, 20))
# 
# # Add the density line
# lines(density(rowSums(df_included_cpgs)), col = "red", lwd = 2)
# 
# 
# 
# # # Plotting heatmaps
# # 
# # # Define the cassettes to include
# # selected_cassettes <- 1:7
# # file_names <- "n_60_beta_10_resampling_3.rds"
# # 
# # # Extract CpGs belonging to each cassette
# # cpg_list <- lapply(selected_cassettes, function(cassette) {
# #   names(list_of_cassettes[[file_names]][list_of_cassettes[[file_names]] == cassette])
# # })
# # 
# # names(cpg_list) <- selected_cassettes  # Assign cassette names
# # 
# # # Flatten the list to get selected CpGs
# # selected_CpGs <- as.character(unlist(cpg_list))
# # 
# # # Subset the adjusted beta values
# # beta_subset <- betaAdj[selected_CpGs, list_of_samples[[file_names]]]
# # 
# # # Create a cassette grouping factor
# # cassette_factor <- factor(list_of_cassettes[[file_names]][selected_CpGs], levels = selected_cassettes)
# # 
# # # Compute CpG counts per cassette
# # cassette_counts <- sapply(cpg_list, length)
# # row_labels <- paste0(selected_cassettes)
# # 
# # # Determine column annotation
# # col_annotation <- HeatmapAnnotation(df=data.frame("Epitype" = orig_epi[list_of_samples[[file_names]], "NMF_atacDistal"]),
# #                                     col=list("Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
# #                                                 "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")))
# # 
# # Heatmap(beta_subset, 
# #         cluster_rows = FALSE, 
# #         cluster_columns = TRUE, 
# #         show_row_names = FALSE, 
# #         show_column_names = FALSE, 
# #         row_split = cassette_factor, 
# #         row_title = row_labels,  
# #         column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
# #         top_annotation = col_annotation,
# #         use_raster = FALSE)
# 
# 
# # SAVING CASSETTE ASSIGNMENTS TO CSV
# 
# list_of_cassettes <- list()
# 
# # Getting cassetttes for n_resapmlings of interest
# for (file in my_files_100) {
#   
#   # Appending cassettes 
#   list_of_cassettes[[file]] <- readRDS(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/",
#                                               file))$colors
# }
# 
# # Getting all CpGs used in every resampling
# all_cpgs <- unique(unlist(lapply(list_of_cassettes, names)))
# 
# # Generating df to store the data
# df_of_cassette_assignments <- as.data.frame(matrix(NA, nrow = length(all_cpgs), ncol = length(list_of_cassettes)))
# rownames(df_of_cassette_assignments) <- all_cpgs
# 
# # Populating dataframe
# for (n_col in seq_along(list_of_cassettes)) {
#   cassette <- list_of_cassettes[[n_col]]
#   common_cpgs <- intersect(all_cpgs, names(cassette))
#   df_of_cassette_assignments[common_cpgs, n_col] <- cassette[common_cpgs]
# }
# 
# # Saving file as csv
# write.csv(df_of_cassette_assignments,
#           "/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/consensus/n_res100_cassettes.csv")
# 
# 
# ### SAVING TOP VARIANCE 5000 CpGs FOR ANALYSIS IN PYTHON
# 
# 
# # Using only distal cpgs
# distal_cpgs <- annoObj$illuminaID[which(annoObj$featureClass == "distal" | annoObj$featureClass == "distal body")]
# 
# filtered_betas <- betaAdj[rownames(betaAdj) %in% distal_cpgs, ]
# 
# # Calculate variance
# variace_of_cpgs <- apply(filtered_betas, MARGIN = 1, FUN = var)
# 
# # Getting first 5000 and 10000 cpgs based on variance
# top_var_5000 <- rownames(filtered_betas)[order(variace_of_cpgs)[1:5000]]
# top_var_10000 <- rownames(filtered_betas)[order(variace_of_cpgs)[1:10000]]
# 
# # Saving vectors as txt files
# write.table(top_var_5000, file = "/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/consensus/distal_top_var_5000.txt", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# 
# write.table(top_var_10000, file = "/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/consensus/distal_top_var_10000.txt", 
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
