#! usr/bin/Rscript

library(ComplexHeatmap)
library(ggplot2)
library(stats)
library(tidyr)

#
# USER DEFINED FUNCTIONS
#

# Function to summarize cassettes using the first principal component
summarize_cassettes <- function(files, output_dir) {
  for (file in files) {
    # Getting CpG cassette
    cpgs_to_cassette <- readRDS(file)$colors
    
    # Initialize a list to store summary dataframes
    summary_list <- list()
    
    # Iterate through cassettes
    for (cassette in unique(cpgs_to_cassette)) {
      # Select betas for the current cassette
      betas <- betaAdj[names(cpgs_to_cassette[cpgs_to_cassette == cassette]),]
      
      # Perform PCA and extract the first principal component
      pca <- prcomp(t(betas), center = TRUE, scale. = TRUE)
      first_pc <- pca$x[, 1]
      
      # Create a summary dataframe
      summary_df <- data.frame(Cassette = cassette, t(first_pc))
      
      # Append to the summary list
      summary_list[[as.character(cassette)]] <- summary_df
    }
    
    # Combine all summary dataframes into one
    final_summary_df <- do.call(rbind, summary_list)
    final_summary_df <- final_summary_df[order(final_summary_df$Cassette), ]
    
    # Define the output file name
    beta <- as.numeric(sub(".*cassettes_beta_(\\d+)\\.rds", "\\1", file))
    output_file <- paste0(output_dir, "/summary_beta_", beta, ".csv")
    
    # Save the summary dataframe to a CSV file
    write.csv(final_summary_df, output_file, row.names = FALSE)
  }
}

#
# LOADING DSATA
#

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")
 
# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Loading gene expression
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")

# Epitype annotations
my_annotations <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID


# Loading validation betas. Variable is called beta.adjusted
load("../../Volumes/Data/Project_3/validation_cohort/PurBeta_adjustedTumor_betaMatrix_V1_V2_reduced_717459commonCpGs_TNBCs_n136.RData")

#
# DISTAL CASSETTES
#

# List all files
distal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/distal/", full.names = TRUE)
#distal_files <- distal_files[!grepl("only", distal_files)]
distal_files <- distal_files[grepl("only_atac_unadjusted", distal_files)]


# Initialize an empty data frame
summary_df <- data.frame(beta = numeric(), num_cassettes = numeric(), mean_cassette_length = numeric())

for (file in distal_files) {
  
  # Getting beta
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Analysing cassettes
  my_data <- readRDS(file)$colors
  
  # Exclude cassette 0
  my_data <- my_data[my_data != 0]
  
  num_cassettes <- length(unique(unname(my_data)))
  mean_cassette_length <- mean(table(my_data))
  
  # Appending to df
  summary_df <- rbind(summary_df, data.frame(beta = beta, num_cassettes = num_cassettes, mean_cassette_length = mean_cassette_length))
}

# Convert beta to a factor with levels in the desired order
summary_df$beta <- factor(summary_df$beta, levels = c(5, 8, 10, 15, 20, 25))

# Create the plot
ggplot(summary_df, aes(x = beta)) +
  geom_bar(aes(y = num_cassettes), stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(aes(y = mean_cassette_length * 25), color = "red", size = 3) +  # Adjust scaling factor as needed
  geom_line(aes(y = mean_cassette_length * 25, group = 1), color = "red", size = 1) +  # Adjust scaling factor as needed
  labs(title = "Distal cassettes (Var > 0.1)",
       x = "Beta",
       y = "Number of Cassettes") +
  scale_y_continuous(limits = c(0, 1050), sec.axis = sec_axis(~ . / 25, name = "Mean Cassette Length")) +  # Adjust scaling factor as needed
  theme_classic()


# Plotting first cassettes with annotatios. 
# List all distal files
distal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/distal/", full.names = TRUE)
#distal_files <- distal_files[!grepl("only", distal_files)]
distal_files <- distal_files[grepl("only_atac_unadjusted", distal_files)]

# Loop through each file
for (file in distal_files) {
  # Extract the beta value from the filename
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Load the corresponding data
  distal <- readRDS(file)
  
  # Define the cassettes to include
  selected_cassettes <- 1:7
  
  # Extract CpGs belonging to each cassette
  cpg_list <- lapply(selected_cassettes, function(cassette) {
    names(distal$colors[distal$colors == cassette])
  })
  names(cpg_list) <- selected_cassettes  # Assign cassette names
  
  # Flatten the list to get selected CpGs
  selected_CpGs <- as.character(unlist(cpg_list))
  
  # Subset the adjusted beta values
  beta_subset <- betaNew[selected_CpGs, ]
  
  # Create a cassette grouping factor
  cassette_factor <- factor(distal$colors[selected_CpGs], levels = selected_cassettes)
  
  # Compute CpG counts per cassette
  cassette_counts <- sapply(cpg_list, length)
  row_labels <- paste0(selected_cassettes)
  
  # Subset row annotation matrix for selected transcription factors
  tf_annotation <- tfMat[selected_CpGs, c("SUZ12", "EZH2", "FOS", "STAT3", "ESR1", "GATA3", "FOXA1"), drop = FALSE]
  
  # Convert to a factor to ensure proper categorical annotation
  tf_annotation <- as.data.frame(lapply(tf_annotation, factor, levels = c(0, 1)))
  
  # Define color mapping for all TFs (0 = white, 1 = black)
  tf_colors <- list(
    SUZ12 = c("0" = "white", "1" = "black"),
    EZH2  = c("0" = "white", "1" = "black"),
    FOS   = c("0" = "white", "1" = "black"),
    STAT3 = c("0" = "white", "1" = "black"),
    ESR1  = c("0" = "white", "1" = "black"),
    GATA3 = c("0" = "white", "1" = "black"),
    FOXA1 = c("0" = "white", "1" = "black")
  )
  
  #Define ATAC annotation
  atac_annotation <- rowAnnotation(df=data.frame("ATAC"=as.factor(annoObj$hasAtacOverlap[annoObj$illuminaID %in% selected_CpGs])),
                                   col = list(ATAC = c("0" = "white", "1" = "black")))
  
  # Create row annotation object
  row_annotation <- rowAnnotation(df = tf_annotation, col = tf_colors, annotation_name_side = "top")
  
  # Epitype annotations
  pam50_annotations <- my_annotations[colnames(beta_subset), "PAM50"]
  tnbc_annotation <- my_annotations[colnames(beta_subset), "TNBC"]
  HRD_annotation <- my_annotations[colnames(beta_subset), "HRD"]
  epi_annotation <- my_annotations[colnames(beta_subset), "NMF_atacDistal"]
  im_annotation <- my_annotations[colnames(beta_subset), "IM"]
  tils_annotation <- as.numeric(x[colnames(beta_subset), "TILs"])
  
  # Create column annotation object
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                         TNBC = tnbc_annotation,
                                         HRD = HRD_annotation,
                                         Epitype = epi_annotation,
                                         IM = im_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                           "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                           "IM"=c("Negative"="grey", "Positive"="black"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate heatmap with row annotation
  heatmap <- Heatmap(beta_subset, 
                     cluster_rows = FALSE, 
                     cluster_columns = TRUE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     row_split = cassette_factor, 
                     row_title = row_labels,  
                     column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
                     top_annotation = column_annotation,
                     left_annotation = row_annotation,
                     right_annotation = atac_annotation,
                     use_raster = FALSE)
  
  # Save the heatmap to a file with double size
  pdf(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/diff_betas/heatmap_beta_", beta, "_atac_unadjusted.pdf"), width = 14, height = 10)  # Adjust width and height as needed
  draw(heatmap)
  dev.off()
}


### ANALYSISNG CASSETTES IN VALIDATION COHORT

# Plotting first cassettes with annotatios. 
# List all distal files
distal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/distal/", full.names = TRUE)
distal_files <- distal_files[grepl("only_atac", distal_files)]
distal_files <- distal_files[!grepl("unadjusted", distal_files)]
distal_files <- distal_files[!grepl("basal", distal_files)]
distal_files <- distal_files[!grepl("nonBasal", distal_files)]

# Loop through each file
for (file in distal_files) {
  # Extract the beta value from the filename
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Load the corresponding data
  distal <- readRDS(file)
  
  # Define the cassettes to include
  selected_cassettes <- 1:7
  
  # Extract CpGs belonging to each cassette
  cpg_list <- lapply(selected_cassettes, function(cassette) {
    names(distal$colors[distal$colors == cassette])
  })
  names(cpg_list) <- selected_cassettes  # Assign cassette names
  
  # Flatten the list to get selected CpGs
  selected_CpGs <- as.character(unlist(cpg_list))
  
  # Subset the adjusted beta values
  beta_subset <- beta.adjusted[selected_CpGs, ]
  
  # Create a cassette grouping factor
  cassette_factor <- factor(distal$colors[selected_CpGs], levels = selected_cassettes)
  
  # Compute CpG counts per cassette
  cassette_counts <- sapply(cpg_list, length)
  row_labels <- paste0(selected_cassettes)
  
  # Subset row annotation matrix for selected transcription factors
  tf_annotation <- tfMat[selected_CpGs, c("SUZ12", "EZH2", "FOS", "STAT3", "ESR1", "GATA3", "FOXA1"), drop = FALSE]
  
  # Convert to a factor to ensure proper categorical annotation
  tf_annotation <- as.data.frame(lapply(tf_annotation, factor, levels = c(0, 1)))
  
  # Define color mapping for all TFs (0 = white, 1 = black)
  tf_colors <- list(
    SUZ12 = c("0" = "white", "1" = "black"),
    EZH2  = c("0" = "white", "1" = "black"),
    FOS   = c("0" = "white", "1" = "black"),
    STAT3 = c("0" = "white", "1" = "black"),
    ESR1  = c("0" = "white", "1" = "black"),
    GATA3 = c("0" = "white", "1" = "black"),
    FOXA1 = c("0" = "white", "1" = "black")
  )
  
  #Define ATAC annotation
  atac_annotation <- rowAnnotation(df=data.frame("ATAC"=as.factor(annoObj$hasAtacOverlap[annoObj$illuminaID %in% selected_CpGs])),
                                   col = list(ATAC = c("0" = "white", "1" = "black")))
  
  # Create row annotation object
  row_annotation <- rowAnnotation(df = tf_annotation, col = tf_colors, annotation_name_side = "top")
  
  # Epitype annotations
  pam50_annotations <- my_annotations[colnames(beta_subset), "PAM50"]
  tnbc_annotation <- my_annotations[colnames(beta_subset), "TNBC"]
  HRD_annotation <- my_annotations[colnames(beta_subset), "HRD"]
  epi_annotation <- my_annotations[colnames(beta_subset), "NMF_atacDistal"]
  im_annotation <- my_annotations[colnames(beta_subset), "IM"]
  tils_annotation <- as.numeric(x[colnames(beta_subset), "TILs"])
  
  # Create column annotation object
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                         TNBC = tnbc_annotation,
                                         HRD = HRD_annotation,
                                         Epitype = epi_annotation,
                                         IM = im_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                           "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                           "IM"=c("Negative"="grey", "Positive"="black"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate heatmap with row annotation
  heatmap <- Heatmap(beta_subset, 
                     cluster_rows = FALSE, 
                     cluster_columns = TRUE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     row_split = cassette_factor, 
                     row_title = row_labels,  
                     column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
                     top_annotation = column_annotation,
                     left_annotation = row_annotation,
                     right_annotation = atac_annotation,
                     use_raster = FALSE)
  
  # Save the heatmap to a file with double size
  pdf(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/diff_betas/heatmap_beta_", beta, "_atac_VALIDATION.pdf"), width = 14, height = 10)  # Adjust width and height as needed
  draw(heatmap)
  dev.off()
}


#
# PROMOTER CASSETTES
#

# List all files
promoter_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/promoter/", full.names = TRUE)
promoter_files <- promoter_files[!grepl("only", promoter_files)]
promoter_files <- promoter_files[!grepl("not_purity_adjusted", promoter_files)]

# Initialize an empty data frame
summary_df <- data.frame(beta = numeric(), num_cassettes = numeric(), mean_cassette_length = numeric())

for (file in promoter_files) {
  
  # Getting beta
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Analysing cassettes
  my_data <- readRDS(file)$colors
  
  # Exclude cassette 0
  my_data <- my_data[my_data != 0]
  
  num_cassettes <- length(unique(unname(my_data)))
  mean_cassette_length <- mean(table(my_data))
  
  # Appending to df
  summary_df <- rbind(summary_df, data.frame(beta = beta, num_cassettes = num_cassettes, mean_cassette_length = mean_cassette_length))
}

# Convert beta to a factor with levels in the desired order
summary_df$beta <- factor(summary_df$beta, levels = c(5, 8, 10, 15, 20, 25))

# Create the plot
ggplot(summary_df, aes(x = beta)) +
  geom_bar(aes(y = num_cassettes), stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(aes(y = mean_cassette_length * 25), color = "red", size = 3) +  # Adjust scaling factor as needed
  geom_line(aes(y = mean_cassette_length * 25, group = 1), color = "red", size = 1) +  # Adjust scaling factor as needed
  labs(title = "Promoter cassettes (Var > 0.05)",
       x = "Beta",
       y = "Number of Cassettes") +
  scale_y_continuous(limits = c(0, 1300), sec.axis = sec_axis(~ . / 25, name = "Mean Cassette Length")) +  # Adjust scaling factor as needed
  theme_classic()


# Plotting first cassettes with annotatios.
# List all promoter files
promoter_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/promoter/", full.names = TRUE)
promoter_files <- promoter_files[!grepl("only", promoter_files)]
promoter_files <- promoter_files[grepl("not_purity_adjusted", promoter_files)]


# Loop through each file
for (file in promoter_files) {
  # Extract the beta value from the filename
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Load the corresponding data
  promoter <- readRDS(file)
  
  # Define the cassettes to include
  selected_cassettes <- 1:7
  
  # Extract CpGs belonging to each cassette
  cpg_list <- lapply(selected_cassettes, function(cassette) {
    names(promoter$colors[promoter$colors == cassette])
  })
  names(cpg_list) <- selected_cassettes  # Assign cassette names
  
  # Flatten the list to get selected CpGs
  selected_CpGs <- unlist(cpg_list)
  
  # Subset the adjusted beta values
  beta_subset <- betaNew[selected_CpGs, ]
  
  # Create a cassette grouping factor
  cassette_factor <- factor(promoter$colors[selected_CpGs], levels = selected_cassettes)
  
  # Compute CpG counts per cassette
  cassette_counts <- sapply(cpg_list, length)
  row_labels <- paste0(selected_cassettes)
  
  # Subset row annotation matrix for selected transcription factors
  tf_annotation <- tfMat[selected_CpGs, c("SUZ12", "EZH2", "FOS", "STAT3", "ESR1", "GATA3", "FOXA1"), drop = FALSE]
  
  # Convert to a factor to ensure proper categorical annotation
  tf_annotation <- as.data.frame(lapply(tf_annotation, factor, levels = c(0, 1)))
  
  # Define color mapping for all TFs (0 = white, 1 = black)
  tf_colors <- list(
    SUZ12 = c("0" = "white", "1" = "black"),
    EZH2  = c("0" = "white", "1" = "black"),
    FOS   = c("0" = "white", "1" = "black"),
    STAT3 = c("0" = "white", "1" = "black"),
    ESR1  = c("0" = "white", "1" = "black"),
    GATA3 = c("0" = "white", "1" = "black"),
    FOXA1 = c("0" = "white", "1" = "black")
  )
  
  # Create row annotation object
  row_annotation <- rowAnnotation(df = tf_annotation, col = tf_colors, annotation_name_side = "top")
  
  # Epitype annotations
  pam50_annotations <- my_annotations[colnames(beta_subset), "PAM50"]
  tnbc_annotation <- my_annotations[colnames(beta_subset), "TNBC"]
  HRD_annotation <- my_annotations[colnames(beta_subset), "HRD"]
  epi_annotation <- my_annotations[colnames(beta_subset), "NMF_atacDistal"]
  im_annotation <- my_annotations[colnames(beta_subset), "IM"]
  tils_annotation <- as.numeric(x[colnames(beta_subset), "TILs"])
  
  # Create column annotation object
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                         TNBC = tnbc_annotation,
                                         HRD = HRD_annotation,
                                         Epitype = epi_annotation,
                                         IM = im_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                           "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                           "IM"=c("Negative"="grey", "Positive"="black"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate heatmap with row annotation
  heatmap <- Heatmap(beta_subset, 
                     cluster_rows = FALSE, 
                     cluster_columns = TRUE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     row_split = cassette_factor, 
                     row_title = row_labels,  
                     column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
                     top_annotation = column_annotation,
                     #left_annotation = row_annotation,
                     use_raster = FALSE)
  
  # Save the heatmap to a file with double size
  pdf(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/diff_betas/heatmap_beta_", beta, "_unadjusted_from_unadjusted.pdf"), width = 14, height = 10)  # Adjust width and height as needed
  draw(heatmap)
  dev.off()
}



#
# PROXIMAL CASSETTES
#

# List all files
proximal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/proximal/", full.names = TRUE)
proximal_files <- proximal_files[!grepl("only", proximal_files)]
proximal_files <- proximal_files[!grepl("not_purity_adjusted", proximal_files)]


# Initialize an empty data frame
summary_df <- data.frame(beta = numeric(), num_cassettes = numeric(), mean_cassette_length = numeric())

for (file in proximal_files) {
  
  # Getting beta
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Analysing cassettes
  my_data <- readRDS(file)$colors
  
  # Exclude cassette 0
  my_data <- my_data[my_data != 0]
  
  num_cassettes <- length(unique(unname(my_data)))
  mean_cassette_length <- mean(table(my_data))
  
  # Appending to df
  summary_df <- rbind(summary_df, data.frame(beta = beta, num_cassettes = num_cassettes, mean_cassette_length = mean_cassette_length))
}

# Convert beta to a factor with levels in the desired order
summary_df$beta <- factor(summary_df$beta, levels = c(5, 8, 10, 15, 20, 25))

# Create the plot
ggplot(summary_df, aes(x = beta)) +
  geom_bar(aes(y = num_cassettes), stat = "identity", fill = "blue", alpha = 0.6) +
  geom_point(aes(y = mean_cassette_length * 25), color = "red", size = 3) +  # Adjust scaling factor as needed
  geom_line(aes(y = mean_cassette_length * 25, group = 1), color = "red", size = 1) +  # Adjust scaling factor as needed
  labs(title = "proximal cassettes (Var > 0.05)",
       x = "Beta",
       y = "Number of Cassettes") +
  theme_classic()


# Plotting first cassettes with annotatios. 
# List all proximal files
proximal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/proximal/", full.names = TRUE)
proximal_files <- proximal_files[!grepl("only", proximal_files)]
proximal_files <- proximal_files[grepl("not_purity_adjusted", proximal_files)]

# Loop through each file
for (file in proximal_files) {
  # Extract the beta value from the filename
  beta <- as.numeric(sub(".*cassettes_beta_(\\d+).*\\.rds", "\\1", file))
  
  # Load the corresponding data
  proximal <- readRDS(file)
  
  # Define the cassettes to include
  selected_cassettes <- 1:7
  
  # Extract CpGs belonging to each cassette
  cpg_list <- lapply(selected_cassettes, function(cassette) {
    names(proximal$colors[proximal$colors == cassette])
  })
  names(cpg_list) <- selected_cassettes  # Assign cassette names
  
  # Flatten the list to get selected CpGs
  selected_CpGs <- unlist(cpg_list)
  
  # Subset the adjusted beta values
  beta_subset <- betaNew[selected_CpGs, ]
  
  # Create a cassette grouping factor
  cassette_factor <- factor(proximal$colors[selected_CpGs], levels = selected_cassettes)
  
  # Compute CpG counts per cassette
  cassette_counts <- sapply(cpg_list, length)
  row_labels <- paste0(selected_cassettes)
  
  # Subset row annotation matrix for selected transcription factors
  tf_annotation <- tfMat[selected_CpGs, c("SUZ12", "EZH2", "FOS", "STAT3", "ESR1", "GATA3", "FOXA1"), drop = FALSE]
  
  # Convert to a factor to ensure proper categorical annotation
  tf_annotation <- as.data.frame(lapply(tf_annotation, factor, levels = c(0, 1)))
  
  # Define color mapping for all TFs (0 = white, 1 = black)
  tf_colors <- list(
    SUZ12 = c("0" = "white", "1" = "black"),
    EZH2  = c("0" = "white", "1" = "black"),
    FOS   = c("0" = "white", "1" = "black"),
    STAT3 = c("0" = "white", "1" = "black"),
    ESR1  = c("0" = "white", "1" = "black"),
    GATA3 = c("0" = "white", "1" = "black"),
    FOXA1 = c("0" = "white", "1" = "black")
  )
  
  # Create row annotation object
  row_annotation <- rowAnnotation(df = tf_annotation, col = tf_colors, annotation_name_side = "top")
  
  # Epitype annotations
  pam50_annotations <- my_annotations[colnames(beta_subset), "PAM50"]
  tnbc_annotation <- my_annotations[colnames(beta_subset), "TNBC"]
  HRD_annotation <- my_annotations[colnames(beta_subset), "HRD"]
  epi_annotation <- my_annotations[colnames(beta_subset), "NMF_atacDistal"]
  im_annotation <- my_annotations[colnames(beta_subset), "IM"]
  tils_annotation <- as.numeric(x[colnames(beta_subset), "TILs"])
  
  # Create column annotation object
  column_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                         TNBC = tnbc_annotation,
                                         HRD = HRD_annotation,
                                         Epitype = epi_annotation,
                                         IM = im_annotation,
                                         col = list(
                                           "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                           "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                           "IM"=c("Negative"="grey", "Positive"="black"),
                                           "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                           "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                       "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue"))
  )
  
  # Generate heatmap with row annotation
  heatmap <- Heatmap(beta_subset, 
                     cluster_rows = FALSE, 
                     cluster_columns = TRUE, 
                     show_row_names = FALSE, 
                     show_column_names = FALSE, 
                     row_split = cassette_factor, 
                     row_title = row_labels,  
                     column_title = paste("CpG Methylation Heatmap by Cassettes (Beta =", beta, ")"),
                     top_annotation = column_annotation,
                     #left_annotation = row_annotation,
                     use_raster = FALSE)
  
  # Save the heatmap to a file with double size
  pdf(paste0("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/diff_betas/heatmap_beta_", beta, "_unadjusted_from_unadjusted.pdf"), width = 14, height = 10)  # Adjust width and height as needed
  draw(heatmap)
  dev.off()
}


#
# SUMMARIZING CASSETTES
#

# Distal cassettes
distal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/distal/", full.names = TRUE)
summarize_cassettes(distal_files, "/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes")

# Promoter cassettes
promoter_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/promoter/", full.names = TRUE)
summarize_cassettes(promoter_files, "/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes")

# Proximal cassettes
proximal_files <- list.files("/Volumes/Data/Project_3/detected_cassettes/proximal/", full.names = TRUE)
summarize_cassettes(proximal_files, "/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes")
