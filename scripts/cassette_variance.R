#! usr/bin/Rscript

#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Example cassette file
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_15.rds")

#
# COMPUTE VARIANCE
#

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Only promoter CpGs
promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass=="promoter")]

variance_of_betas_promoter <- apply(betaAdj[promoter_cpgs,], MARGIN = 1, FUN=var)
cassettes <- sort(unique(promoter_15$colors))

plot(density(variance_of_betas_promoter[variance_of_betas_promoter > 0.05]))

# Initialize a dataframe to store results
cassette_stats <- data.frame(
  cassette = character(),
  mean_var = numeric(),
  sd_var = numeric(),
  stringsAsFactors = FALSE
)

for (cassette in cassettes) {
  # Subset rows belonging to the current cassette
  cassette_rows <- names(which(promoter_15$colors == cassette))
  
  if (length(cassette_rows) > 0) {
    # Compute variance of betas for the cassette
    variance_of_betas <- apply(betaAdj[cassette_rows, , drop = FALSE], MARGIN = 1, FUN = var)
    
    # Compute mean and standard deviation of variance
    mean_var <- mean(variance_of_betas, na.rm = TRUE)
    sd_var <- sd(variance_of_betas, na.rm = TRUE)
    
    # Append to the dataframe
    cassette_stats <- rbind(cassette_stats, data.frame(
      cassette = cassette,
      mean_var = mean_var,
      sd_var = sd_var
    ))
  }
}

# View the results
plot(density(cassette_stats$mean_var))

