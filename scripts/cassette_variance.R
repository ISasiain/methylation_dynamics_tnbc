#! usr/bin/Rscript

#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID


#
# ANALYSING PROMOTER
#

# Example cassette file
promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_5.rds")

# COMPUTE VARIANCE

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Only promoter CpGs
promoter_cpgs <- annoObj$illuminaID[which(annoObj$featureClass=="promoter")]

variance_of_betas_promoter <- apply(betaAdj[promoter_cpgs,], MARGIN = 1, FUN=var)
cassettes <- sort(unique(promoter_15$colors))

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

# PLOT CASSETTE VARIANCE

# Plot the results
par(mfrow=c(2,1), mar=c(3,3,3,3))

# Define cassette to plot
cas_num <- 1

# Compute density of variances in every cassette
density_data <- density(cassette_stats$mean_var)
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "Mean cassette variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size


# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)


lines(c(rep(cassette_stats[cas_num + 1,"mean_var"], 2)), c(0, density_data$y[which.min(abs(density_data$x - cassette_stats[cas_num + 1,"mean_var"]))]), col = "red", lwd = 2, lty = 1)



# Compute density for variance_of_betas_promoter values > 0.05
density_data <- density(variance_of_betas_promoter[variance_of_betas_promoter > 0.05])
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "CpG variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size

# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)

# Identify CpG names for cassette 6
cpg_sites <- names(promoter_15$colors)[promoter_15$colors == cas_num]

if (length(cpg_sites) > 20) {
  
  # Add density line
  density_data <- density(variance_of_betas_promoter[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)

  lines(density_data,lwd = 0.7, col = "red", )
  
} else {
  
  # Add density line
  density_data <- density(variance_of_betas_promoter[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)
  
  lines(density_data,lwd = 1.2, col = "red",)
  
  # Add vertical lines only under the polygon
  for (cpg in cpg_sites) {
    lines(c(rep(var(betaAdj[cpg, ]), 2)), c(0, density_data$y[which.min(abs(density_data$x - var(betaAdj[cpg, ])))]), col = "red", lwd = 0.25, lty = 1)
  } 
}


# ANALYSING LOCATION OF THE TOP VARIABLE CPGS IN CASSETTES

# Getting top 1000, 5000, 10000 and 20000 most variable CpGs
top_1000 <- names(sort(variance_of_betas_promoter, decreasing = T)[1:1000])
top_5000 <- names(sort(variance_of_betas_promoter, decreasing = T)[1:5000])
top_10000 <- names(sort(variance_of_betas_promoter, decreasing = T)[1:10000])
top_17725 <- names(sort(variance_of_betas_promoter, decreasing = T)[1:17725])

# Create a combined data frame for top 1000, 5000, 10000, and 17725
df_combined <- data.frame()

# Loop through different sets of CpGs (top 1000, 5000, etc.)
for (top_n in list(top_1000, top_5000, top_10000, top_17725)) {
  # Count frequencies for each level based on the current top CpG set
  level_counts <- table(as.factor(promoter_15$colors[top_n]))
  
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

# Define the order for the top CpG sets and the cassettes
top_cpg_order <- c("Top 1000", "Top 5000", "Top 10000", "Top 17725")

# Plot with specific order for bars and fill
ggplot(df_combined, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), 
                        y = Frequency, fill = factor(Cassette))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cassette", y = "Frequency", fill="Cassettes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_classic()


#
# ANALYSING DISTAL
#

# Example cassette file
distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")

# COMPUTE VARIANCE

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Only distal CpGs
distal_cpgs <- annoObj$illuminaID[which( ( (annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body") ) )]

variance_of_betas_distal <- apply(betaAdj[distal_cpgs,], MARGIN = 1, FUN=var)
cassettes <- sort(unique(distal_15$colors))

# Initialize a dataframe to store results
cassette_stats <- data.frame(
  cassette = character(),
  mean_var = numeric(),
  sd_var = numeric(),
  stringsAsFactors = FALSE
)

for (cassette in cassettes) {
  # Subset rows belonging to the current cassette
  cassette_rows <- names(which(distal_15$colors == cassette))
  
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

# PLOT CASSETTE VARIANCE

# Plot the results
par(mfrow=c(2,1), mar=c(3,3,3,3))

# Define cassette to plot
cas_num <- 40

# Compute density of variances in every cassette
density_data <- density(cassette_stats$mean_var)
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "Mean cassette variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size


# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)


lines(c(rep(cassette_stats[cas_num + 1,"mean_var"], 2)), c(0, density_data$y[which.min(abs(density_data$x - cassette_stats[cas_num + 1,"mean_var"]))]), col = "red", lwd = 2, lty = 1)



# Compute density for variance_of_betas_distal values > 0.05
density_data <- density(variance_of_betas_distal[variance_of_betas_distal > 0.1])
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "CpG variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size

# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)

# Identify CpG names for cassette 6
cpg_sites <- names(distal_15$colors)[distal_15$colors == cas_num]

if (length(cpg_sites) > 20) {
  
  # Add density line
  density_data <- density(variance_of_betas_distal[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)
  
  lines(density_data,lwd = 0.7, col = "red", )
  
} else {
  
  # Add density line
  density_data <- density(variance_of_betas_distal[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)
  
  lines(density_data,lwd = 1.2, col = "red",)
  
  # Add vertical lines only under the polygon
  for (cpg in cpg_sites) {
    lines(c(rep(var(betaAdj[cpg, ]), 2)), c(0, density_data$y[which.min(abs(density_data$x - var(betaAdj[cpg, ])))]), col = "red", lwd = 0.25, lty = 1)
  } 
}


# ANALYSING LOCATION OF THE TOP VARIABLE CPGS IN CASSETTES

# Getting top 1000, 5000, 10000 and 20000 most variable CpGs
top_1000 <- names(sort(variance_of_betas_distal, decreasing = T)[1:1000])
top_5000 <- names(sort(variance_of_betas_distal, decreasing = T)[1:5000])
top_10000 <- names(sort(variance_of_betas_distal, decreasing = T)[1:10000])
top_20000 <- names(sort(variance_of_betas_distal, decreasing = T)[1:20000])
top_34481 <- names(sort(variance_of_betas_distal, decreasing = T)[1:34481])

# Create a combined data frame for top 1000, 5000, 10000, and 17725
df_combined <- data.frame()

# Loop through different sets of CpGs (top 1000, 5000, etc.)
for (top_n in list(top_1000, top_5000, top_10000, top_20000, top_34481)) {
  # Count frequencies for each level based on the current top CpG set
  level_counts <- table(as.factor(distal_15$colors[top_n]))
  
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

# Define the order for the top CpG sets and the cassettes
top_cpg_order <- c("Top 1000", "Top 5000", "Top 10000", "Top 20000","Top 34481")

# Plot with specific order for bars and fill
ggplot(df_combined, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), 
                        y = Frequency, fill = factor(Cassette))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cassette", y = "Frequency", fill="Cassettes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_classic()


#
# ANALYSING PROXIMAL
#


# Example cassette file
proximal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_5.rds")

# COMPUTE VARIANCE

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Only proximal CpGs
proximal_cpgs <- annoObj$illuminaID[which( ( (annoObj$featureClass=="proximal up") | (annoObj$featureClass=="proximal dn") ) )]

variance_of_betas_proximal <- apply(betaAdj[proximal_cpgs,], MARGIN = 1, FUN=var)
cassettes <- sort(unique(proximal_15$colors))

# Initialize a dataframe to store results
cassette_stats <- data.frame(
  cassette = character(),
  mean_var = numeric(),
  sd_var = numeric(),
  stringsAsFactors = FALSE
)

for (cassette in cassettes) {
  # Subset rows belonging to the current cassette
  cassette_rows <- names(which(proximal_15$colors == cassette))
  
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


# PLOT CASSETTE VARIANCE

# Plot the results
par(mfrow=c(2,1), mar=c(3,3,3,3))

# Define cassette to plot
cas_num <- 6

# Compute density of variances in every cassette
density_data <- density(cassette_stats$mean_var)
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "Mean cassette variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size


# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)


lines(c(rep(cassette_stats[cas_num + 1,"mean_var"], 2)), c(0, density_data$y[which.min(abs(density_data$x - cassette_stats[cas_num + 1,"mean_var"]))]), col = "red", lwd = 2, lty = 1)



# Compute density for variance_of_betas_proximal values > 0.05
density_data <- density(variance_of_betas_proximal[variance_of_betas_proximal > 0.05])
density_data$y = density_data$y/max(density_data$y)

# Plot density curve with improved aesthetics
plot(density_data, 
     main = "CpG variance", 
     xlab = "Variance of Betas", 
     ylab = "Density", 
     lwd = 2, 
     col = "blue", 
     bty = "n",  # Remove box around the plot
     las = 1,    # Rotate axis labels for readability
     cex.axis = 1, # Increase axis label size
     cex.lab = 1,  # Increase axis title size
     cex.main = 1)  # Increase title size

# Fill area under the density curve with a semi-transparent blue
polygon(density_data, col = rgb(0, 0, 1, 0.3), border = NA)

# Identify CpG names for cassette 6
cpg_sites <- names(proximal_15$colors)[proximal_15$colors == cas_num]

if (length(cpg_sites) > 20) {
  
  # Add density line
  density_data <- density(variance_of_betas_proximal[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)
  
  lines(density_data,lwd = 0.7, col = "red", )
  
} else {
  
  # Add density line
  density_data <- density(variance_of_betas_proximal[cpg_sites])
  density_data$y = density_data$y/max(density_data$y)
  
  lines(density_data,lwd = 1.2, col = "red",)
  
  # Add vertical lines only under the polygon
  for (cpg in cpg_sites) {
    lines(c(rep(var(betaAdj[cpg, ]), 2)), c(0, density_data$y[which.min(abs(density_data$x - var(betaAdj[cpg, ])))]), col = "red", lwd = 0.25, lty = 1)
  } 
}


# ANALYSING LOCATION OF THE TOP VARIABLE CPGS IN CASSETTES

# Getting top 1000, 5000, 10000 and 20000 most variable CpGs
top_1000 <- names(sort(variance_of_betas_proximal, decreasing = T)[1:1000])
top_5000 <- names(sort(variance_of_betas_proximal, decreasing = T)[1:5000])
top_10000 <- names(sort(variance_of_betas_proximal, decreasing = T)[1:10000])
top_20000 <- names(sort(variance_of_betas_proximal, decreasing = T)[1:20000])
top_33536 <- names(sort(variance_of_betas_proximal, decreasing = T)[1:33536])

# Create a combined data frame for top 1000, 5000, 10000, and 17725
df_combined <- data.frame()

# Loop through different sets of CpGs (top 1000, 5000, etc.)
for (top_n in list(top_1000, top_5000, top_10000, top_20000, top_33536)) {
  # Count frequencies for each level based on the current top CpG set
  level_counts <- table(as.factor(proximal_15$colors[top_n]))
  
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

# Define the order for the top CpG sets and the cassettes
top_cpg_order <- c("Top 1000", "Top 5000", "Top 10000", "Top 20000","Top 33536")

# Plot with specific order for bars and fill
ggplot(df_combined, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), 
                        y = Frequency, fill = factor(Cassette))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cassette", y = "Frequency", fill="Cassettes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_classic()


#
# ANALYSING ALL CPGS REGARDLESS OF CONTEXT
#

# Total beta variance
variance_of_betas <- apply(betaAdj, MARGIN = 1, FUN=var)

# Cassettes labels
proximal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_8.rds")$colors
promoter_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_8.rds")$colors
distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_8.rds")$colors

# Merging labels
cpg_labels <- c(sapply(proximal_5, function(x) {paste0("proximal_", x)}),
                sapply(distal_5, function(x) {paste0("distal_", x)}),
                sapply(promoter_5, function(x) {paste0("promoter_", x)}))

# Getting most variant CpGs
top_1000 <- names(sort(variance_of_betas, decreasing = T)[1:1000])
top_5000 <- names(sort(variance_of_betas, decreasing = T)[1:5000])
top_10000 <- names(sort(variance_of_betas, decreasing = T)[1:10000])
top_20000 <- names(sort(variance_of_betas, decreasing = T)[1:20000])
top_30000 <- names(sort(variance_of_betas, decreasing = T)[1:30000])
top_50000 <- names(sort(variance_of_betas, decreasing = T)[1:50000])


# Create a combined data frame for top 1000, 5000, 10000, and 17725
df_combined <- data.frame()

# Loop through different sets of CpGs (top 1000, 5000, etc.)
for (top_n in list(top_1000, top_5000, top_10000, top_20000, top_30000, top_50000)) {
  # Count frequencies for each level based on the current top CpG set
  level_counts <- table(as.factor(cpg_labels[top_n]))
  
  # Identify levels with frequency less than 0.5 %
  infrequent_levels <- names(level_counts[level_counts < length(top_n) * 0.01])
  others <- sum(level_counts[infrequent_levels])
  
  # Combine infrequent levels into 'Others'
  level_counts <- level_counts[level_counts >= length(top_n) * 0.015]
  level_counts["Others"] <- others
  
  # Add data for current top_n to the combined data frame
  df_temp <- data.frame(
    Cassette = names(level_counts),
    Frequency = as.vector(level_counts),
    Top_CPG_Set = paste("Top", length(top_n))
  )
  
  df_combined <- rbind(df_combined, df_temp)
}

# Define the order for the top CpG sets and the cassettes
top_cpg_order <- c("Top 1000", "Top 5000", "Top 10000", "Top 20000","Top 30000","Top 50000")

# Plot with specific order for bars and fill
ggplot(df_combined, aes(x = factor(Top_CPG_Set, levels = top_cpg_order), 
                        y = Frequency, fill = factor(Cassette))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cassette", y = "Frequency", fill="Cassettes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_classic()


Heatmap(betaAdj[names(cpg_labels[top_5000][!cpg_labels[top_5000] %in% c("distal_0", "promoter_0", "proximal_0", "distal_1")]),], 
        use_raster=F)

Heatmap(betaAdj[top_5000 ,], 
        use_raster=F,
        )

