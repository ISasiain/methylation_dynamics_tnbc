#! usr/bin/Rscript

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

#
# USER DEFINED FUNCTIONS
#

# Function to calculate CpG density
calculate_cpg_density <- function(cpg_data, sequence_length=1000) {
  # Sort data by Chromosome ID and Position to ensure correct ordering
  cpg_data <- cpg_data[order(cpg_data$chromosome, cpg_data$position), ]
  
  # Initialize a vector to store the CpG density for each CpG site
  cpg_density <- numeric(nrow(cpg_data))
  
  # Loop through each CpG to calculate CpG density in the surrounding sequence length
  for (i in 1:nrow(cpg_data)) {
    # Get the current CpG position and chromosome
    chrom <- cpg_data$chromosome[i]
    position <- cpg_data$position[i]
    
    # Define the window for the current CpG (sequence_length before and after)
    window_start <- position - sequence_length / 2
    window_end <- position + sequence_length / 2
    
    # Subset data to get all CpGs in the same chromosome within the window range
    window_cpgs <- subset(cpg_data, chromosome == chrom & 
                            position >= window_start & position <= window_end)
    
    # Calculate the density as the number of CpGs in the window divided by the sequence length
    cpg_density[i] <- nrow(window_cpgs) / sequence_length
  }
  
  # Add CpG density as a new column to the data
  cpg_data$cpg_density <- cpg_density
  
  return(cpg_data)
}

# Formatting data to be used
cpg_data <- data.frame(
  chromosome = annoObj$chr,
  position = annoObj$start,
  cpg_id = annoObj$illuminaID
)

# Calculating CpG density taking different sequence lengths
lengths <- c(500, 1000, 2000)

for (l in lengths) {
  
  # Call the function with sequence_length of 200
  result <- calculate_cpg_density(cpg_data, sequence_length = l)
  
  filename <- paste0("./PhD/Projects/project_3/analysis/cpg_desnity_", l, "_bp.rds")
  
  saveRDS(result, file="./PhD/Projects/project_3/analysis/cpg_desnity.rds")
}



