#! usr/bin/Rscript

library(WGCNA)


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

# Define number for resampling
samples_n <- c(60, 100, 140)

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
    
    # Calculating variance of CpGs
    variance_dis <- sapply(1:nrow(sampled_betas), FUN = function(row) {var(sampled_betas[row,])})
    
    # Filtering CpGs with low variance
    sampled_betas <-  t(sampled_betas[variance_dis > 0.1,])
    
    # Defining betas to use 
    betas <- c(5, 10, 15)
    
    # Running WGCNA
    for (beta in betas) {
      
      netwk <- blockwiseModules(sampled_betas,               
                                corrType="bicor", # Using biweight midcorrelation 
                                nThreads = 10,
                                
                                # == Adjacency Function ==
                                power = beta,             
                                networkType = "unsigned", #Choosing unsigned to account for positive and negative correlations
                                
                                # == Tree and Block Options ==
                                deepSplit = 2,
                                pamRespectsDendro = F,
                                # detectCutHeight = 0.75,
                                minModuleSize = 3,
                                maxBlockSize = 6000,
                                
                                # == Module Adjustments ==
                                reassignThreshold = 0,
                                mergeCutHeight = 0.25,
                                
                                # == Output Options
                                numericLabels = T,
                                verbose = 3)
      
      # Saving network
      my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/cassettes_resampling/all_distal/", "n_", sample_n, "_beta_", beta, "_resampling_", run, ".rds" )
      saveRDS(netwk, file = my_filename)
      
    }
  }
}
