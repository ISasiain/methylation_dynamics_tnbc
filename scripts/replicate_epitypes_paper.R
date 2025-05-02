#! usr/bin/Rscript

library(ComplexHeatmap)
library(NMF)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(WGCNA)

#
# LOADING DATA
#

# Loading corrected betas and annotations
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load epitype annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID
orig_epi <- read.table("/Users/isasiain//PhD/Projects/immune_spatial/ecosystem_analysis/data/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt", sep = "\t")


# Loading CpG cassettes (Beta = 5)
distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5.rds")
nonbasal_distal_5_all <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_nonBasal_cassettes_beta_5.rds")
basal_distal_5_all <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_basal_cassettes_beta_5.rds")

distal_5_atac <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5_only_atac.rds")
nonbasal_distal_5_atac <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_nonBasal_cassettes_beta_5_only_atac.rds")
basal_distal_5_atac <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_basal_cassettes_beta_5_only_atac.rds")


# Filtering Non Distal CpGs
my_cpgs <- annoObj$illuminaID[which( ((annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body")))]
my_betas <- betaAdj[rownames(betaAdj) %in% my_cpgs, ]

# Filtering Non Distal CpGs and Cassettes 0 and 1
my_cpgs_filtered <- my_cpgs[!distal_5$colors[my_cpgs] %in% c(0, 1)]
my_betas_filtered <- betaAdj[rownames(betaAdj) %in% my_cpgs_filtered, ]

# Filtering Non Distal CpGs and Cassettes not enriched for TFs
my_cpgs_TFfiltered <- my_cpgs[distal_5$colors[my_cpgs] %in% c(2, 3, 6)]
my_betas_TFfiltered <- betaAdj[rownames(betaAdj) %in% my_cpgs_TFfiltered, ]

# Filtering Non Distal Non ATAC CpGs
my_cpgs_atac <- distal_cpgs <- annoObj$illuminaID[which((annoObj$hasAtacOverlap == 1) & 
                                                          (annoObj$featureClass == "distal" | annoObj$featureClass == "distal body"))]
my_betas_atac <- betaAdj[rownames(betaAdj) %in% my_cpgs_atac, ]

# 
# NMF CLUSTERING 1
#

# Filtering data based on variance

  # ALL

# Calculating variance of CpGs
variance_vector <- apply(X=my_betas, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_all <- betaAdj[top_variance_cpgs,]

  # FILTERING 0 and 1 CASSETTES

# Calculating variance of CpGs
variance_vector <- apply(X=my_betas_filtered, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]

# Getting beta matrix
betas_to_analyse1_filtered <- betaAdj[top_variance_cpgs,]


  # FILTERING NON TF ENRICHED CASSETTES

# Calculating variance of CpGs
variance_vector <- apply(X=my_betas_TFfiltered, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]

# Getting beta matrix
betas_to_analyse1_TFfiltered <- betaAdj[top_variance_cpgs,]


  # ONLY ATAC

# Calculating variance of CpGs
variance_vector <- apply(X=my_betas_atac, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. 
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]

# Getting beta matrix
betas_to_analyse1_atac <- betaAdj[top_variance_cpgs,]

# Running NMF clustering

# 5000 most variable distal
r.nmf_all <- nmf(betas_to_analyse1_all,
             rank = 2,
             nrun = 100,
             method = 'brunet',
             seed = 221027)

# 5000 most variable distal minus cassettes 0 and 1
r.nmf_filtered <-  nmf(betas_to_analyse1_filtered,
                       rank = 2,
                       nrun = 100,
                       method = 'brunet',
                       seed = 221027)

# 5000 most variable distal in TF enriched cassettes
r.nmf_TFfiltered <-  nmf(betas_to_analyse1_TFfiltered,
                       rank = 2,
                       nrun = 100,
                       method = 'brunet',
                       seed = 221027)

# 5000 most variable atac distal
r.nmf_atac <-  nmf(betas_to_analyse1_atac,
                   rank = 2,
                   nrun = 100,
                   method = 'brunet',
                   seed = 221027)

# Getting clusters based on NMF
clusters.1.all <- predict(r.nmf_all, "chc")
clusters.1.filtered <- predict(r.nmf_filtered, "chc")
clusters.1.TFfiltered <- predict(r.nmf_TFfiltered, "chc")
clusters.1.atac <- predict(r.nmf_atac, "chc")

# Heatmap of clusters and CpG betas
subtypes_1 <- data.frame(
  
  "orig_epi"=as.factor(orig_epi[names(clusters.1.all), "NMF_atacDistal"]),
  "PAM50"= as.factor(orig_epi[names(clusters.1.all), "PAM50"]),
  "new_epi_all"=as.factor(clusters.1.all),
  "new_epi_filtered"=as.factor(clusters.1.filtered),
  "new_epi_TFfiltered"=as.factor(clusters.1.TFfiltered),
  "new_epi_atac"= as.factor(clusters.1.atac)
  
)

# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters.1.TFfiltered

original_epitypes <- orig_epi[names(new_epitypes), "NMF_atacDistal"]
original_epitypes <- substr(original_epitypes, 1, nchar(original_epitypes) - 1)
original_epitypes <- paste0(original_epitypes, "_epi")

pam50 <- x[names(new_epitypes), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(new_epitypes),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")


#
# WGCNA FOR DETECTED NMF MAIN SUBTYPES
#

# Assign the detected subtypes as Basal and nonBasal. Highest number -> Basal
for (nmf_sub_1 in colnames(subtypes_1)) {
  
  if(!nmf_sub_1 %in% c("orig_epi", "PAM50")) {
    
    # Getting frequency of assignments
    table_of_clusters <- table(subtypes_1[, nmf_sub_1])
    
    # Renaming dataframe
    subtypes_1[, nmf_sub_1] <- sapply(subtypes_1[, nmf_sub_1], function(subtype) {
      
      if (table_of_clusters[subtype] >= length(subtypes_1[, nmf_sub_1]) / 2) {
        "Basal"
      } else {"nonBasal"}
      
    })
    
  }
  
}



# Looping over main subtype assignments and running WGCNA for basal and nonBasal samples
for (nmf_sub_1 in colnames(subtypes_1)) {
  
  if(!nmf_sub_1 %in% c("orig_epi", "PAM50")) {
    
    # Running in basal and saving
    betas <- c(5, 10, 15)
    
    # Iterate through betas
    for (beta in betas) {
      
      if (nmf_sub_1 == "new_epi_atac") {
        
        dis_to_analyse <- my_betas_atac[,rownames(subtypes_1)[subtypes_1[, nmf_sub_1] == "Basal"]]
        
      } else {
        
        dis_to_analyse <- my_betas[,rownames(subtypes_1)[subtypes_1[, nmf_sub_1] == "Basal"]]
        
      }
      
      # Running WGCNA
      netwk_basal <- blockwiseModules(dis_to_analyse,               
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
                                      
                                      # == TOM == Archive the run results in TOM file (saves time)
                                      saveTOMs = F,
                                      saveTOMFileBase = "ER",
                                      
                                      # == Output Options
                                      numericLabels = T,
                                      verbose = 3)
      
      # Saving
      my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/wgcna_nmf_subtypes/Basal_cassettes_beta_", beta, "_", nmf_sub_1, ".rds" )
      saveRDS(netwk_basal, file = my_filename)
      
      
    }
    
    
    # Running in nonBasal and saving
    # Running in basal and saving
    betas <- c(5, 10, 15)
    
    # Iterate through betas
    for (beta in betas) {
      
      if (nmf_sub_1 == "new_epi_atac") {
        
        dis_to_analyse <- my_betas_atac[,rownames(subtypes_1)[subtypes_1[, nmf_sub_1] == "nonBasal"]]
        
      } else {
        
        dis_to_analyse <- my_betas[,rownames(subtypes_1)[subtypes_1[, nmf_sub_1] == "nonBasal"]]
        
      }
      
      
      # Running WGCNA
      netwk_nonbasal <- blockwiseModules(dis_to_analyse,               
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
                                      
                                      # == TOM == Archive the run results in TOM file (saves time)
                                      saveTOMs = F,
                                      saveTOMFileBase = "ER",
                                      
                                      # == Output Options
                                      numericLabels = T,
                                      verbose = 3)
      
      # Saving
      my_filename <- paste0("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/wgcna_nmf_subtypes/nonBasal_cassettes_beta_", beta, "_", nmf_sub_1, ".rds" )
      saveRDS(netwk_nonbasal, file = my_filename)
      
      
    }
  }
  
}



# 
# NMF CLUSTERING 2.1 BASAL
#


# Filtering Non Distal CpGs
my_cpgs <- annoObj$illuminaID[which( ((annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body")))]
my_betas <- betaAdj[rownames(betaAdj) %in% my_cpgs, ]

# Filtering Non Distal CpGs and Cassettes 0 and 1
my_cpgs_filtered <- my_cpgs[!basal_distal_5$colors[my_cpgs] %in% c(0, 1)]
my_betas_filtered <- betaAdj[rownames(betaAdj) %in% my_cpgs_filtered, ]

# Filtering Non Distal CpGs and Cassettes not enriched fo TFs
my_cpgs_TFfiltered <- my_cpgs[basal_distal_5$colors[my_cpgs] %in% c(3, 5, 7, 8, 9, 10, 13)]
my_betas_TFfiltered <- betaAdj[rownames(betaAdj) %in% my_cpgs_TFfiltered, ]

# Filtering Non Distal Non ATAC CpGs
my_cpgs_atac <- distal_cpgs <- annoObj$illuminaID[which((annoObj$hasAtacOverlap == 1) & 
                                                          (annoObj$featureClass == "distal" | annoObj$featureClass == "distal body"))]
my_betas_atac <- betaAdj[rownames(betaAdj) %in% my_cpgs_atac, ]


#
# Filtering based on cassettes from all samples
#


## ALL
distal_all_Basal <- my_betas[,names(clusters.1.all)[clusters.1.all == 1]]

# Calculating variance of CpGs
variance_vector_Basal <- apply(X=distal_all_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_Basal <- names(sort(variance_vector_Basal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_Basal_all <- distal_all_Basal[top_variance_cpgs_Basal,]


## FILTERING CASSETTE 0 AND 1
distal_filtered_Basal <- my_betas_filtered[,names(clusters.1.filtered)[clusters.1.filtered == 1]]

# Calculating variance of CpGs
variance_vector_Basal <- apply(X=distal_filtered_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_Basal <- names(sort(variance_vector_Basal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_Basal_filtered <- distal_filtered_Basal[top_variance_cpgs_Basal,]


## FILTERING NON TF ENRICHED CASSETTES
distal_TFfiltered_Basal <- my_betas_atac[,names(clusters.1.TFfiltered)[clusters.1.TFfiltered == 2]]

# Calculating variance of CpGs
variance_vector_Basal <- apply(X=distal_TFfiltered_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_Basal <- names(sort(variance_vector_Basal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_Basal_TFfiltered <- distal_TFfiltered_Basal[top_variance_cpgs_Basal,]


## ONLY ATAC
distal_atac_Basal <- my_betas_atac[,names(clusters.1.atac)[clusters.1.atac == 2]]

# Calculating variance of CpGs
variance_vector_Basal <- apply(X=distal_atac_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_Basal <- names(sort(variance_vector_Basal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_Basal_atac <- distal_atac_Basal[top_variance_cpgs_Basal,]




# Running NMF

# 5000 most variable distal Basal
r.nmf_all_basal <- nmf(betas_to_analyse1_Basal_all,
                 rank = 3,
                 nrun = 100,
                 method = 'brunet',
                 seed = 221027)

# 5000 most variable distal minus cassettes 0 and 1
r.nmf_filtered_basal <-  nmf(betas_to_analyse1_Basal_filtered,
                       rank = 3,
                       nrun = 100,
                       method = 'brunet',
                       seed = 221027)

# 5000 most variable distal in TF enriched cassettes
r.nmf_TFfiltered_basal <-  nmf(betas_to_analyse1_Basal_TFfiltered,
                         rank = 3,
                         nrun = 100,
                         method = 'brunet',
                         seed = 221027)

# 5000 most variable atac distal
r.nmf_atac_basal <-  nmf(betas_to_analyse1_Basal_atac,
                   rank = 3,
                   nrun = 100,
                   method = 'brunet',
                   seed = 221027)


# Getting clusters based on NMF
clusters.1.all_basal <- predict(r.nmf_all_basal, "chc")
clusters.1.filtered_basal <- predict(r.nmf_filtered_basal, "chc")
clusters.1.TFfiltered_basal <- predict(r.nmf_TFfiltered_basal, "chc")
clusters.1.atac_basal <- predict(r.nmf_atac_basal, "chc")


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters.1.atac_basal

original_epitypes <- orig_epi[names(new_epitypes), "NMF_atacDistal"]

pam50 <- x[names(new_epitypes), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(new_epitypes),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")




# 
# NMF CLUSTERING 2.1 NON-BASAL
#

## ALL
distal_all_NonBasal <- my_betas[,names(clusters.1.all)[clusters.1.all == 2]]

# Calculating variance of CpGs
variance_vector_NonBasal <- apply(X=distal_all_NonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_NonBasal <- names(sort(variance_vector_NonBasal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_NonBasal_all <- distal_all_NonBasal[top_variance_cpgs_NonBasal,]


## FILTERING CASSETTE 0 AND 1
distal_filtered_NonBasal <- my_betas_filtered[,names(clusters.1.filtered)[clusters.1.filtered == 2]]

# Calculating variance of CpGs
variance_vector_NonBasal <- apply(X=distal_filtered_NonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_NonBasal <- names(sort(variance_vector_NonBasal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_NonBasal_filtered <- distal_filtered_NonBasal[top_variance_cpgs_NonBasal,]


## FILTERING NON TF ENRICHED CASSETTES
distal_TFfiltered_NonBasal <- my_betas_atac[,names(clusters.1.TFfiltered)[clusters.1.TFfiltered == 1]]

# Calculating variance of CpGs
variance_vector_NonBasal <- apply(X=distal_TFfiltered_NonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_NonBasal <- names(sort(variance_vector_NonBasal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_NonBasal_TFfiltered <- distal_TFfiltered_NonBasal[top_variance_cpgs_NonBasal,]


## ONLY ATAC
distal_atac_NonBasal <- my_betas_atac[,names(clusters.1.atac)[clusters.1.atac == 1]]

# Calculating variance of CpGs
variance_vector_NonBasal <- apply(X=distal_atac_NonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance
top_variance_cpgs_NonBasal <- names(sort(variance_vector_NonBasal, decreasing = T))[1:5000]

# Getting beta matriX
betas_to_analyse1_NonBasal_atac <- distal_atac_NonBasal[top_variance_cpgs_NonBasal,]




# Running NMF

# 5000 most variable distal NonBasal
r.nmf_all_NonBasal <- nmf(betas_to_analyse1_NonBasal_all,
                       rank = 2,
                       nrun = 100,
                       method = 'brunet',
                       seed = 221027)

# 5000 most variable distal minus cassettes 0 and 1
r.nmf_filtered_NonBasal <-  nmf(betas_to_analyse1_NonBasal_filtered,
                             rank = 2,
                             nrun = 100,
                             method = 'brunet',
                             seed = 221027)

# 5000 most variable distal in TF enriched cassettes
r.nmf_TFfiltered_NonBasal <-  nmf(betas_to_analyse1_NonBasal_TFfiltered,
                               rank = 2,
                               nrun = 100,
                               method = 'brunet',
                               seed = 221027)

# 5000 most variable atac distal
r.nmf_atac_NonBasal <-  nmf(betas_to_analyse1_NonBasal_atac,
                         rank = 2,
                         nrun = 100,
                         method = 'brunet',
                         seed = 221027)


# Getting clusters based on NMF
clusters.1.all_NonBasal <- predict(r.nmf_all_NonBasal, "chc")
clusters.1.filtered_NonBasal <- predict(r.nmf_filtered_NonBasal, "chc")
clusters.1.TFfiltered_NonBasal <- predict(r.nmf_TFfiltered_NonBasal, "chc")
clusters.1.atac_NonBasal <- predict(r.nmf_atac_NonBasal, "chc")


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters.1.TFfiltered_NonBasal

original_epitypes <- orig_epi[names(new_epitypes), "NMF_atacDistal"]

pam50 <- x[names(new_epitypes), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(new_epitypes),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")

