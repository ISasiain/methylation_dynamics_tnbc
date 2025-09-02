#! usr/bin/Rscript

#
# LOADING DATA
#

# CpG cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")
all_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/all/cassettes_beta_10.rds")


# Loading EPIC methylation matrix and annotations
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

rownames(annoObj) <- annoObj$illuminaID

#
# CONVERT TO CSV
#

# Convert to data frame
df_cassettes_promoter <- data.frame(
  CpG = names(promoter_10$colors),
  Cassette = as.factor(promoter_10$colors),
  row.names = NULL
)

df_cassettes_proximal <- data.frame(
  CpG = names(proximal_10$colors),
  Cassette = as.factor(proximal_10$colors),
  row.names = NULL
)

df_cassettes_distal <- data.frame(
  CpG = names(distal_10$colors),
  Cassette = as.factor(distal_10$colors),
  row.names = NULL
)

df_cassettes_all <- data.frame(
  CpG = names(all_10$colors),
  Cassette = as.factor(all_10$colors),
  row.names = NULL
)

# Adding CpG context to data frame
df_cassettes_promoter$Context <- annoObj[df_cassettes_promoter$CpG, "featureClass"]
df_cassettes_proximal$Context <- annoObj[df_cassettes_proximal$CpG, "featureClass"]
df_cassettes_distal$Context <- annoObj[df_cassettes_distal$CpG, "featureClass"]
df_cassettes_all$Context <- annoObj[df_cassettes_all$CpG, "featureClass"]


# Adding adjusted beta values to dataframe
df_cassettes_promoter <- cbind(df_cassettes_promoter,betaAdj[df_cassettes_promoter$CpG,])
df_cassettes_proximal <- cbind(df_cassettes_proximal,betaAdj[df_cassettes_proximal$CpG,])
df_cassettes_distal <- cbind(df_cassettes_distal,betaAdj[df_cassettes_distal$CpG,])
df_cassettes_all <- cbind(df_cassettes_all,betaAdj[df_cassettes_all$CpG,])


# Saving as csv files
write.csv(df_cassettes_promoter, "/Users/isasiain/PhD/Projects/project_3/data/supp_data/prom_cassettes_10.csv")
write.csv(df_cassettes_proximal, "/Users/isasiain/PhD/Projects/project_3/data/supp_data/prox_cassettes_10.csv")
write.csv(df_cassettes_distal, "/Users/isasiain/PhD/Projects/project_3/data/supp_data/dis_cassettes_10.csv")
write.csv(df_cassettes_all, "/Users/isasiain/PhD/Projects/project_3/data/supp_data/all_cassettes_10.csv")

