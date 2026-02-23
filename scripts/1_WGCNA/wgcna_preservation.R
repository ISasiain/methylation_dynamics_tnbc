library(WGCNA)
library(patchwork)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
#
# LOADING DATA
#

# Loading dicovery betas
#load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

load("/Volumes/Data/Project_3/validation_cohort/Combined_annotations_rel4SCANB_deNovoMainNMF_distalAtac5000_supervisedSubNMF.RData")
load("/Volumes/Data/Project_3/validation_cohort/FPKM_rel4TNBC_validationCohort_n136.RData")
load("/Volumes/Data/Project_3/validation_cohort/PurBeta_adjustedTumor_betaMatrix_V1_V2_reduced_717459commonCpGs_TNBCs_n136.RData")

# Defining gene-cpg dictionary
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")

# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID


# Object name: SCANBrel4_rdata
load("/Users/in2245sa/PhD/Projects/project_3/data/SCANBrel4valcohort_annotations.RData") 
rownames(SCANBrel4_rdata) <- SCANBrel4_rdata$id

# Reading cassettes detected in discovery cohort
promoter_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_7.rds")
distal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_7.rds")
proximal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_7.rds")

atac_distal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/atac_cassettes_beta_7.rds")

#
# PERFORMING MODULE PRESERVATION ANALYSIS
#

# PROMOTER

# Getting shared CpGs for module preservation
shared_cpgs <- rownames(beta.adjusted)[rownames(beta.adjusted) %in% rownames(betaAdj)]
shared_cpgs_promoter <- shared_cpgs[shared_cpgs %in% names(promoter_7$colors)]

# Prepare input data for module preservation

# Subset and transpose (samples x CpGs)
disc <- t(betaAdj[shared_cpgs_promoter, ])
val  <- t(beta.adjusted[shared_cpgs_promoter, ])

# M-value transform
eps <- 0.01
disc <- pmin(pmax(disc, eps), 1 - eps)
val  <- pmin(pmax(val,  eps), 1 - eps)

disc <- log2(disc / (1 - disc))
val  <- log2(val / (1 - val))

# Set a tiny threshold
var_cutoff <- 1e-6

disc <- disc[, apply(disc, 2, var) > var_cutoff]
val  <- val[,  apply(val,  2, var) > var_cutoff]

# Re-align columns
common_cpgs <- intersect(colnames(disc), colnames(val))
disc <- disc[, common_cpgs]
val  <- val[,  common_cpgs]

multiColor <- list(
  Discovery = promoter_7$colors[common_cpgs]
)

multiData <- list(
  Discovery = list(data = disc),  
  Validation = list(data = val) 
)

preservation_results <- modulePreservation(
  multiData = multiData,         
  multiColor = multiColor,
  networkType = "unsigned",
  referenceNetworks = 1,         
  nPermutations = 200,      
  randomSeed = 12345,             
  verbose = 3,                      
)

saveRDS(preservation_results, file = "/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_promoter_7.rds")

# PROXIMAL

#Getting shared CpGs for module preservation
shared_cpgs <- rownames(beta.adjusted)[rownames(beta.adjusted) %in% rownames(betaAdj)]
shared_cpgs_proximal <- shared_cpgs[shared_cpgs %in% names(proximal_7$colors)]

# Prepare input data for module preservation

# Subset and transpose (samples x CpGs)
disc <- t(betaAdj[shared_cpgs_proximal, ])
val  <- t(beta.adjusted[shared_cpgs_proximal, ])

# M-value transform
eps <- 0.01
disc <- pmin(pmax(disc, eps), 1 - eps)
val  <- pmin(pmax(val,  eps), 1 - eps)

disc <- log2(disc / (1 - disc))
val  <- log2(val / (1 - val))

# Set a tiny threshold
var_cutoff <- 1e-6

disc <- disc[, apply(disc, 2, var) > var_cutoff]
val  <- val[,  apply(val,  2, var) > var_cutoff]

# Re-align columns
common_cpgs <- intersect(colnames(disc), colnames(val))
disc <- disc[, common_cpgs]
val  <- val[,  common_cpgs]

multiColor <- list(
  Discovery = proximal_7$colors[common_cpgs]
)

multiData <- list(
  Discovery = list(data = disc),  
  Validation = list(data = val) 
)

preservation_results <- modulePreservation(
  multiData = multiData,         
  multiColor = multiColor,
  networkType = "unsigned",
  referenceNetworks = 1,         
  nPermutations = 200,      
  randomSeed = 12345,             
  verbose = 3                      
)

saveRDS(preservation_results, file = "/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_proximal_7.rds")


# DISTAL

#Getting shared CpGs for module preservation
shared_cpgs <- rownames(beta.adjusted)[rownames(beta.adjusted) %in% rownames(betaAdj)]
shared_cpgs_distal <- shared_cpgs[shared_cpgs %in% names(distal_7$colors)]

# Prepare input data for module preservation

# Subset and transpose (samples x CpGs)
disc <- t(betaAdj[shared_cpgs_distal, ])
val  <- t(beta.adjusted[shared_cpgs_distal, ])

# M-value transform
eps <- 0.01
disc <- pmin(pmax(disc, eps), 1 - eps)
val  <- pmin(pmax(val,  eps), 1 - eps)

disc <- log2(disc / (1 - disc))
val  <- log2(val / (1 - val))

# Set a tiny threshold
var_cutoff <- 1e-6

disc <- disc[, apply(disc, 2, var) > var_cutoff]
val  <- val[,  apply(val,  2, var) > var_cutoff]

# Re-align columns
common_cpgs <- intersect(colnames(disc), colnames(val))
disc <- disc[, common_cpgs]
val  <- val[,  common_cpgs]

multiColor <- list(
  Discovery = distal_7$colors[common_cpgs]
)

multiData <- list(
  Discovery = list(data = disc),  
  Validation = list(data = val) 
)

preservation_results <- modulePreservation(
  multiData = multiData,         
  multiColor = multiColor,
  networkType = "unsigned",
  referenceNetworks = 1,         
  nPermutations = 200,      
  randomSeed = 12345,             
  verbose = 3                      
)

saveRDS(preservation_results, file = "/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_distal_7.rds")


#
# PLOTTING PRESERVATION METRICS
#

preservation_proximal <- readRDS("/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_proximal_7.rds")
preservation_distal <- readRDS("/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_distal_7.rds")
preservation_promoter <- readRDS("/Users/in2245sa/PhD/Projects/project_3/analysis/cassette_preservation/preservation_promoter_7.rds")

# PROMOTER

# Preservation stats
pres <- preservation_promoter$preservation

# Zsummary
Zsummary <- pres$Z$ref.Discovery$inColumnsAlsoPresentIn.Validation$Zsummary.pres
p_adjust.zsum <- pres$log.pBonf$ref.Discovery$inColumnsAlsoPresentIn.Validation$log.p.Bonfsummary.pres

df_promoter <- data.frame(
  module = rownames(pres$observed$ref.Discovery$inColumnsAlsoPresentIn.Validation),
  Zsummary = Zsummary,
  pZsummary = p_adjust.zsum,
  Context = "Promoter"
)

# Remove uncalssified and gold modules
df_promoter <- df_promoter[!df_promoter$module %in% c("0", "0.1"),]

# PROXIMAL

# Preservation stats
pres <- preservation_proximal$preservation

# Zsummary
Zsummary <- pres$Z$ref.Discovery$inColumnsAlsoPresentIn.Validation$Zsummary.pres
p_adjust.zsum <- pres$log.pBonf$ref.Discovery$inColumnsAlsoPresentIn.Validation$log.p.Bonfsummary.pres

df_proximal <- data.frame(
  module = rownames(pres$observed$ref.Discovery$inColumnsAlsoPresentIn.Validation),
  Zsummary = Zsummary,
  pZsummary = p_adjust.zsum,
  Context = "Proximal"
)

# Remove uncalssified and gold modules
df_proximal <- df_proximal[!df_proximal$module %in% c("0", "0.1"),]

# DISTAL

# Preservation stats
pres <- preservation_distal$preservation

# Zsummary
Zsummary <- pres$Z$ref.Discovery$inColumnsAlsoPresentIn.Validation$Zsummary.pres
p_adjust.zsum <- pres$log.pBonf$ref.Discovery$inColumnsAlsoPresentIn.Validation$log.p.Bonfsummary.pres

df_distal <- data.frame(
  module = rownames(pres$observed$ref.Discovery$inColumnsAlsoPresentIn.Validation),
  Zsummary = Zsummary,
  pZsummary = p_adjust.zsum,
  Context = "Distal"
)

# Remove uncalssified and gold modules
df_distal <- df_distal[!df_distal$module %in% c("0", "0.1"),]

# Merge and plor
df_preservation_all <- rbind(df_promoter, df_proximal, df_distal)

ggplot(df_preservation_all, aes(x = Context, y = Zsummary)) +
  geom_beeswarm(size = 1, aes(colour = ifelse(10^(pZsummary) < 0.05, "signif", "ns"))) +
  geom_hline(yintercept = 2, linetype = "dashed") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  geom_text_repel(
    aes(label = module),
    size = 3,
    max.overlaps = 20
  ) +
  scale_colour_manual(
    name = "p value",
    values = c("signif" = "black", "ns" = "red")
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Context",
    y = "Zsummary"
  )
