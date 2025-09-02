#! usr/bin/Rscript

library(dplyr)
library(ggplot2)
library(patchwork)

# Loading betas
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

# Loading promoter cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_5.rds")


#
# PERFORMING DIFFERENTIAL METHYLATION ANALYSIS: Basal/nonBasal
#

# Detect differential methylated cpgs from adjusted data
p_vals <- apply(betaAdj, MARGIN = 1,
      FUN = function(my_betas) {
        wilcox.test(my_betas ~ x[colnames(betaAdj), "PAM50_Basal_NCN"])$p.value
      })

# Adjust for mulyiple testing. Bonferroni
bonferroni_p_vals <- p.adjust(p_vals, method = "bonferroni")

# Getting CpG ids 
sign_cpgs <- names(na.omit(bonferroni_p_vals)[na.omit(bonferroni_p_vals) <= 0.05])

dis_sign_cpgs <- sign_cpgs[annoObj[sign_cpgs,"CpG_context"] == "Distal"]

# Function to group small assettes into "Others"
collapse_counts <- function(tbl, min_size = 10) {
  df <- as.data.frame(tbl)
  colnames(df) <- c("Cassette", "Count")
  df$Cassette <- as.character(df$Cassette)
  
  # Collapse small groups
  df <- df %>%
    mutate(Cassette = ifelse(Count < min_size, "Others", Cassette)) %>%
    group_by(Cassette) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  return(df)
}

# Getting distal cassettes
dist_tbl <- table(distal_10$colors[dis_sign_cpgs])
dist_df  <- collapse_counts(dist_tbl, min_size = sum(dist_tbl) / 10) %>%
  mutate(Context = "Distal")

dist_df$Cassette <- factor(dist_df$Cassette, 
                           levels = c("1","2","3","4","5","6","7","8","9","10","11","12","23","28","30","33","46","Others"))

# Barplot
ggplot(dist_df, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(x = "Distal cassettes", y = "Number of significant CpGs")

#
# PERFORMING DIFFERENTIAL METHYLATION ANALYSIS: Lehmann4
#

# Helper function to group cassettes with low number of CpGs
collapse_counts <- function(tbl, min_size = sum(tbl) / 20) {
  df <- as.data.frame(tbl)
  colnames(df) <- c("Cassette", "Count")
  df$Cassette <- as.character(df$Cassette)
  
  df <- df %>%
    mutate(Cassette = ifelse(Count < min_size, "Others", Cassette)) %>%
    group_by(Cassette) %>%
    summarise(Count = sum(Count), .groups = "drop")
  
  return(df)
}

# Function to compare lehmann groups
compare_groups <- function(g1, g2, context = "Distal") {
  
  # Select samples for the two groups
  my_samples <- colnames(betaAdj)[x[colnames(betaAdj), "TNBCtype4_n235_notPreCentered"] %in% c(g1, g2)]
  group_labels <- x[my_samples, "TNBCtype4_n235_notPreCentered"]
  
  # Differential test (Wilcoxon)
  p_vals <- apply(betaAdj[, my_samples], 1, function(my_betas) {
    wilcox.test(my_betas ~ group_labels)$p.value
  })
  
  # Multiple testing
  bonferroni_p_vals <- p.adjust(p_vals, method = "bonferroni")
  
  # Significant CpGs for given context
  sign_cpgs <- names(na.omit(bonferroni_p_vals)[na.omit(bonferroni_p_vals) <= 0.05])
  ctx_cpgs  <- sign_cpgs[annoObj[sign_cpgs, "CpG_context"] == context]
  
  # Get cassette counts
  if (context == "Promoter") {
    tbl <- table(promoter_10$colors[ctx_cpgs])
  } else if (context == "Proximal") {
    tbl <- table(proximal_10$colors[ctx_cpgs])
  } else {
    tbl <- table(distal_10$colors[ctx_cpgs])
  }
  
  df <- collapse_counts(tbl) %>%
    mutate(Context = context,
           Comparison = paste(g1, "vs", g2))
  
  return(df)
}

# Run for all pairwise comparisons
groups <- c("BL1", "BL2", "M", "LAR")
contexts <- c("Distal")

all_results <- list()

for (ctx in contexts) {
  combs <- combn(groups, 2, simplify = FALSE)
  for (pair in combs) {
    df <- compare_groups(pair[1], pair[2], context = ctx)
    all_results[[paste(pair, collapse = "_") %>% paste0("_", ctx)]] <- df
  }
}

  
p1 <- ggplot(all_results$BL1_BL2_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "BL1 vs BL2", x = "Distal cassettes", y = "Number of significant CpGs")

p2 <- ggplot(all_results$BL1_M_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "BL1 vs M", x = "Distal cassettes", y = "Number of significant CpGs")

p3 <- ggplot(all_results$BL1_LAR_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "BL1 vs LAR", x = "Distal cassettes", y = "Number of significant CpGs")

p4 <- ggplot(all_results$BL2_M_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "BL2 vs M", x = "Distal cassettes", y = "Number of significant CpGs")

p5 <- ggplot(all_results$BL2_LAR_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "BL2 vs LAR", x = "Distal cassettes", y = "Number of significant CpGs")

p6 <- ggplot(all_results$M_LAR_Distal, aes(x = Cassette, y = Count)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  labs(title = "M vs LAR", x = "Distal cassettes", y = "Number of significant CpGs")

final_plot <- (p1 | p2 | p3) /
  (p4 | p5 | p6)

final_plot

