#! usr/bin/Rscript

library(ggplot2)

#
# LOADING DATA 
#

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

prom_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_promoter_beta_15.csv")
rownames(prom_cassettes) <- prom_cassettes$Cassette
prom_cassettes$Cassette <- NULL

# Association of cassettes with TILs

# DGetting TILs
main_var <- as.numeric(x[colnames(prom_cassettes),"TILs"])

# Define the variables to test
variables <- rownames(prom_cassettes)

# Compute Kendall’s Tau and p-values
results <- data.frame(
  Cassette = variables,
  Tau = sapply(variables, function(v) cor.test(main_var, as.numeric(prom_cassettes[as.character(v),]), method = "kendall")$estimate),
  P_value = sapply(variables, function(v) cor.test(main_var, as.numeric(prom_cassettes[as.character(v),]), method = "kendall")$p.value)
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$P_value, method = "fdr")  # Recommended for multiple comparisons

#Change rownames
rownames(results) <- results$Cassette

# PLOTTING. VOLCANO PLOT

# Convert p-values to -log10 scale
results$logP <- -log10(results$P_adj_Bonf)

# Define significance threshold
threshold <- 0.05

# Create the volcano plot
ggplot(results, aes(x = Tau, y = logP, label = Cassette)) +
  geom_point(aes(color = P_adj_Bonf < threshold), size = 3, alpha = 0.8) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(vjust = 1.5, size = 4, check_overlap = TRUE) +  # Add labels
  theme_minimal() +
  labs(
    title = "Volcano Plot of Kendall’s Tau Correlations",
    x = "Kendall’s Tau",
    y = "-log10(Bonf. P-value)",
    color = "Significant (P < 0.05)"
  ) +
  theme(
    legend.position = "top",
    text = element_text(size = 14)
  ) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "blue")


# Assuming your dataframe is named df
results_sorted <- results[order(-abs(results$Tau)), ]

promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")
data <- promoter_15$colors
cpgs <- names(data)[data == 336]
genes[cpgs]


# Association of cassettes with HRD

# Define the main categorical variable (binary for Mann-Whitney)
main_var <- x[colnames(prom_cassettes), "HRD.2.status"]

# Define the variables to test
variables <- rownames(prom_cassettes)

# Function to compute Mann-Whitney U test p-value
mann_whitney_p <- function(v) {
  wilcox.test(as.numeric(prom_cassettes[as.character(v),]) ~ main_var)$p.value
}

# Function to compute eta-squared
eta_squared_value <- function(v) {
  model <- aov(as.numeric(prom_cassettes[as.character(v),]) ~ main_var)
  ss_total <- sum(model$residuals^2) + sum((model$fitted.values - mean(model$fitted.values))^2)
  ss_between <- sum((model$fitted.values - mean(model$fitted.values))^2)
  return(ss_between / ss_total)  # Eta-squared
}

# Compute results
results <- data.frame(
  Cassette = variables,
  Eta_Squared = sapply(variables, eta_squared_value),
  P_value = sapply(variables, mann_whitney_p)
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$P_value, method = "fdr")  # Recommended for multiple comparisons

# View results
head(results)

#Change rownames
rownames(results) <- results$Cassette

# Convert p-values to -log10 scale
results$logP <- -log10(results$P_adj_Bonf)

# Define significance threshold
threshold <- 0.05

# Create the volcano plot
ggplot(results, aes(x = Eta_Squared, y = logP, label = Cassette)) +
  geom_point(aes(color = P_adj_Bonf < threshold), size = 3, alpha = 0.8) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text(vjust = 1.5, size = 4, check_overlap = TRUE) +  # Add labels
  theme_minimal() +
  labs(
    title = "Volcano Plot of Eta_Squared",
    x = "Eta_Squared",
    y = "-log10(Bonf. P-value)",
    color = "Significant (P < 0.05)"
  ) +
  theme(
    legend.position = "top",
    text = element_text(size = 14)
  ) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "blue")


# Assuming your dataframe is named df
results_sorted <- results[order(-abs(results$Eta_Squared)), ]

promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
data <- promoter_15$colors
cpgs <- names(data)[data == 71]
genes[cpgs]

cor(as.numeric(prom_cassettes["675",]),
    as.numeric(fpkm_data["IRX1",colnames(prom_cassettes)]))




results_sorted
