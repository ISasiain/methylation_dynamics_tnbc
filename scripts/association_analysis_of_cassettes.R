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

head(results_sorted)

promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_15.rds")
data <- promoter_15$colors
cpgs <- names(data)[data == 250]
genes[cpgs]

# Analyse genes with the highest associations


# Analyse genes with the highest associations

# GBP4

# Define variables
x_values <- as.numeric(prom_cassettes["141",])  # X-axis variable
c <- as.numeric(fpkm_data["GBP4", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "GBP4 Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()


# ZBP1

# Define variables
x_values <- as.numeric(prom_cassettes["547",])  # X-axis variable
c <- as.numeric(fpkm_data["ZBP1", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "ZBP1 Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()

# NOSTRIN

# Define variables
x_values <- as.numeric(prom_cassettes["453",])  # X-axis variable
c <- as.numeric(fpkm_data["NOSTRIN", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "NOSTRIN Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()
