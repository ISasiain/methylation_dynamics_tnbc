#! usr/bin/Rscript

library(ggplot2)
library(survival)
library(survminer)
library(ggrepel)


#
# LOADING DATA 
#

# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

prom_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_15.csv")
rownames(prom_cassettes) <- prom_cassettes$Cassette
prom_cassettes$Cassette <- NULL

dis_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_15.csv")
rownames(dis_cassettes) <- dis_cassettes$Cassette
dis_cassettes$Cassette <- NULL

prox_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_15.csv")
rownames(prox_cassettes) <- prox_cassettes$Cassette
prox_cassettes$Cassette <- NULL



# Generate obkects for genes linked to CpGs
genes <- annoObj$nameUCSCknownGeneOverlap <- sapply(annoObj$nameUCSCknownGeneOverlap, function(x) {
  if (grepl("\\|", x)) {
    strsplit(x, "\\|")[[1]][2]
  } else {
    x
  }
})

names(genes) <- annoObj$illuminaID


# Loading gene expression
fpkm_data <- read.table("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/gexFPKM_unscaled.csv")


#
# ASSOCIATION OF CASSETTES WITH TILS
#

# Getting TILs
main_var <- as.numeric(x[colnames(prom_cassettes),"TILs"])

# PROMOTER

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

# Plotting results. Volcano plot

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

head(results_sorted, 30)

promoter_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_8.rds")
data <- promoter_15$colors
cpgs <- names(data)[data == 8]
genes[cpgs]


# DISTAL

# Define the variables to test
variables <- rownames(dis_cassettes)

# Compute Kendall’s Tau and p-values
results <- data.frame(
  Cassette = variables,
  Tau = sapply(variables, function(v) cor.test(main_var, as.numeric(dis_cassettes[as.character(v),]), method = "kendall")$estimate),
  P_value = sapply(variables, function(v) cor.test(main_var, as.numeric(dis_cassettes[as.character(v),]), method = "kendall")$p.value)
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$P_value, method = "fdr")  # Recommended for multiple comparisons

#Change rownames
rownames(results) <- results$Cassette

# Plotting results. Volcano plot

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

distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
data <- distal_15$colors
cpgs <- names(data)[data == 40]
genes[cpgs]


# PROXIMAL

# Define the variables to test
variables <- rownames(prox_cassettes)

# Compute Kendall’s Tau and p-values
results <- data.frame(
  Cassette = variables,
  Tau = sapply(variables, function(v) cor.test(main_var, as.numeric(prox_cassettes[as.character(v),]), method = "kendall")$estimate),
  P_value = sapply(variables, function(v) cor.test(main_var, as.numeric(prox_cassettes[as.character(v),]), method = "kendall")$p.value)
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$P_value, method = "fdr")  # Recommended for multiple comparisons

#Change rownames
rownames(results) <- results$Cassette

# Plotting results. Volcano plot

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

head(results_sorted, 20)

proximal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
data <- proximal_15$colors
cpgs <- names(data)[data == 6]
genes[cpgs]


# Analyse genes with the highest associations


# GBP4

# Define variables
x_values <- as.numeric(prom_cassettes["141",])  # X-axis variable
y_values <- as.numeric(fpkm_data["GBP4", colnames(prom_cassettes)])  # Y-axis variable
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
y_values <- as.numeric(fpkm_data["ZBP1", colnames(prom_cassettes)])  # Y-axis variable
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

# CELF2

# Define variables
x_values <- as.numeric(prox_cassettes["6",])  # X-axis variable
y_values <- as.numeric(fpkm_data["CELF2", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "CELF2 Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()

# PPP1R36

# Define variables
x_values <- as.numeric(prom_cassettes["250",])  # X-axis variable
y_values <- as.numeric(fpkm_data["PPP1R36", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "PPP1R36 Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()

# PCDHGA2

# Define variables
x_values <- as.numeric(prom_cassettes["183",])  # X-axis variable
y_values <- as.numeric(fpkm_data["PCDHGA2", colnames(prom_cassettes)])  # Y-axis variable
tils_values <- as.numeric(main_var) / 100  

# Create a dataframe for ggplot
plot_data <- data.frame(X = x_values, Y = y_values, TILs = tils_values)

# Create scatterplot with both color and size representing TILs
ggplot(plot_data, aes(x = X, y = Y)) +
  geom_point(aes(color = TILs, size = TILs), alpha = 0.7) +
  scale_color_gradient(low = "black", high = "red") +  # Black-to-Red gradient
  scale_size(range = c(0.3, 1.5)) +  # Adjust point sizes
  labs(x = "Cassette PC1", 
       y = "PCDHGA2 Expression (FPKM)", 
       size = "TILs (%)") +  # Keep only the size label
  theme_classic()


#
# IMPACT IN OUTCOME
#

IDFS <- x[colnames(prom_cassettes),"IDFS"]
IDFS_bin <- x[colnames(prom_cassettes),"IDFSbin"]

# Use scaled values
prom_cassettes <- scale(prom_cassettes)
dis_cassettes <- scale(dis_cassettes)
prox_cassettes <- scale(prox_cassettes)

# PROMOTER

# Define all cassette genes
all_cassettes <- rownames(prom_cassettes)
all_cassettes_names <- all_cassettes # Assuming gene names are the same as rownames

# Initialize a data frame to store results
results <- data.frame(Cassette = character(), HR = numeric(), Lower = numeric(), Upper = numeric(), pvalue = numeric())

# Loop through all genes
for (i in 1:length(all_cassettes)) {
  predictor_values <- as.numeric(prom_cassettes[all_cassettes[i], ])
  
  cox_model <- coxph(Surv(IDFS, IDFS_bin) ~ predictor_values)
  model_summary <- summary(cox_model)
  
  # Extract HR, confidence intervals, and p-value
  hr <- model_summary$coefficients[,"exp(coef)"]
  lower_ci <- model_summary$conf.int[,"lower .95"]
  upper_ci <- model_summary$conf.int[,"upper .95"]
  pvalue <- model_summary$coefficients[,"Pr(>|z|)"]
  
  # Store results
  results <- rbind(results, data.frame(Cassette = all_cassettes_names[i], HR = hr, Lower = lower_ci, Upper = upper_ci, pvalue = pvalue))
}

# FDR correction
results$pvalue_adj <- p.adjust(results$pvalue, method = "fdr")

# Create volcano plot
ggplot(results, aes(x = log2(HR), y = -log10(pvalue_adj), label = Cassette)) +
  geom_point(aes(color = pvalue_adj < 0.05), size = 3) +
  geom_text_repel() +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  xlab("Log2 Hazard Ratio") + 
  ylab("-Log10 Adjusted P-value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  ggtitle("Volcano Plot: Association of All Cassettes PC1 to Outcome")

my_genes <- genes[names(promoter_15$colors)[promoter_15$colors==250]]
results[548,]

# Plotting metastasis type vs cassette 183

metastasis_type <- as.character(x[colnames(prom_cassettes), "Metastasis_type"])
metastasis_type[is.na(metastasis_type)] <- "NA"  # Convert NA to a string
boxplot(as.numeric(prom_cassettes["183",]) ~ as.factor(metastasis_type),
        xlab="Metastasis type",
        ylab="Cassette PC1",
        main=paste0("KW p=", round(kruskal.test(as.numeric(prom_cassettes["183", ]) ~ as.factor(metastasis_type))$p.value, 5)),
        border = "black",  # Dark borders for better contrast
        las = 1,  # Horizontal axis labels for readability
        notch = F,  # Add notches for confidence intervals
        cex.axis = 1.2,  # Increase axis text size
        cex.lab = 1.4,  # Increase label text size
        frame = FALSE,  # Remove default box around plot
        main = "PC1 Cassette 1 across PAM50 Subtypes",  # Add a clear title
        outpch = 16,  # Use solid circles for outliers
        outcol = "red"  # Highlight outliers in red
)



# DISTAL

# Define all cassette genes
all_cassettes <- rownames(dis_cassettes)
all_cassettes_names <- all_cassettes # Assuming gene names are the same as rownames

# Initialize a data frame to store results
results <- data.frame(Cassette = character(), HR = numeric(), Lower = numeric(), Upper = numeric(), pvalue = numeric())

# Loop through all genes
for (i in 1:length(all_cassettes)) {
  predictor_values <- as.numeric(dis_cassettes[all_cassettes[i], ])
  
  cox_model <- coxph(Surv(IDFS, IDFS_bin) ~ predictor_values)
  model_summary <- summary(cox_model)
  
  # Extract HR, confidence intervals, and p-value
  hr <- model_summary$coefficients[,"exp(coef)"]
  lower_ci <- model_summary$conf.int[,"lower .95"]
  upper_ci <- model_summary$conf.int[,"upper .95"]
  pvalue <- model_summary$coefficients[,"Pr(>|z|)"]
  
  # Store results
  results <- rbind(results, data.frame(Cassette = all_cassettes_names[i], HR = hr, Lower = lower_ci, Upper = upper_ci, pvalue = pvalue))
}

# FDR correction
results$pvalue_adj <- p.adjust(results$pvalue, method = "fdr")

# Create volcano plot
ggplot(results, aes(x = log2(HR), y = -log10(pvalue_adj), label = Cassette)) +
  geom_point(aes(color = pvalue_adj < 0.05), size = 3) +
  geom_text_repel() +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  xlab("Log2 Hazard Ratio") + 
  ylab("-Log10 Adjusted P-value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  ggtitle("Volcano Plot: Association of All Cassettes PC1 to Outcome")

genes[names(distal_15$colors)[distal_15$colors==59]]
results[60,]


# PROXIMAL


# Define all cassette genes
all_cassettes <- rownames(prox_cassettes)
all_cassettes_names <- all_cassettes # Assuming gene names are the same as rownames

# Initialize a data frame to store results
results <- data.frame(Cassette = character(), HR = numeric(), Lower = numeric(), Upper = numeric(), pvalue = numeric())

# Loop through all genes
for (i in 1:length(all_cassettes)) {
  predictor_values <- as.numeric(prox_cassettes[all_cassettes[i], ])
  
  cox_model <- coxph(Surv(IDFS, IDFS_bin) ~ predictor_values)
  model_summary <- summary(cox_model)
  
  # Extract HR, confidence intervals, and p-value
  hr <- model_summary$coefficients[,"exp(coef)"]
  lower_ci <- model_summary$conf.int[,"lower .95"]
  upper_ci <- model_summary$conf.int[,"upper .95"]
  pvalue <- model_summary$coefficients[,"Pr(>|z|)"]
  
  # Store results
  results <- rbind(results, data.frame(Cassette = all_cassettes_names[i], HR = hr, Lower = lower_ci, Upper = upper_ci, pvalue = pvalue))
}

# FDR correction
results$pvalue_adj <- p.adjust(results$pvalue, method = "fdr")

# Create volcano plot
ggplot(results, aes(x = log2(HR), y = -log10(pvalue_adj), label = Cassette)) +
  geom_point(aes(color = pvalue_adj < 0.05), size = 3) +
  geom_text_repel() +
  scale_color_manual(values = c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  xlab("Log2 Hazard Ratio") + 
  ylab("-Log10 Adjusted P-value") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  ggtitle("Volcano Plot: Association of All Cassettes PC1 to Outcome")

genes[names(proximal_15$colors)[proximal_15$colors==40]]


H