#! usr/bin/Rscript

library(ggplot2)
library(survival)
library(survminer)
library(ggrepel)
library(energy)
library(mclust)
library(reshape2)


#
# LOADING DATA 
#
s
# Loading EPIC methylation matrix
load("/Volumes/Data/Project_3/TNBC_epigenetics/workspace_full_trim235_updatedSampleAnno_withNmfClusters.RData")

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

prom_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(prom_cassettes) <- prom_cassettes$Cassette
prom_cassettes$Cassette <- NULL

dis_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/distal_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(dis_cassettes) <- dis_cassettes$Cassette
dis_cassettes$Cassette <- NULL

prox_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/proximal_cassettes/summary_cassettes/summary_beta_10.csv")
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

# CpG cassettes
promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")
distal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_10.rds")



#
# ASSOCIATION OF CASSETTES WITH TILS
#

# Getting TILs
main_var <- as.numeric(x[colnames(prom_cassettes),"TILs"])

# PROMOTER

# Define the variables to test
variables <- rownames(prom_cassettes)

# Compute correlations and wilocoxon p values. Kendall + Distance
results <- data.frame(
  Cassette = variables,
  Wilcoxon_P_value = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    if (length(unique(y)) < 2) return(NA)
    
    # K-means clustering into 2 groups
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Wilcoxon test between groups for TIL levels
    wilcox.test(x ~ cluster)$p.value
  }),
  
  Kendall_Tau = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    # Kendall's Tau correlation
    cor.test(x, y, method = "kendall")$estimate
  }),
  
  Distance_Correlation = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    if (all(is.na(y))) return(NA)
    
    # Keep only non-NA pairs
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    # Distance Correlation
    dcor(x, y)
  }),
  
  Median_TIL_Hypo = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Selecting hypermethylated cluster
    if (mean(y[cluster == 1]) > mean(y[cluster == 2])) {
      mean(x[cluster == 1])
    } else {
      median(x[cluster == 2])
    }
    
    
  }),
  
  Median_TIL_Hyper = sapply(variables, function(v) {
    y <- as.numeric(prom_cassettes[as.character(v), ])
    keep <- complete.cases(main_var, y)
    x <- main_var[keep]
    y <- y[keep]
    if (length(unique(y)) < 2) return(NA)
    cluster <- kmeans(y, centers = 2)$cluster
    
    # Selecting hypermethylated cluster
    if (mean(y[cluster == 2]) > mean(y[cluster == 1])) {
      mean(x[cluster == 1])
    } else {
      median(x[cluster == 2])
    }
  }),
  
  gene = sapply(variables, function(v) {
    
    data <- promoter_10$colors
    cpgs <- names(data)[data == v]
    paste0(unique(unname(genes[cpgs])), collapse = ",")
  }) 
  
)

# Apply multiple testing correction
results$P_adj_Bonf <- p.adjust(results$Wilcoxon_P_value, method = "bonferroni")
results$P_adj_FDR <- p.adjust(results$Wilcoxon_P_value, method = "fdr")  # Recommended for multiple comparisons

# Convert p-values to -log10 scale
results$logP_FDR <- -log10(results$P_adj_FDR)
results$logP_bonferroni <- -log10(results$P_adj_Bonf)

#Change rownames
rownames(results) <- results$Cassette

# Set significance threshold and number of top hits to label
threshold <- 0.05
top_n <- 30

# Create a subset with only the top significant results to label
top_labels <- results %>%
  filter(P_adj_FDR < threshold) %>%
  arrange(desc(logP_FDR)) %>%
  slice_head(n = top_n)

# Set threshold and number of top hits to label
threshold <- 0.05
top_n <- 30

# Subset the top significant hits to label
top_labels <- results %>%
  filter(P_adj_FDR < threshold) %>%
  arrange(desc(logP_FDR)) %>%
  slice_head(n = top_n)

# Volcano plot
ggplot(results, aes(x = Kendall_Tau, y = logP_FDR)) +
  geom_point(aes(color = P_adj_FDR < threshold), size = 3, alpha = 0.8) + 
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(
    data = top_labels,
    aes(label = Cassette),
    size = 4,
    max.overlaps = Inf,
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", color = "blue") +
  theme_bw() +
  labs(
    x = "Kendall’s Tau",
    y = "-log10(FDR. P-value)"
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 14)
  )


results_sorted <- results[order(-abs(results$Kendall_Tau)), ]

head(results_sorted, 30)


Heatmap(betaNew[cpgs,])

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


results_sorted <- results[order(-abs(results$Kendall_Tau)), ]

head(results_sorted)

distal_15 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/cassettes_beta_15.rds")
data <- distal_15$colors
cpgs <- names(data)[data == 867]
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


results_sorted <- results[order(-abs(results$Tau)), ]

head(results_sorted, 20)

proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")
data <- proximal_10$colors
cpgs <- names(data)[data == 6]
genes[cpgs]


#
# IMPACT IN OUTCOME
#

# Getting data
IDFS <- x[colnames(prom_cassettes),"DRFI"]
IDFS_bin <- x[colnames(prom_cassettes), "DRFIbin"]

chemo <- unname(sapply(x[colnames(prom_cassettes),"Chemotherapy"], function(x) {
  if (is.na(x)) {return(NA)}
  else if (x %in% c("Metastatic from beginnig", "Palliative", "FEC100X1 (STOPPED BY PATIENT)")) {return(NA)}
  else if (x == "None") {return("No_chemo")}
  else {return("Chemo")}
}))


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
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  xlab("Log2 Hazard Ratio") + 
  ylab("-Log10 Adjusted P-value") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14))


#
# PLOTTING CPGS AFFECTING 
#

data <- promoter_10$colors
cpgs <- names(data)[data == 10]
genes[cpgs]

pam50_annotations <- my_annotations[colnames(betaAdj), "PAM50"]
tnbc_annotation <- my_annotations[colnames(betaAdj), "TNBC"]
HRD_annotation <- my_annotations[colnames(betaAdj), "HRD"]
epi_annotation <- my_annotations[colnames(betaAdj), "NMF_atacDistal"]
im_annotation <- my_annotations[colnames(betaAdj), "IM"]
tils_annotation <- as.numeric(x[colnames(betaAdj), "TILs"])

my_fpkm_data <- as.numeric(fpkm_data["EPO", colnames(betaAdj)])

# Create top anotation
top_annotation <- HeatmapAnnotation(PAM50 = pam50_annotations,
                                    TNBC = tnbc_annotation,
                                    HRD = HRD_annotation,
                                    IM = im_annotation,
                                    Epitype = epi_annotation,
                                    TILs = anno_points(tils_annotation,
                                                       ylim=c(0,100),
                                                       size=unit(0.75, "mm"),
                                                       axis_param = list(
                                                         side="left",
                                                         at=c(0,25,50,75,100),
                                                         labels=c("0","25","50","75","100")
                                                       )),
                                    col = list(
                                      "PAM50"=c("Basal"="indianred1", "Her2"="pink", "LumA"="darkblue", "LumB"="lightblue", "Normal"="darkgreen", "Uncl."="grey"),
                                      "HRD"=c("High"="darkred", "Low/Inter"="lightcoral"),
                                      "IM"=c("Negative"="grey", "Positive"="black"),
                                      "TNBC"=c("BL1"="red", "BL2"="blue", "LAR"="green", "M"="grey"),
                                      "Epitype"=c("Basal1" = "tomato4", "Basal2" = "salmon2", "Basal3" = "red2", 
                                                  "nonBasal1" = "cadetblue1", "nonBasal2" = "dodgerblue")
                                    )
)


# Updated left_annotation with color scale
right_annotation <- rowAnnotation(
  "Normal beta" = rowMeans(betaNorm[cpgs,]),
  "ATAC" = annoObj$hasAtacOverlap[annoObj$illuminaID %in% cpgs],
  col = list("Normal beta" = colorRamp2(c(0, 0.5, 1), c("darkblue", "white", "darkred")),
             "ATAC" = c("0" = "white", "1"= "black"))
)

# CpG context annotation
left_annotation <- rowAnnotation("Context"= annoObj$featureClass[annoObj$illuminaID %in% cpgs]
)

# Expression annotation
bottom_annotation <- HeatmapAnnotation(
  "FPKM" = anno_barplot(my_fpkm_data)
)

# Cluster based on methylation
cluster_promoter <- kmeans(t(betaAdj[cpgs,]), centers = 2)

# Determine hypo and hypermethylated cluster
promoter_state <- if (mean(betaAdj[cpgs,cluster_promoter$cluster==1]) >
                           mean(betaAdj[cpgs,cluster_promoter$cluster==2])) {
  
  as.factor(ifelse(cluster_promoter$cluster == 1, "Hypermethylated", "Hypomethylated"))
  
} else {
  
  as.factor(ifelse(cluster_promoter$cluster == 2, "Hypermethylated", "Hypomethylated"))
  
}


# Heatmap of genes
Heatmap(
  betaAdj[cpgs,],
  cluster_rows = FALSE,
  row_order = order(annoObj$start[annoObj$illuminaID %in% cpgs]),
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = FALSE,
  column_split = promoter_state,
  top_annotation = top_annotation,
  right_annotation = right_annotation,
  bottom_annotation = bottom_annotation,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  name = "Tumor beta"
)

# Kaplan meier plot of based on methylation state

# Create the survival object
surv_obj <- Surv(time = IDFS, event = IDFS_bin)

# Fit Kaplan-Meier survival curves stratified by promoter_state
fit <- survfit(surv_obj ~ promoter_state)

# Plot Kaplan-Meier curves
ggsurvplot(
  fit,
  data = data.frame(IDFS = IDFS, IDFS_bin = IDFS_bin, promoter_state = promoter_state),
  pval = TRUE,                  # Add p-value from log-rank test
  conf.int = TRUE,              # Add confidence interval
  risk.table = TRUE,            # Show number at risk table
  legend.title = "Promoter State",
  xlab = "Time (IDFS)",
  ylab = "Survival Probability",
  palette = "Dark2"             
)

# Defining 4 groups, Split based on chemo
groups <- interaction(promoter_state, chemo)

# Fit Kaplan-Meier survival curves stratified by promoter_state
fit <- survfit(surv_obj ~ groups)

# Plot Kaplan-Meier curves
ggsurvplot(
  fit,
  data = data.frame(IDFS = IDFS, IDFS_bin = IDFS_bin, promoter_state = promoter_state),
  pval = TRUE,                 
  conf.int = FALSE,              
  risk.table = TRUE,
  legend.labs = NULL,
  legend.title = "Promoter State",
  xlab = "Time",
  ylab = "IDFS Probability",
  palette = "Dark2"             
)

# Plotting metastasis type vs cassette 233
metastasis_type <- as.character(x[colnames(prom_cassettes), "Metastasis_type"])
metastasis_type[is.na(metastasis_type)] <- "NA"  # Convert NA to a string
boxplot(as.numeric(prom_cassettes["233",]) ~ as.factor(metastasis_type),
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
