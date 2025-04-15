#! usr/bin/Rscript

library(ComplexHeatmap)
library(NMF)
library(ggplot2)
library(dplyr)
library(ggalluvial)

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
nonbasal_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_nonBasal_cassettes_beta_5.rds")
basal_distal_5 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/distal/only_basal_cassettes_beta_5.rds")


# Filtering Distal ATAC CpGs
distal_atac_cpgs <- annoObj$illuminaID[which( ((annoObj$featureClass=="distal") | (annoObj$featureClass=="distal body")) & annoObj$hasAtacOverlap == 1)]
distal_atac_betas <- betaAdj[rownames(betaAdj) %in% distal_atac_cpgs, ]


# 
# NMF CLUSTERING 1
#

# Filtering data

# Calculating variance of CpGs
variance_vector <- apply(X=distal_atac_betas, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. Exclude cassette 1
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]
filtered_top_variance_cpgs <- top_variance_cpgs[!distal_5$colors[top_variance_cpgs] %in% c(0, 1)]

# Getting beta matrix
betas_to_analyse1 <- betaAdj[filtered_top_variance_cpgs,]
Heatmap(betas_to_analyse1, use_raster = F)
  

r.nmf <- nmf(betas_to_analyse1,
             rank = 2,
             nrun = 100,
             method = 'brunet',
             seed = 221027)



W <- basis(r.nmf)
H <- coef(r.nmf)

clusters.1 <- predict(r.nmf, "chc")

# Heatmap of clusters and CpG betas
annotations <- HeatmapAnnotation(
  
  "orig_epi"=orig_epi[names(clusters.1), "NMF_atacDistal"],
  "new_epi"=clusters.1
  
)

my_heatmap <- Heatmap(betas_to_analyse1, use_raster=F,
                      top_annotation = annotations,
                      column_order = order(orig_epi[names(clusters.1), "NMF_atacDistal"]))

# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters.1

original_epitypes <- orig_epi[names(clusters.1), "NMF_atacDistal"]
original_epitypes <- substr(original_epitypes, 1, nchar(original_epitypes) - 1)
original_epitypes <- paste0(original_epitypes, "_epi")

pam50 <- x[names(clusters.1), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(clusters.1),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
clusters_plot <- ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")

# Saving plots

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/all_filtered_with_all_heatmap.pdf")
print(my_heatmap)  
dev.off()

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/all_filtered_with_all_clustering.pdf")
print(clusters_plot)  
dev.off()



# 
# NMF CLUSTERING 2.1 BASAL
#

#
# Filtering based on cassettes from all samples
#

distal_atac_Basal <- distal_atac_betas[,names(clusters.1)[clusters.1 == 2]]

# Calculating variance of CpGs
variance_vector <- apply(X=distal_atac_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. Exclude cassette 1
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]
filtered_top_variance_cpgs <- top_variance_cpgs[!distal_5$colors[top_variance_cpgs] %in% c(0, 1)]

# Getting beta matrix
betas_to_analyse2.1 <- distal_atac_Basal[filtered_top_variance_cpgs,]
Heatmap(betas_to_analyse2.1, use_raster = F)

# Running NMF
r.nmf <- nmf(betas_to_analyse2.1,
              rank = 3,
              nrun = 100,
              method = 'brunet',
              seed = 221027)

W <- basis(r.nmf)
H <- coef(r.nmf)

clusters <- predict(r.nmf,"chc")

# Heatmap of clusters and CpG betas
annotations <- HeatmapAnnotation(
  
  "orig_epi"=orig_epi[names(clusters), "NMF_atacDistal"],
  "new_epi"=clusters
  
)

my_heatmap <- Heatmap(betas_to_analyse2.1, use_raster=F,
                      top_annotation = annotations,
                      column_order = order(orig_epi[names(clusters), "NMF_atacDistal"]))


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters

original_epitypes <- orig_epi[names(clusters), "NMF_atacDistal"]

pam50 <- x[names(clusters), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(clusters),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
clusters_plot <- ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")

# Saving plots

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/basal_filtered_with_all_heatmap.pdf")
print(my_heatmap)  
dev.off()

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/basal_filtered_with_all_clustering.pdf")
print(clusters_plot)  
dev.off()


#
# Filtering based on cassettes from Basal samples
#

distal_atac_Basal <- distal_atac_betas[,names(clusters.1)[clusters.1 == 2]]

# Calculating variance of CpGs
variance_vector <- apply(X=distal_atac_Basal, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. Exclude cassette 1
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]
filtered_top_variance_cpgs <- top_variance_cpgs[!basal_distal_5$colors[top_variance_cpgs] %in% c(0, 1)]

# Getting beta matrix
betas_to_analyse2.1 <- distal_atac_Basal[filtered_top_variance_cpgs,]
Heatmap(betas_to_analyse2.1, use_raster = F)

# Running NMF
r.nmf <- nmf(betas_to_analyse2.1,
                rank = 3,
                nrun = 100,
                method = 'brunet',
                seed = 221027)

W <- basis(r.nmf)
H <- coef(r.nmf)


clusters <- predict(r.nmf,"sample")

# Heatmap of clusters and CpG betas
annotations <- HeatmapAnnotation(
  
  "orig_epi"=orig_epi[names(clusters), "NMF_atacDistal"],
  "new_epi"=clusters
  
)

my_heatmap <- Heatmap(betas_to_analyse2.1, use_raster=F,
                      top_annotation = annotations,
                      column_order = order(orig_epi[names(clusters), "NMF_atacDistal"]))


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters

original_epitypes <- orig_epi[names(clusters), "NMF_atacDistal"]

pam50 <- x[names(clusters), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(clusters),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
clusters_plot <- ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")

# Saving plots

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/basal_filtered_with_basal_heatmap.pdf")
print(my_heatmap)  
dev.off()

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/basal_filtered_with_basal_clustering.pdf")
print(clusters_plot)  
dev.off()


# 
# NMF CLUSTERING 2.1 NON-BASAL
#

#
# Filtering based on cassettes from all samples
#

distal_atac_nonBasal <- distal_atac_betas[,names(clusters.1)[clusters.1 == 1]]

# Calculating variance of CpGs
variance_vector <- apply(X=distal_atac_nonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. Exclude cassette 1
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]
filtered_top_variance_cpgs <- top_variance_cpgs[!distal_5$colors[top_variance_cpgs] %in% c(0, 1)]

# Getting beta matrix
betas_to_analyse2.2 <- distal_atac_nonBasal[filtered_top_variance_cpgs,]
Heatmap(betas_to_analyse2.2, use_raster = F)

# Running NMF
r.nmf <- nmf(betas_to_analyse2.2,
             rank = 2,
             nrun = 100,
             method = 'brunet',
             seed = 221027)

W <- basis(r.nmf)
H <- coef(r.nmf)

clusters <- predict(r.nmf,"chc")

# Heatmap of clusters and CpG betas
annotations <- HeatmapAnnotation(
  
  "orig_epi"=orig_epi[names(clusters), "NMF_atacDistal"],
  "new_epi"=clusters
  
)

my_heatmap <- Heatmap(betas_to_analyse2.2, use_raster=F,
                      top_annotation = annotations,
                      column_order = order(orig_epi[names(clusters), "NMF_atacDistal"]))


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters

original_epitypes <- orig_epi[names(clusters), "NMF_atacDistal"]

pam50 <- x[names(clusters), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(clusters),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
clusters_plot <- ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")


# Saving plots

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/non_basal_filtered_with_all_heatmap.pdf")
print(my_heatmap)  
dev.off()

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/non_basal_filtered_with_all_clustering.pdf")
print(clusters_plot)  
dev.off()



#
# Filtering based on cassettes from nonBasal samples
#

distal_atac_nonBasal <- distal_atac_betas[,names(clusters.1)[clusters.1 == 2]]

# Calculating variance of CpGs
variance_vector <- apply(X=distal_atac_nonBasal, MARGIN = 1, FUN=var)

# Filtering based on variance and based on detected cassettes. Exclude cassette 1
top_variance_cpgs <- names(sort(variance_vector, decreasing = T))[1:5000]
filtered_top_variance_cpgs <- top_variance_cpgs[!nonbasal_distal_5 $colors[top_variance_cpgs] %in% c(0, 1)]

# Getting beta matrix
betas_to_analyse2.1 <- distal_atac_nonBasal[filtered_top_variance_cpgs,]
Heatmap(betas_to_analyse2.1, use_raster = F)

# Running NMF<<
r.nmf <- nmf(betas_to_analyse2.1,
             rank = 2,
             nrun = 100,
             method = 'brunet',
             seed = 221027)

W <- basis(r.nmf)
H <- coef(r.nmf)


clusters <- predict(r.nmf,"chc")

# Heatmap of clusters and CpG betas
annotations <- HeatmapAnnotation(
  
  "orig_epi"=orig_epi[names(clusters), "NMF_atacDistal"],
  "new_epi"=clusters
  
)

my_heatmap <- Heatmap(betas_to_analyse2.1, use_raster=F,
        top_annotation = annotations,
        column_order = order(orig_epi[names(clusters), "NMF_atacDistal"]))


# Visualize clusterings

# Getting sample groups / clusters
new_epitypes <- clusters

original_epitypes <- orig_epi[names(clusters), "NMF_atacDistal"]

pam50 <- x[names(clusters), "PAM50_Basal_NCN"]

# Create a data frame of transitions
df <- data.frame(
  Sample = names(clusters),
  New_epitypes = new_epitypes,
  Original_epitypes = original_epitypes,
  PAM50 = pam50
)

# Plot Alluvial Diagram
clusters_plot <- ggplot(df, aes(axis1 = Original_epitypes, axis2 = new_epitypes, axis3 = PAM50)) +
  geom_alluvium(aes(fill = as.factor(new_epitypes))) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void() +
  labs(x = "Clusterings", y = "Count", title = "Comparison of Cluster Assignments")

# Saving plots

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/non_basal_filtered_with_nonbasal_heatmap.pdf")
print(my_heatmap)  
dev.off()

pdf("/Users/isasiain/PhD/Projects/project_3/analysis/nmf_comparison/non_basal_filtered_with_nonbasal_clustering.pdf")
print(clusters_plot)  
dev.off()


