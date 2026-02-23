#!usr/bin/Rscript

library(missMethyl)
library(ggplot2)
library(patchwork)

#
# LOADING DATA
#

promoter_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_7.rds")
proximal_7 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_7.rds")

hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))

# Load annotation file
load("/Users/in2245sa/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID

# Create a new grouped feature class
annoObj$CpG_context <- feature_class_grouped <- dplyr::case_when(
  annoObj$featureClass %in% c("distal", "distal body") ~ "Distal",
  annoObj$featureClass %in% c("promoter") ~ "Promoter",
  annoObj$featureClass %in% c("proximal dn", "proximal up") ~ "Proximal",
  TRUE ~ as.character(annoObj$featureClass) 
)

#
# PERFORMING GO ENRICHMENT ANALYSIS PER CASSETTE
#

# PROMOTER
list_of_go <- list()

# Performing GO enrichment per cassette CpGs. Using first 60 cassettes
for (cassette in 1:60) {

  # Getting CpGs from cassette
  cpg_subset <- names(promoter_7$colors)[promoter_7$colors == cassette]

  # Performing GO enrichment
  list_of_go[[cassette]] <- gometh(
    cpg_subset,
    all.cpg = names(promoter_7$colors),
    collection = "go",
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_go, "PhD/Projects/project_3/analysis/GO_enrichment_promoter_7.rds")

# PROXIMAL

list_of_go <- list()

# Performing GO enrichment per cassette CpGs. 
for (cassette in 1:60) {
  
  # Getting CpGs from cassette
  cpg_subset <- names(proximal_7$colors)[proximal_7$colors == cassette]
  
  # Performing GO enrichment
  list_of_go[[cassette]] <- gometh(
    cpg_subset,
    all.cpg = names(proximal_7$colors),
    collection = "go",
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_go, "PhD/Projects/project_3/analysis/GO_enrichment_proximal_7.rds")


#
# PERFORMING HALLMARK ENRICHMENT ANALYSIS PER CASSETTE
#

# PROMOTER


list_of_hallmarks <- list()

for (cassette in 1:60) {
  
  # Getting CpGs from cassette
  cpg_subset <- names(promoter_7$colors)[promoter_7$colors == cassette]
  
  # Performing GO enrichment
  list_of_hallmarks[[cassette]] <- gsameth(
    cpg_subset,
    all.cpg = names(promoter_7$colors),
    collection = hallmark,
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_hallmarks, "PhD/Projects/project_3/analysis/hallmark_enrichment_promoter_7.rds")


# PROXIMAL

list_of_hallmarks <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:60) {

  # Getting CpGs from cassette
  cpg_subset <- names(proximal_7$colors)[proximal_7$colors == cassette]

  # Performing GO enrichment
  list_of_hallmarks[[cassette]] <- gsameth(
    cpg_subset,
    all.cpg = names(proximal_7$colors),
    collection = hallmark,
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_hallmarks, "PhD/Projects/project_3/analysis/hallmark_enrichment_proximal_7.rds")


#
# PERFORMING KEGG ENRICHMENT ANALYSIS PER CASSETTE
#

# PROMOTER
list_of_kegg <- list()

# Performing enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:60) {
  
  # Getting CpGs from cassette
  cpg_subset <- names(promoter_7$colors)[promoter_7$colors == cassette]
  
  # Performing GO enrichment
  list_of_kegg[[cassette]] <- gometh(
    cpg_subset,
    all.cpg = names(promoter_7$colors),
    collection = "kegg",
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_kegg, "PhD/Projects/project_3/analysis/kegg_enrichment_promoter_7.rds")


# PROXIMAL

list_of_kegg <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:3) {

  # Getting CpGs from cassette
  cpg_subset <- names(proximal_7$colors)[proximal_7$colors == cassette]

  # Performing GO enrichment
  list_of_kegg[[cassette]] <- gometh(
    cpg_subset,
    all.cpg = names(proximal_7$colors),
    collection = "kegg",
    array.type = "EPIC",
    plot.bias = FALSE,
    prior.prob = TRUE,
    equiv.cpg = TRUE,
    fract.counts = TRUE,
    genomic.features = "ALL",
    sig.genes = FALSE
  )
}

# Saving output
saveRDS(list_of_kegg, "PhD/Projects/project_3/analysis/kegg_enrichment_proximal_7.rds")


#
# READING FILES
#

list_of_kegg_promoter <- readRDS("PhD/Projects/project_3/analysis/kegg_enrichment_promoter_7.rds")
list_of_kegg_proximal <- readRDS("PhD/Projects/project_3/analysis/kegg_enrichment_proximal_7.rds")

list_of_hallmarks_promoter <- readRDS("PhD/Projects/project_3/analysis/hallmark_enrichment_promoter_7.rds")
list_of_hallmarks_proximal <- readRDS("PhD/Projects/project_3/analysis/hallmark_enrichment_proximal_7.rds")

list_of_go_promoter <- readRDS("PhD/Projects/project_3/analysis/GO_enrichment_promoter_7.rds")
list_of_go_proximal <- readRDS("PhD/Projects/project_3/analysis/GO_enrichment_proximal_7.rds")

#
# PLOTTING PROMOTER
#


# GO

# Cassette 1
go_df <- list_of_go_promoter[[1]] 
go_df <- go_df[go_df$FDR <= 0.05,]

go_df <- go_df[order(go_df$FDR, decreasing = TRUE),]
go_df$TERM <- factor(go_df$TERM, levels = go_df$TERM)


go_cas_1 <- ggplot(go_df, aes(x = -log10(FDR), y = TERM)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 10) +
  scale_size_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

go_cas_1

# KEGG

# Cassette 1
kegg_df <- list_of_kegg_promoter[[1]] 
kegg_df <- kegg_df[kegg_df$FDR <= 0.05,]

kegg_df <- kegg_df[order(kegg_df$FDR, decreasing = TRUE),]
kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)


kegg_cas_1 <- ggplot(kegg_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 20) +
  scale_size_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

# HALLMARKS

# Cassette 1
hallmarks_df <- list_of_hallmarks_promoter[[1]]
hallmarks_df <- hallmarks_df[hallmarks_df$FDR <= 0.05,]

formatted_rownames <- tolower(gsub("_", " ", gsub("^HALLMARK_", "", rownames(hallmarks_df))))
formatted_rownames <- sub("^([a-z])", "\\U\\1", formatted_rownames, perl = TRUE)

hallmarks_df$Description <- formatted_rownames

hallmarks_df <- hallmarks_df[order(hallmarks_df$FDR, decreasing = TRUE),]
hallmarks_df$Description <- factor(hallmarks_df$Description, levels = hallmarks_df$Description)


hallmarks_cas_1 <- ggplot(hallmarks_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 20) +
  scale_size_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "Hallmarks",
    size = "Gene Ratio"
  )

kegg_cas_1 / hallmarks_cas_1 +
  plot_layout(heights = c(10, 2), guides = "collect") &
  theme(legend.position = "bottom")


#
# PLOTTING PROXIMAL
#

# GO

go_df <- list_of_go_proximal[[1]] 
go_df <- go_df[go_df$FDR <= 0.05,]

go_df <- go_df[order(go_df$FDR, decreasing = TRUE),]
go_df$TERM <- factor(go_df$TERM, levels = go_df$TERM)


go_cas_1 <- ggplot(go_df, aes(x = -log10(FDR), y = TERM)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 12) +
  scale_size_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

go_cas_1

# KEGG

# Cassette 1
kegg_df <- list_of_kegg_proximal[[1]] 
kegg_df <- kegg_df[kegg_df$FDR <= 0.05,]

kegg_df <- kegg_df[order(kegg_df$FDR, decreasing = TRUE),]
kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)


kegg_cas_1 <- ggplot(kegg_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 20) +
  scale_size_continuous(limits = c(0, 0.6), breaks = c(0, 0.2, 0.4, 0.6)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )

kegg_cas_1

# HALLMARKS

# Cassette 1
hallmarks_df <- list_of_hallmarks_proximal[[1]]
hallmarks_df <- hallmarks_df[hallmarks_df$FDR <= 0.05,]

formatted_rownames <- tolower(gsub("_", " ", gsub("^HALLMARK_", "", rownames(hallmarks_df))))
formatted_rownames <- sub("^([a-z])", "\\U\\1", formatted_rownames, perl = TRUE)

hallmarks_df$Description <- formatted_rownames

hallmarks_df <- hallmarks_df[order(hallmarks_df$FDR, decreasing = TRUE),]
hallmarks_df$Description <- factor(hallmarks_df$Description, levels = hallmarks_df$Description)


hallmarks_cas_1 <- ggplot(hallmarks_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 20) +
  scale_size_continuous(limits = c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "Hallmarks",
    size = "Gene Ratio"
  )

kegg_cas_1 / hallmarks_cas_1 +
  plot_layout(heights = c(10, 3), guides = "collect") &
  theme(legend.position = "bottom")


#
# GO ENRICHMENT OF IMMUNE CASSETTE
#

# Cassette 1
go_df <- list_of_go_promoter[[9]]
go_df <- go_df[go_df$FDR <= 0.05,]

go_df <- go_df[order(go_df$FDR, decreasing = TRUE),]
go_df$Description <- factor(go_df$TERM, levels = go_df$TERM)


ggplot(go_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 3) +
  scale_size_continuous(limits = c(0, 0.2), breaks = c(0, 0.1, 0.2)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "GO terms",
    size = "Gene Ratio"
  )

go_df <- list_of_go_promoter[[1]]
go_df <- go_df[go_df$FDR <= 0.05,]
nrow(go_df)

go_df <- go_df[order(go_df$FDR, decreasing = F),]

write_csv(go_df, "/Users/in2245sa/PhD/Projects/project_3/analysis/go_enrichment_promoter_cassette_1.csv")

go_df <- list_of_go_proximal[[1]]
go_df <- go_df[go_df$FDR <= 0.05,]
nrow(go_df)

go_df <- go_df[order(go_df$FDR, decreasing = F),]

write_csv(go_df, "/Users/in2245sa/PhD/Projects/project_3/analysis/go_enrichment_proximal_cassette_1.csv")
