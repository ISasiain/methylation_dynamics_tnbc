#! usr/bin/Rscript

library(missMethyl)

#
# LOADING DATA
#

promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")
proximal_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/proximal/cassettes_beta_10.rds")

hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))

#
# PERFORMING GO ENRICHMENT ANALYSIS PER CASSETTE
#

list_of_go <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:60) {

  # Getting CpGs from cassette
  cpg_subset <- names(proximal_10$colors)[proximal_10$colors == cassette]

  # Performing GO enrichment
  list_of_go[[cassette]] <- gometh(
    cpg_subset,
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
saveRDS(list_of_go, "PhD/Projects/project_3/analysis/GO_enrichment_proximal_10.rds")


#
# PERFORMING HALLMARK ENRICHMENT ANALYSIS PER CASSETTE
#

list_of_hallmarks <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:60) {

  # Getting CpGs from cassette
  cpg_subset <- names(proximal_10$colors)[proximal_10$colors == cassette]

  # Performing GO enrichment
  list_of_hallmarks[[cassette]] <- gsameth(
    cpg_subset,
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
saveRDS(list_of_hallmarks, "PhD/Projects/project_3/analysis/hallmark_enrichment_proximal_10.rds")


#
# PERFORMING KEGG ENRICHMENT ANALYSIS PER CASSETTE
#

list_of_kegg <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:60) {

  # Getting CpGs from cassette
  cpg_subset <- names(proximal_10$colors)[proximal_10$colors %in% c(1,3)]

  # Performing GO enrichment
  list_of_kegg[[cassette]] <- gometh(
    cpg_subset,
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
saveRDS(list_of_kegg, "PhD/Projects/project_3/analysis/kegg_enrichment_proximal_10.rds")


#
# READING FILES
#

list_of_kegg_promoter <- readRDS("PhD/Projects/project_3/analysis/kegg_enrichment_promoter_10.rds")
list_of_kegg_proximal <- readRDS("PhD/Projects/project_3/analysis/kegg_enrichment_proximal_10.rds")

list_of_hallmarks_promoter <- readRDS("PhD/Projects/project_3/analysis/hallmark_enrichment_promoter_10.rds")
list_of_hallmarks_proximal <- readRDS("PhD/Projects/project_3/analysis/hallmark_enrichment_proximal_10.rds")

#
# PLOTTING
#


# KEGG

# Cassette 1
kegg_df <- list_of_kegg_proximal[[3]] 
kegg_df <- kegg_df[kegg_df$FDR <= 0.05,]

kegg_df <- kegg_df[order(kegg_df$FDR, decreasing = TRUE),]
kegg_df$Description <- factor(kegg_df$Description, levels = kegg_df$Description)


kegg_cas_1 <- ggplot(kegg_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 6) +
  scale_size_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "KEGG Pathways",
    size = "Gene Ratio"
  )



# HALLMARKS

# Cassette 1
hallmarks_df <- list_of_hallmarks_proximal[[3]]
hallmarks_df <- hallmarks_df[hallmarks_df$FDR <= 0.05,]

formatted_rownames <- tolower(gsub("_", " ", gsub("^HALLMARK_", "", rownames(hallmarks_df))))
formatted_rownames <- sub("^([a-z])", "\\U\\1", formatted_rownames, perl = TRUE)

hallmarks_df$Description <- formatted_rownames

hallmarks_df <- hallmarks_df[order(hallmarks_df$FDR, decreasing = TRUE),]
hallmarks_df$Description <- factor(hallmarks_df$Description, levels = hallmarks_df$Description)


hallmarks_cas_1 <- ggplot(hallmarks_df, aes(x = -log10(FDR), y = Description)) +
  geom_point(aes(size = DE/N), color = "black") +
  xlim(1, 6) +
  scale_size_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "Hallmarks",
    size = "Gene Ratio"
  )

kegg_cas_1 / hallmarks_cas_1 +
  plot_layout(heights = c(10, 3), guides = "collect") &
  theme(legend.position = "bottom")


# GO

# Cassette 1
GO_df <- list_of_go[[3]]
GO_df <- GO_df[GO_df$FDR <= 0.05,]

# Capitalize first letter
GO_df$TERM <- sub("^([a-z])", "\\U\\1", GO_df$TERM, perl = TRUE)

# Sorting
GO_df <- GO_df[order(GO_df$FDR, decreasing = TRUE),]
GO_df$TERM <- factor(GO_df$TERM, levels = GO_df$TERM)


ggplot(GO_df, aes(x = -log10(FDR), y = TERM)) +
  geom_point(aes(size = DE/N), color = "black") +
  scale_size_continuous(limits = c(0, 0.4), breaks = c(0, 0.1, 0.2, 0.3, 0.4)) +
  theme_bw() +
  labs(
    x = "-log10(FDR p value)",
    y = "GO terms",
    size = "Gene Ratio"
  )

