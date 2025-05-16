#! usr/bin/Rscript

library(missMethyl)

#
# LOADING DATA
#

promoter_10 <- readRDS("/Volumes/Data/Project_3/detected_cassettes/promoter/cassettes_beta_10.rds")

#
# PERFORMING GO ENRICHMENT ANALYSIS PER CASSETTE
#

list_of_go <- list()

# Performing GO enrichment per cassette CpGs. Using first 140 CpGs (More than 10 CpGs)
for (cassette in 1:140) {
  
  # Getting CpGs from cassette
  cpg_subset <- names(promoter_10$colors)[promoter_10$colors == cassette]
  
  # Performing GO enrichment
  list_of_go[[cassette]] <- gometh(
    cpg_subset,
    collection = "GO",
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
saveRDS(list_of_go, "PhD/Projects/project_3/analysis/GO_enrichment_promoter_10.rds")


# Example for cassette 10
go_df <- list_of_go[[10]] %>%
  filter(FDR < 0.05) %>%
  mutate(term = reorder(TERM, -log10(FDR))) %>%
  slice_max(order_by = -log10(FDR), n = 15)

ggplot(go_df, aes(x = -log10(FDR), y = term)) +
  geom_point(aes(size = DE/N), color = "steelblue") +
  theme_minimal() +
  labs(
    x = "-log10(FDR)",
    y = "GO Term",
    size = "Gene Ratio"
  )

