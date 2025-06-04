#! usr/bin/Rscript

library(Boruta)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)


#
# LOADING DATA 
#

# Load annotation file
load("/Users/isasiain/PhD/Projects/immune_spatial/ecosystem_analysis/data/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
rownames(x) <- x$PD_ID


# Load cassettes
prom_cassettes <- read.csv("/Users/isasiain/PhD/Projects/project_3/analysis/promoter_cassettes/summary_of_cassettes/summary_beta_10.csv")
rownames(prom_cassettes) <- prom_cassettes$Cassette
prom_cassettes$Cassette <- NULL
prom_cassettes <- prom_cassettes[, x$PD_ID]


#
# DEFINING IMPORTANT CASSETTES USING BORUTA ALGORITHM
#

# Getting data to run boruta
X <- t(prom_cassettes)
tils <- x$TILs

# Remove NAs (from TILs)
to_remove <- is.na(tils)

X <- X[!to_remove,]
tils <- tils[!to_remove]

# Combine into data frame
df_boruta <- data.frame(X)
df_boruta$TILs <- tils

# Run Boruta for feture selection
set.seed(123)
boruta_result <- Boruta(TILs ~ ., data = df_boruta, doTrace = 2, maxRuns = 500)

# Check results
print(boruta_result)

# Plot boruta results
plot(boruta_result, las = 2, cex.axis = 0.7)

# Get confirmed attributes
confirmed_features <- getSelectedAttributes(boruta_result, withTentative = FALSE)
cat("Confirmed important features:\n")
print(confirmed_features)

# Get feature importance scores
feature_importance <- attStats(boruta_result)


# Plot importance of each selected variable

# 1. Extract variable importance and convert to a tidy data frame
importance_df <- data.frame("Median_Importance" = feature_importance[confirmed_features, "medianImp"],
                            "Cassette" = rownames(feature_importance[confirmed_features, ]))

importance_df$Cassette <- substring(importance_df$Cassette, 2)

importance_df <- tibble(
  variable = importance_df$Cassette,
  importance = as.numeric(importance_df$Median_Importance)
)

# 2. Extract Kendall Tau named vector and convert to data frame
kendall_tau_vec <- results_sorted[, "Kendall_Tau"] 
names(kendall_tau_vec) <- rownames(results_sorted)
kendall_df <- tibble(
  variable = names(kendall_tau_vec),
  kendall_tau = as.numeric(kendall_tau_vec)
)

# 3. Join importance and Kendall Tau
merged_df <- importance_df %>%
  left_join(kendall_df, by = "variable")

# 4. Sort op importance
top_vars <- merged_df %>%
  arrange(desc(importance)) %>%
  mutate(variable = factor(variable, levels = rev(variable)))

# 5. Dot plot (left)
p1 <- ggplot(top_vars, aes(x = kendall_tau, y = variable)) +
  geom_point(aes(
    size = abs(kendall_tau),
    color = "red"
  )) +
  scale_size(range = c(0.5, 3)) +
  scale_color_identity() +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(-0.4, 0.4) +
  xlab("Kendall Tau")

# 6. Bar plot (right)
p2 <- ggplot(top_vars, aes(x = importance, y = variable)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  theme(
    axis.title.y = element_blank()
  ) +
  xlab("Boruta Median Importance")

# 7. Combine plots
p1 + p2 + plot_layout(ncol = 2, widths = c(1, 2))
