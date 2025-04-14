###############################################################################
# Fisher's Exact Test for Genetic Variant Comparisons
# Author: Thais Tavares
# Date: 09/04/2025
# Description: Performs Fisher's exact tests comparing variant sets with adjusted
#              probability expectations (25% equal, 75% different alleles)
###############################################################################

# --------------------------
# 1. Package Installation
# --------------------------

# Check and install required packages
if (!requireNamespace("stats", quietly = TRUE)) {
  install.packages("stats")
}

# --------------------------
# 2. Data Preparation
# --------------------------

# Define variant set comparisons with counts:
# - set1_total: Total variants in set 1
# - set2_total: Total variants in set 2
# - intersection: Variants present in both sets
# - observed_equal: Variants with identical alleles
comparisons <- list(
  "SNPs vs COSMIC" = list(
    set1_total = 16040,
    set2_total = 3662887,
    intersection = 9575,
    observed_equal = 9358
  ),
  "SNPs vs Benign" = list(
    set1_total = 16040,
    set2_total = 56794,
    intersection = 2728,
    observed_equal = 2728
  ),
  "Rares vs COSMIC" = list(
    set1_total = 496791,
    set2_total = 3662887,
    intersection = 126723,
    observed_equal = 97619
  ),
  "Rares vs Benign" = list(
    set1_total = 496791,
    set2_total = 56794,
    intersection = 17777,
    observed_equal = 17454
  ),
  "Benign vs Pathogenic" = list(
    set1_total = 56794,
    set2_total = 13098,
    intersection = 192,
    observed_equal = 0
  )
)

# --------------------------
# 3. Analysis Function
# --------------------------

#' Perform Fisher's exact test with adjusted probability expectations
#' 
#' @param data List containing:
#'   - intersection: Total overlapping variants
#'   - observed_equal: Variants with identical alleles
#' @return List with test results including:
#'   - Observed and expected counts
#'   - Odds ratio
#'   - p-value
calculate_fisher_with_probabilities <- function(data) {
  # Calculate observed different alleles
  observed_equal <- data$observed_equal
  observed_different <- data$intersection - observed_equal
  
  # Define expected probabilities (25% equal, 75% different)
  prob_equal <- 0.25
  prob_different <- 0.75
  
  # Calculate expected counts
  expected_equal <- data$intersection * prob_equal
  expected_different <- data$intersection * prob_different
  
  # Perform Fisher's exact test
  contingency_table <- matrix(
    c(observed_equal, observed_different,
      expected_equal, expected_different),
    nrow = 2
  )
  
  fisher_result <- fisher.test(contingency_table)
  
  # Return comprehensive results
  return(list(
    observed_equal = observed_equal,
    observed_different = observed_different,
    expected_equal = round(expected_equal, 2),
    expected_different = round(expected_different, 2),
    odds_ratio = round(fisher_result$estimate, 2),
    p_value = format.pval(fisher_result$p.value, digits = 5)
  ))
}

# --------------------------
# 4. Run Analyses
# --------------------------

# Apply analysis to all comparisons
results <- lapply(comparisons, calculate_fisher_with_probabilities)

# Convert results to data frame
results_df <- do.call(rbind, lapply(names(results), function(name) {
  data.frame(
    Comparison = name,
    Observed_Equal = results[[name]]$observed_equal,
    Observed_Different = results[[name]]$observed_different,
    Expected_Equal = results[[name]]$expected_equal,
    Expected_Different = results[[name]]$expected_different,
    Odds_Ratio = results[[name]]$odds_ratio,
    P_Value = results[[name]]$p_value,
    stringsAsFactors = FALSE
  )
}))

# --------------------------
# 5. Output Results
# --------------------------

# Print formatted results to console
cat("Fisher's Exact Test Results with Adjusted Probabilities:\n")
print(results_df, row.names = FALSE)

# Save results to CSV
write.csv(results_df, 
          "fisher_test_results_with_probabilities.csv", 
          row.names = FALSE)

# --------------------------
# 6. Session Information
# --------------------------

# Save session info for reproducibility
sink("session_info.txt")
sessionInfo()
sink()