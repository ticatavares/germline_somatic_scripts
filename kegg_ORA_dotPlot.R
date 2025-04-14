###############################################################################
# GO/KEGG Enrichment Dot Plot Visualization
# Author: Thais Tavares
# Date: 11/04/2025
# Description: Processes enrichment results and creates publication-ready dot plots
###############################################################################

# --------------------------
# 1. Setup and Configuration
# --------------------------

# Load required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(ggplot2)
library(dplyr)

# --------------------------
# 2. Data Processing Functions
# --------------------------

#' Process enrichment result files for dot plot visualization
#'
#' @param file_path Path to enrichment result file
#' @return Dataframe with extracted enrichment metrics
process_enrichment_file <- function(file_path) {
  # Read header information
  header_info <- readLines(file_path, n = 2)
  
  # Extract GO/KEGG term information
  term_id <- sub('.*ID: (.*); Name:.*', '\\1', header_info[1])
  term_name <- sub('.*Name: (.*)', '\\1', header_info[1])
  
  # Extract enrichment statistics
  enrichment_ratio <- as.numeric(
    sub('.*enrichmentRatio=([0-9\\.eE\\-]+);.*', '\\1', header_info[2])
  )
  fdr <- as.numeric(
    sub('.*FDR=([0-9\\.eE\\-]+).*', '\\1', header_info[2])
  )
  
  # Count genes in the pathway
  gene_data <- read.csv(file_path, skip = 2)
  gene_count <- nrow(gene_data)
  
  # Return structured data
  return(data.frame(
    Term = term_name,
    Term_ID = term_id,
    Enrichment_Ratio = enrichment_ratio,
    FDR = fdr,
    Gene_Count = gene_count,
    stringsAsFactors = FALSE
  ))
}

# --------------------------
# 3. Data Loading and Processing
# --------------------------

# Process all enrichment files (modify path as needed)
enrichment_files <- list.files(
  path = "data/enrichment_results/ORA_analysis",
  pattern = "*.csv",
  full.names = TRUE
)

# Combine all enrichment results
enrichment_data <- do.call(
  rbind,
  lapply(enrichment_files, process_enrichment_file)
)

# Filter and sort results (optional)
enrichment_data <- enrichment_data %>%
  filter(FDR < 0.05) %>%          # Keep significant results
  arrange(desc(Enrichment_Ratio)) # Sort by enrichment

# --------------------------
# 4. Visualization
# --------------------------

# Create dot plot visualization
enrichment_dotplot <- ggplot(
  enrichment_data,
  aes(
    x = Enrichment_Ratio,
    y = reorder(Term, Enrichment_Ratio),
    size = Gene_Count,
    color = -log10(FDR)
  )
) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(
    low = "blue",
    high = "red",
    name = "-log10(FDR)"
  ) +
  scale_size_continuous(
    range = c(3, 8),
    name = "Gene Count"
  ) +
  labs(
    x = "Enrichment Ratio",
    y = NULL,
    title = "Pathway Enrichment Analysis"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

# --------------------------
# 5. Save Outputs
# --------------------------

# Create output directory
output_dir <- "results/enrichment_analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save plots
ggsave(
  filename = paste0(output_dir, "enrichment_dotplot.pdf"),
  plot = enrichment_dotplot,
  width = 8,
  height = 6
)

ggsave(
  filename = paste0(output_dir, "enrichment_dotplot.png"),
  plot = enrichment_dotplot,
  width = 8,
  height = 6,
  dpi = 300
)

# Save processed data
write.table(
  enrichment_data,
  paste0(output_dir, "enrichment_results.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# --------------------------
# 6. Session Information
# --------------------------

sink(paste0(output_dir, "session_info.txt"))
sessionInfo()
sink()
