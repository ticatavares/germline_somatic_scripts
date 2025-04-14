###############################################################################
# AlphaMissense Pathogenicity Score Analysis
# Author: Thais Tavares
# Date: 11/04/2025
# Description: Analyzes and visualizes AlphaMissense pathogenicity scores
#              across multiple variant datasets (Cosmic, SNP, Rares, Benign, Pathogenic)
###############################################################################

# --------------------------
# 1. Setup and Configuration
# --------------------------

# Set working directory (modify as needed)
setwd('/path/to/your/data')

# Load required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggridges", quietly = TRUE)) {
  install.packages("ggridges")
}
library(ggplot2)
library(ggridges)

# --------------------------
# 2. Data Loading
# --------------------------

# Load AlphaMissense pathogenicity scores for each dataset
cosmic <- read.table(
  "cosmic_CDS_am_pathog.bed.gz",
  header = FALSE,
  col.names = c("Gene", "V2", "V3", "V4", "pathogenity_score")
)

benign <- read.table(
  "benign_CDS_am_pathog.bed.gz",
  header = FALSE,
  col.names = c("Gene", "V2", "V3", "V4", "pathogenity_score")
)

pathog <- read.table(
  "patho_CDS_am_pathog.bed.gz",
  header = FALSE,
  col.names = c("Gene", "V2", "V3", "V4", "pathogenity_score")
)

snp <- read.table(
  "SNP1_subt_rares2_amPatho.bed.gz",
  header = FALSE,
  col.names = c("Gene", "V2", "V3", "V4", "pathogenity_score")
)

rares <- read.table(
  "rares_total_CDS_am_pathog.bed.gz",
  header = FALSE,
  col.names = c("Gene", "V2", "V3", "V4", "pathogenity_score")
)

# --------------------------
# 3. Data Preparation
# --------------------------

# Combine all datasets with appropriate labels
combined_data <- rbind(
  transform(cosmic, Dataset = "Cosmic"),
  transform(benign, Dataset = "Benign"),
  transform(pathog, Dataset = "Pathogenic"),
  transform(snp, Dataset = "SNP"),
  transform(rares, Dataset = "Rares")
)

# Set factor levels for consistent ordering
combined_data$Dataset <- factor(
  combined_data$Dataset,
  levels = c("Cosmic", "SNP", "Rares", "Benign", "Pathogenic")
)

# Define color palette
color_palette <- c(
  "#42454C",  # Cosmic
  "#00B4F5",  # SNP
  "#d45500ff", # Rares
  "#74A02C",  # Benign
  "#8D3B72"   # Pathogenic
)

# --------------------------
# 4. Visualization
# --------------------------

# Create ridgeline density plot
pathogenicity_plot <- ggplot(
  combined_data,
  aes(
    x = pathogenity_score,
    y = Dataset,
    fill = Dataset,
    color = Dataset
  )
) +
  geom_density_ridges(
    alpha = 0.5,
    rel_min_height = 0.001  # Adjusts curve smoothness
  ) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  labs(
    x = "AlphaMissense Pathogenicity Score",
    y = "Dataset",
    title = "Pathogenicity Score Distribution by Dataset"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    legend.position = "none",  # Remove legend (redundant with y-axis)
    panel.grid.major.y = element_line(
      colour = "#707070",
      size = 0.15
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(
      size = 13,
      margin = margin(t = 0, r = 5, b = 0, l = 0)
    ),
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 10, r = 0, b = 0, l = 0)
    ),
    plot.margin = margin(0.8, 1, 0.8, 1, "cm")
  )

# --------------------------
# 5. Save Outputs
# --------------------------

# Create output directory if it doesn't exist
output_dir <- "results/pathogenicity_analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save plots
ggsave(
  filename = paste0(output_dir, "pathogenicity_ridgeline.pdf"),
  plot = pathogenicity_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = paste0(output_dir, "pathogenicity_ridgeline.png"),
  plot = pathogenicity_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# --------------------------
# 6. Session Information
# --------------------------

sink(paste0(output_dir, "session_info.txt"))
sessionInfo()
sink()
  
