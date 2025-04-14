###############################################################################
# Genomic Consequences Analysis of Variant Datasets
# Author: Thais Tavares
# Date: 11/04/2025
# Description: Analyzes and visualizes variant consequences across multiple 
#              genomic datasets (SNPs, rare variants, COSMIC, benign, pathogenic)
###############################################################################

# --------------------------
# 1. Setup and Configuration
# --------------------------

# Set working directory (modify paths as needed)
# setwd("/path/to/raw/data")

# Load required packages
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --------------------------
# 2. Define Genomic Categories
# --------------------------

# Categorize VEP consequences into functional groups
consequence_categories <- list(
  Splicing_site = c(
    'splice_donor_variant', 'splice_acceptor_variant', 
    'splice_region_variant', 'splice_donor_5th_base',
    'splice_donor_region', 'splice_polypyrimidine_tract',
    'splice_donor_5th_base_variant', 'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant'
  ),
  Intron = c('intron_variant'),
  Intergenic = c(
    'regulatory_region', 'TF_binding_site',
    'TF_binding_site_variant', 'intergenic_variant',
    'upstream_gene_variant', 'downstream_gene_variant',
    'regulatory_region_variant'
  ),
  Exon = c(
    'synonymous_variant', 'missense_variant',
    'inframe_insertion', 'inframe_deletion',
    'stop_gained', 'frameshift_variant',
    'coding_sequence_variant', 'start_retained_variant',
    'start_lost', 'stop_lost',
    'stop_retained_variant', 'incomplete_terminal_codon_variant',
    'non_coding_transcript_exon_variant'
  ),
  UTR = c(
    '5_prime_UTR_variant',
    '3_prime_UTR_variant'
  ),
  Others = c(
    'transcript_ablation', 'transcript_amplification',
    'protein_altering_variant', 'mature_miRNA_variant',
    'NMD_transcript_variant', 'non_coding',
    'transcript_exon_variant', 'TFBS_ablation/amplification',
    'feature_elongation', 'feature_truncation',
    'coding_transcript', 'variant',
    'sequence_variant'
  )
)

# --------------------------
# 3. Data Processing Functions
# --------------------------

#' Aggregate variant counts by consequence categories
#'
#' @param df Input dataframe with Count and Consequence_type columns
#' @param categories List of consequence categories
#' @return Aggregated dataframe with summed counts by category
aggregate_consequences <- function(df, categories) {
  result <- data.frame(
    Consequence_type = character(),
    Count = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (category_name in names(categories)) {
    category_terms <- categories[[category_name]]
    subset_df <- df %>% filter(Consequence_type %in% category_terms)
    category_count <- sum(subset_df$Count)
    
    result <- rbind(result, data.frame(
      Consequence_type = category_name,
      Count = category_count
    ))
  }
  
  return(result)
}

# --------------------------
# 4. Load and Process Datasets
# --------------------------

# 4.1 Load intersection data
intersect_data <- read.table(
  "SNP1_subtRares2_SNP2inters_WG_severe.txt",
  header = FALSE,
  col.names = c("Count", "Consequence_type")
)

# Calculate percentages
total_intersect <- sum(intersect_data$Count)
intersect_data$Percent <- (intersect_data$Count / total_intersect) * 100

# 4.2 Process each variant dataset
process_dataset <- function(filename, dataset_name, categories) {
  df <- read.table(
    filename,
    header = FALSE,
    col.names = c("Count", "Consequence_type")
  )
  
  df <- aggregate_consequences(df, categories)
  df$Dataset <- dataset_name
  df$Percent <- (df$Count / sum(df$Count)) * 100
  return(df)
}

# Process all datasets
rare_variants <- process_dataset("rares_total_CDS_severe.txt", "Rares", consequence_categories)
snp_variants <- process_dataset("SNP1_CDS_severe.txt", "SNP", consequence_categories)
cosmic_variants <- process_dataset("cosmic_CDS_severe.txt", "Cosmic", consequence_categories)
benign_variants <- process_dataset("benign_CDS_severe.txt", "Benign", consequence_categories)
pathogenic_variants <- process_dataset("patho_CDS_severe.txt", "Pathogenic", consequence_categories)

# 4.3 Combine all datasets
combined_data <- bind_rows(
  rare_variants,
  snp_variants,
  cosmic_variants,
  benign_variants,
  pathogenic_variants
)

# Set factor levels for consistent ordering
combined_data$Consequence_type <- factor(
  combined_data$Consequence_type,
  levels = names(consequence_categories)
)

combined_data$Dataset <- factor(
  combined_data$Dataset,
  levels = c("Cosmic", "SNP", "Rares", "Benign", "Pathogenic")
)

# --------------------------
# 5. Visualization
# --------------------------

# 5.1 Create stacked bar plot
consequence_plot <- ggplot(
  combined_data,
  aes(x = Dataset, y = Percent, fill = Consequence_type)
) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +
  coord_flip() +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(
    values = c(
      "#42454C",  # Splicing
      "#00B4F5",  # Intron
      "#d45500ff", # Intergenic
      "#74A02C",  # Exon
      "#8D3B72",  # UTR
      "#2935A3"   # Others
    ),
    name = "Genomic Region"
  ) +
  labs(
    title = "Variant Consequences by Dataset - CDS",
    x = "Dataset",
    y = "Percentage of Variants"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    axis.text.y = element_text(size = 11),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(size = 18, face = "bold"),
    panel.grid.major.y = element_line(colour = "#707070", size = 0.15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 13, margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.title.x = element_text(size = 13, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    plot.margin = margin(0.8, 1, 0.8, 1, "cm")
  )

# 5.2 Save plots
output_dir <- "/home/thais/Documents/UBUNTU/phd/final_datasets/bed_Ensembl/reduced_CDS/images_CDS/imagens_v7/"

ggsave(
  filename = paste0(output_dir, "genomicRegions_CDS.pdf"),
  plot = consequence_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = paste0(output_dir, "genomicRegions_CDS.png"),
  plot = consequence_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# --------------------------
# 6. Save Processed Data
# --------------------------

write.csv(
  combined_data,
  paste0(output_dir, "CDS_mostSevereVEP.csv"),
  row.names = FALSE
)

write.table(
  intersect_data,
  paste0(output_dir, "SNPintersCOSMIC_mostSevereVEP.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# --------------------------
# 7. Session Information
# --------------------------

sink(paste0(output_dir, "session_info.txt"))
sessionInfo()
sink()