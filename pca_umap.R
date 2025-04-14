###############################################################################
# UMAP Analysis of Population Genetic Structure
# Author: Thais Tavares
# Date: 10/04/2025
# Description: Performs UMAP dimensionality reduction on PCA results (PCA-UMAP) from PLINK
#              to visualize population genetic structure in the HGDP dataset
###############################################################################

# --------------------------
# 1. Setup and Configuration
# --------------------------

# Set working directory (user should modify)
setwd("/path/to/your/data")

# Load required packages
required_packages <- c("umap", "ggplot2", "dplyr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# --------------------------
# 2. Data Loading and Preparation
# --------------------------

# 2.1 Load PCA results from PLINK (.eigenvec file)
eigenvec <- read.table(
  "PCA_combined_hgdp_biall_pass_noTrans_geno_mInd_LDin_pruned.eigenvec",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Rename columns appropriately
colnames(eigenvec) <- c("FID", "Sample_ID", paste0("PC", 1:(ncol(eigenvec)-2)))

# Set sample IDs as row names and remove unnecessary columns
rownames(eigenvec) <- eigenvec$Sample_ID
eigenvec <- eigenvec[, -c(1, 2)]  # Remove FID and Sample_ID columns

# 2.2 Load eigenvalues
eigenval <- scan(
  "PCA_combined_hgdp_biall_pass_noTrans_geno_mInd_LDin_pruned.eigenval",
  quiet = TRUE
)

# --------------------------
# 3. UMAP Dimensionality Reduction
# --------------------------

# 3.1 Prepare PCA matrix (PLINK already normalizes the scores)
pca_matrix <- as.matrix(eigenvec)

# 3.2 Perform UMAP with reproducible parameters
set.seed(20)  # For reproducibility
umap_result <- umap(
  pca_matrix,
  n_neighbors = 50,    # Balances local vs global structure
  min_dist = 0.3,      # Controls cluster tightness
  spread = 2,          # Controls separation between clusters
  metric = "euclidean" # Standard distance metric for genetic data
)

# 3.3 Create UMAP coordinates dataframe
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# --------------------------
# 4. Population Metadata Integration
# --------------------------

# 4.1 Load population information
sample_to_pop <- read.csv(
  "samplesPopulatins_hgdp_limpo.csv",
  sep = "\t",
  stringsAsFactors = FALSE
)

# 4.2 Match population labels to samples
umap_df$Population <- sample_to_pop$Superpopulation.name[
  match(rownames(umap_df), sample_to_pop$Sample.name)
]

# 4.3 Remove samples with missing population information
umap_df <- na.omit(umap_df)

# --------------------------
# 5. Visualization
# --------------------------

# 5.1 Create UMAP plot
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Population)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2", name = "Superpopulation") +
  labs(
    title = "UMAP Visualization of Population Genetic Structure",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

# 5.2 Save plot to file
ggsave(
  "population_structure_umap.pdf",
  plot = umap_plot,
  width = 8,
  height = 6,
  device = "pdf"
)

# 5.3 Display plot
print(umap_plot)

# --------------------------
# 6. Save Results
# --------------------------

# Save UMAP coordinates with population information
write.csv(
  umap_df,
  "umap_coordinates_with_population.csv",
  row.names = TRUE
)

# Save session information for reproducibility
sink("umap_analysis_session_info.txt")
sessionInfo()
sink()