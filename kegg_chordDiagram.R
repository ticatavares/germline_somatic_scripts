###############################################################################
# KEGG Pathway Chord Diagram Visualization
# Author: Thais Tavares
# Date: 14/04/2025
# Description: Creates a chord diagram showing gene-KEGG pathway associations
###############################################################################
# 1. Load Packages ------------------------------------------------------------
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")
if (!requireNamespace("ggsci", quietly = TRUE)) install.packages("ggsci")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")

library(dplyr)
library(circlize)
library(ggsci)
library(viridis)

# 2. Data Loading ------------------------------------------------------------
kegg_dir <- "yourpath/"
kegg_files <- list.files(kegg_dir, pattern = "*.csv", full.names = TRUE)

# 3. Data Processing ---------------------------------------------------------
process_kegg_file <- function(file_path) {
  header <- readLines(file_path, n = 2)
  pathway_id <- sub('.*ID: (.*); Name:.*', '\\1', header[1])
  pathway_name <- sub('.*Name: (.*)', '\\1', header[1])
  
  genes <- read.csv(file_path, skip = 2) %>% 
    select(Gene.Symbol) %>%
    mutate(Pathway_Term = pathway_name,
           Pathway_ID = pathway_id)
  
  return(genes)
}

# Create connection_data properly
kegg_data <- lapply(kegg_files, process_kegg_file) %>% 
  bind_rows() %>%
  filter(!is.na(Gene.Symbol))  # Remove NA genes

connection_data <- kegg_data %>%
  group_by(Gene.Symbol) %>%
  filter(n() > 2) %>%          # Genes in >2 pathways
  ungroup() %>%
  select(Pathway_Term, Gene.Symbol) %>%
  distinct()                   # Remove duplicates

# 4. Visualization Setup -----------------------------------------------------
pathway_colors <- pal_igv()(length(unique(connection_data$Pathway_Term)))
names(pathway_colors) <- unique(connection_data$Pathway_Term)

gene_colors <- viridis(length(unique(connection_data$Gene.Symbol)), option = "inferno")
names(gene_colors) <- unique(connection_data$Gene.Symbol)

all_colors <- c(pathway_colors, gene_colors)

# 5. Create Chord Diagram ----------------------------------------------------
output_dir <- "results/kegg_analysis/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#pdf(paste0(output_dir, "kegg_chord.pdf"), width = 10, height = 8)
layout(matrix(c(1, 2), ncol = 2), widths = c(8, 2))

# Main plot
par(mar = c(0, 0, 0, 0))
circos.clear()

# Check connection_data exists
if (exists("connection_data") && nrow(connection_data) > 0) {
  chordDiagram(connection_data,
               grid.col = all_colors,
               transparency = 0.3,
               annotationTrack = "grid",
               preAllocateTracks = 1,
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow")
  
  # Add gene labels
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.index <- CELL_META$sector.index
      if (sector.index %in% connection_data$Gene.Symbol) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], 
                    sector.index, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
      }
    },
    bg.border = NA
  )
} else {
  plot.new()
  text(0.5, 0.5, "No significant gene-pathway connections found", cex = 1.2)
}

# Legend
par(mar = c(0, 2, 0, 0))
plot.new()
legend("center", legend = names(pathway_colors), 
       fill = pathway_colors, title = "KEGG Pathways",
       cex = 0.8, ncol = 1, bty = "n")

#dev.off()

# 6. Save Data --------------------------------------------------------------
write.table(kegg_data, paste0(output_dir, "kegg_connections.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# 7. Session Info -----------------------------------------------------------
sink(paste0(output_dir, "session_info.txt"))
sessionInfo()
sink()