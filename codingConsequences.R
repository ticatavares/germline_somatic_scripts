# Set working directory
#setwd("/home/thais/link_chrom_file/paper_figures_dez2023/arquivosanotacaovep")
setwd("/home/thais/Documents/UBUNTU/phd/final_datasets/bed_Ensembl/reduced_CDS")

# Import packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)


#FOR CDS
categorias <- list(
  SG = 'stop_gained',
  MISS = 'missense_variant',
  SYN = 'synonymous_variant',
  LOST = c('start_lost','stop_lost'),
  RET = 'stop_retained_variant',
  SPL = 'splice_region_variant'
  
)

# Função para agrupar e somar as categorias
agrupar_categorias <- function(df, categorias) {
  resultado <- data.frame(Consequence_type = character(), Count = numeric(), stringsAsFactors = FALSE)
  
  for (nome_categoria in names(categorias)) {
    categoria <- categorias[[nome_categoria]]
    subset_df <- subset(df, Consequence_type %in% categoria)
    count_sum <- sum(subset_df$Count)
    
    nova_linha <- data.frame(Consequence_type = nome_categoria, Count = count_sum)
    resultado <- rbind(resultado, nova_linha)
  }
  
  return(resultado)
}


# CDS:
df_raras <- read.table("rares_total_CDS_severe.txt", header = FALSE, col.names = c("Count","Consequence_type"))
# Aplicando a função ao DataFrame_1
df_raras <- agrupar_categorias(df_raras, categorias)
# Calcular o total de variantes
df_raras$Dataset <- 'Rares'
total_variantes1 <- sum(df_raras$Count)
# Calcular a porcentagem de cada linha em relação ao total_variantes
df_raras$Percent <- (df_raras$Count / total_variantes1) * 100

# CDS
df_SNP <- read.table("SNP1_CDS_severe.txt", header = FALSE, col.names = c("Count","Consequence_type"))
df_SNP <- agrupar_categorias(df_SNP, categorias)
df_SNP$Dataset <- 'SNP'
total_variantes1 <- sum(df_SNP$Count)
df_SNP$Percent <- (df_SNP$Count / total_variantes1) * 100

# CDS
df_cosmic <- read.table("cosmic_CDS_severe.txt", header = FALSE, col.names = c("Count","Consequence_type"))
df_cosmic <- agrupar_categorias(df_cosmic, categorias)
df_cosmic$Dataset <- 'Cosmic'
total_variantes1 <- sum(df_cosmic$Count)
df_cosmic$Percent <- (df_cosmic$Count / total_variantes1) * 100

# CDS
df_benign <- read.table("benign_CDS_severe.txt", header = FALSE, col.names = c("Count","Consequence_type"))
df_benign <- agrupar_categorias(df_benign, categorias)
df_benign$Dataset <- 'Benign'
total_variantes1 <- sum(df_benign$Count)
df_benign$Percent <- (df_benign$Count / total_variantes1) * 100

# CDS
df_pathog <- read.table("patho_CDS_severe.txt", header = FALSE, col.names = c("Count","Consequence_type"))
df_pathog <- agrupar_categorias(df_pathog, categorias)
df_pathog$Dataset <- 'Pathogenic'
total_variantes1 <- sum(df_pathog$Count)
df_pathog$Percent <- (df_pathog$Count / total_variantes1) * 100


# Combinar todos os DataFrames em um único DataFrame
all_data <- rbind(df_raras, df_SNP, df_cosmic, df_benign, df_pathog)

# Garantir que os tipos de dados estão corretos para plotagem
all_data$Consequence_type <- factor(all_data$Consequence_type, levels = names(categorias))
# Reordenar o fator Dataset de acordo com a ordem desejada
all_data$Dataset <- factor(all_data$Dataset, levels = c("Cosmic", "SNP", "Rares", "Benign","Pathogenic"))


# Criar o gráfico de barras horizontais
pdf("/home/thais/Documents/UBUNTU/phd/final_datasets/bed_Ensembl/reduced_CDS/images_CDS/imagens_v7/codingConsequence_CDS.pdf", height = 5, width = 8)
plt1 = ggplot(all_data, aes(x = Consequence_type, 
                     y = Percent, 
                     fill = Dataset)) +
  geom_bar(stat = "identity", 
           position = "dodge",  
           width = 0.9) +  
  facet_wrap(~ Consequence_type, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c("#42454C", "#00B4F5", "#d45500ff", "#74A02C", "#8D3B72", "#2935A3")) +
  labs(y = "% of SNVs",
       x = NULL,
       fill = "Dataset") +
  theme_minimal() +
  ylim(0, 100) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12, ),
        axis.text.y = element_text(size = 11),
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
dev.off()




ggsave(filename = "/home/thais/Documents/UBUNTU/phd/final_datasets/bed_Ensembl/reduced_CDS/images_CDS/imagens_v7/codingconsequence_CDS.png", 
       plot = plt1, 
       width = 8, 
       height = 5, 
       dpi = 300)

write.csv(all_data, "/home/thais/Documents/UBUNTU/phd/final_datasets/bed_Ensembl/reduced_CDS/images_CDS/imagens_v7/codingconsequence_CDS.csv", row.names = FALSE)
