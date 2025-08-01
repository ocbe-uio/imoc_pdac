# 0. Load libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)

# 1. Load data ----

patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

cna_tcga <- read_csv("data/data_omics/TCGA/raw data/cancer_data_PAAD_GISTIC_Peaks-20160128.csv")

cna_tcga_ID <- cna_tcga[,c(3,13:ncol(cna_tcga))]

cna_tcga_ID <- as.data.frame(t(cna_tcga_ID))

colnames(cna_tcga_ID) <- cna_tcga_ID[1,]
cna_tcga_ID <- cna_tcga_ID[-1,]

# Retrieve patients names

patient <- rownames(cna_tcga_ID)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

cna_tcga_ID <- apply(cna_tcga_ID, 2, as.numeric)

cna_tcga_ID <- as.data.frame(cna_tcga_ID)

cna_tcga_ID$patient <- patient

#rownames(met_tcga) <- met_tcga$patient

## Assign cluster IDs

cna_tcga_ID <- left_join(cna_tcga_ID, patient_cluster[,c(1,3)], by = join_by(patient == patient))

cna_tcga_ID$cluster <- as.factor(cna_tcga_ID$cluster)

cna_tcga_ID_long <- cna_tcga_ID %>% pivot_longer(cols = 1:52, names_to = "Descriptor", values_to = "CNA")
cna_tcga_ID_long$CNA[cna_tcga_ID_long$CNA == 0] <- NA

openxlsx2::write_xlsx(cna_tcga_ID_long, file = "data/data_omics/cna/cna_tcga_ID_long.xlsx")

head(cna_tcga_ID_long,20)

df <- cna_tcga_ID_long %>% group_by(patient, cluster) %>% summarise(counts = mean(!is.na(CNA)))

df %>% ggplot(aes(x = counts)) + geom_density(aes(colour = cluster))

# 2.  Copy number alteration by patient ----

df %>% ggplot(aes(x = cluster, y = counts)) + 
  geom_violin(aes(color = cluster, fill = cluster), alpha = 0.15, trim = F,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(aes(color = cluster), shape = 1, width = 0.15) +
  ylab("Number of CNA by Patient") + xlab("") +
  theme_cowplot() + labs(title = "Mean Copy Number Alterations by Patient",
                         caption = "p-value < 2.2e-16 Wilcoxon rank sum test")
wilcox.test(counts ~ cluster, data = df)

# 3. Contingency table ----

contin_cna <- cna_tcga_ID_long %>% group_by(Descriptor, cluster) %>% summarise(Altered = sum(!is.na(CNA)))

table(cna_tcga_ID$cluster == "Cluster_2")

contin_cna <- contin_cna %>% mutate(Normal = c(32,121) - Altered)


# Create list of contingency tables
contingency_tables <- split(contin_cna, contin_cna$Descriptor) %>%
  lapply(function(sub_df) {
    table_matrix <- as.matrix(sub_df[, c("Altered", "Normal")])
    rownames(table_matrix) <- sub_df$cluster
    return(table_matrix)
  })

contingency_tables

cna_fisher <- lapply(contingency_tables, FUN = fisher.test)

pvals <- cna_fisher %>% map_dfr(\(x){
  tibble(pval = x %>% pluck("p.value"))
  })

pvals$Descriptor <- names(contingency_tables)

pvals$fisher_adjpval <- p.adjust(pvals$pval, method = "BH")


OddsRatio <- function(mat) {
  OR <- ((mat[2]+0.5)/(mat[4]+0.5))/((mat[1]+0.5)/(mat[3]+0.5))
  return(OR)
}

OR <- lapply(contingency_tables, OddsRatio)
OR <- data.frame(OR = unlist(OR))
OR$Descriptor <- rownames(OR)

cna_tcga <- left_join(cna_tcga, pvals[,c("Descriptor", "fisher_adjpval")], by = join_by(Descriptor == Descriptor))
cna_tcga <- left_join(cna_tcga, OR, by = join_by(Descriptor == Descriptor))

openxlsx2::write_xlsx(cna_tcga, "data/data_omics/cna/cna_tcga.xlsx")

cna_tcga %>% filter(fisher_adjpval < 0.001) %>% arrange(OR) %>% 
  mutate(Descriptor = factor(Descriptor, levels = Descriptor)) %>% ggplot(aes(x = OR, y = Descriptor)) +
  geom_segment(aes(xend = 0, yend = Descriptor), color = "gray80") +
  geom_point(aes(fill = -log10(fisher_adjpval), shape = type), size = 4) +
  scale_shape_manual(values = c(21,23), name = "Type of CNA") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 Fisher's\ntest FDR", limits = c(3,16)) +
  theme_cowplot() +
  labs(title = "Fisher Exact Test Results CNA Cluster 2 vs Cluster 1", subtitle = "Fisher's test FDR < 0.001") +
  xlab("Odds Ratio") + ylab("Cytoband")

sel_cna <- c("17p12", "9p21.3", "21q11.2", "18q21.2")
cytobands <- cna_tcga %>% filter(Descriptor %in% sel_cna)

genes_amp <- read.delim("data/data_omics/TCGA/raw data/data_gistic_genes_amp.txt")
genes_del <- read.delim("data/data_omics/TCGA/raw data/data_gistic_genes_del.txt")

genes <- rbind(genes_amp, genes_del)

cytobands <- left_join(cytobands, genes[,c("cytoband","genes_in_peak","genes_in_region")], by = join_by(Descriptor == cytoband))
save(cytobands, file = "data/data_omics/cna/cytobands.RData")

genes_in_peak <- cytobands$genes_in_peak %>% str_split(pattern = "[,|]")
names(genes_in_peak) <- cytobands$Descriptor

genes_in_peak <- lapply(genes_in_peak, str_remove_all, pattern = "[\\[\\]]")
genes_in_peak <- lapply(genes_in_peak, function(x){x[x!=""]})

genes_in_region <- cytobands$genes_in_region %>% str_split(pattern = "[,|]")
names(genes_in_region) <- cytobands$Descriptor

genes_in_region <- lapply(genes_in_region, str_remove_all, pattern = "[\\[\\]]")
genes_in_region <- lapply(genes_in_region, function(x){x[x!=""]})


cytobands_select <- cytobands %>% dplyr::select(Descriptor, type, fisher_adjpval, OR)
cytobands_select$genes_in_peak <- genes_in_peak
cytobands_select$genes_in_region <- genes_in_region

save(cytobands_select, file = "data/data_omics/cna/cytobands_select.Rdata")

cytobands_select <- subset(cytobands_select, select = -genes_in_region)

cytobands_select <- cytobands_select %>% unnest_longer(col = genes_in_peak)
openxlsx::write.xlsx(cytobands_select, "data/data_omics/cna/cytobands_select.xlsx")

cytobands_select %>% arrange(OR) %>% 
  mutate(Descriptor = factor(Descriptor, levels = c(unique(Descriptor)))) %>% ggplot(aes(x = OR, y = Descriptor)) +
  geom_segment(aes(xend = 0, yend = Descriptor), color = "gray80") +
  geom_point(aes(fill = -log10(fisher_adjpval), shape = type), size = 5) +
  geom_text_repel(aes(label = genes_in_peak, color = type), size = 4, max.overlaps = Inf, force = 1, segment.linetype =2, segment.alpha = 0.5) +  # Add gene names
  scale_shape_manual(values = c(23), name = "Type of CNA") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 Fisher's\ntest FDR", limits = c(3,16)) +
  scale_color_manual(values = c("firebrick", "navy"), name = "Type of CNA")+
  theme_cowplot() +
  labs(title = "Genes in Significant Fisher Test Cytobands") +
  xlab("Odds Ratio") + ylab("Cytoband")


