# 0. Load libraries ----

library(tidyverse)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)

# 1. Load data ----

patient_cluster <- read_csv("results/cluster_analysis/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

cna_tcga <- read_csv("data/TCGA/omics_data/raw/cancer_data_PAAD_GISTIC_Peaks-20160128.csv")

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

## Assign cluster IDs

cna_tcga_ID <- left_join(cna_tcga_ID, patient_cluster[,c(1,3)], by = join_by(patient == patient))

cna_tcga_ID$cluster <- as.factor(cna_tcga_ID$cluster)

cna_tcga_ID_long <- cna_tcga_ID %>% pivot_longer(cols = 1:52, names_to = "Descriptor", values_to = "CNA")
cna_tcga_ID_long$CNA[cna_tcga_ID_long$CNA == 0] <- NA

openxlsx2::write_xlsx(cna_tcga_ID_long, file = "results/omics_analysis/CNA/cna_tcga_ID_long.xlsx")

# 2. Contingency table ----

contin_cna <- cna_tcga_ID_long %>% group_by(Descriptor, cluster) %>% summarise(Altered = sum(!is.na(CNA)))

contin_cna <- contin_cna %>% mutate(Normal = c(32,121) - Altered)

# Create list of contingency tables
contingency_tables <- split(contin_cna, contin_cna$Descriptor) %>%
  lapply(function(sub_df) {
    table_matrix <- as.matrix(sub_df[, c("Altered", "Normal")])
    rownames(table_matrix) <- sub_df$cluster
    return(table_matrix)
  })

cna_fisher <- lapply(contingency_tables, FUN = fisher.test)

pvals <- cna_fisher %>% map_dfr(\(x){
  tibble(pval = x %>% pluck("p.value"))
  })

pvals$Descriptor <- names(contingency_tables)

pvals$fisher_adjpval <- p.adjust(pvals$pval, method = "BH")


OddsRatio <- function(mat) {
  a <- mat[2] + 0.5
  b <- mat[4] + 0.5 
  c <- mat[1] + 0.5
  d <- mat[3] + 0.5
  OR <- ((a)/(b))/((c)/(d))
  SElogOR <- sqrt(1/a+1/b+1/c+1/d)
  lwr <- exp(log(OR) - 1.96*SElogOR)
  upr <- exp(log(OR) + 1.96*SElogOR)
  OR <- data.frame(OR = OR, lwr = lwr, upr = upr)
  return(OR)
}



OR <- lapply(contingency_tables, OddsRatio)
OR <- imap_dfr(OR, ~mutate(.x, Descriptor = .y))

cna_tcga <- left_join(cna_tcga, pvals[,c("Descriptor", "fisher_adjpval")], by = join_by(Descriptor == Descriptor))
cna_tcga <- left_join(cna_tcga, OR, by = join_by(Descriptor == Descriptor))

openxlsx2::write_xlsx(cna_tcga, "results/omics_analysis/CNA/cna_tcga.xlsx")

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

## Assign genes to selected cytobands ----

# Loading Ensembl database
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Checking names of filters and attributes
filters <- listFilters(ensembl)
attributes = listAttributes(ensembl)

# Preparing the list of coordinates

WidePeakLimits <- cytobands$Wide.Peak.Limits %>% str_remove("chr") %>% str_remove("\\(.*") %>% str_replace_all(pattern = "-", replacement = ":")
PeakLimits <- cytobands$Peak.Limits %>% str_remove("chr") %>% str_remove("\\(.*") %>% str_replace_all(pattern = "-", replacement = ":")

# Wide Peak Limits
genesWidePeak <- lapply(WidePeakLimits, function(r) {
  coords <- strsplit(r, ":")[[1]]
  getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand","description"),
        filters = c("chromosome_name", "start", "end"),
        values = list(coords[1], coords[2], coords[3]),
        mart = ensembl)
})

names(genesWidePeak) <- cytobands$Descriptor

# Flatten genes without extra commas
flat_genes <- lapply(genesWidePeak, function(x){
  symbols <- x$hgnc_symbol[x$hgnc_symbol != ""]
  str_flatten_comma(symbols)
})

# Convert to a dataframe
flat_genes <- data.frame(
  Descriptor = names(flat_genes),
  Genes.Wide.Peak = unlist(flat_genes),
  stringsAsFactors = FALSE
)

cytobands <- cytobands %>% left_join(flat_genes, by = "Descriptor")

# Peak Limits

genesPeakLimits <- lapply(PeakLimits, function(r) {
  coords <- strsplit(r, ":")[[1]]
  getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand","description"),
        filters = c("chromosome_name", "start", "end"),
        values = list(coords[1], coords[2], coords[3]),
        mart = ensembl)
})

names(genesPeakLimits) <- cytobands$Descriptor

# Flatten genes without extra commas
flat_genes <- lapply(genesPeakLimits, function(x){
  symbols <- x$hgnc_symbol[x$hgnc_symbol != ""]
  str_flatten_comma(symbols)
})

# Convert to a dataframe
flat_genes <- data.frame(
  Descriptor = names(flat_genes),
  Genes.Peak.Limits = unlist(flat_genes),
  stringsAsFactors = FALSE
)

cytobands <- cytobands %>% left_join(flat_genes, by = "Descriptor")

save(cytobands, file = "results/omics_analysis/CNA/cytobands.RData")

genes_in_peak <- cytobands$Genes.Peak.Limits %>% str_split(pattern = ", ")
names(genes_in_peak) <- cytobands$Descriptor

genes_in_wide_peak <- cytobands$Genes.Peak.Limits %>% str_split(pattern = ", ")
names(genes_in_wide_peak) <- cytobands$Descriptor

cytobands_select <- cytobands %>% dplyr::select(Descriptor, type, fisher_adjpval, OR)
cytobands_select$genes_in_peak <- genes_in_peak
cytobands_select$genes_in_wide_peak <- genes_in_wide_peak

cytobands_select <- subset(cytobands_select, select = -genes_in_wide_peak)

cytobands_select <- cytobands_select %>% unnest_longer(col = genes_in_peak)
openxlsx::write.xlsx(cytobands_select, "results/omics_analysis/CNA/cytobands_select.xlsx")

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

## Proportion of homozygous vs hemizygous deletions ----

counts_df <- cna_tcga_ID_long
counts_df$CNA[is.na(counts_df$CNA)] <- 0
counts_df$CNA <- factor(counts_df$CNA, levels = c("0", "1", "2"))

counts_df <- counts_df %>%
  group_by(cluster, CNA) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(cluster, CNA)

counts_df <- counts_df %>%
  tidyr::complete(cluster, CNA, fill = list(n = 0))


counts_df %>% ggplot(aes(x = CNA, y = n)) + 
  geom_bar(aes( color = cluster, fill = cluster), alpha = 0.4,stat = "identity", position = position_dodge()) + 
  labs(title = "Alleles affected by copy number alterations") + xlab("Alleles affected") + ylab("Count") + 
  ylim(0,5000) + theme_cowplot()

cytobands_long <- cna_tcga_ID_long %>% filter(Descriptor %in% cytobands$Descriptor)
cytobands_long$CNA[is.na(cytobands_long$CNA)] <- 0
cytobands_long$CNA <- factor(cytobands_long$CNA, levels = c("0", "1", "2"))

cytobands_long_count <- cytobands_long %>%
  group_by(cluster, CNA) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(cluster, CNA)

cytobands_long_count <- cytobands_long_count %>%
  tidyr::complete(cluster, CNA, fill = list(n = 0))


cytobands_long_count %>% ggplot(aes(x = CNA, y = n)) + 
  geom_bar(aes( color = cluster, fill = cluster), alpha = 0.4,stat = "identity", position = position_dodge()) + 
  labs(title = "Alleles affected by copy number alterations\nin selected cytobands") + xlab("Alleles affected") + ylab("Count") + 
  ylim(0,400) +
  theme_cowplot()

