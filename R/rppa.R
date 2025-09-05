# 0. Load libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(limma)
library(edgeR)

# 1. Load data ----

patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

rppa_tgca <- read_csv("data/data_omics/TCGA/processed_data/RPPA_processed.csv")

rppa_tgca$patient <- rppa_tgca$...1
rppa_tgca <- rppa_tgca[,-1]

## Assign cluster IDs

rppa_tgca <- left_join(rppa_tgca, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rppa_tgca$cluster <- as.factor(rppa_tgca$cluster)


# 2. Differential Expression ----

rppa <- as.data.frame(t(subset(rppa_tgca, select = -c(patient, cluster))))
colnames(rppa) <- rppa_tgca$patient %>% str_replace_all(pattern = "-", replacement = "_")

# The clusters are our factor of interest
cluster <- rppa_tgca$cluster

# use the above to create a design matrix
design <- model.matrix(~0+cluster)
colnames(design) <- c(levels(cluster))

# fit the linear model 
fit <- lmFit(object = rppa, design = design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(Cluster_2 - Cluster_1, levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2, trend = TRUE)

rppa_ann <- read_delim("data/data_omics/rppa/TCGA_antibodies_descriptions.gencode.v36.tsv", delim = "\t")

rppa_toptable <- topTable(fit2, number = Inf)
peptide_target <- rownames(rppa_toptable)
rppa_toptable$peptide_target <- rownames(rppa_toptable)
rppa_toptable$peptide_target <- rppa_toptable$peptide_target %>% str_to_lower() %>% str_replace_all(c("_" = "", "-" = ""))

rppa_ann$peptide_target_ed <- rppa_ann$peptide_target %>% str_to_lower() %>% str_replace_all(c("_" = "", "-" = ""))

rppa_toptable <- left_join(rppa_toptable, subset(rppa_ann, select = c("peptide_target_ed", "gene_name", "gene_id", "source", "catalog_number")), 
                           by = join_by(peptide_target == peptide_target_ed))

rppa_toptable$peptide_target <- peptide_target
rppa_toptable[is.na(rppa_toptable$gene_id),]

write_csv(rppa_toptable, file = "data/data_omics/rppa/rppa_toptable.csv")

### Volcano Plot ----

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  rppa_toptable$logFC < 0 & rppa_toptable$P.Value < 0.05, 'slateblue4',
  ifelse(rppa_toptable$logFC > 0 & rppa_toptable$P.Value < 0.05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'


EnhancedVolcano(toptable = rppa_toptable, lab = rppa_toptable$peptide_target, x = "logFC", y = "P.Value", FCcutoff = 0,
                               xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                               colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 1 vs Cluster 0",
                               title = "Volcano Plot of Differentially Expressed/Phosphorylated Proteins",
                               boxedLabels = T, labCol = keyvals[names(keyvals) != "NS"], labSize = 4, ylim = c(0,5), xlim = c(-1,1),
                               drawConnectors = T, max.overlaps = Inf, legendPosition = "top", gridlines.minor = F, gridlines.major = F)

# 3. ClusterProfiler ----

library(clusterProfiler)
library(org.Hs.eg.db)

rppa_toptable_sig <- rppa_toptable %>% filter(P.Value < 0.05)

# Convert IDs
gene.df <- bitr(rppa_toptable_sig$gene_name, fromType = "ALIAS", toType = c("ENTREZID", "ENSEMBL"),
                OrgDb = org.Hs.eg.db)

### KEGG pathway over-representation analysis ----

kk <- enrichKEGG(
  gene = gene.df$ENTREZID,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.01,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)


kk_results <- kk@result

GeneID <- kk_results$geneID

GeneID <- GeneID %>% str_split(pattern = "/")

mapENTREZ <- function(id, df){
  ALIAS <- df[(df$ENTREZID %in% id), "ALIAS"]
  return(ALIAS)
}

GeneID <- lapply(GeneID, FUN = mapENTREZ, df = gene.df)

kk_results$ALIAS <- GeneID
openxlsx2::write_xlsx(kk_results, file = "data/data_omics/rppa/KEGG_ORA.xlsx")

kk_results_long <- kk_results %>% unnest(ALIAS)

sel_DEP <- rppa_toptable_sig[!duplicated(rppa_toptable_sig$gene_name),]

kk_results_long <- left_join(kk_results_long, sel_DEP[,c("logFC", "gene_name")], by = join_by(ALIAS == gene_name))
kk_results_long <- distinct(kk_results_long)

openxlsx2::write_xlsx(kk_results_long, file = "data/data_omics/rppa/KEGG_ORA_long.xlsx")

kk_results_long %>% filter(p.adjust < 0.001) %>%
  arrange(RichFactor) %>% mutate(Description = factor(Description, levels = unique(Description))) %>% 
  ggplot(aes(x = RichFactor, y = Description)) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray50") +
  geom_text_repel(aes(label = ALIAS, color = logFC), max.overlaps = Inf, fontface = "bold") +
  geom_point(aes(fill = -log10(p.adjust)), shape = 21, size = 8, color = "gray50") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 FDR") +
  scale_color_gradient(low = "blue", high= "red", limits = c(-0.5,0.5)) +
  theme_cowplot() + labs(title = "KEGG Pathways Over-Representation Analysis",
                         subtitle = "FDR < 0.001")

