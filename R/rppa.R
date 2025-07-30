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

# 2. Quick EDA ----

# PCA and HPCP

rppa_mat <- as.matrix(subset(rppa_tgca, select = -c(cluster,patient)))

rownames(rppa_mat) <- rppa_tgca$patient

str(rppa_mat)


rppa_pca <- PCA(rppa_mat)

p2 <- fviz_pca_ind(rppa_pca, habillage = rppa_tgca$cluster, repel = TRUE, addEllipses = T, label = "none",
                   mean.point = FALSE) + theme_cowplot()

rppa_hcpc <- HCPC(rppa_pca, graph = FALSE)

p1 <- fviz_dend(rppa_hcpc, rect = T, rect_fill = T, cex = 0.2, ggtheme = cowplot::theme_cowplot())

p3 <- fviz_cluster(rppa_hcpc,
                   repel = TRUE,            # Avoid label overlapping
                   show.clust.cent = FALSE, # Show cluster centers
                   #palette = "jco",         # Color palette see ?ggpubr::ggpar
                   ggtheme = cowplot::theme_cowplot(),
                   main = "Factor map", geom = "point") + geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2)

plot_grid(p1, p2, p3)

rppa_tgca_long <- pivot_longer(subset(rppa_tgca, select = -patient), cols = -cluster, names_to = "patient", values_to = "beta")

rppa_tgca_long %>% ggplot(aes(x = beta)) + geom_density(aes(color = cluster, fill = cluster), alpha = 0.1, lwd = 1) + theme_cowplot() +
  xlab("Normalized Expression") + ylab("Density") + labs(title = "Normalized Expression Values")

# 3. Differential Expression ----

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

deg_count <- summary(decideTests(fit2, adjust.method = "none"))


barplot(height = c(deg_count["Up",], deg_count["Down",], deg_count["NotSig",]), names.arg = c("Up", "Down", "NotSig"), 
        col = c("pink", "lightblue", "gray80"), border = c("firebrick3", "slateblue4", "gray50"),
        ylab = "Count", main = "Differentially Expressed/Phosphorylated\nProteins at p-value < 0.05", ylim = c(0,200))

top5 <- rppa_toptable[1:5,]

top5 <- rppa_tgca %>% select(c(top5$peptide_target, cluster))

top5 %>% pivot_longer(cols = -cluster) %>% ggplot(aes(x = name, y = value)) +
  geom_boxplot(aes(fill = cluster), alpha = 0.5, outliers = F) +
  geom_dotplot(aes(fill = cluster), binaxis='y', stackdir='center',
                 position=position_jitterdodge(0.2), dotsize = 0.2) + theme_cowplot()

write_csv(rppa_toptable, file = "data/data_omics/rppa/rppa_toptable.csv")

## 3.1 Volcano Plot ----

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

## 4 ClusterProfiler ----

### GO over-representation analysis ----

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

browseKEGG(kk, pathID = "hsa04066")

library(pathview)

table(sel_DEP$gene_name == unique(gene.df$ALIAS))

df <- sel_DEP$logFC
names(df) <- gene.df[!duplicated(gene.df$ALIAS),"ENTREZID"]

hsa04066 <- pathview(gene.data = df,
                     pathway.id = "hsa04066",
                     species = "hsa", low = "slateblue2", high = "firebrick2")



cowplot::plot_grid(DEP_volcano, kegg_ora)
