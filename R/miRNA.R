# 0. Load libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)

library(clusterProfiler)


## 1. Load data ----

patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

mirna_tgca <- read_csv("data/imo_clustering/tests_jupyter/TCGA/raw data/cancer_data_PAAD_miRNASeqGene-20160128.csv")

mirna_tgca <- as.data.frame(t(mirna_tgca))

mirna_tgca[1,] <- mirna_tgca[1,] %>% str_replace(pattern = "mir", replacement = "miR")

colnames(mirna_tgca) <- mirna_tgca[1,]
mirna_tgca <- mirna_tgca[-1,]



# Retrieve patients names

patient <- rownames(mirna_tgca)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with 0 values in more than 20% of patient

mirna_tgca <- apply(mirna_tgca, 2, as.numeric)

table(is.na(mirna_tgca))


mirna_tgca <- mirna_tgca[, which(!colMeans((mirna_tgca == 0)) > 0.2)]

mirna_tgca <- as.data.frame(mirna_tgca)

mirna_tgca$patient <- patient

#rownames(met_tgca) <- met_tgca$patient

## Assign cluster IDs

mirna_tgca <- left_join(mirna_tgca, patient_cluster[,c(1,3)], by = join_by(patient == patient))

mirna_tgca$cluster <- as.factor(mirna_tgca$cluster)

## 2. Quick EDA ----

# PCA and HPCP

mirna_mat <- as.matrix(subset(mirna_tgca, select = -c(cluster,patient)))

rownames(mirna_mat) <- mirna_tgca$patient

png(filename = "results/miRNA/hist_reads.png", width = 12, height = 8, units = "in", res = 300)
barplot(rowSums(mirna_mat), col = "firebrick2", 
        main = "miRNASeqGene", sub = "Histogram of Reads by Patient")
dev.off()

png(filename = "results/miRNA/plotRLEmiRNA.png", 
    width = 12, height = 5, res = 300, units = "in")
EDASeq::plotRLE(t(mirna_mat), outline = FALSE,
                col = as.numeric(mirna_tgca$cluster),
                main = "TCGA miRNASeqGene")
dev.off()


# Distribution of counts

mirna_tgca_long <- pivot_longer(subset(mirna_tgca, select = -patient), cols = -cluster, names_to = "patient", values_to = "Reads")

mirna_tgca_long %>% ggplot(aes(x = log(Reads))) + geom_density()

## 3. Differential Expression Integrated with mRNA ----

# https://www.bioconductor.org/packages/release/bioc/vignettes/MIRit/inst/doc/MIRit.html 

library(MIRit)

# miRNA data
mirnaCounts <- t(mirna_mat)

# mRNA data

rna_tgca <- read_csv("data/imo_clustering/tests_jupyter/TCGA/raw data/cancer_data_PAAD_RNASeq2Gene-20160128.csv")
rna_tgca <- as.data.frame(rna_tgca)
rownames(rna_tgca) <- rna_tgca[,1]
rna_tgca <- rna_tgca[,-1]

# Retrieve patients names

patient <- colnames(rna_tgca)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with 0 values in more than 20% of patient

rna_tgca <- rna_tgca[which(rowMeans(!(rna_tgca == 0)) > 0.2),]

geneCounts <- round(as.matrix(rna_tgca),0)
colnames(geneCounts) <- patient

# Define sample metadata

patient_cluster <- patient_cluster %>% filter(patient %in% colnames(mirnaCounts))

## create a data.frame containing sample metadata
meta <- data.frame(primary = patient_cluster$patient,
                   mirnaCol = patient_cluster$patient,
                   geneCol = patient_cluster$patient,
                   cluster = patient_cluster$cluster)

# Create a MirnaExperiment object

## create the MirnaExperiment object
experiment <- MirnaExperiment(mirnaExpr = mirnaCounts,
                              geneExpr = geneCounts,
                              samplesMetadata = meta,
                              pairedSamples = TRUE)

# Differential expression analysis

# Visualize expression variability

geneMDS <- plotDimensions(experiment,
                          assay = "genes",
                          condition = "cluster",
                          title = "MDS plot for genes")

mirnaMDS <- plotDimensions(experiment,
                           assay = "microRNA",
                           condition = "cluster",
                           title = "MDS plot for miRNAs")

ggpubr::ggarrange(geneMDS, mirnaMDS,
                  nrow = 1, labels = "AUTO", common.legend = TRUE)

# Perform miRNA and gene differential expression

# Model design

model <- ~ cluster

## perform differential expression for genes
experiment <- performGeneDE(experiment,
                            method = "voom",
                            group = "cluster",
                            contrast = "Cluster_1-Cluster_0",
                            design = model,
                            pCutoff = 0.01, 
                            logFC = 1)

experiment <- performMirnaDE(experiment,
                             method = "voom",
                             group = "cluster",
                             contrast = "Cluster_1-Cluster_0",
                             design = model,
                             pCutoff = 0.01, 
                             logFC = 0)

# Visualize DEFs

## access DE results for genes
deGenes <- geneDE(experiment)
openxlsx::write.xlsx(deGenes, "results/miRNA/deGenes.xlsx")

## access DE results for miRNAs
deMirnas <- mirnaDE(experiment)
openxlsx::write.xlsx(deMirnas, "results/miRNA/deMirnas.xlsx")

# Create a volcano plot for miRNAs and genes

## create a volcano plot for genes
geneVolcano <- plotVolcano(experiment,
                           assay = "genes",
                           title = "Gene differential expression") +
  theme(title = element_text(face = "bold"))

## create a volcano plot for miRNAs
mirnaVolcano <- plotVolcano(experiment,
                            assay = "microRNA",
                            title = "miRNA differential expression") +
  theme(title = element_text(face = "bold"))

## plot graphs side by side
ggpubr::ggarrange(geneVolcano, mirnaVolcano,
                  nrow = 1, labels = "AUTO", common.legend = TRUE)
ggsave("results/miRNA/DEGs_RNA_miRNA_volcano.png", width = 12, height = 8, dpi = 300)




## 4. Functional Enrichment Analysis ----

# ORA GO BP 
ora_go <- enrichGenes(experiment,
                      method = "ORA",
                      database = "GO",
                      category = "bp",
                      organism = "Homo sapiens",
                      pCutoff = 0.05,
                      pAdjustment = "fdr")


## create a dot plot for ORA
enrichmentDotplot(ora_go$downregulated, title = "GO BP Downregulated", showTerms = 20)
enrichmentDotplot(ora_go$upregulated, title = "GO BP Upregulated", showTerms = 20)

# ORA KEGG Pathways
ora_kegg <- enrichGenes(experiment,
                      method = "ORA",
                      database = "KEGG",
                      category = "pathway",
                      organism = "Homo sapiens",
                      pCutoff = 0.05,
                      pAdjustment = "fdr")


## create a dot plot for ORA
enrichmentDotplot(ora_kegg$downregulated, title = "KEGG Pathways Downregulated")
enrichmentDotplot(ora_kegg$upregulated, title = "KEGG Pathways Upregulated")

# GSEA Hallmarks
set.seed(123)

gsea_hallmarks <- enrichGenes(experiment,
                              method = "GSEA",
                              database = "MsigDB",
                              category = "H",
                              organism = "Homo sapiens",
                              pCutoff = 0.05,
                              pAdjustment = "fdr")


## create a dot plot for GSEA
enrichmentDotplot(gsea_hallmarks, title = "Top 20 Hallmarks", showTerms = 20,
                  showTermsParam = "ratio") + theme(title = element_text(face = "bold"))
ggsave("results/miRNA/enrichment_gsea_hallmarks_dotplot.pdf", width = 12, height = 10)


enrichmentResults(gsea_hallmarks)
openxlsx::write.xlsx(enrichmentResults(gsea_hallmarks), "results/miRNA/enrichment_gsea_hallmarks.xlsx")


# GSEA C6 Cancer

set.seed(123)
gsea_c6 <- enrichGenes(experiment,
                              method = "GSEA",
                              database = "MsigDB",
                              category = "C6",
                              organism = "Homo sapiens",
                              pCutoff = 0.05,
                              pAdjustment = "fdr")


## create a dot plot for GSEA
enrichmentDotplot(gsea_c6, title = "Top 20 C6 Oncogenic Terms", showTerms = 20) +
  theme(title = element_text(face = "bold"))
ggsave("results/miRNA/enrichment_gsea_C6_dotplot.pdf", width = 12, height = 10)


gsea_c6_results <- enrichmentResults(gsea_c6)
openxlsx::write.xlsx(enrichmentResults(gsea_c6), "results/miRNA/enrichment_gsea_C6.xlsx")

# ORA Network of Cancer Genes

ora_ncg <- enrichGenes(experiment,
                       method = "ORA",
                       database = "NCG",
                       category = "v7",
                       organism = "Homo sapiens",
                       pCutoff = 0.25,
                       pAdjustment = "fdr")


## create a dot plot for ORA

#enrichmentDotplot(ora_ncg$upregulated, title = "Upregulated NCG", showTerms = 20)
enrichmentDotplot(ora_ncg$downregulated, title = "Downregulated NCG", showTerms = 20)
ora_ncg_down <- enrichmentResults(ora_ncg$downregulated)


# GSEA Network of Cancer Genes
set.seed(123)
gsea_ncg <- enrichGenes(experiment,
                        method = "GSEA",
                        database = "NCG",
                        category = "v7",
                        organism = "Homo sapiens",
                        pCutoff = 0.25,
                        pAdjustment = "fdr")


## create a dot plot for GSEA
enrichmentDotplot(gsea_ncg, title = "Top NCG", showTerms = 20)

## 5. Retrieve miRNA targets ----

# problems accessing to https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/miRTarBase_MTI.csv

# Download table manually
miRTarBase <- read.csv("data/miRTarBase_MTI.csv")

library(BiocParallel)

experiment <- getTargets(experiment,
                         organism = "Homo sapiens",
                         score = "Very High",
                         includeValidated = T,
                         evidence = "all",
                         local = miRTarBase) 

# perform a correlation analysis

# Load modified correlation function (one in package produces errors)

source("R/correlateMirnaTargets.R")

experiment <- correlateMirnaTargets(
  mirnaObj = experiment,
  corMethod = "spearman",
  corCutoff = 0.5,
  pCutoff = 0.05,
  pAdjustment = "fdr",
  BPPARAM = bpparam("MulticoreParam")
)

integrationResults <- integration(experiment)

openxlsx::write.xlsx(integrationResults, "results/miRNA/integrationResults.xlsx")

leadingEdge <- gsea_c6_results$leadingEdge
names(leadingEdge) <- gsea_c6_results$pathway

leadingEdge <- plyr::ldply(leadingEdge, rbind)

leadingEdge <- as.data.frame(t(leadingEdge))

colnames(leadingEdge) <- leadingEdge[1,]
leadingEdge <- leadingEdge[-1,]

leadingEdge <- leadingEdge %>% pivot_longer(cols = everything(), names_to = "Gene_Sets", values_to = "Genes")

miRGene_Sets <- integrationResults %>% filter(Target %in% leadingEdge$Genes)

miRGene_Sets <- left_join(miRGene_Sets, leadingEdge, by = join_by(Target == Genes))

miRGene_Sets <- left_join(miRGene_Sets, gsea_c6_results[,c("pathway", "padj", "NES")], by = join_by(Gene_Sets == pathway))

miRGene_Sets$logpadj <- -log10(miRGene_Sets$padj)

miRGene_Sets <- miRGene_Sets %>%
  mutate(NES_sign = case_when(
    NES < 0 ~ "Enriched in Cluster 1",
    NES > 0 ~ "Enriched in Cluster 0",
    TRUE ~ NA_character_
  ))


library(ggsankeyfier)

miRGene_long <- pivot_stages_longer(
  data = miRGene_Sets,
  stages_from = c("microRNA", "Target", "Gene_Sets"),
  values_from = "Corr.Coefficient",
  additional_aes_from = c("NES_sign", "padj", "microRNA")
  )


pos <- position_sankey(v_space = "auto", align = "center",order = "ascending")
pos_text <- position_sankey(v_space = "auto", align = "center", order = "ascending", nudge_x = 0.1)

ggplot(
  miRGene_long,
  aes(
    x = stage,
    y = Corr.Coefficient,
    group = node,
    connector = connector,
    edge_id = edge_id
    )) +
  geom_sankeynode(aes(color = NES_sign), position = pos) +
  geom_sankeyedge(aes(fill = microRNA, alpha = Corr.Coefficient),position = pos) +
  geom_text(aes(label = node), stat = "sankeynode", position = pos_text, hjust = 0, cex = 4) +
  scale_x_discrete(expand = expansion(add = c(0.2, .6))) + 
  scale_alpha_continuous(name = "miRNA-Target\nCorrelation Coefficient",limits = c(-0.8, -0.5), range = c(1,0.6)) +
  theme_void() + labs(title = "miRNA-Target-Oncogenic Gene Set Relationship") +
  theme(title = element_text(face = "bold"))
ggsave("results/miRNA/miRNA_Target_GeneSet_sankey.pdf", width = 14, height = 10)


cor1 <- plotCorrelation(experiment,
                        mirna = "hsa-miR-200b",
                        gene = "BCL2",
                        condition = "cluster")

cor2 <- plotCorrelation(experiment,
                        mirna = "hsa-miR-200b",
                        gene = "SESN1",
                        condition = "cluster")

cor3 <- plotCorrelation(experiment,
                        mirna = "hsa-miR-96",
                        gene = "BCL2",
                        condition = "cluster")

cor4 <- plotCorrelation(experiment,
                        mirna = "hsa-miR-135b",
                        gene = "GIMAP6",
                        condition = "cluster")

## plot graphs side by side
cowplot::plot_grid(cor1, cor3, cor4, ncol = 2)
ggsave("results/miRNA/miRNA_targets_cor.pdf", width = 10, height = 8)

# Functional enrichment of integrated target genes

oraTarg_hallmarks <- enrichTargets(
  experiment,
  database = "MsigDB",
  category = "H",
  organism = "Homo sapiens",
  pCutoff = 0.05,
  pAdjustment = "none"
)



enrichmentResults(oraTarg_hallmarks$upregulated)
openxlsx::write.xlsx(enrichmentResults(oraTarg_hallmarks$upregulated), "results/miRNA/oraTarg_hallmarks_up.xlsx")
enrichmentResults(oraTarg_hallmarks$downregulated)
openxlsx::write.xlsx(enrichmentResults(oraTarg_hallmarks$downregulated), "results/miRNA/oraTarg_hallmarks_down.xlsx")


enrichmentDotplot(oraTarg_hallmarks$upregulated)
enrichmentDotplot(oraTarg_hallmarks$downregulated)


oraTarg_C6 <- enrichTargets(
  experiment,
  database = "MsigDB",
  category = "C6",
  organism = "Homo sapiens",
  pCutoff = 0.05,
  pAdjustment = "none"
)

oraTarg_C6_up <- enrichmentResults(oraTarg_C6$upregulated)
openxlsx::write.xlsx(enrichmentResults(oraTarg_C6$upregulated), "results/miRNA/oraTarg_C6_up.xlsx")
oraTarg_C6_down <- enrichmentResults(oraTarg_C6$downregulated)
openxlsx::write.xlsx(enrichmentResults(oraTarg_C6$downregulated), "results/miRNA/oraTarg_C6_down.xlsx")

oraTarg_C6_up$Enrich <- "Upregulated"
oraTarg_C6_down$Enrich <- "Downregulated"

oraTarg_C6_up_down <- rbind(oraTarg_C6_up, oraTarg_C6_down)

oraTarg_C6_up_down %>% arrange(foldEnrichment) %>%
  mutate(pathway = factor(pathway, levels = pathway)) %>%
  ggplot(aes(x = foldEnrichment, y = pathway)) +
  geom_segment(aes(xend = 0, yend = pathway, color = Enrich)) +
  geom_point(aes(fill = padj, shape = Enrich, size = overlap)) +
  scale_fill_gradient(name = "Adj p-value", low = "red", high = "blue") +
  scale_size_continuous(name = "Count", range = c(5,10), breaks = c(2,4,6)) +
  geom_text_repel(aes(label = overlapGenes, color = Enrich)) +
  scale_color_manual(name = "", values = c("navy", "firebrick")) +
  scale_shape_manual(name = "Cluster 1 vs Cluster 0", values = c(21,23)) + theme_cowplot() +
  labs(title = "Oncogenic Signature ORA of integrated miRNA targets")
ggsave("results/miRNA/oraTarg_Oncogenic.pdf", width = 12, height = 12)


enrichmentDotplot(oraTarg_C6$upregulated)
enrichmentDotplot(oraTarg_C6$downregulated)

## 6. Topology-Aware Integrative Pathway Analysis (TAIPA) ----

## create miRNA-augmented networks using KEGG pathways
networks <- preparePathways(experiment,
                            database = "KEGG",
                            organism = "Homo sapiens",
                            minPc = 50, 
                            BPPARAM = bpparam("MulticoreParam"))


## set seed for reproducible results
set.seed(123)
## perform TAIPA with 1000 permutations
taipa <- topologicalAnalysis(experiment,
                             pathways = networks,
                             progress = T,
                             tasks = 20,
                             nPerm = 1000,
                             pCutoff = 0.05,
                             pAdjustment = "none")

save(taipa, file = "results/miRNA/taipa.Rdata")

load("results/miRNA/taipa.Rdata")

## produce a dotplot that shows the most affected networks
integrationDotplot(taipa, showTerms = 20)

## extract the results of TAIPA
perturbedNetworks <- integratedPathways(taipa)
openxlsx::write.xlsx(perturbedNetworks, "results/miRNA/perturbedNetworks.xlsx")


## plot the impaired network responsible for reduced TG synthesis

visualizeNetwork(taipa, "Glycerolipid metabolism", fontsize = 100)

