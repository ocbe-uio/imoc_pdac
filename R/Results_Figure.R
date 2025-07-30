# 0. Load Libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(limma)
library(EnhancedVolcano)
library(missMDA)
library(minfi)
library(openxlsx2)

# 1. Load Data ----

## Methylomics ----

# Selected features
met_features_genes <- openxlsx2::read_xlsx("data/data_omics/methylomics/met_features_genes.xlsx")

met_tcga <- read_csv("data/data_omics/methylomics/met_tcga.csv")

# Differential Methylation Analysis (all CpGs)
DMPs <- openxlsx::read.xlsx("data/data_omics/methylomics/DMPs.xlsx")

# Selected CpGs (abs(logFC) > 2 & adj.P.Val < 1e-05) with annotations
DCpGs <- openxlsx::read.xlsx("data/data_omics/methylomics/DCpGs.xlsx")

## CNA ----

cna_tcga_ID_long <- openxlsx2::read_xlsx("data/data_omics/cna/cna_tcga_ID_long.xlsx")
cna_tcga <- read_xlsx("data/data_omics/cna/cna_tcga.xlsx")
cytobands_select <- read_xlsx("data/data_omics/cna/cytobands_select.xlsx")

sel_cna <- c("17p12", "9p21.3", "21q11.2", "18q21.2")

## RNA ----

rna_mat_norm <- read.csv("data/data_omics/rna_seq/rna_mat_norm.csv", row.names = 1)
rna_tcga <- read.csv("data/data_omics/rna_seq/rna_tcga.csv", row.names = 1)
limma_results <- read.csv("data/data_omics/rna_seq/limma_results.csv", row.names = 1)

## RPPA ----

rppa_toptable <- read_csv("data/data_omics/rppa/rppa_toptable.csv")
kk_results_long <- openxlsx2::read_xlsx("data/data_omics/rppa/KEGG_ORA_long.xlsx")

# 2. Plots ----

### CNA by patient ----

df <- cna_tcga_ID_long %>% group_by(patient, cluster) %>% summarise(counts = sum(!is.na(CNA)))

cna_violin <- df %>% ggplot(aes(x = cluster, y = counts)) + 
  geom_violin(aes(color = cluster, fill = cluster), alpha = 0.15, trim = T,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(aes(color = cluster), shape = 1, width = 0.15) +
  ylab("Number of CNA by Patient") + xlab("") +
  theme_cowplot() + labs(title = "Total Copy Number Alterations by Patient") +
  annotate(geom = "text", x = 1, y = 30, label = "p-value < 2.2e-16 Wilcoxon rank sum test")
cna_violin

### Odds Ratio CNA ----

OR_cna <- cna_tcga %>% filter(fisher_adjpval < 0.001) %>% arrange(OR) %>% 
  mutate(Descriptor = factor(Descriptor, levels = Descriptor)) %>% ggplot(aes(x = OR, y = Descriptor)) +
  geom_segment(aes(xend = 0, yend = Descriptor), color = "gray80") +
  geom_point(aes(fill = -log10(fisher_adjpval), shape = type), size = 6) +
  scale_shape_manual(values = c(21,23), name = "Type of CNA") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 Fisher's\ntest FDR", limits = c(3,16)) +
  theme_cowplot() +
  labs(title = "Odds Ratio CNA Cluster 1 vs Cluster 0", subtitle = "Fisher's test FDR < 0.001") +
  xlab("Odds Ratio") + ylab("Cytoband")
OR_cna

### Odss Ratio with Genes ----

OR_cna_genes <- cytobands_select %>% arrange(OR) %>% 
  mutate(Descriptor = factor(Descriptor, levels = c(unique(Descriptor)))) %>% ggplot(aes(x = OR, y = Descriptor)) +
  geom_segment(aes(xend = 0, yend = Descriptor), color = "gray80") +
  geom_point(aes(fill = -log10(fisher_adjpval), shape = type), size = 5) +
  geom_text_repel(aes(label = genes_in_peak, color = type), size = 4, max.overlaps = Inf, force = 1, segment.linetype =2, segment.alpha = 0.5) +  # Add gene names
  scale_shape_manual(values = c(23), name = "Type of CNA") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 Fisher's\ntest FDR", limits = c(3,16)) +
  scale_color_manual(values = c("firebrick", "navy"), name = "Type of CNA")+
  theme_cowplot() +
  labs(title = "Genes in Selected Cytobands") +
  xlab("Odds Ratio") + ylab("Cytoband")
OR_cna_genes

### Volcano plot Methylation ----

keyvals <- ifelse(
  DMPs$logFC < -2 & DMPs$adj.P.Val < 1e-05, 'slateblue4',
  ifelse(DMPs$logFC > 2 & DMPs$adj.P.Val < 1e-05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'

rownames(DMPs) <- DMPs$cpg_name

dmp_volcano <- EnhancedVolcano(toptable = DMPs, lab = rownames(DMPs), x = "logFC", y = "adj.P.Val", FCcutoff = 2,
                               xlab = "Log2 fold-change M value",
                               colCustom = keyvals, ylab = "-log10 Adjusted -value", subtitle = "Cluster 2 vs Cluster 1",
                               title = "Volcano Plot of Differentially Methylated CpG sites",
                               selectLab = met_features_genes$cpg_name,
                               boxedLabels = T, xlim = c(-4,4),
                               drawConnectors = TRUE, max.overlaps = Inf, legendPosition = "top", gridlines.minor = F, gridlines.major = F,
                               raster = TRUE)
dmp_volcano

### Volcano plot RNA ----

top.table <- limma_results

keyvals <- ifelse(
  top.table$logFC < -0.25 & top.table$adj.P.Val < 0.05, 'slateblue4',
  ifelse(top.table$logFC > 0.25 & top.table$adj.P.Val < 0.05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'

met_features_genes$UCSC_RefGene_Name <- met_features_genes$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")
rna_volcano <- EnhancedVolcano(toptable = top.table, lab = top.table$Gene, x = "logFC",
                               y = "adj.P.Val", FCcutoff = 0.25,
                               xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                               colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 2 vs Cluster 1",
                               title = "Volcano Plot of Differentially Expressed Genes", xlim = c(-4,4),
                               boxedLabels = T, labSize = 4, selectLab = c(met_features_genes$UCSC_RefGene_Name, cytobands_select$genes_in_peak),
                               drawConnectors = T, max.overlaps = Inf, legendPosition = "top", gridlines.minor = F, gridlines.major = F,
                               raster = TRUE)
rna_volcano

## Correlation plots ----

rna_features <- as.data.frame(t(rna_mat_norm))
colnames(rna_features) <- rna_tcga$patient
rna_features$gene <- rownames(rna_features)

rna_features <- rna_features %>% filter(gene %in% met_features_genes$UCSC_RefGene_Name)
rna_features <- subset(rna_features, select = -gene)
rna_features <- as.data.frame(t(rna_features))
rna_features$patien <- rownames(rna_features) %>% str_replace_all(pattern = "-", "_")

Bval <- as.data.frame(t(subset(met_tcga, select = -c(patient, cluster))))
colnames(Bval) <- met_tcga$patient %>% str_replace_all(pattern = "-", replacement = "_")

Bval$cpg <- rownames(Bval)

Bval_features <- Bval %>% filter(cpg %in% met_features_genes$cpg_name)
Bval_features <- subset(Bval_features, select = -cpg)
Bval_features <- as.data.frame(t(Bval_features))
Bval_features$patien <- rownames(Bval_features)

Bval_features <- Bval_features %>% filter(patien %in% rna_features$patien)

table(Bval_features$patien == rna_features$patien)

shapiro.test(Bval_features[,"cg07095230"])
ggpubr::ggqqplot(Bval_features[,"cg07095230"])
shapiro.test(rna_features[,"TBX2"])
ggpubr::ggqqplot(rna_features[,"TBX2"])

cor.test(x = Bval_features[,"cg07095230"], y = rna_features[,"TBX2"], method = "kendall")
plot(x = Bval_features[,"cg07095230"], y = rna_features[,"TBX2"])

TBX2 <- data.frame(Bval = Bval_features[,"cg07095230"], RSEM_TPM = rna_features[,"TBX2"])
TBX2_cor <- ggpubr::ggscatter(data = TBX2, x = "Bval", y = "RSEM_TPM", 
                              add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "slateblue",
                              xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "Correlation RNA - Methylation TBX2") + theme_cowplot()

 "#db1456"

MEOX2 <- data.frame(Bval = Bval_features[,"cg00839579"], RSEM_TPM = rna_features[,"MEOX2"])
MEOX2_cor <- ggpubr::ggscatter(data = MEOX2, x = "Bval", y = "RSEM_TPM", 
                               add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "slateblue",
                               xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "Correlation RNA - Methylation MEOX2") + theme_cowplot()

### C6 GSEA Methylomics ----

library(missMethyl)

OncoSig <- GSA::GSA.read.gmt("data/data_omics/c6.all.v2025.1.Hs.entrez.gmt")
OncoSig.genesets <- OncoSig %>% pluck("genesets")
names(OncoSig.genesets) <- OncoSig$geneset.names
gsa_OncoSig<- gsameth(
  sig.cpg = DCpGs$cpg_name,
  all.cpg = DMPs$cpg_name,
  collection = OncoSig.genesets,
  array.type = "450K",
  plot.bias = TRUE,
  sig.genes = TRUE
)
topOncoSig<- topGSA(gsa_OncoSig)
topOncoSig$logFDR <- -log10(topOncoSig$FDR)
topOncoSig$propDE <- topOncoSig$DE/topOncoSig$N
topOncoSig$OncoSig <- rownames(topOncoSig)
SigGenesInSet <- strsplit(topOncoSig$SigGenesInSet, ",")
names(SigGenesInSet) <- topOncoSig$OncoSig
topOncoSig$SigGenesInSet <- SigGenesInSet
topOncoSig_unnest <- topOncoSig %>% unnest_longer(SigGenesInSet)

c6_methyl <- topOncoSig %>% filter(FDR < 0.1) %>% arrange(propDE) %>% 
  mutate(OncoSig = factor(OncoSig, levels = OncoSig)) %>% 
  ggplot(aes(x = propDE, y = OncoSig)) + 
  geom_segment(aes(xend = 0, yend = OncoSig)) +
  geom_point(aes(fill = FDR, size = DE/N), color = "gray50", shape = 21) +
  scale_fill_gradient(name = "FDR", high = "white", low = "#fc0352") +
  labs(title = "Significantly Enriched Oncogenic Gene Sets - Methylomics",
       subtitle = "FDR < 0.1") + ylab("") + theme_cowplot(font_size = 14) +
  scale_size_continuous(name = "Gene Ratio") +
  xlab("Gene Ratio")
c6_methyl

topOncoSig_unnest %>% filter(FDR < 0.05) %>% arrange(propDE) %>% 
  mutate(OncoSig = factor(OncoSig, levels = unique(OncoSig))) %>% 
  ggplot(aes(x = propDE, y = OncoSig)) + 
  geom_text_repel(aes(label = SigGenesInSet), max.overlaps = Inf, segment.colour = "gray", segment.alpha = 0.5) +
  geom_point(aes(fill = FDR, size = DE/N), color = "gray50", shape = 21) +
  scale_fill_gradient(name = "FDR",low = "white", high = "#fc0352") +
  labs(title = "Significantly Enriched Oncogenic Gene Sets",
       subtitle = "FDR < 0.05") + ylab("") + theme_cowplot(font_size = 14) +
  scale_size_continuous(name = "Gene Ratio") +
  xlab("Gene Ratio")

### C6 GSEA RNA-Seq ----

enrich_cluster2 <- read_delim("data/data_omics/rna_seq/GSEA/C6_Oncogenic_Norm.Gsea.1753876686244/gsea_report_for_Cluster_2_1753876686244.tsv",
                              delim = "\t")
enrich_cluster1 <- read_delim("data/data_omics/rna_seq/GSEA/C6_Oncogenic_Norm.Gsea.1753876686244/gsea_report_for_Cluster_1_1753876686244.tsv",
                              delim = "\t")

enrich_gsea_c6 <- rbind(enrich_cluster2, enrich_cluster1)

enrich_gsea_c6$`LEADING EDGE` <- enrich_gsea_c6$`LEADING EDGE` %>% str_remove(pattern = "%.*") %>% str_remove("tags=")
enrich_gsea_c6$`LEADING EDGE` <- as.numeric(enrich_gsea_c6$`LEADING EDGE`)

c6_rna <- enrich_gsea_c6 %>% filter(`FDR q-val`< 0.25) %>%
  dplyr::arrange(NES) %>% mutate(NAME = factor(NAME, levels = NAME)) %>% ggplot(aes(x = NES, y = NAME)) +
  geom_segment(aes(xend = 0, yend = NAME), color = "gray50") +
  geom_point(aes(fill = `FDR q-val`, size = `LEADING EDGE`/100), shape = 21, color = "gray50") +
  geom_vline(xintercept = 0, color = "gray") +
  scale_fill_gradient(name = "FDR", high = "white", low = "#fc0352") +
  scale_size_area(name = "Gene Ratio") +
  labs(title = "Significantly Enriched Oncogenic Gene Sets - Transcriptomics",
       subtitle = "FDR < 0.25") + ylab("") + theme_cowplot()
c6_rna

## RPPA ----

keyvals <- ifelse(
  rppa_toptable$logFC < 0 & rppa_toptable$P.Value < 0.05, 'slateblue4',
  ifelse(rppa_toptable$logFC > 0 & rppa_toptable$P.Value < 0.05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'

DEP_volcano <- EnhancedVolcano(toptable = rppa_toptable, lab = rppa_toptable$peptide_target, x = "logFC", y = "P.Value", FCcutoff = 0,
                               xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                               colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 2 vs Cluster 1",
                               title = "Volcano Plot of Differentially Expressed/Phosphorylated Proteins",
                               boxedLabels = T, labCol = keyvals[names(keyvals) != "NS"], labSize = 4, ylim = c(0,5), xlim = c(-1,1),
                               drawConnectors = T,maxoverlapsConnectors = Inf, typeConnectors = "open",max.overlaps = Inf, legendPosition = "top", gridlines.minor = F, gridlines.major = F,
                               raster = TRUE)

kegg_ora <- kk_results_long %>% filter(p.adjust < 0.001) %>%
  arrange(RichFactor) %>% mutate(Description = factor(Description, levels = unique(Description))) %>% 
  ggplot(aes(x = RichFactor, y = Description)) +
  geom_segment(aes(xend = 0, yend = Description), color = "gray50") +
  geom_text_repel(aes(label = ALIAS, color = logFC), max.overlaps = Inf, fontface = "bold") +
  geom_point(aes(fill = -log10(p.adjust)), shape = 21, size = 8, color = "gray50") +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 FDR") +
  scale_color_gradient(low = "blue", high= "red", limits = c(-0.6,0.6)) +
  theme_cowplot() + labs(title = "KEGG Pathways Over-Representation Analysis",
                         subtitle = "FDR < 0.001")

# FIGURE ----

## Features Information Tables ----

# Load libraries for annotation of CpGs
library(annotate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


cpg_sites <- met_features_genes$cpg_name
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
names(annot)
head(annot)

cpg_info <- annot[cpg_sites, c("Name","chr", "pos", "strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
cpg_info <- as.data.frame(cpg_info)
openxlsx2::write_xlsx(cpg_info, "data/data_omics/methylomics/Fig_Table_A.xlsx")

load("data/data_omics/cna/cytobands.RData")
cytobands <- as.data.frame(cytobands)

tblB <- cytobands %>% dplyr::select(Descriptor, rowRanges, Unique.Name, type, genes_in_peak)
openxlsx2::write_xlsx(tblB, "data/data_omics/cna/Fig_Table_B.xlsx")

## Figures ----

cd <- cowplot::plot_grid(cna_violin, OR_cna, ncol = 2, labels = c("C", "D"), rel_widths = c(0.4, 0.6))
ef <- cowplot::plot_grid(dmp_volcano, rna_volcano, labels = c("E", "F"))
gh <- cowplot::plot_grid(c6_methyl, c6_rna, labels = c("G", "H"))

fig <- cowplot::plot_grid(NULL, NULL, cd,ef,gh, ncol = 1, rel_heights = c(0.1,0.1, 0.25, 0.25, 0.25))
ggsave(plot = fig, filename = "figures/omics_analysis/fig.pdf", height = 3*11.69, width = 3*8.27)

# SUPP FIGURE ----

abc <- cowplot::plot_grid(OR_cna_genes, MEOX2_cor, TBX2_cor, nrow = 1, labels = "AUTO")
de <- cowplot::plot_grid(DEP_volcano, kegg_ora, labels = c("D", "E"))
supp_fig <- plot_grid(abc, de, ncol = 1, rel_heights = c(0.4, 0.6))
ggsave("figures/omics_analysis/supp_fig.pdf", plot = supp_fig, height = 2*11.69, width = 3*8.27)

