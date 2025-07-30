# 0. Load libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(clusterProfiler)

## Methylation features ----

met_features <- openxlsx2::read_xlsx("data/data_omics/methylomics/met_features_genes.xlsx")

met_features$gene <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")

# Normalized Data ----

## 1. Load data ----

patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

rna_tcga <- read_csv("data/data_omics/TCGA/raw data/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv")

rna_tcga <- as.data.frame(t(rna_tcga))

colnames(rna_tcga) <- rna_tcga[1,]
rna_tcga <- rna_tcga[-1,]

# Retrieve patients names

patient <- rownames(rna_tcga)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with 0 values in more than 20% of patient

rna_tcga <- apply(rna_tcga, 2, as.numeric)


rna_tcga <- rna_tcga[, which(!colMeans((rna_tcga == 0)) > 0.2)]

rna_tcga <- as.data.frame(rna_tcga)

rna_tcga$patient <- patient

#rownames(met_tcga) <- met_tcga$patient

## Assign cluster IDs

rna_tcga <- left_join(rna_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rna_tcga$cluster <- as.factor(rna_tcga$cluster)

## 2. Quick EDA ----

# PCA and HPCP

rna_mat <- as.matrix(subset(rna_tcga, select = -c(cluster,patient)))
write.csv(rna_mat, "data/data_omics/rna_seq/rna_mat_norm.csv")
rna_mat_norm <- rna_mat
rownames(rna_mat) <- rna_tcga$patient

barplot(rowSums(rna_mat), col = "firebrick2", 
        main = "RNASeg2GeneNorm", sub = "Histogram of Counts by Patient")

EDASeq::plotRLE(t(rna_mat), outline = FALSE,
                col = as.numeric(rna_tcga$cluster),
                main = "TCGA RNASeq2GeneNorm")


rna_pca <- PCA(rna_mat)

p2 <- fviz_pca_ind(rna_pca, habillage = rna_tcga$cluster, repel = TRUE, addEllipses = T, label = "none",
                   mean.point = FALSE) + theme_cowplot()

rna_hcpc <- HCPC(rna_pca, graph = FALSE)


p1 <- fviz_dend(rna_hcpc, rect = T, rect_fill = T, cex = 0.5, ggtheme = cowplot::theme_cowplot())

p3 <- fviz_cluster(rna_hcpc,
                   repel = TRUE,            # Avoid label overlapping
                   show.clust.cent = FALSE, # Show cluster centers
                   #palette = "jco",         # Color palette see ?ggpubr::ggpar
                   ggtheme = cowplot::theme_cowplot(),
                   main = "Factor map", geom = "point")

plot_grid(p1, p2, p3)

# Distribution of counts

rna_tcga_long <- pivot_longer(subset(rna_tcga, select = -patient), cols = -cluster, names_to = "patient", values_to = "TPM")

rna_tcga_long %>% ggplot(aes(x = log(TPM))) + geom_density()


# Heatmap of most variable genes

# TPM values
tpm <- as.data.frame(t(subset(rna_tcga, select = -c(patient, cluster))))
colnames(tpm) <- patient

#compute the variance of each gene across samples
V <- apply(tpm, 1, var)
#sort the results by variance in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:1000])

annotation_col <- data.frame(Cluster = rna_tcga$cluster)
rownames(annotation_col) <- colnames(tpm)


pheatmap::pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE, annotation_col = annotation_col, treeheight_row = 0,
         main = "Heatmap of 1000 most variable genes", treeheight_col = 10, color = hcl.colors(50, palette = "RdBu", rev = T),
         width = 10, height = 10, fontsize_col = 5)

## 3. Differential Expression Analysis ----

# As these are normalized TPM values, we cannot perform differential expression
# analysis using limma or DESeq2. But, we can do GSEA analysis using permutation
# of samples: https://docs.gsea-msigdb.org/#GSEA/GSEA_and_RNA-Seq/

# Files for GSEA

countsNormalized <- t(rna_mat)

gsea_df <- add_column(as.data.frame(countsNormalized), NAME = rownames(countsNormalized), DESCRIPTION = NA, .before = 1)

write_delim(gsea_df, file = "data/data_omics/rna_seq/GSEA/NormCounts.txt", delim = "\t")

phenotype <- data.frame(patient = colnames(countsNormalized))
phenotype <- left_join(phenotype, subset(rna_tcga, select = c("patient", "cluster")),
                       by = join_by(patient == patient))

c1 <- phenotype %>% filter(cluster == "Cluster_1")
c2 <- phenotype %>% filter(cluster == "Cluster_2")

table(c2$patient %in% c1$patient)


# Use phenotype table to create phenotypes in GSEA "on the fly"

write_csv(phenotype, file = "data/data_omics/rna_seq/GSEA/NormCounts_phenotype.csv")


# "Raw" Counts ----

## 0. Download data ----

#BiocManager::install("TCGAutils")
#BiocManager::install("curatedTCGAData")

library(TCGAutils)
library(curatedTCGAData)

# choose the code
TCGA_CODE <- "PAAD"
CANCER_CODE <- "PDAC"
PLATFORM_CODE <- "TCGA"

# check experiments
metadata <- curatedTCGAData(diseaseCode = TCGA_CODE, version = "2.0.1")

# filter by omics
omics <- c("RNASeqGene")
curatedTCGAData(diseaseCode = TCGA_CODE, assays = "*", version = "2.1.1")

# download data
cancer_data <- curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, dry.run = FALSE, version = "2.0.1")

# filter by primary tumor
sampleTables(cancer_data)
tums <- TCGAsampleSelect(barcodes = colnames(cancer_data), sampleCodes = "01")
cancer_data <- cancer_data[,tums,]

# visualize data
upsetSamples(cancer_data)

# clinical data
clinical_data <- colData(cancer_data)[, c('patient.patient_id', 'histological_type')]
dim(clinical_data)

# remove histological nan
table(colData(cancer_data)$histological_type)
clinical_data <- clinical_data[!is.na(clinical_data$histological_type),]

# select histology
histological_type <- names(sort(table(clinical_data$histological_type), decreasing = T))[1]

# select patients
patients <- rownames(clinical_data[clinical_data$histological_type == histological_type,])

# filter patients by histology
cancer_data <- cancer_data[,colnames(cancer_data)[substr(colnames(cancer_data), 1, 12) %in% patients],]

exportClass(cancer_data, dir = "./data/data_omics/TCGA/TCGA/raw data/", 
            fmt = "csv", ext = ".csv")

## 1. Load data ----

patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

c1 <- patient_cluster %>% filter(cluster == "Cluster_1")
c2 <- patient_cluster %>% filter(cluster == "Cluster_2")

table(c2$patient %in% c1$patient)

rna_tcga <- read_csv("data/data_omics/TCGA/raw data/cancer_data_PAAD_RNASeq2Gene-20160128.csv")

rna_tcga <- as.data.frame(t(rna_tcga))

colnames(rna_tcga) <- rna_tcga[1,]
rna_tcga <- rna_tcga[-1,]

# Retrieve patients names

patient <- rownames(rna_tcga)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with 0 values in more than 20% of patient

rna_tcga <- apply(rna_tcga, 2, as.numeric)


rna_tcga <- rna_tcga[, which(!colMeans((rna_tcga == 0)) > 0.2)]

rna_tcga <- as.data.frame(rna_tcga)

rna_tcga$patient <- patient

#rownames(met_tcga) <- met_tcga$patient

## Assign cluster IDs

rna_tcga <- left_join(rna_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rna_tcga$cluster <- as.factor(rna_tcga$cluster)

write.csv(rna_tcga, "data/data_omics/rna_seq/rna_tcga.csv")

## 2. Quick EDA ----

# PCA and HPCP

rna_mat <- as.matrix(subset(rna_tcga, select = -c(cluster,patient)))

rownames(rna_mat) <- rna_tcga$patient

barplot(rowSums(rna_mat), col = "slateblue2", 
        main = "RNASeg2Gene", sub = "Histogram of Counts by Patient")

range(rna_mat)

str(rna_mat)

rna_pca <- PCA(rna_mat)

p2 <- fviz_pca_ind(rna_pca, habillage = rna_tcga$cluster, repel = TRUE, addEllipses = T, label = "none",
                   mean.point = FALSE) + theme_cowplot()

rna_hcpc <- HCPC(rna_pca, graph = FALSE)


p1 <- fviz_dend(rna_hcpc, rect = T, rect_fill = T, cex = 0.5, ggtheme = cowplot::theme_cowplot())

p3 <- fviz_cluster(rna_hcpc,
                   repel = TRUE,            # Avoid label overlapping
                   show.clust.cent = FALSE, # Show cluster centers
                   #palette = "jco",         # Color palette see ?ggpubr::ggpar
                   ggtheme = cowplot::theme_cowplot(),
                   main = "Factor map", geom = "point")

plot_grid(p1, p2, p3)

# Distribution of counts

rna_tcga_long <- pivot_longer(subset(rna_tcga, select = -patient), cols = -cluster, names_to = "Gene", values_to = "Counts")

rna_tcga_long %>% ggplot(aes(x = log(Counts))) + geom_density()

# Heatmap of most variable genes

# TPM values
tpm <- as.data.frame(t(subset(rna_tcga, select = -c(patient, cluster))))
colnames(tpm) <- patient

#compute the variance of each gene across samples
V <- apply(tpm, 1, var)
#sort the results by variance in decreasing order 
#and select the top 100 genes 
selectedGenes <- names(V[order(V, decreasing = T)][1:1000])

annotation_col <- data.frame(Cluster = rna_tcga$cluster)
rownames(annotation_col) <- colnames(tpm)


pheatmap::pheatmap(tpm[selectedGenes,], scale = 'row', show_rownames = FALSE, annotation_col = annotation_col, treeheight_row = 0,
         main = "Heatmap of 1000 most variable genes", treeheight_col = 10, color = hcl.colors(50, palette = "RdBu", rev = T),
         fontsize_col = 5)

## 3. Differential Expression Analysis ----

### 3.2 Limma ----

# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html


#define the read count table (rows = features, columns = patients)
countData <- t(rna_mat)

#define the experimental setup 
colData <- rna_tcga %>% dplyr::select(cluster)
rownames(colData) <- rna_tcga$patient

# Create DGEList object

d0 <- DGEList(countData)

# Calculate normalization factors

d0 <- calcNormFactors(d0)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)

# Multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(colData$cluster))

#Specify the model to be fitted. We do this before using voom since voom uses 
#variances of the model residuals (observed - fitted)

mm <- model.matrix(~0 + colData$cluster)
colnames(mm) <- levels(rna_tcga$cluster)

# Voom

y <- voom(d, mm, plot = T)

tmp <- voom(d0, mm, plot = T)

#  Fitting linear models in limma

fit <- lmFit(y, mm)
head(coef(fit))

# contrast matrix

contMatrix <- makeContrasts(Cluster_2 - Cluster_1, levels=levels(rna_tcga$cluster))
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# DEGs
top.table <- topTable(fit2, sort.by = "P", n = Inf)

top.table$Gene <- rownames(top.table)

write.csv(top.table, "data/data_omics/rna_seq/limma_results.csv")

## Volcano Plot Methylation ----

keyvals <- ifelse(
  top.table$logFC < -0.25 & top.table$adj.P.Val < 0.05, 'slateblue4',
  ifelse(top.table$logFC > 0.25 & top.table$adj.P.Val < 0.05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'

met_features <- openxlsx2::read_xlsx("data/data_omics/methylomics/met_features_genes.xlsx")
met_features$UCSC_RefGene_Name <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")
EnhancedVolcano(toptable = top.table, lab = top.table$Gene, x = "logFC",
                y = "adj.P.Val", FCcutoff = 0.25,
                xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 1 vs Cluster 0",
                title = "Volcano Plot of Differentially Expressed Genes",
                boxedLabels = T, labSize = 4, selectLab = met_features$UCSC_RefGene_Name,
                drawConnectors = T, max.overlaps = Inf, legendPosition = "right", gridlines.minor = F, gridlines.major = F)

## Volcano Plot CNA

cna_features <- openxlsx2::read_xlsx("data/data_omics/cna/cytobands_select.xlsx")
met_features$UCSC_RefGene_Name <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")
EnhancedVolcano(toptable = top.table, lab = top.table$Gene, x = "logFC",
                y = "adj.P.Val", FCcutoff = 0.25,
                xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 1 vs Cluster 0",
                title = "Volcano Plot of Differentially Expressed Genes",
                boxedLabels = T, labSize = 4, selectLab = c(met_features$UCSC_RefGene_Name, cna_features$genes_in_peak),
                drawConnectors = T, max.overlaps = Inf, legendPosition = "right", gridlines.minor = F, gridlines.major = F)


## 4. Correlation Methylation - Expression ----

rna_features <- as.data.frame(t(rna_mat_norm))
colnames(rna_features) <- rna_tcga$patient
rna_features$gene <- rownames(rna_features)

rna_features <- rna_features %>% filter(gene %in% met_features$gene)
rna_features <- subset(rna_features, select = -gene)
rna_features <- as.data.frame(t(rna_features))
rna_features$patien <- rownames(rna_features) %>% str_replace_all(pattern = "-", "_")

### Methylation Data ----
patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")


met_tcga <- read_csv("data/imo_clustering/tests_jupyter/TCGA/raw data/cancer_data_PAAD_Methylation-20160128.csv")

met_tcga <- as.data.frame(t(met_tcga))

colnames(met_tcga) <- met_tcga[1,]
met_tcga <- met_tcga[-1,]

# Retrieve patients names

patient <- rownames(met_tcga)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with more NAs in more than 20% of patient

met_tcga <- met_tcga[, which(!colMeans(is.na(met_tcga)) > 0.2)]

met_tcga <- apply(met_tcga, 2, as.numeric)
met_tcga <- as.data.frame(met_tcga)

str(met_tcga)

met_tcga$patient <- patient

## Assign cluster IDs

met_tcga <- left_join(met_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

met_tcga$cluster <- as.factor(met_tcga$cluster)
met_mat <- as.matrix(subset(met_tcga, select = -c(cluster,patient)))

Bval <- as.data.frame(t(subset(met_tcga, select = -c(patient, cluster))))
colnames(Bval) <- met_tcga$patient %>% str_replace_all(pattern = "-", replacement = "_")

# Convert the beta values in M values

Mval <- log2(Bval/(1-Bval))

Bval$cpg <- rownames(Bval)
Mval$cpg <- rownames(Mval)

Bval_features <- Bval %>% filter(cpg %in% met_features$cpg_name)
Mval_features <- Mval %>% filter(cpg %in% met_features$cpg_name)

Bval_features <- Bval %>% filter(cpg %in% met_features$cpg_name)
Bval_features <- subset(Bval_features, select = -cpg)
Bval_features <- as.data.frame(t(Bval_features))
Bval_features$patien <- rownames(Bval_features)

Bval_features <- Bval_features %>% filter(patien %in% rna_features$patien)

Mval_features <- Mval %>% filter(cpg %in% met_features$cpg_name)
Mval_features <- subset(Mval_features, select = -cpg)
Mval_features <- as.data.frame(t(Mval_features))
Mval_features$patien <- rownames(Mval_features)

Mval_features <- Mval_features %>% filter(patien %in% rna_features$patien)

table(Bval_features$patien == rna_features$patien)

shapiro.test(Bval_features[,"cg07095230"])
ggpubr::ggqqplot(Bval_features[,"cg07095230"])
shapiro.test(rna_features[,"TBX2"])
ggpubr::ggqqplot(rna_features[,"TBX2"])

cor.test(x = Bval_features[,"cg07095230"], y = rna_features[,"TBX2"], method = "kendall")
plot(x = Bval_features[,"cg07095230"], y = rna_features[,"TBX2"])



shapiro.test(Bval_features[,"cg00839579"])
ggpubr::ggqqplot(Bval_features[,"cg00839579"])
shapiro.test(rna_features[,"MEOX2"])
ggpubr::ggqqplot(rna_features[,"MEOX2"])

cor.test(x = Bval_features[,"cg00839579"], y = rna_features[,"MEOX2"], method = "kendall")
plot(x = Bval_features[,"cg00839579"], y = rna_features[,"MEOX2"])

shapiro.test(Mval_features[,"cg07095230"])
ggpubr::ggqqplot(Mval_features[,"cg07095230"])
shapiro.test(rna_features[,"TBX2"])
ggpubr::ggqqplot(rna_features[,"TBX2"])

cor.test(x = Mval_features[,"cg07095230"], y = rna_features[,"TBX2"], method = "kendall")
plot(x = Mval_features[,"cg07095230"], y = rna_features[,"TBX2"])





shapiro.test(Mval_features[,"cg00839579"])
ggpubr::ggqqplot(Mval_features[,"cg00839579"])
shapiro.test(rna_features[,"MEOX2"])
ggpubr::ggqqplot(rna_features[,"MEOX2"])

cor.test(x = Mval_features[,"cg00839579"], y = rna_features[,"MEOX2"], method = "kendall")
plot(x = Mval_features[,"cg00839579"], y = rna_features[,"MEOX2"])


TBX2 <- data.frame(Bval = Bval_features[,"cg07095230"], Mval = Mval_features[,"cg07095230"], RSEM_TPM = rna_features[,"TBX2"])
ggpubr::ggscatter(data = TBX2, x = "Bval", y = "RSEM_TPM", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "slateblue",
                  xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "TBX2") + theme_cowplot()

ggpubr::ggscatter(data = TBX2, x = "Mval", y = "RSEM_TPM", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "slateblue",
                  xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "TBX2") + theme_cowplot()



MEOX2 <- data.frame(Bval = Bval_features[,"cg00839579"], Mval = Mval_features[,"cg00839579"], RSEM_TPM = rna_features[,"MEOX2"])
ggpubr::ggscatter(data = MEOX2, x = "Bval", y = "RSEM_TPM", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "#db1456",
                  xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "MEOX2") + theme_cowplot()
ggpubr::ggscatter(data = MEOX2, x = "Mval", y = "RSEM_TPM", 
                  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "kendall", color = "#db1456",
                  xlab = "Methylation B values", ylab = "RSEM TPM Normalized Counts") + 
  labs(title = "MEOX2") + theme_cowplot()
