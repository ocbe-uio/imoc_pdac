# 0. Load libraries ----

library(tidyverse)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(limma)
library(edgeR)
library(clusterProfiler)

## Methylation features ----

met_features <- openxlsx2::read_xlsx("results/omics_analysis/Methylation/met_features_genes.xlsx")

met_features$gene <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")

# Normalized Data ----

## 1. Load data ----

patient_cluster <- read_csv("results/cluster_analysis/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

rna_tcga <- read_csv("data/TCGA/omics_data/raw/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv")

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

## Assign cluster IDs

rna_tcga <- left_join(rna_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rna_tcga$cluster <- as.factor(rna_tcga$cluster)

rna_mat <- as.matrix(subset(rna_tcga, select = -c(cluster,patient)))

rownames(rna_mat) <- rna_tcga$patient

write.csv(rna_mat, "results/omics_analysis/RNA_Seq/rna_mat_norm.csv", row.names = T)


## 2. Files for GSEA ----

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


# Raw Counts ----

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
metadata <- curatedTCGAData(diseaseCode = TCGA_CODE, version = "2.1.1")

# filter by omics
omics <- c("RNASeq2Gene")
curatedTCGAData(diseaseCode = TCGA_CODE, assays = "*", version = "2.1.1")

# download data
cancer_data <- curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, dry.run = FALSE, version = "2.1.1")

# filter by primary tumor
sampleTables(cancer_data)
tums <- TCGAsampleSelect(barcodes = colnames(cancer_data), sampleCodes = "01")
cancer_data <- cancer_data[,tums,]

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

exportClass(cancer_data, dir = ".data/TCGA/omics_data/raw/", 
            fmt = "csv", ext = ".csv")

## 1. Load data ----

patient_cluster <- read_csv("results/cluster_analysis/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

c1 <- patient_cluster %>% filter(cluster == "Cluster_1")
c2 <- patient_cluster %>% filter(cluster == "Cluster_2")

table(c2$patient %in% c1$patient)

rna_tcga <- read_csv("data/TCGA/omics_data/raw/cancer_data_PAAD_RNASeq2Gene-20160128.csv")

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


## Assign cluster IDs

rna_tcga <- left_join(rna_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rna_tcga$cluster <- as.factor(rna_tcga$cluster)

write.csv(rna_tcga, "results/omics_analysis/RNA_Seq/rna_tcga.csv")

rna_mat <- as.matrix(subset(rna_tcga, select = -c(cluster,patient)))

rownames(rna_mat) <- rna_tcga$patient


## 2. Differential Expression Analysis ----

### 2.1 Limma ----

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

write.csv(top.table, "results/omics_analysis/RNA_Seq/limma_results.csv")

## Volcano Plot Methylation ----

keyvals <- ifelse(
  top.table$logFC < -0.25 & top.table$adj.P.Val < 0.05, 'slateblue4',
  ifelse(top.table$logFC > 0.25 & top.table$adj.P.Val < 0.05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'

met_features <- openxlsx2::read_xlsx("results/omics_analysis/Methylation/met_features_genes.xlsx")
met_features$UCSC_RefGene_Name <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")
EnhancedVolcano(toptable = top.table, lab = top.table$Gene, x = "logFC",
                y = "adj.P.Val", FCcutoff = 0.25,
                xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 1 vs Cluster 0",
                title = "Volcano Plot of Differentially Expressed Genes",
                boxedLabels = T, labSize = 4, selectLab = met_features$UCSC_RefGene_Name,
                drawConnectors = T, max.overlaps = Inf, legendPosition = "right", gridlines.minor = F, gridlines.major = F)

## Volcano Plot CNA

cna_features <- openxlsx2::read_xlsx("results/omics_analysis/CNA/cytobands_select.xlsx")
met_features$UCSC_RefGene_Name <- met_features$UCSC_RefGene_Name %>% str_remove(pattern = ";.*")
EnhancedVolcano(toptable = top.table, lab = top.table$Gene, x = "logFC",
                y = "adj.P.Val", FCcutoff = 0.25,
                xlab = "Log2 fold-change Normalized Expression", pCutoff = 0.05,
                colCustom = keyvals, ylab = "-log10 p-value", subtitle = "Cluster 1 vs Cluster 0",
                title = "Volcano Plot of Differentially Expressed Genes",
                boxedLabels = T, labSize = 4, selectLab = c(met_features$UCSC_RefGene_Name, cna_features$genes_in_peak),
                drawConnectors = T, max.overlaps = Inf, legendPosition = "right", gridlines.minor = F, gridlines.major = F)


