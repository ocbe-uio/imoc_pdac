# 0. Load libraries ----

library(tidyverse)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(limma)
library(EnhancedVolcano)
library(missMDA)
library(minfi)
library(openxlsx)

# Refs
# https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html#introduction
# https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#introduction

# 1. Load data ----
patient_cluster <- read_csv("data/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")

met_tcga <- read_csv("data/data_omics/TCGA/raw data/cancer_data_PAAD_Methylation-20160128.csv")

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

# Assign patients
met_tcga$patient <- patient

#rownames(met_tcga) <- met_tcga$patient

## Assign cluster IDs

met_tcga <- left_join(met_tcga, patient_cluster[,c(1,3)], by = join_by(patient == patient))

met_tcga$cluster <- as.factor(met_tcga$cluster)

write_csv(met_tcga, file = "data/data_omics/methylomics/met_tcga.csv")

## 1.2 Methylation features ----

met_features <- c("cg06785999",	"cg07095230",	"cg00839579")

# 2. Quick EDA ----

# PCA and HPCP

met_mat <- as.matrix(subset(met_tcga, select = -c(cluster,patient)))

rownames(met_mat) <- met_tcga$patient

str(met_mat)

#nb <- estim_ncpPCA(met_mat, ncp.max = 5, verbose = TRUE)
res.comp <- imputePCA(met_mat,ncp=15)
met_pca <- PCA(res.comp$completeObs)

#met_pca <- PCA(met_mat)

p2 <- fviz_pca_ind(met_pca, habillage = met_tcga$cluster, repel = TRUE, addEllipses = T, label = "none",
             mean.point = FALSE) + theme_cowplot()

met_hcpc <- HCPC(met_pca, graph = FALSE)


p1 <- fviz_dend(met_hcpc, rect = T, rect_fill = T, cex = 0.5, ggtheme = cowplot::theme_cowplot())

p3 <- fviz_cluster(met_hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = FALSE, # Show cluster centers
             #palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = cowplot::theme_cowplot(),
             main = "Factor map", geom = "point")

plot_grid(p1, p2, p3)

# Distribution of beta values

met_tcga_long <- pivot_longer(subset(met_tcga, select = -patient), cols = -cluster, names_to = "patient", values_to = "beta")

met_tcga_long %>% ggplot(aes(x = beta)) + geom_density(aes(color = cluster, fill = cluster), alpha = 0.1, lwd = 1) + theme_cowplot() +
  xlab("Beta Values") + ylab("Density") + labs(title = "Distribution of Beta Values")

# 2. Differential Methylation Analysis ----

# For the DMA we need to have the data as a matrix-like data object containing beta values, with rows corresponding to cg IDs and columns to samples.

Bval <- as.data.frame(t(subset(met_tcga, select = -c(patient, cluster))))
colnames(Bval) <- met_tcga$patient %>% str_replace_all(pattern = "-", replacement = "_")

# Convert the beta values in M values

Mval <- log2(Bval/(1-Bval))

# The clusters are our factor of interest
cluster <- met_tcga$cluster

# use the above to create a design matrix
design <- model.matrix(~0+cluster, data = met_tcga)
colnames(design) <- c(levels(cluster))

# fit the linear model 
fit <- lmFit(object = Mval, design = design)

# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(Cluster_2 - Cluster_1, levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

barplot(height = c(68492, 57253, 270320), names.arg = c("Up", "Down", "NotSig"), 
        col = c("pink", "lightblue", "gray80"), border = c("firebrick3", "slateblue4", "gray50"),
        ylab = "Differentially Methylated CpGs Count", main = "DM CpGs at FDR < 0.05", ylim = c(0,300000))


# get the table of results for the first contrast (naive - rTreg)
DMPs <- topTable(fit2, num=Inf)
head(DMPs)

DMPs$cpg_name <- rownames(DMPs)

openxlsx::write.xlsx(DMPs, file = "data/data_omics/methylomics/DMPs.xlsx")

## DMPs ----

DMPs <- openxlsx::read.xlsx("data/data_omics/methylomics/DMPs.xlsx")

# Plot most significant differentially methylated CpG

par(mfrow = c(2,2))
minfi::plotCpg(Bval, cpg = rownames(DMPs)[1:4], pheno=cluster, ylab = "Beta values")
par(mfrow = c(1,1))

## Volcano Plot ----

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  DMPs$logFC < -2 & DMPs$adj.P.Val < 1e-05, 'slateblue4',
  ifelse(DMPs$logFC > 2 & DMPs$adj.P.Val < 1e-05, 'firebrick3',
         'gray50'))
keyvals[is.na(keyvals)] <- 'gray50'
names(keyvals)[keyvals == 'firebrick3'] <- 'Up'
names(keyvals)[keyvals == 'slateblue4'] <- 'Down'
names(keyvals)[keyvals == 'gray50'] <- 'NS'


rownames(DMPs) <- DMPs$cpg_name

EnhancedVolcano(toptable = DMPs, lab = rownames(DMPs), x = "logFC", y = "adj.P.Val", FCcutoff = 2,
                xlab = "Log2 fold-change M value",
                colCustom = keyvals, ylab = "-log10 adjusted p-value", subtitle = "Cluster 1 vs Cluster 0",
                title = "Volcano Plot of Differentially Methylated CpG sites",
                selectLab = met_features,
                boxedLabels = T,
                drawConnectors = TRUE, max.overlaps = Inf, legendPosition = "right", gridlines.minor = F, gridlines.major = F)


# 3. Gene Enrichment ----

# https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#additional-analyses 

library(missMethyl)

DCpGs <- DMPs %>% filter(abs(logFC) > 2 & adj.P.Val < 1e-05)

# Load libraries for annotation of DCpGs
library(annotate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Convert CpG IDs to genomic coordinates
cpg_sites <- DCpGs$cpg_name
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
names(annot)

gene_info <- annot[cpg_sites, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_RefGene_Accession")]

gene_info <- as.data.frame(gene_info)

## DCpGs ----

DCpGs <- left_join(DCpGs, gene_info, by = join_by(cpg_name == Name))
openxlsx::write.xlsx(DCpGs, "data/data_omics/methylomics/DCpGs.xlsx")
DCpGs <- openxlsx::read.xlsx("results/methylomics/DCpGs.xlsx")

met_features_genes <- DCpGs %>% filter(cpg_name %in% met_features)
openxlsx2::write_xlsx(met_features_genes, file = "data/data_omics/methylomics/met_features_genes.xlsx")


