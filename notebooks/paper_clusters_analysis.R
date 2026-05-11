# 0. Load libraries ----

library(tidyverse)
library(cowplot)
library(ggsankey)
library(EnhancedVolcano)
library(limma)
library(edgeR)
library(pheatmap)
library(ggpubr)
library(GSVA)
library(ComplexUpset)
library(openxlsx2)

library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(survival)
library(survRM2)

# To run this script, open the imoc_pdac.Rproj file in the main folder (imoc_pdac)

# 1. Load data ----

papers_cluster_df <- read_csv("data/TCGA/comparison_papers/clusters_papers.csv")
papers_cluster_df$`My Clusters`[papers_cluster_df$`My Clusters` == 1] <- "Cluster 1"
papers_cluster_df$`My Clusters`[papers_cluster_df$`My Clusters` == 2] <- "Cluster 2"
papers_cluster_df <- papers_cluster_df %>% rename(IMOC = `My Clusters`)

patient_cluster <- read_csv("results/cluster_analysis/patient_clusters.csv")

patient_cluster <- patient_cluster %>% pivot_longer(cols = 2:ncol(patient_cluster), , values_to = "patient", names_to = "ID")

colnames(patient_cluster) <- c("cluster", "ID", "patient")

patient_cluster <- patient_cluster %>% filter(!is.na(patient))

patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "1", "2")
patient_cluster$cluster <- patient_cluster$cluster %>% str_replace(pattern = "0", "1")


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

## Clinical data for survival analysis

clinical_data_file <- read_delim("data/TCGA/omics_data/raw/cancer_data_PAAD_clinical_data.tsv", delim = "\t")


## Function for DGE Analysis

limma_contr <- function(data, cluster_data, patientCol, clusterCol, pairwise = TRUE, ref_level = NA, cutoff = 2) {
  
  # Prepare data
  data <- left_join(data, cluster_data[,c(patientCol,clusterCol)], by = patientCol) # Assign cluster IDs
  data <- data[!is.na(data[[clusterCol]]),]
  
  if (pairwise == TRUE) {
    data[[clusterCol]] <- as.factor(data[[clusterCol]]) # Convert to factor
  } else {
    if (is.na(ref_level)) {
      stop("'ref_level' must be defined when 'pairwise = FALSE'.")
    } 
    data[[clusterCol]] <- ifelse(data[[clusterCol]] == ref_level, ref_level, "Other") %>% factor(levels = c(ref_level, "Other"))
  } 
  
  # Prepare expression matrix
  rna_mat <- as.matrix(data[, !(names(data) %in% c(patientCol, clusterCol)), drop = FALSE]) # Convert to matrix; remove ID columns
  rownames(rna_mat) <- data[[patientCol]] # Add rownames
  
  # --Differential expression analysis--
  
  countData <- t(rna_mat) #define the read count table (rows = features, columns = patients)
  
  #define the experimental setup 
  colData <- data[[clusterCol]]
  names(colData) <- data[[patientCol]]
  
  # Create DGEList object
  
  
  
  d0 <- edgeR::DGEList(countData)
  
  # Calculate normalization factors
  
  d0 <- edgeR::calcNormFactors(d0)
  cutoff <- cutoff
  drop <- which(apply(edgeR::cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  cat(paste("Cutoff filtering \n","Genes:", nrow(d), "\nSamples:", ncol(d)))
  
  
  #Specify the model to be fitted. We do this before using voom since voom uses 
  #variances of the model residuals (observed - fitted)
  
  mm <- model.matrix(~ 0 + colData)
  colnames(mm) <- levels(data[[clusterCol]])
  
  # Voom
  
  
  tmp <- limma::voom(d0, mm, plot = T)
  title("\n\n All counts")
  
  y <- limma::voom(d, mm, plot = T)
  title(paste("\n\n Cutoff:", 2))
  
  #  Fitting linear models in limma
  fit <- limma::lmFit(y, mm)
  plotSA(fit, main="Final model: Mean-variance trend")
  
  # contrast matrix
  lvls <- colnames(coef(fit))   
  cmb <- combn(lvls, 2)
  x <- x <- setNames(
    apply(cmb, 2, paste, collapse = "-"),   # values: "A-B", "A-C", ...
    apply(cmb, 2, paste0, collapse = "vs")   # names:  "AvsB", "AvsC", ...
  )
  
  contMatrix <- makeContrasts(contrasts = x, levels=lvls)
  # fit the contrasts
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  # DGE analysis
  
  coefs <- colnames(contMatrix)
  DEGs <- list()
  
  for (i in 1:length(coefs)){
    top.table <- topTable(fit2, sort.by = "P", n = Inf, coef = coefs[i])
    top.table$Gene <- rownames(top.table)
    DEGs[[i]] <- top.table
  }
  
  names(DEGs) <- coefs
  
  return(DEGs)
}


# 2. Sankey plots ----

sankey_bailey <- papers_cluster_df %>% make_long(IMOC, bailey)

p_bailey <- ggplot(sankey_bailey, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", rep("#5ab4ac",4))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Bailey subtypes")

sankey_collisson <- papers_cluster_df %>% make_long(IMOC, collisson)

p_collisson <- ggplot(sankey_collisson, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4",rep("#d8b365",3))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Collisson subtypes")

sankey_moffitt <- papers_cluster_df %>% make_long(IMOC, moffitt)

p_moffitt <- ggplot(sankey_moffitt, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.33) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4", rep("#8c510a",2))) +
  geom_sankey_label(size = 4, color = "white", alpha = 0.75) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5)) +
  ggtitle("Moffitt subtypes")

FigXA <- plot_grid(p_bailey, p_collisson, p_moffitt, nrow = 1)

bailey_tbl <- table(papers_cluster_df$IMOC, papers_cluster_df$bailey)

rcompanion::pairwiseNominalIndependence(
  bailey_tbl,
  compare = "column",
  fisher = T,
  gtest = F,
  chisq = F,
  method = "fdr",
  cramer = TRUE
)

collisson_tbl <- table(papers_cluster_df$IMOC, papers_cluster_df$collisson)

rcompanion::pairwiseNominalIndependence(
  collisson_tbl,
  compare = "column",
  fisher = T,
  gtest = F,
  chisq = F,
  method = "fdr",
  cramer = TRUE
)

moffitt_tbl <- table(papers_cluster_df$IMOC, papers_cluster_df$moffitt)

rcompanion::pairwiseNominalIndependence(
  moffitt_tbl,
  compare = "column",
  fisher = T,
  gtest = F,
  chisq = F,
  method = "fdr",
  cramer = TRUE
)


# 3. Differential Gene Expression Analysis ----

## 3.1 IMOC Clusters ----

cluster_DEGs <- limma_contr(
  data = rna_tcga,
  cluster_data = patient_cluster,
  patientCol = "patient",
  clusterCol = "cluster",
  cutoff = 2
)

## 3.2 Bailey Classification ----

rna_tcga <- rename(rna_tcga, PatientID = patient)

bailey_DEGs <- list()

for (l in unique(papers_cluster_df$bailey)){
  print(l)
  DEGs <- limma_contr(
    data = rna_tcga,
    cluster_data = papers_cluster_df,
    patientCol = "PatientID",
    clusterCol = "bailey",
    pairwise = FALSE, 
    ref_level = l,
    cutoff = 2
  )
  bailey_DEGs <- c(bailey_DEGs, DEGs)
  rm(DEGs)
}


## 3.3 Collison Classification ----

collisson_DEGs <- list()

for (l in unique(papers_cluster_df$collisson)){
  print(l)
  DEGs <- limma_contr(
    data = rna_tcga,
    cluster_data = papers_cluster_df,
    patientCol = "PatientID",
    clusterCol = "collisson",
    pairwise = FALSE, 
    ref_level = l,
    cutoff = 2
  )
  collisson_DEGs <- c(collisson_DEGs, DEGs)
  rm(DEGs)
}


## 3.4 Moffitt Classification ----

papers_cluster_df$moffitt <- papers_cluster_df$moffitt %>% str_replace(pattern = "-", "_")

moffitt_DEGs <- limma_contr(
  data = rna_tcga,
  cluster_data = papers_cluster_df,
  patientCol = "PatientID",
  clusterCol = "moffitt",
  cutoff = 2
)

# 4. Common DEGs ----

## 4.1 Extract DEGs ----

Clust1vsClust2 <- cluster_DEGs$Cluster_1vsCluster_2 %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)

ADEX <- bailey_DEGs$ADEXvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
Squamous <- bailey_DEGs$SquamousvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
Progenitor <- bailey_DEGs$ProgenitorvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
Immunogenic <- bailey_DEGs$ImmunogenicvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)


Classical_collison <- collisson_DEGs$ClassicalvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
QM <- collisson_DEGs$QMvsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)
Exocrine <- collisson_DEGs$ExocrinevsOther %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)


Basal_likevsClassical <- moffitt_DEGs$Basal_likevsClassical %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.25)

bailey_DEGs$intersect_DEGs <- Reduce(intersect,
                                     list(
                                       ADEX$Gene,
                                       Squamous$Gene,
                                       Progenitor$Gene,
                                       Immunogenic$Gene
                                     ))

bailey_DEGs$union_DEGs <- Reduce(union,
                                 list(
                                   ADEX$Gene,
                                   Squamous$Gene,
                                   Progenitor$Gene,
                                   Immunogenic$Gene
                                 ))

collisson_DEGs$intersect_DEGs <- Reduce(intersect,
                                        list(
                                          Classical_collison$Gene,
                                          QM$Gene,
                                          Exocrine$Gene
                                        ))

collisson_DEGs$union_DEGs <- Reduce(union,
                                    list(
                                      Classical_collison$Gene,
                                      QM$Gene,
                                      Exocrine$Gene
                                    ))



## 4.2 Upset plot ----

sets <- list(
  "Cluster 1" = Clust1vsClust2$Gene,
  Squamous = Squamous$Gene,
  Progenitor = Progenitor$Gene,
  ADEX = ADEX$Gene, 
  Immunogenic = Immunogenic$Gene,
  QM = QM$Gene,
  Classical = Classical_collison$Gene,
  Exocrine = Exocrine$Gene,
  "Basal-Like" = Basal_likevsClassical$Gene
)

all_genes <- sort(unique(unlist(sets)))

# incidence (wide) table: one row per gene, TRUE/FALSE per comparison
incidence <- data.frame(
  gene = all_genes,
  lapply(sets, function(x) all_genes %in% x),
  check.names = FALSE
)

subtype_metadata <- data.frame(
  set = c("Cluster 1", "Squamous", "Progenitor", "ADEX", "Immunogenic", "QM", "Classical", "Exocrine", "Basal-Like"),
  Classification = c("IMOC", "Bailey", "Bailey", "Bailey", "Bailey", "Collisson", "Collisson", "Collisson", "Moffitt")
)

FigYA <- upset(
  incidence, 
  intersect = names(sets),
  name = "DEG set intersections", 
  base_annotations = list(
    "Intersection size (DEGs)" = intersection_size()
  ),
  set_sizes = upset_set_size(), 
  width_ratio = 0.1, height_ratio = 1, n_intersections = 15,
  stripes=upset_stripes(
    mapping=aes(color=Classification),
    colors=c(
      "IMOC"="#01665e",
      "Bailey"="#5ab4ac",
      "Collisson"="#d8b365",
      "Moffitt" = "#8c510a"
    ),
    data=subtype_metadata
  )
) + labs(title = "First 15 DEG intersections")


# 5. Normalized Data ----

## 5.1 Heatmaps ----

rna_tcga_norm <- read_csv("data/TCGA/omics_data/raw/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv")

rna_tcga_norm <- as.data.frame(t(rna_tcga_norm))

colnames(rna_tcga_norm) <- rna_tcga_norm[1,]
rna_tcga_norm <- rna_tcga_norm[-1,]

# Retrieve patients names

patient <- rownames(rna_tcga_norm)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

# Filter to remove features with 0 values in more than 20% of patient

rna_tcga_norm <- apply(rna_tcga_norm, 2, as.numeric)

rna_tcga_norm <- rna_tcga_norm[, which(!colMeans((rna_tcga_norm == 0)) > 0.2)]

rna_tcga_norm <- as.data.frame(rna_tcga_norm)

rna_tcga_norm$patient <- patient


## Assign cluster IDs

rna_tcga_norm <- left_join(rna_tcga_norm, patient_cluster[,c(1,3)], by = join_by(patient == patient))

rna_tcga_norm$cluster <- as.factor(rna_tcga_norm$cluster)

rna_tcga_norm <- left_join(rna_tcga_norm, papers_cluster_df[,c("PatientID", "bailey", "collisson", "moffitt")], by = join_by(patient == PatientID))


rna_tcga_norm_long <- rna_tcga_norm %>% pivot_longer(cols = 1:17123, names_to = "Gene", values_to = "Expression")
# Mean expression accross subtypes/clusters

rna_avg_clusters <- rna_tcga_norm_long %>% 
  group_by(cluster, Gene) %>% 
  summarise(avg = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = avg)

rna_avg_bailey <- rna_tcga_norm_long %>% filter(!is.na(bailey)) %>%
  group_by(bailey, Gene) %>% 
  summarise(avg = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bailey, values_from = avg)

rna_avg_collisson <- rna_tcga_norm_long %>% filter(!is.na(collisson)) %>%
  group_by(collisson, Gene) %>% 
  summarise(avg = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = collisson, values_from = avg)

rna_avg_collisson <- rename(rna_avg_collisson, Classical_Collisson = Classical)


rna_avg_moffitt <- rna_tcga_norm_long %>% filter(!is.na(moffitt)) %>%
  group_by(moffitt, Gene) %>% 
  summarise(avg = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = moffitt, values_from = avg)
rna_avg_moffitt <- rename(rna_avg_moffitt, Classical_Moffitt = Classical)


rna_avg <- rna_avg_clusters %>% 
  left_join(rna_avg_bailey, by = "Gene") %>%
  left_join(rna_avg_collisson, by = "Gene") %>%
  left_join(rna_avg_moffitt, by = "Gene")

rna_avg <- column_to_rownames(rna_avg, var = "Gene")
rna_avg <- as.matrix(rna_avg)

annotation_col <- data.frame(
  Classification = factor(c(rep("IMOC", 2), rep("Bailey", 4), rep("Collisson", 3), rep("Moffitt", 2))),
  row.names = colnames(rna_avg)
)

annotation_colors <- list(
  Classification = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")
)
names(annotation_colors) <- colnames(annotation_col)


union_genes <-  Reduce(union,
                       list(
                         Clustering = Clust1vsClust2$Gene,
                         Squamous = Squamous$Gene,
                         Progenitor = Progenitor$Gene,
                         ADEX = ADEX$Gene, 
                         Immunogenic = Immunogenic$Gene,
                         QM = QM$Gene,
                         Classical_Collisson = Classical_collison$Gene,
                         Exocrine = Exocrine$Gene,
                         Basal_Like = Basal_likevsClassical$Gene
                       ))

rna_avg_union <-  rna_avg[rownames(rna_avg) %in% union_genes,]
dim(rna_avg_union)



(FigXC <- pheatmap(t(rna_avg),
                   scale = "column", 
                   #cluster_rows = F,
                   color =  colorRampPalette(c("navy",  "white", "firebrick2")) (1000), border_color = "white", 
                   annotation_row = annotation_col, angle_col = 45,
                   annotation_colors = annotation_colors, 
                   cutree_rows = 4, 
                   clustering_method = "average", 
                   show_colnames = F,
                   treeheight_col = 0,
                   main = "Centroid gene expression heatmap (whole transcriptome)",
                   fontsize = 16
))

FigYB <- pheatmap(t(rna_avg_union),
                   scale = "column",
                   #cluster_rows = F,
                   color =  colorRampPalette(c("navy",  "white", "firebrick2")) (1000), border_color = "white",
                   annotation_row = annotation_col, 
                   annotation_colors = annotation_colors, 
                   cutree_rows = 4, 
                   clustering_method = "average", 
                   show_colnames = F,
                   treeheight_col = 0, angle_col = 45,
                   main = "Centroid gene expression heatmap (union of DEGs)",
                   fontsize = 16
                  )


## 5.2 Correlation to Centroids ----

# Fig. XD
rna_avg_z <- t(scale(t(rna_avg)))
cor_mat_z <- cor(rna_avg_z, method = "spearman", use = "pairwise.complete.obs")

col <- colorRampPalette(c("#4450AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

testRes = corrplot::cor.mtest(cor_mat_z, conf.level = 0.95)

# We need to save the file in order to use it with `cowplot::plot_grid()`

svg(file = "figures/omics_analysis/CorPlotAllGenes.svg", width = 10, height = 10)
p1 <- corrplot::corrplot(cor_mat_z, method = "circle", type = "upper", col = col(200),
                         tl.srt = 45, tl.col = "black",
                         diag = FALSE, title = "Centroid correlation (whole transcriptome)", 
                         p.mat = testRes$p, insig='blank',
                         number.cex = 0.8, 
                         mar = c(1,1,4,1),
                         pch.cex = 3, tl.cex = 1.5, cl.cex = 1.5)$corrPos
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()

#. Fig. YD

rna_avg_z <- t(scale(t(rna_avg_union)))
cor_mat_z <- cor(rna_avg_z, method = "spearman", use = "pairwise.complete.obs")

col <- colorRampPalette(c("#4450AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

testRes = corrplot::cor.mtest(cor_mat_z, conf.level = 0.95)

svg(file = "figures/omics_analysis/CorPlotUDEGs.svg", width = 10, height = 10)
p1 <- corrplot::corrplot(cor_mat_z, method = "circle", type = "upper", col = col(200),
                         tl.srt = 45, tl.col = "black",
                         diag = FALSE, title = "Centroid correlation (union of DEGs)",
                         p.mat = testRes$p, insig='blank',
                         number.cex = 0.8,
                         mar = c(1,1,4,1),
                         pch.cex = 3, tl.cex = 1.5, cl.cex = 1.5)$corrPos
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()

## 5.3 MDS ----

# Calculate distances from correlations
rna_avg_z <- t(scale(t(rna_avg)))
cor_mat_z <- cor(rna_avg_z, method = "spearman", use = "pairwise.complete.obs")

X <- 1-cor_mat_z
# D <- as.dist(X)
D <- as.dist(sqrt(2 * X)) # Euclidian distance

(mds <- cmdscale(D, k = 2, eig = TRUE, add = TRUE))

var_exp <- mds$eig / sum(pmax(mds$eig, 0))
var_exp <- round(100 * var_exp[1:2], 1)

class_map <- c(
  Cluster_1="IMOC", Cluster_2="IMOC",
  Squamous="Bailey", Progenitor="Bailey", ADEX="Bailey", Immunogenic="Bailey",
  Classical_Collisson="Collisson", Exocrine="Collisson", QM="Collisson",
  Classical_Moffitt="Moffitt", Basal_like="Moffitt"
)

class_map <- factor(class_map, levels = c("IMOC", "Bailey", "Collisson", "Moffitt"))

mds <- as_tibble(mds$points, .name_repair = "minimal") %>%
  setNames(c("Dim.1","Dim.2")) %>%
  mutate(Subtypes = rownames(mds$points),
         Classification = factor(class_map[Subtypes]))

FigXE <- mds %>% ggplot(aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(fill = Classification, colour = Classification, shape = Classification), size = 5, alpha = 0.5) +
    geom_label_repel(aes(label = Subtypes, fill = Classification), color = "white", show.legend = F, fontface = "bold", vjust = 1) +
    scale_shape_manual(values = c(21:24)) +
    scale_fill_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    scale_color_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    theme_cowplot(font_size = 18) +
    theme(legend.position = "top") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.2) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
    xlab(paste("MDS Dim.1 (", var_exp[1], "%)", sep = "")) +
    ylab(paste("MDS Dim.2 (", var_exp[2], "%)", sep = "")) +
    labs(title = "Centroid gene expression MDS (whole transcriptome)", 
         subtitle = paste("MDS computed on n=", nrow(rna_avg), " genes", sep = "")) 

## Union of DEGs

rna_avg_z <- t(scale(t(rna_avg_union)))
cor_mat_z <- cor(rna_avg_z, method = "spearman", use = "pairwise.complete.obs")

X <- 1-cor_mat_z
# D <- as.dist(X)
D <- as.dist(sqrt(2 * X)) # Euclidian distance

(mds <- cmdscale(D, k = 2, eig = TRUE, add = TRUE))

var_exp <- mds$eig / sum(pmax(mds$eig, 0))
var_exp <- round(100 * var_exp[1:2], 1)

class_map <- c(
  Cluster_1="IMOC", Cluster_2="IMOC",
  Squamous="Bailey", Progenitor="Bailey", ADEX="Bailey", Immunogenic="Bailey",
  Classical_Collisson="Collisson", Exocrine="Collisson", QM="Collisson",
  Classical_Moffitt="Moffitt", Basal_like="Moffitt"
)

class_map <- factor(class_map, levels = c("IMOC", "Bailey", "Collisson", "Moffitt"))

mds <- as_tibble(mds$points, .name_repair = "minimal") %>%
  setNames(c("Dim.1","Dim.2")) %>%
  mutate(Subtypes = rownames(mds$points),
         Classification = factor(class_map[Subtypes]))

FigYD <- mds %>% ggplot(aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(fill = Classification, colour = Classification, shape = Classification), size = 5, alpha = 0.5) +
    geom_label_repel(aes(label = Subtypes, fill = Classification), color = "white", show.legend = F, fontface = "bold", vjust = 1) +
    scale_shape_manual(values = c(21:24)) +
    scale_fill_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    scale_color_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    theme_cowplot(font_size = 18) +
    theme(legend.position = "top") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.2) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
    xlab(paste("MDS Dim.1 (", var_exp[1], "%)", sep = "")) +
    ylab(paste("MDS Dim.2 (", var_exp[2], "%)", sep = "")) +
    labs(title = "Centroid gene expression MDS (union of DEGs)", 
         subtitle = paste("MDS computed on n=", nrow(rna_avg_union), " genes (DEG union)", sep = "")) 


## 5.4 GSEA ----

# Prepare phenotype files for GSEA

phenotypes <- papers_cluster_df

phenotypes$collisson <- phenotypes$collisson %>% str_replace(pattern = "Classical", replacement = "Classical_Collisson")

clusterCol <- c("bailey", "collisson", "moffitt")

for (c in clusterCol) {
  lvls <- unique(phenotypes[[c]])
  for (l in lvls) {
    phenotypes[[l]] <- ifelse(phenotypes[[c]] == l, l, "Other") %>% factor(levels = c(l, "Other"))
  }
}

openxlsx2::write_xlsx(phenotypes, "data/data_omics/rna_seq/GSEA/NormCounts_phenotype_subtypes.xlsx")

## Use "data/data_omics/rna_seq/GSEA/NormCounts.txt" file to run GSEA in GSEA software

# Load results from GSEA

folders <- list.dirs("results/omics_analysis/RNA_Seq/subtypes_GSEA", recursive = FALSE, full.names = TRUE)

# Add the path to Clustering GSEA results

folders <- c(folders, "results/omics_analysis/RNA_Seq")

folders <- data.frame(dir = folders, subtype = c("ADEX", "Immunogenic", "Progenitor", "Squamous", "Classical_Collison", "Exocrine", "QM", "Basal_Like", "IMOC"))

gsea_reports <- list()

for (i in 1:nrow(folders)) {
  
  found_files <- list.files(
    path = folders$dir[i],
    pattern = paste0("^", "gsea_report", ".*.tsv"), # Matches the exact start (^) and end ($) of the filename
    full.names = TRUE, # Returns the full path
    recursive = FALSE # Set to TRUE if searching in subdirectories
  )
  
  enrich_cluster2 <- read_delim(found_files[1],
                                delim = "\t")
  enrich_cluster1 <- read_delim(found_files[2],
                                delim = "\t")
  
  enrich_gsea_c6 <- rbind(enrich_cluster2, enrich_cluster1)
  
  enrich_gsea_c6$`LEADING EDGE` <- enrich_gsea_c6$`LEADING EDGE` %>% str_remove(pattern = "%.*") %>% str_remove("tags=")
  enrich_gsea_c6$`LEADING EDGE` <- as.numeric(enrich_gsea_c6$`LEADING EDGE`)
  
  gsea_reports[[i]] <- enrich_gsea_c6
  names(gsea_reports)[i] <- folders$subtype[i]
}

# Keep only significant (FDR < 0.25) pathways

gsea_reports <- map(
  gsea_reports,
  ~ filter(.x, `FDR q-val` < 0.25)
)

df <- sapply(gsea_reports, nrow)
df <- data.frame(GSEA_Pathways = df, Subtype = factor(names(df)), 
                 Classification = factor(c(rep("Bailey", 4), rep("Collisson", 3), "Moffitt", "IMOC")))

FigYE <- df %>% mutate(Subtype = fct_reorder(Subtype, GSEA_Pathways, .desc = T)) %>% 
    ggplot(aes(x = Subtype, y = GSEA_Pathways)) +
    geom_bar(aes(color = Classification, fill = Classification), stat = "identity", alpha = 0.5) +
    scale_fill_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    scale_color_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    theme_cowplot(font_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
    ylab("Number of enriched oncogenic gene sets (FDR < 0.25)") +
    xlab("Subtype contrast") +
    labs(title = "Subtype GSEA enriched oncogenic pathways")


c6_fdr <- map(gsea_reports, "NAME") %>% unlist(use.names = FALSE) %>% unique()

## 5.5 GSVA ----

# Prepare gene-sample matrix

rna_mat <- read.csv("data/TCGA/omics_data/raw/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv",
                    row.names = 1,
                    check.names = F)

colnames(rna_mat) <-  colnames(rna_mat)  %>% str_replace(pattern = "-01A.*", replacement = "")
rna_mat <- as.matrix(rna_mat)

C6_msigdb <- jsonlite::fromJSON("data/data_omics/rna_seq/GSEA/c6.all.v2025.1.Hs.txt")

msig <- tibble()

for (i in 1:length(C6_msigdb)){
  df <- C6_msigdb[[i]]
  df <- tibble(
    geneSymbols = df$geneSymbols,
    msig = names(C6_msigdb)[i]
  )
  msig <- bind_rows(msig, df)
}

genesets <- split(msig$geneSymbols, msig$msig)

gsvaPar <- gsvaParam(rna_mat, geneSets = genesets, kcdf = "auto")

gsva.es <- gsva(gsvaPar)

### 5.5.1 GSVA Heatmaps ----

## Assign cluster IDs

gsva.es_df <- as.data.frame(t(gsva.es))

gsva.es_df$patient <- rownames(gsva.es_df)

gsva.es_df <- left_join(gsva.es_df, patient_cluster[,c(1,3)], by = "patient")

gsva.es_df$cluster <- as.factor(gsva.es_df$cluster)

gsva.es_df <- left_join(gsva.es_df, papers_cluster_df[,c("PatientID", "bailey", "collisson", "moffitt")], by = join_by(patient == PatientID))

gsva.es_df_long <- gsva.es_df %>% pivot_longer(cols = 1:189, names_to = "Pathway", values_to = "GSVA_Score")

# Mean expression accross subtypes/clusters

gsva_avg_clusters <- gsva.es_df_long %>% 
  group_by(cluster, Pathway) %>% 
  summarise(avg = mean(GSVA_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = avg)

gsva_avg_bailey <- gsva.es_df_long %>% filter(!is.na(bailey)) %>%
  group_by(bailey, Pathway) %>% 
  summarise(avg = mean(GSVA_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = bailey, values_from = avg)

gsva_avg_collisson <- gsva.es_df_long %>% filter(!is.na(collisson)) %>%
  group_by(collisson, Pathway) %>% 
  summarise(avg = mean(GSVA_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = collisson, values_from = avg)
gsva_avg_collisson <- rename(gsva_avg_collisson, Classical_Collisson = Classical)

gsva_avg_moffitt <- gsva.es_df_long %>% filter(!is.na(moffitt)) %>%
  group_by(moffitt, Pathway) %>% 
  summarise(avg = mean(GSVA_Score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = moffitt, values_from = avg)
gsva_avg_moffitt <- rename(gsva_avg_moffitt, Classical_Moffitt = Classical)

gsva_avg <- gsva_avg_clusters %>% 
  left_join(gsva_avg_bailey, by = "Pathway") %>%
  left_join(gsva_avg_collisson, by = "Pathway") %>%
  left_join(gsva_avg_moffitt, by = "Pathway")

gsva_avg <- column_to_rownames(gsva_avg, var = "Pathway")
gsva_avg <- as.matrix(gsva_avg)

annotation_col <- data.frame(
  Classification = factor(c(rep("IMOC", 2), rep("Bailey", 4), rep("Collisson", 3), rep("Moffitt", 2))),
  row.names = colnames(rna_avg)
)

annotation_colors <- list(
  Classification = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")
)
names(annotation_colors) <- colnames(annotation_col)

FigXF <- pheatmap::pheatmap(t(gsva_avg),
                             scale = "column",
                             #cluster_rows = F,
                             color =  colorRampPalette(c("navy",  "white", "firebrick2")) (1000),
                             annotation_row = annotation_col, border_color = "white", 
                             annotation_colors = annotation_colors, 
                             cutree_rows = 4, 
                             clustering_method = "average", 
                             show_colnames = F, fontsize_col = 8, angle_col = 45,
                             treeheight_col = 0,
                             main = "Centroid GSVA heatmap on oncogenic pathways",
                             fontsize = 16
)

FigYF <- pheatmap::pheatmap(t(gsva_avg[rownames(gsva_avg) %in% c6_fdr, ]),
                   scale = "column", border_color = "white", 
                   #cluster_rows = F,
                   color =  colorRampPalette(c("navy",  "white", "firebrick2")) (1000),
                   annotation_row = annotation_col, 
                   annotation_colors = annotation_colors, 
                   cutree_rows = 4, 
                   clustering_method = "average", 
                   show_colnames = T, fontsize_col = 8, angle_col = 45,
                   treeheight_col = 0,
                   main = "Centroid GSVA heatmap on enriched oncogenic pathways",
                   fontsize = 16
)

### 5.5.2 Correlation to Centroids ----


# Fig XG
gsva_avg_z <- t(scale(t(gsva_avg)))
cor_gsva_z <- cor(gsva_avg_z, method = "spearman", use = "pairwise.complete.obs")

col <- colorRampPalette(c("#4450AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

testRes = corrplot::cor.mtest(cor_gsva_z, conf.level = 0.95)

svg(file = "figures/omics_analysis/CorPlotAllPathways.svg", width = 10, height = 10)
p1 <- corrplot::corrplot(cor_gsva_z, method = "circle", type = "upper", col = col(200),
                         tl.srt = 45, tl.col = "black",
                         diag = FALSE, title = "Centroid oncogenic pathways GSVA correlation", 
                         p.mat = testRes$p, insig='blank',
                         number.cex = 0.8,
                         mar = c(1,1,4,1),
                         pch.cex = 3, tl.cex = 1.5, cl.cex = 1.5)$corrPos
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()

# Fig YG 

gsva_avg_z <- t(scale(t(gsva_avg[rownames(gsva_avg) %in% c6_fdr, ])))
cor_gsva_z <- cor(gsva_avg_z, method = "spearman", use = "pairwise.complete.obs")

col <- colorRampPalette(c("#4450AA", "#77AADD", "#FFFFFF", "#EE9988", "#BB4444"))

testRes = corrplot::cor.mtest(cor_gsva_z, conf.level = 0.95)


svg(file = "figures/omics_analysis/CorPlotSigPathways.svg", width = 10, height = 10)
p1 <- corrplot::corrplot(cor_gsva_z, method = "circle", type = "upper", col = col(200),
                         tl.srt = 45, tl.col = "black",
                         diag = FALSE, title = "Centroid enriched oncogenic pathways GSVA correlation", 
                         p.mat = testRes$p, insig='blank',
                         number.cex = 0.8,
                         mar = c(1,1,4,1),
                         pch.cex = 3, tl.cex = 1.5, cl.cex = 1.5)$corrPos
text(p1$x, p1$y, round(p1$corr, 2))
dev.off()

### 5.5.3 MDS ----

# Calculate distances from correlations

gsva_avg_z <- t(scale(t(gsva_avg)))
cor_gsva_z <- cor(gsva_avg_z, method = "spearman", use = "pairwise.complete.obs")

X <- 1-cor_gsva_z
# D <- as.dist(X)
D <- as.dist(sqrt(2 * X)) # Euclidian distance

(mds <- cmdscale(D, k = 2, eig = TRUE, add = TRUE))

var_exp <- mds$eig / sum(pmax(mds$eig, 0))
var_exp <- round(100 * var_exp[1:2], 1)

class_map <- c(
  Cluster_1="IMOC", Cluster_2="IMOC",
  Squamous="Bailey", Progenitor="Bailey", ADEX="Bailey", Immunogenic="Bailey",
  Classical_Collisson="Collisson", Exocrine="Collisson", QM="Collisson",
  Classical_Moffitt="Moffitt", Basal_like="Moffitt"
)

class_map <- factor(class_map, levels = c("IMOC", "Bailey", "Collisson", "Moffitt"))

mds <- as_tibble(mds$points, .name_repair = "minimal") %>%
  setNames(c("Dim.1","Dim.2")) %>%
  mutate(Subtypes = rownames(mds$points),
         Classification = factor(class_map[Subtypes]))


FigXH <- mds %>% ggplot(aes(x = Dim.1, y = Dim.2)) +
    geom_point(aes(fill = Classification, colour = Classification, shape = Classification), size = 5, alpha = 0.5) +
    geom_label_repel(aes(label = Subtypes, fill = Classification), color = "white", show.legend = F, fontface = "bold", vjust = 1) +
    scale_shape_manual(values = c(21:24)) +
    scale_fill_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    scale_color_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
    theme_cowplot(font_size = 18) +
    theme(legend.position = "top") +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.2) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
    xlab(paste("MDS Dim.1 (", var_exp[1], "%)", sep = "")) +
    ylab(paste("MDS Dim.2 (", var_exp[2], "%)", sep = "")) +
    labs(title = "Centroid GSVA MDS on oncogenic pathways", 
         subtitle = paste("MDS computed on n=", nrow(gsva_avg_z), " oncogenic gene sets", sep = "")) 



gsva_avg_z <- t(scale(t(gsva_avg[rownames(gsva_avg) %in% c6_fdr, ])))
cor_gsva_z <- cor(gsva_avg_z, method = "spearman", use = "pairwise.complete.obs")

X <- 1-cor_gsva_z
# D <- as.dist(X)
D <- as.dist(sqrt(2 * X)) # Euclidian distance

(mds <- cmdscale(D, k = 2, eig = TRUE, add = TRUE))

var_exp <- mds$eig / sum(pmax(mds$eig, 0))
var_exp <- round(100 * var_exp[1:2], 1)


class_map <- c(
  Cluster_1="IMOC", Cluster_2="IMOC",
  Squamous="Bailey", Progenitor="Bailey", ADEX="Bailey", Immunogenic="Bailey",
  Classical_Collisson="Collisson", Exocrine="Collisson", QM="Collisson",
  Classical_Moffitt="Moffitt", Basal_like="Moffitt"
)

class_map <- factor(class_map, levels = c("IMOC", "Bailey", "Collisson", "Moffitt"))

mds <- as_tibble(mds$points, .name_repair = "minimal") %>%
  setNames(c("Dim.1","Dim.2")) %>%
  mutate(Subtypes = rownames(mds$points),
         Classification = factor(class_map[Subtypes]))


FigYH <- mds %>% ggplot(aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(fill = Classification, colour = Classification, shape = Classification), size = 5, alpha = 0.5) +
  geom_label_repel(aes(label = Subtypes, fill = Classification), color = "white", show.legend = F, fontface = "bold", vjust = 1) +
  scale_shape_manual(values = c(21:24)) +
  scale_fill_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
  scale_color_manual(values = c("IMOC" = "#01665e", "Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
  theme_cowplot(font_size = 18) +
  theme(legend.position = "top") +
  geom_vline(xintercept = 0, lty = 2, lwd = 0.2) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
  xlab(paste("MDS Dim.1 (", var_exp[1], "%)", sep = "")) +
  ylab(paste("MDS Dim.2 (", var_exp[2], "%)", sep = "")) +
  labs(title = "Centroid GSVA MDS on enriched (FDR < 0.25) oncogenic pathways", 
       subtitle = paste("MDS computed on n=", nrow(gsva_avg_z), " enriched oncogenic gene sets", sep = ""))

# 6. Survival Comparison ----


survival_data <- clinical_data_file %>% select('Patient ID', 
                                               'Disease Free Status',
                                               'Overall Survival (Months)', 
                                               'Overall Survival Status', 
                                               'Disease Free (Months)')


survival_data <- survival_data %>% 
  mutate(
    `Overall Survival Status` = recode_values(`Overall Survival Status`, "0:LIVING" ~ 0, "1:DECEASED" ~ 1),
    `Disease Free Status` = recode_values(`Disease Free Status`, "0:DiseaseFree" ~ 0, "1:Recurred/Progressed" ~ 1)
  )

survival_data <- survival_data %>% left_join(papers_cluster_df %>% select(PatientID, IMOC, bailey, collisson, moffitt),
                                             by = c("Patient ID" = "PatientID")) %>% filter(!is.na(IMOC))

survival_data <- survival_data %>% rename(
  SurvivalTime = `Overall Survival (Months)`,
  SurvivalStatus = `Overall Survival Status`,
  DFStime = `Disease Free (Months)`,
  DFSstatus = `Disease Free Status`
) %>% mutate(
  IMOC = as.factor(IMOC),
  bailey = as.factor(bailey),
  collisson = as.factor(collisson),
  moffitt = as.factor(moffitt)
)


### Overal survival ----

IMOC_surv <- 
  survfit2(Surv(SurvivalTime, SurvivalStatus) ~ IMOC, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "IMOC",
    x = "Time (months)",
    y = "Survival probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  theme_ggsurvfit_KMunicate()  

bailey_surv <- 
  survfit2(Surv(SurvivalTime, SurvivalStatus) ~ bailey, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Bailey",
    x = "Time (months)",
    y = "Survival probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme_ggsurvfit_KMunicate()  


collisson_surv <- 
  survfit2(Surv(SurvivalTime, SurvivalStatus) ~ collisson, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Collisson",
    x = "Time (months)",
    y = "Survival probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_fill_brewer(type = "qual", palette = 6) +
  theme_ggsurvfit_KMunicate()  
collisson_surv

moffitt_surv <- 
  survfit2(Surv(SurvivalTime, SurvivalStatus) ~ moffitt, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Moffitt",
    x = "Time (months)",
    y = "Survival probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 7) +
  scale_fill_brewer(type = "qual", palette = 7) +
  theme_ggsurvfit_KMunicate()  
moffitt_surv


### Disease Free Progression ----

IMOC_dfs <- 
  survfit2(Surv(DFStime, DFSstatus) ~ IMOC, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "IMOC",
    x = "Time (months)",
    y = "Disease free probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  theme_ggsurvfit_KMunicate()  


bailey_dfs <- 
  survfit2(Surv(DFStime, DFSstatus) ~ bailey, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Bailey",
    x = "Time (months)",
    y = "Disease free probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 2) +
  scale_fill_brewer(type = "qual", palette = 2) +
  theme_ggsurvfit_KMunicate()  


collisson_dsf <- 
  survfit2(Surv(DFStime, DFSstatus) ~ collisson, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Collisson",
    x = "Time (months)",
    y = "Disease free probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_fill_brewer(type = "qual", palette = 6) +
  theme_ggsurvfit_KMunicate()  
collisson_dsf

moffitt_dsf <- 
  survfit2(Surv(DFStime, DFSstatus) ~ moffitt, data = survival_data) %>%
  ggsurvfit(theme = theme_cowplot()) +
  labs(
    title = "Moffitt",
    x = "Time (months)",
    y = "Desease free probability"
  ) + add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
  scale_ggsurvfit(x_scales = list(breaks = seq(0,78,6))) +
  add_censor_mark() +
  scale_color_brewer(type = "qual", palette = 7) +
  scale_fill_brewer(type = "qual", palette = 7) +
  theme_ggsurvfit_KMunicate()  
moffitt_dsf


### RMST ----

# Helper: identify worst- and best-prognosis subgroups for a scheme

identify_extreme_groups <- function(data, time_var, event_var, scheme,
                                    tau_for_ranking) {
  # Ranks subgroups by their KM-estimated survival probability at tau_for_ranking.
  # Returns the subgroup names with the lowest (worst) and highest (best) survival.
  
  fml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", scheme))
  km  <- survfit(fml, data = data)
  
  # Get survival at tau_for_ranking for each stratum
  sm <- summary(km, times = tau_for_ranking, extend = TRUE)
  surv_by_group <- data.frame(
    group = sub(paste0("^", scheme, "="), "", sm$strata),
    surv  = sm$surv,
    stringsAsFactors = FALSE
  )
  
  worst_group <- surv_by_group$group[which.min(surv_by_group$surv)]
  best_group  <- surv_by_group$group[which.max(surv_by_group$surv)]
  
  list(worst = worst_group, best = best_group, ranking = surv_by_group)
}

# Core function: compute RMST difference at a single tau
rmst_difference_at_tau <- function(data, time_var, event_var, scheme,
                                   tau, worst_group = NULL, best_group = NULL,
                                   ranking_tau = NULL) {
  # Subset to complete cases
  d <- data[complete.cases(data[, c(time_var, event_var, scheme)]), ]
  d[[scheme]] <- droplevels(factor(d[[scheme]]))
  
  # Auto-identify extreme groups if not provided
  if (is.null(worst_group) || is.null(best_group)) {
    rt <- if (is.null(ranking_tau)) tau else ranking_tau
    extremes <- identify_extreme_groups(d, time_var, event_var, scheme, rt)
    worst_group <- extremes$worst
    best_group  <- extremes$best
  }
  
  # Restrict to the two extreme subgroups
  d_pair <- d[d[[scheme]] %in% c(worst_group, best_group), ]
  d_pair[[scheme]] <- droplevels(factor(d_pair[[scheme]]))
  
  # survRM2 expects arm coded as 0/1, with arm=1 being the "treatment" group.
  # We code arm=1 = best (reference), arm=0 = worst, so the RMST difference
  # arm1 - arm0 = best - worst, which will be POSITIVE for a useful predictor.
  arm <- ifelse(d_pair[[scheme]] == best_group, 1, 0)
  
  # survRM2 requirement: tau must be <= min of the largest observed time per arm
  max_time_arm0 <- max(d_pair[[time_var]][arm == 0])
  max_time_arm1 <- max(d_pair[[time_var]][arm == 1])
  tau_max <- min(max_time_arm0, max_time_arm1)
  tau_used <- min(tau, tau_max)
  
  # Need at least one event per arm for stable estimation
  events_arm0 <- sum(d_pair[[event_var]][arm == 0])
  events_arm1 <- sum(d_pair[[event_var]][arm == 1])
  
  if (events_arm0 < 1 || events_arm1 < 1 ||
      length(unique(arm)) < 2 || tau_used <= 0) {
    return(tibble(
      scheme       = scheme,
      tau          = tau,
      tau_used     = NA_real_,
      worst_group  = worst_group,
      best_group   = best_group,
      n_worst      = sum(arm == 0),
      n_best       = sum(arm == 1),
      events_worst = events_arm0,
      events_best  = events_arm1,
      rmst_worst   = NA_real_,
      rmst_best    = NA_real_,
      rmst_diff    = NA_real_,
      ci_lower     = NA_real_,
      ci_upper     = NA_real_,
      p_value      = NA_real_,
      note         = "insufficient events or follow-up"
    ))
  }
  
  fit <- tryCatch(
    rmst2(time   = d_pair[[time_var]],
          status = d_pair[[event_var]],
          arm    = arm,
          tau    = tau_used),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      rmst2(time   = d_pair[[time_var]],
            status = d_pair[[event_var]],
            arm    = arm,
            tau    = tau_used)
    )
  )
  
  if (is.null(fit)) {
    return(tibble(
      scheme       = scheme, tau = tau, tau_used = tau_used,
      worst_group  = worst_group, best_group = best_group,
      n_worst      = sum(arm == 0), n_best = sum(arm == 1),
      events_worst = events_arm0,   events_best = events_arm1,
      rmst_worst   = NA_real_, rmst_best = NA_real_,
      rmst_diff    = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      p_value      = NA_real_,
      note         = "rmst2 failed"
    ))
  }
  
  # Extract per-arm RMST and the unadjusted difference (best - worst)
  rmst_worst <- fit$RMST.arm0$rmst["Est."]
  rmst_best  <- fit$RMST.arm1$rmst["Est."]
  diff_row   <- fit$unadjusted.result["RMST (arm=1)-(arm=0)", ]
  
  tibble(
    scheme       = scheme,
    tau          = tau,
    tau_used     = tau_used,
    worst_group  = worst_group,
    best_group   = best_group,
    n_worst      = sum(arm == 0),
    n_best       = sum(arm == 1),
    events_worst = events_arm0,
    events_best  = events_arm1,
    rmst_worst   = unname(rmst_worst),
    rmst_best    = unname(rmst_best),
    rmst_diff    = unname(diff_row["Est."]),
    ci_lower     = unname(diff_row["lower .95"]),
    ci_upper     = unname(diff_row["upper .95"]),
    p_value      = unname(diff_row["p"]),
    note         = if (tau_used < tau) "tau capped to follow-up" else NA_character_
  )
}

# Wrapper: run RMST analysis for one scheme across multiple tau values
rmst_across_taus <- function(data, time_var, event_var, scheme,
                             taus, ranking_tau = NULL,
                             worst_group = NULL, best_group = NULL) {
  # If ranking_tau not specified, pick a tau in the middle of the range for ranking
  # (so all tau values use the same worst/best assignment, ensuring consistency)
  if (is.null(ranking_tau) && (is.null(worst_group) || is.null(best_group))) {
    d <- data[complete.cases(data[, c(time_var, event_var, scheme)]), ]
    # Use a mid-range tau, capped by available follow-up
    ranking_tau <- min(24, quantile(d[[time_var]], 0.75, na.rm = TRUE))
    extremes <- identify_extreme_groups(d, time_var, event_var, scheme, ranking_tau)
    worst_group <- extremes$worst
    best_group  <- extremes$best
    cat(sprintf("  [%s] worst group = '%s', best group = '%s' (ranked at tau=%.1f)\n",
                scheme, worst_group, best_group, ranking_tau))
    # If a ranking_tau is specified without specifying the groups
  } else if (!is.null(ranking_tau) && (is.null(worst_group) || is.null(best_group))) {
    d <- data[complete.cases(data[, c(time_var, event_var, scheme)]), ]
    extremes <- identify_extreme_groups(d, time_var, event_var, scheme, ranking_tau)
    worst_group <- extremes$worst
    best_group  <- extremes$best
    cat(sprintf("  [%s] worst group = '%s', best group = '%s' (ranked at tau=%.1f)\n",
                scheme, worst_group, best_group, ranking_tau))
  }
  
  bind_rows(lapply(taus, function(tau) {
    rmst_difference_at_tau(
      data        = data,
      time_var    = time_var,
      event_var   = event_var,
      scheme      = scheme,
      tau         = tau,
      worst_group = worst_group,
      best_group  = best_group
    )
  }))
}

# Run the analysis

schemes <- c("IMOC", "bailey", "collisson", "moffitt")

# Identify worst groups at 12 months
for (s in schemes) {
  rank <- identify_extreme_groups(
    survival_data,
    time_var = "SurvivalTime",
    event_var = "SurvivalStatus",
    scheme = s,
    tau_for_ranking = 12
  )
  cat(sprintf("Scheme '%s':\nWorst group: '%s' / Best group: '%s'\n", s,rank$worst, rank$best))
}

for (s in schemes) {
  rank <- identify_extreme_groups(
    survival_data,
    time_var = "DFStime",
    event_var = "DFSstatus",
    scheme = s,
    tau_for_ranking = 12
  )
  cat(sprintf("Scheme '%s':\nWorst group: '%s' / Best group: '%s'\n", s,rank$worst, rank$best))
}


taus    <- c(6, 12, 18, 24, 36, 60)  # "Overall" added below

# --- Overall Survival ---
# "Overall" tau = min of largest observed time across all schemes' extreme pairs
tau_overall_os <- floor(min(survival_data$SurvivalTime[survival_data$SurvivalStatus == 1] %>%
                              max(),
                            max(survival_data$SurvivalTime)))
taus_os <- c(taus, tau_overall_os)
cat(sprintf("OS: tau values = %s\n", paste(taus_os, collapse = ", ")))


## IMOC
results_os <- 
  rmst_across_taus(
    data      = survival_data,
    time_var  = "SurvivalTime",
    event_var = "SurvivalStatus",
    scheme    = "IMOC",
    ranking_tau = 12,
    taus      = taus_os
  ) 

## Moffitt

results_os <- results_os %>% bind_rows(
  rmst_across_taus(
    data      = survival_data,
    time_var  = "SurvivalTime",
    event_var = "SurvivalStatus",
    scheme    = "moffitt",
    ranking_tau = 12,
    taus      = taus_os
  )
)

## Bailey

types <- c("Immunogenic", "Progenitor", "ADEX")

for (t in types) {
  results_os <- results_os %>% bind_rows(
    rmst_across_taus(
      data      = survival_data,
      time_var  = "SurvivalTime",
      event_var = "SurvivalStatus",
      scheme    = "bailey", 
      worst_group = "Squamous", 
      best_group = t,
      taus      = taus_os
    )
  )
}

## Collisson

types <- c("Classical", "Exocrine")

for (t in types) {
  results_os <- results_os %>% bind_rows(
    rmst_across_taus(
      data      = survival_data,
      time_var  = "SurvivalTime",
      event_var = "SurvivalStatus",
      scheme    = "collisson", 
      worst_group = "QM", 
      best_group = t,
      taus      = taus_os
    )
  )
}

results_os <- results_os %>% mutate(
  endpoint = "Overall survival", .before = 1
)

# --- Disease-Free Survival ---
tau_overall_dfs <- floor(min(survival_data$DFStime[survival_data$DFSstatus == 1] %>%
                               max(),
                             max(survival_data$DFStime, na.rm = TRUE), na.rm = TRUE))
taus_dfs <- c(taus, tau_overall_dfs)
cat(sprintf("\nDFS: tau values = %s\n", paste(taus_dfs, collapse = ", ")))


## IMOC
results_dfs <- 
  rmst_across_taus(
    data      = survival_data,
    time_var  = "DFStime",
    event_var = "DFSstatus",
    scheme    = "IMOC",
    ranking_tau = 12,
    taus      = taus_os
  ) 

## Moffitt

results_dfs <- results_dfs %>% bind_rows(
  rmst_across_taus(
    data      = survival_data,
    time_var  = "DFStime",
    event_var = "DFSstatus",
    scheme    = "moffitt",
    ranking_tau = 12,
    taus      = taus_os
  )
)

## Bailey

types <- c("Immunogenic", "Progenitor", "ADEX")

for (t in types) {
  results_dfs <- results_dfs %>% bind_rows(
    rmst_across_taus(
      data      = survival_data,
      time_var  = "DFStime",
      event_var = "DFSstatus",,
      scheme    = "bailey", 
      worst_group = "Squamous", 
      best_group = t,
      taus      = taus_os
    )
  )
}

## Collisson

types <- c("Classical", "Exocrine")

for (t in types) {
  results_dfs <- results_dfs %>% bind_rows(
    rmst_across_taus(
      data      = survival_data,
      time_var  = "DFStime",
      event_var = "DFSstatus",
      scheme    = "collisson", 
      worst_group = "QM", 
      best_group = t,
      taus      = taus_os
    )
  )
}

results_dfs <- results_dfs %>% mutate(
  endpoint = "Disease free progression", .before = 1
)



# --- Combine and apply BH correction (within endpoint × scheme family) ---
results_all <- bind_rows(results_os, results_dfs) %>%
  group_by(endpoint, scheme) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# --- View results ---
print(results_all, n = Inf)

results_all <- results_all %>% mutate(
  scheme = str_to_title(scheme) %>% str_replace("Imoc", "IMOC")
) %>% mutate(
  Comparison = paste0(scheme, " (", best_group, " vs. ", worst_group,")")
) %>% mutate(
  endpoint = factor(endpoint, levels = c("Overall survival", "Disease free progression")),
  scheme = factor(scheme, levels = c("IMOC", "Bailey", "Collisson", "Moffitt")),
  Comparison = factor(Comparison, 
                      levels = c("IMOC (Cluster 1 vs. Cluster 2)", 
                                 "Bailey (ADEX vs. Squamous)",
                                 "Bailey (Immunogenic vs. Squamous)",
                                 "Bailey (Progenitor vs. Squamous)",
                                 "Collisson (Classical vs. QM)",
                                 "Collisson (Exocrine vs. QM)",
                                 "Moffitt (Classical vs. Basal-like)"))
)


# Compute z-statistic from RMST difference and CI

# survRM2 returns CI as: estimate ± 1.96 × SE
# So SE = (upper - lower) / (2 × 1.96), and z = estimate / SE.
# This works for any analysis where the CI was constructed via normal
# approximation (which survRM2 uses).

add_z_statistic <- function(results) {
  results %>%
    mutate(
      se     = (ci_upper - ci_lower) / (2 * 1.96),
      z_stat = rmst_diff / se,
      # Two-sided p-value from z (sanity check vs survRM2's reported p)
      p_from_z = 2 * pnorm(-abs(z_stat))
    )
}

# Apply to combined Test 1 results

results_with_z <- results_all %>%    
  add_z_statistic()

# Quick sanity check: p_from_z should match p_value to within rounding
results_with_z %>%
  select(endpoint, scheme, tau, rmst_diff, se, z_stat, p_value, p_from_z) %>%
  mutate(p_diff = abs(p_value - p_from_z)) %>%
  arrange(desc(p_diff)) %>%
  head(10)
# If p_diff is consistently < 0.001, the z derivation is correct.
# If not, survRM2 may have returned a non-symmetric CI and you'll need
# to extract SE more carefully (see note at end).

# Summary table

z_summary <- results_with_z %>%
  select(endpoint, scheme, Comparison, tau,
         rmst_diff, se, z_stat, p_value, p_adj_bh) %>%
  mutate(across(c(rmst_diff, se, z_stat), \(x) round(x, 3))) %>%
  arrange(endpoint, tau, desc(z_stat))

print(z_summary, n = Inf)


### RMST difference within stratum ----


# Compute IMOC RMST difference within a single stratum

rmst_imoc_within_stratum <- function(data, time_var, event_var, tau,
                                     stratum_var, stratum_level,
                                     imoc_var = "IMOC",
                                     imoc_worst = "Cluster 2",
                                     imoc_best  = "Cluster 1") {
  # Subset to this stratum
  d <- data[data[[stratum_var]] == stratum_level, ]
  d <- d[complete.cases(d[, c(time_var, event_var, imoc_var)]), ]
  d[[imoc_var]] <- droplevels(factor(d[[imoc_var]]))
  
  # Code arm: 1 = best (Cluster 1), 0 = worst (Cluster 2) — matches Test 1 convention
  if (!all(c(imoc_worst, imoc_best) %in% levels(d[[imoc_var]]))) {
    return(tibble(
      stratum = stratum_level, tau = tau, tau_used = NA_real_,
      n_worst = sum(d[[imoc_var]] == imoc_worst, na.rm = TRUE),
      n_best  = sum(d[[imoc_var]] == imoc_best,  na.rm = TRUE),
      events_worst = NA_integer_, events_best = NA_integer_,
      rmst_diff = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
      note = "missing IMOC level in stratum"
    ))
  }
  
  arm <- ifelse(d[[imoc_var]] == imoc_best, 1, 0)
  
  # Tau capping
  max_t0 <- max(d[[time_var]][arm == 0])
  max_t1 <- max(d[[time_var]][arm == 1])
  tau_used <- min(tau, max_t0, max_t1)
  
  events_worst <- sum(d[[event_var]][arm == 0])
  events_best  <- sum(d[[event_var]][arm == 1])
  n_worst <- sum(arm == 0)
  n_best  <- sum(arm == 1)
  
  if (is.na(tau_used) || events_worst < 1 || events_best < 1 || n_worst < 2 || n_best < 2 ||
      tau_used <= 0) {
    return(tibble(
      stratum = stratum_level, tau = tau, tau_used = tau_used,
      n_worst = n_worst, n_best = n_best,
      events_worst = events_worst, events_best = events_best,
      rmst_diff = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
      note = "insufficient events or sample size"
    ))
  }
  
  fit <- tryCatch(
    suppressWarnings(rmst2(time = d[[time_var]],
                           status = d[[event_var]],
                           arm = arm, tau = tau_used)),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      stratum = stratum_level, tau = tau, tau_used = tau_used,
      n_worst = n_worst, n_best = n_best,
      events_worst = events_worst, events_best = events_best,
      rmst_diff = NA_real_, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
      note = "rmst2 failed"
    ))
  }
  
  diff_row <- fit$unadjusted.result["RMST (arm=1)-(arm=0)", ]
  est <- unname(diff_row["Est."])
  lo  <- unname(diff_row["lower .95"])
  hi  <- unname(diff_row["upper .95"])
  # Recover SE from the 95% CI (survRM2 uses normal-approx CI: est ± 1.96*SE)
  se  <- (hi - lo) / (2 * 1.96)
  
  tibble(
    stratum = stratum_level, tau = tau, tau_used = tau_used,
    n_worst = n_worst, n_best = n_best,
    events_worst = events_worst, events_best = events_best,
    rmst_diff = est, se = se,
    ci_lower = lo, ci_upper = hi,
    p_value = unname(diff_row["p"]),
    note = if (tau_used < tau) "tau capped to follow-up" else NA_character_
  )
}

# Pool within-stratum RMST differences via inverse-variance weighting

pool_rmst_strata <- function(stratum_results) {
  # Standard fixed-effects inverse-variance meta-analytic pooling.
  # Strata with NA RMST/SE are dropped before pooling.
  valid <- stratum_results %>%
    filter(!is.na(rmst_diff), !is.na(se), se > 0)
  
  if (nrow(valid) == 0) {
    return(tibble(
      pooled_rmst_diff = NA_real_, pooled_se = NA_real_,
      pooled_ci_lower = NA_real_, pooled_ci_upper = NA_real_,
      pooled_p_value = NA_real_, n_strata_pooled = 0L,
      n_strata_total = nrow(stratum_results)
    ))
  }
  
  w <- 1 / valid$se^2
  pooled_est <- sum(w * valid$rmst_diff) / sum(w)
  pooled_se  <- sqrt(1 / sum(w))
  z <- pooled_est / pooled_se
  p <- 2 * pnorm(-abs(z))
  
  # Cochran's Q test for heterogeneity
  Q_stat <- sum(w * (valid$rmst_diff - pooled_est)^2)
  Q_df   <- nrow(valid) - 1
  Q_p    <- if (Q_df > 0) pchisq(Q_stat, df = Q_df, lower.tail = FALSE) else NA_real_
  # I-squared
  I2     <- if (Q_df > 0) max(0, (Q_stat - Q_df) / Q_stat) * 100 else NA_real_
  
  tibble(
    pooled_rmst_diff = pooled_est,
    pooled_se        = pooled_se,
    pooled_ci_lower  = pooled_est - 1.96 * pooled_se,
    pooled_ci_upper  = pooled_est + 1.96 * pooled_se,
    pooled_p_value   = p,
    Q_stat           = Q_stat,
    Q_df             = Q_df,
    Q_p              = Q_p,
    I2               = I2,
    n_strata_pooled  = nrow(valid),
    n_strata_total   = nrow(stratum_results)
  )
}


# Run "added value" analysis for IMOC adjusted by one taxonomy at one tau

imoc_added_value <- function(data, time_var, event_var, tau,
                             taxonomy_var,
                             imoc_var = "IMOC",
                             imoc_worst = "Cluster 2",
                             imoc_best  = "Cluster 1") {
  d <- data[complete.cases(data[, c(time_var, event_var, imoc_var, taxonomy_var)]), ]
  d[[taxonomy_var]] <- droplevels(factor(d[[taxonomy_var]]))
  strata <- levels(d[[taxonomy_var]])
  
  per_stratum <- bind_rows(lapply(strata, function(s) {
    rmst_imoc_within_stratum(
      data = d, time_var = time_var, event_var = event_var,
      tau = tau, stratum_var = taxonomy_var, stratum_level = s,
      imoc_var = imoc_var, imoc_worst = imoc_worst, imoc_best = imoc_best
    )
  })) %>%
    mutate(taxonomy = taxonomy_var, .before = 1)
  
  pooled <- pool_rmst_strata(per_stratum) %>%
    mutate(taxonomy = taxonomy_var, tau = tau, .before = 1)
  
  list(per_stratum = per_stratum, pooled = pooled)
}


# Run across taus and taxonomies
run_test2 <- function(data, time_var, event_var, taus,
                      taxonomies = c("bailey", "collisson", "moffitt"),
                      imoc_var = "IMOC",
                      imoc_worst = "Cluster 2",
                      imoc_best  = "Cluster 1",
                      endpoint_label) {
  
  grid <- expand.grid(taxonomy = taxonomies, tau = taus,
                      stringsAsFactors = FALSE)
  
  per_stratum_all <- list()
  pooled_all      <- list()
  
  for (i in seq_len(nrow(grid))) {
    tx <- grid$taxonomy[i]
    tu <- grid$tau[i]
    res <- imoc_added_value(
      data = data, time_var = time_var, event_var = event_var,
      tau = tu, taxonomy_var = tx,
      imoc_var = imoc_var, imoc_worst = imoc_worst, imoc_best = imoc_best
    )
    per_stratum_all[[i]] <- res$per_stratum %>% mutate(tau = tu, .after = taxonomy)
    pooled_all[[i]]      <- res$pooled
  }
  
  list(
    per_stratum = bind_rows(per_stratum_all) %>%
      mutate(endpoint = endpoint_label, .before = 1),
    pooled      = bind_rows(pooled_all) %>%
      mutate(endpoint = endpoint_label, .before = 1) %>%
      group_by(endpoint, taxonomy) %>%
      mutate(p_adj_bh = p.adjust(pooled_p_value, method = "BH")) %>%
      ungroup()
  )
}


# Use the same tau grids as before
tau_overall_os <- floor(min(
  max(survival_data$SurvivalTime[survival_data$SurvivalStatus == 1]),
  max(survival_data$SurvivalTime)
))
taus_os <- c(6, 12, 18, 24, 36, 60, tau_overall_os)

tau_overall_dfs <- floor(min(
  max(survival_data$DFStime[survival_data$DFSstatus == 1], na.rm = TRUE),
  max(survival_data$DFStime, na.rm = TRUE),
  na.rm = TRUE
))
taus_dfs <- c(6, 12, 18, 24, 36, 60, tau_overall_dfs)

# --- Overall Survival ---
test2_os <- run_test2(
  data           = survival_data,
  taxonomies = "moffitt", 
  time_var       = "SurvivalTime",
  event_var      = "SurvivalStatus",
  taus           = taus_os,
  endpoint_label = "Overall survival"
)

# --- Disease-Free Survival ---
test2_dfs <- run_test2(
  data = survival_data,
  taxonomies = "moffitt",
  time_var = "DFStime",
  event_var = "DFSstatus",
  taus = taus_dfs,
  endpoint_label = "Disease-free survival"
)

# --- Combine ---
test2_per_stratum <- bind_rows(test2_os$per_stratum, test2_dfs$per_stratum)
test2_per_stratum <- test2_per_stratum %>% group_by(endpoint, stratum) %>%  
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()


# --- View ---

cat("\n=== PER-STRATUM RESULTS ===\n")
print(test2_per_stratum, n = Inf)


# 7. Figures ----

height <- 11.69
width <- 8.27

## Omics analysis

FigXD <- grImport2::readPicture("figures/omics_analysis/CorPlotAllGenes.svg")
FigXD <- grImport2::pictureGrob(FigXD)
FigXD <- ggdraw() + draw_grob(FigXD)


FigXG <- grImport2::readPicture("figures/omics_analysis/CorPlotAllPathways.svg")
FigXG <- grImport2::pictureGrob(FigXG)
FigXG <- ggdraw() + draw_grob(FigXG)

AMI_ARI_clusters <- tibble(
  AMI = c(0.17, 0.11, -0.006),
  ARI = c(0.08, 0.03, -0.001),
  Classification = c("Bailey", "Collisson", "Moffitt")
) %>% pivot_longer(cols = c("AMI", "ARI"), names_to = "Metric")

FigXB <- AMI_ARI_clusters %>% ggplot(aes(x = Classification, y = value)) +
  geom_segment(aes(x = Classification, y = 0, yend = value, colour = Classification)) +
  geom_hline(yintercept = 0, lty = 2, color = "gray") +
  geom_hline(yintercept = 1, lty = 2, color = "gray") +
  geom_point(aes(x = Classification, y = value, colour = Classification, fill = Classification), size = 5) +
  ylim(c(-0.01, 1)) + 
  ylab("Metric value") +
  scale_color_manual(values = c("Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
  scale_fill_manual(values = c("Bailey" = "#5ab4ac", "Collisson" = "#d8b365", "Moffitt" = "#8c510a")) +
  theme_cowplot() +
  theme(legend.position = "none") +
  facet_wrap(~Metric)

FigXAB <- plot_grid(FigXA, FigXB, nrow = 1, rel_widths = c(0.66, 0.33), labels = "AUTO")
FigXCH <- plot_grid(plotlist = list(FigXC[[4]], FigXD, FigXE, FigXF[[4]], FigXG, FigXH), byrow = F, ncol = 2, labels = c("C", "D", "E", "F", "G", "H"))
FigX <- plot_grid(FigXAB, FigXCH, rel_heights = c(0.2,0.8), ncol = 1)
ggsave(filename = "figures/omics_analysis/FigSubtypes.pdf",plot = FigX, height = 2.5*height, width = 2.5*width)


FigYC <- grImport2::readPicture("figures/omics_analysis/CorPlotUDEGs.svg")
FigYC <- grImport2::pictureGrob(FigYC)
FigYC <- ggdraw() + draw_grob(FigYC)

FigYG <- grImport2::readPicture("figures/omics_analysis/CorPlotSigPathways.svg")
FigYG <- grImport2::pictureGrob(FigYG)
FigYG <- ggdraw() + draw_grob(FigYG)

FigY <- plot_grid(plotlist = list(FigYA, FigYB[[4]], FigYC, FigYD, FigYE, FigYF[[4]], FigYG, FigYH), ncol = 2, byrow = F, labels = c("A","B", "C", "D", "E", "F", "G","H"))

ggsave(plot = FigY, filename = "figures/omics_analysis/FigSubtypesSig.pdf", height = 2.6*height, width = 2.6*width)

## Survival analysis


plot_grid(
  ggsurvfit_build(IMOC_surv),
  ggsurvfit_build(bailey_surv),
  ggsurvfit_build(collisson_surv),
  ggsurvfit_build(moffitt_surv),
  ncol = 2
)
ggsave("figures/OS_curves_classifications.pdf", height = heigth, width = width*2.5)


plot_grid(
  ggsurvfit_build(IMOC_dfs),
  ggsurvfit_build(bailey_dfs),
  ggsurvfit_build(collisson_dsf),
  ggsurvfit_build(moffitt_dsf),
  ncol = 2
)

ggsave("figures/DFS_curves_classifications.pdf", height = heigth, width = width*2.5)


(RMST_Diff <- 
    results_all %>% filter(tau_used <= 36) %>% 
    ggplot(aes(x = tau_used, y = rmst_diff, group = Comparison)) +
    geom_segment(
      aes(y = ci_lower, yend = ci_upper, x = tau_used, color = Comparison),
      arrow = arrow(
        angle = 90,
        length = unit(0.01, "npc"),
        ends = "both"
      ), position = position_dodge(width = 2), alpha = 0.5)  +
    geom_point(aes(shape = scheme, color = Comparison), position = position_dodge(width = 2), size = 3) +
    geom_text(aes(y = ci_upper, label = paste0("p=",round(p_adj_bh,3)), 
                  color = Comparison), show.legend = FALSE, size = 3,
              position = position_dodge(width = 2), angle = 45, hjust = 0, vjust = -0.5) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.25) +
    scale_x_continuous(breaks = round(unique(results_all$tau_used), 0), name = "Time (months)") +
    scale_color_brewer(type = "qual", palette = 2) +
    scale_shape_manual(name = "Classification", values = c(19,15, 18, 17)) + 
    theme_cowplot(font_size = 20) +
    labs(
      title = "RMST differences between worst- and best-prognosis subgroups",
      x = "Time (months)",
      y = expression(Delta~"RMST (\u00B1 95% CI)")
    ) +  facet_wrap( ~ endpoint, ncol = 1) #+ theme(legend.position = "top")
)



(RMST_zval <- 
  ggplot(results_with_z %>% filter(tau <= 36),  # match 36-month window
         aes(x = tau, y = z_stat,
             color = Comparison,
             shape = scheme,
             group = Comparison)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", color = "grey70") +
  geom_line(alpha = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ endpoint, ncol = 1) +
  scale_x_continuous(breaks = c(6, 12, 18, 24, 36)) +
  labs(
    title = "Signal-to-noise ratio by comparison",
    x = "Time (months)",
    y = expression("Z-statistic"~(hat(Delta)/widehat(SE))),
    color = "Comparison",
    shape = "Classification",
    caption = "Dotted lines at z = \u00B1 1.96 mark the threshold for nominal significance (p ~ 0.05).",
  ) +
  scale_shape_manual(values = c(19,15, 18, 17)) +
  scale_color_brewer(type = "qual", palette = 2) +
  theme_cowplot(font_size = 20) #+ theme(legend.position = "top")
)
plot_grid(RMST_Diff, RMST_zval, ncol = 1, labels = "AUTO")
ggsave("figures/RMST_Classification.pdf", height = heigth*2, width = width*2.5)


(rmst_diff_strata <- 
  test2_per_stratum %>% mutate(
    endpoint = factor(endpoint, levels = c("Overall survival", "Disease-free survival"))
  ) %>%
  ggplot(aes(x = tau_used, y = rmst_diff, group = stratum)) +
  geom_segment(
    aes(y = ci_lower, yend = ci_upper, x = tau_used),
    arrow = arrow(
      angle = 90,
      length = unit(0.01, "npc"),
      ends = "both"
    ), color = "gray50", position = position_dodge(width = 1.5))  +
  geom_point(aes(shape = stratum, fill = p_adj_bh),color = "gray50", size = 5,
             position = position_dodge(width = 1.5)) +
  geom_text(aes(x = tau_used, y = ci_upper, label = paste0("p=",round(p_adj_bh, 3)), hjust = 0, vjust = -0.5),
            position = position_dodge(width = 1.5), angle = 45, color = "gray50", size = 5) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.25) +
  scale_x_continuous(breaks = round(unique(test2_per_stratum$tau_used), 0), name = "Time (months)") +
  scale_shape_manual(name = "Moffitt subtype", values = c(22,23)) +
  scale_fill_gradient(name = "Adj. p-value", high = "#fee5d9", low = "firebrick") +
  theme_cowplot(font_size = 20) +
  labs(
    title = "Moffitt IMOC RMST difference within each subtype (Cluster 2 vs. Cluster 1)",
    x = "Time (months)",
    y = expression(Delta~"RMST (\u00B1 95% CI)")
  ) +  facet_wrap( ~ endpoint, ncol = 1)
)

ggsave("figures/RMST_Diff_Moffitt.pdf",plot = rmst_diff_strata, height = heigth, width = width*2)


# 8. Model to test uniqueness ----

meta <- gsva.es_df_long %>% select(-c(Pathway, GSVA_Score)) %>% group_by(patient, cluster, bailey, collisson, moffitt) %>% unique() %>% ungroup()

meta$cluster <- factor(meta$cluster)   
meta$bailey <- factor(meta$bailey)
meta$collisson <- factor(meta$collisson)
meta$moffitt <- factor(meta$moffitt)

limma_gsva <- function(meta, scores, new_class, old_class, ref_level, target_level, adjust.method = "BH") {
  
  meta <- meta %>% filter(!is.na(.data[[old_class]])) # Some samples do not have a classic subtype assigned
  
  # align score columns to meta order
  stopifnot(all(meta$patient %in% colnames(scores)))
  scores <- scores[, meta$patient, drop = FALSE]
  
  meta[[new_class]] <- relevel(factor(meta[[new_class]]), ref = ref_level)
  
  meta[[old_class]] <- factor(meta[[old_class]])
  
  
 
  design <- model.matrix(reformulate(termlabels = c(new_class, old_class), response = NULL), data = meta)
  fit <- lmFit(scores, design)
  fit <- eBayes(fit)
  
  coef_name <- paste0(new_class, target_level)
  stopifnot(coef_name %in% colnames(design))
  
  res <- topTable(fit, coef = coef_name, number = Inf, adjust.method = adjust.method, confint = TRUE) # Use Cluster_2 as coefficient of interest
  return(res)
}


limma_gsva_bailey <- limma_gsva(meta = meta, scores = gsva.es, new_class = "cluster", old_class = "bailey", ref_level = "Cluster_1", target_level = "Cluster_2")
limma_gsva_bailey$Pathway <- rownames(limma_gsva_bailey)
limma_gsva_bailey$Classification <- "Bailey"

limma_gsva_collisson <- limma_gsva(meta = meta, scores = gsva.es, new_class = "cluster", old_class = "collisson", ref_level = "Cluster_1", target_level = "Cluster_2")
limma_gsva_collisson$Pathway <- rownames(limma_gsva_collisson)
limma_gsva_collisson$Classification <- "Collisson"

limma_gsva_moffitt <- limma_gsva(meta = meta, scores = gsva.es, new_class = "cluster", old_class = "moffitt", ref_level = "Cluster_1", target_level = "Cluster_2")
limma_gsva_moffitt$Pathway <- rownames(limma_gsva_moffitt)
limma_gsva_moffitt$Classification <- "Moffitt"


limma_gsva_df <- bind_rows(limma_gsva_bailey, limma_gsva_collisson, limma_gsva_moffitt)
rownames(limma_gsva_df) <- NULL

top10 <- limma_gsva_df %>% group_by(Classification) %>% slice_min(order_by = adj.P.Val, n = 10) %>% ungroup()

top10 %>% mutate(Pathway = fct_reorder(Pathway, logFC, .desc = F)) %>% ggplot(aes(x = logFC, y = Pathway)) + 
  geom_segment(aes(x = CI.L, xend = CI.R, y = Pathway), arrow = arrow(angle = 90, length = unit(0.02, "npc"),ends = "both"), color = "gray15") +
  geom_point(aes(fill = adj.P.Val), size = 5, shape = 21) +
  scale_fill_gradient(low = "firebrick1", high = "white", name = "FDR") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "gray15") +
  facet_wrap(~Classification) +
  theme_bw(base_size = 18) +
  labs(title = "Cluster-associated oncogenic pathway activity after adjustment for classic PDAC subtype labels",
       subtitle = "limma-estimated difference in GSVA scores (Cluster 2 − Cluster 1), adjusted for Bailey/Collisson/Moffitt taxonomy") +
  xlab("Adjusted difference in GSVA score (Cluster 2 − Cluster 1)") +
  ylab("Oncogenic gene set (MSigDB C6)") +
  theme(title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 16))
ggsave("figures/omics_analysis/limma_pathways.pdf", width = 2.5*width, height = height)


# 9. Data S3 ----

DataS3 <- list(GSVA_Score = gsva.es_df)

wb <- wb_workbook()

for(i in names(DataS3)) {
  wb$add_worksheet(i)
  wb$add_data(sheet = i, x = DataS3[[i]])
}

wb$save("results/DataS3.xlsx")

# 10. Data S4 ---- 

ADEX <- bailey_DEGs$ADEXvsOther
Squamous <- bailey_DEGs$SquamousvsOther 
Progenitor <- bailey_DEGs$ProgenitorvsOther
Immunogenic <- bailey_DEGs$ImmunogenicvsOther


Classical_collison <- collisson_DEGs$ClassicalvsOther
QM <- collisson_DEGs$QMvsOther 
Exocrine <- collisson_DEGs$ExocrinevsOther

Basal_likevsClassical <- moffitt_DEGs$Basal_likevsClassical 

DataS4 <- list(ADEX_Bailey = ADEX, Squamous_Bailey = Squamous, Progenitor_Bailey = Progenitor, 
               Immunogenic_Bailey = Immunogenic, Classical_Collisson = Classical_collison, QM_Collisson = QM, 
               Exocrine_Collisson = Exocrine, BasalvsClassical_Moffitt = Basal_likevsClassical)

wb <- wb_workbook()

for(i in names(DataS4)) {
  wb$add_worksheet(i)
  wb$add_data(sheet = i, x = DataS4[[i]])
}


wb$save("results/DataS4.xlsx")

# 11. Data S5 ---- 

folders <- list.dirs("results/omics_analysis/RNA_Seq/subtypes_GSEA", recursive = FALSE, full.names = TRUE)

# Add the path to Clustering GSEA results

folders <- c(folders, "results/omics_analysis/RNA_Seq")

folders <- data.frame(dir = folders, subtype = c("ADEX", "Immunogenic", "Progenitor", "Squamous", "Classical_Collison", "Exocrine", "QM", "Basal_Like", "IMOC"))

gsea_reports <- list()

for (i in 1:nrow(folders)) {
  
  found_files <- list.files(
    path = folders$dir[i],
    pattern = paste0("^", "gsea_report", ".*.tsv"), # Matches the exact start (^) and end ($) of the filename
    full.names = TRUE, # Returns the full path
    recursive = FALSE # Set to TRUE if searching in subdirectories
  )
  
  enrich_cluster2 <- read_delim(found_files[1],
                                delim = "\t")
  enrich_cluster1 <- read_delim(found_files[2],
                                delim = "\t")
  
  enrich_gsea_c6 <- rbind(enrich_cluster2, enrich_cluster1)
  
  enrich_gsea_c6$`LEADING EDGE` <- enrich_gsea_c6$`LEADING EDGE` %>% str_remove(pattern = "%.*") %>% str_remove("tags=")
  enrich_gsea_c6$`LEADING EDGE` <- as.numeric(enrich_gsea_c6$`LEADING EDGE`)
  
  gsea_reports[[i]] <- enrich_gsea_c6
  names(gsea_reports)[i] <- folders$subtype[i]
}


ADEX <- gsea_reports$ADEX
Squamous <- gsea_reports$Squamous
Progenitor <- gsea_reports$Progenitor
Immunogenic <- gsea_reports$Immunogenic


Classical_collisson <- gsea_reports$Classical_Collison
QM <- gsea_reports$QM
Exocrine <- gsea_reports$Exocrine

Basal_likevsClassical <- gsea_reports$Basal_Like


DataS5 <- list(ADEX_Bailey = ADEX, Squamous_Bailey = Squamous, Progenitor_Bailey = Progenitor, 
               Immunogenic_Bailey = Immunogenic, Classical_Collisson = Classical_collisson, QM_Collisson = QM, 
               Exocrine_Collisson = Exocrine, BasalvsClassical_Moffitt = Basal_likevsClassical)

wb <- wb_workbook()

for(i in names(DataS5)) {
  wb$add_worksheet(i)
  wb$add_data(sheet = i, x = DataS5[[i]])
}


wb$save("results/DataS5.xlsx")

# 12. Data S6 ----

limma_gsva_bailey

DataS6 <- list(limma_GSVA_Bailey = limma_gsva_bailey,
               limma_GSVA_Collisson = limma_gsva_collisson,
               limma_GSVA_MOffitt = limma_gsva_moffitt)

wb <- wb_workbook()

for(i in names(DataS6)) {
  wb$add_worksheet(i)
  wb$add_data(sheet = i, x = DataS6[[i]])
}

wb$save("results/DataS6.xlsx")


# 13. Data S7 ----

write_xlsx(results_with_z, "results/DataS7.xlsx")

# 14. Data S8 ----

write_xlsx(test2_per_stratum, "results/DataS8.xlsx")


