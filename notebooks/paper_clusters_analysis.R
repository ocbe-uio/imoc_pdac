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


# 6. Figures ----

height <- 11.69
width <- 8.27

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


# 7. Model to test uniqueness ----

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

# 8. Data S3 ----

DataS3 <- list(GSVA_Score = gsva.es_df)

wb <- wb_workbook()

for(i in names(DataS3)) {
  wb$add_worksheet(i)
  wb$add_data(sheet = i, x = DataS3[[i]])
}

wb$save("results/DataS3.xlsx")

# 9. Data S4 ---- 

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

# 10. Data S5 ---- 

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

# 11. Data S6 ----

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
