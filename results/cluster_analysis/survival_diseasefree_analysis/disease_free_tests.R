
devtools::install_github("Shamir-Lab/Logrank-Inaccuracies/logrankHeinze")
install.packages("survRM2")
library("survRM2")
library("logrankHeinze")
library("dplyr")
library("survival")
library("tidyverse")
library("rms")
library("cowplot")

wd <- ".../imo_clustering/results/cluster_analysis"    # edit to set working directory
setwd(wd)

clusters_diseasefree <- read.csv('results/cluster_analysis/survival_diseasefree_analysis/disease_free_final_clusters.csv')
pdac_diseasefree <- clusters_diseasefree[, c("Patient.ID", "Disease.Free.Status",
                                             "Disease.Free..Months.", "Cluster")]
row.names(pdac_diseasefree) <- pdac_diseasefree$Patient.ID
pdac_diseasefree <- pdac_diseasefree[, -1]
colnames(pdac_diseasefree) <- c("DiseaseFreeStatus", "DiseaseFreeMonths", "Cluster")
pdac_diseasefree$Cluster <- as.numeric(as.factor(pdac_diseasefree$Cluster))



### Log-rank tests, RMST test at different timepoints and checking if PH model can be used

# Complete survival timeline
result <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                        data=pdac_diseasefree, max.num.perms = 1000, seed=42)
print(result$pvalue)

diseasefree_rmst <- rmst2(pdac_diseasefree$DiseaseFreeMonths, pdac_diseasefree$DiseaseFreeStatus, 
                       as.numeric(pdac_diseasefree$Cluster == "1"), tau=max(pdac_diseasefree$DiseaseFreeMonths))
print(diseasefree_rmst)

fit.unstrat <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree)
cox.zph(fit.unstrat)

# First six months
pdac_diseasefree_6m <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 6, 0, DiseaseFreeStatus))
results_6m <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_6m, max.num.perms = 1000, seed=42)
print(results_6m$pvalue)

diseasefree_rmst_6m <- rmst2(pdac_diseasefree_6m$DiseaseFreeMonths, pdac_diseasefree_6m$DiseaseFreeStatus, 
                       as.numeric(pdac_diseasefree_6m$Cluster == "1"), tau=6)
print(diseasefree_rmst_6m)

fit.unstrat_6m <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_6m)
cox.zph(fit.unstrat_6m)


# First year
pdac_diseasefree_y1 <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 12, 0, DiseaseFreeStatus))
results_y1 <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_y1, max.num.perms = 1000, seed=42)
print(results_y1$pvalue)

diseasefree_rmst_y1 <- rmst2(pdac_diseasefree_y1$DiseaseFreeMonths, pdac_diseasefree_y1$DiseaseFreeStatus, 
                             as.numeric(pdac_diseasefree_y1$Cluster == "1"), tau=12)
print(diseasefree_rmst_y1)

fit.unstrat_y1 <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_y1)
cox.zph(fit.unstrat_y1)


# First 18 months (1.5 years)
pdac_diseasefree_18m <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 18, 0, DiseaseFreeStatus))
results_18m <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_18m, max.num.perms = 1000, seed=42)
print(results_18m$pvalue)

diseasefree_rmst_18m <- rmst2(pdac_diseasefree_18m$DiseaseFreeMonths, pdac_diseasefree_18m$DiseaseFreeStatus, 
                             as.numeric(pdac_diseasefree_18m$Cluster == "1"), tau=18)
print(diseasefree_rmst_18m)

fit.unstrat_18m <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_18m)
cox.zph(fit.unstrat_18m)

# First 2 years
pdac_diseasefree_y2 <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 24, 0, DiseaseFreeStatus))
results_y2 <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_y2, max.num.perms = 1000, seed=42)
print(results_y2$pvalue)

diseasefree_rmst_y2 <- rmst2(pdac_diseasefree_y2$DiseaseFreeMonths, pdac_diseasefree_y2$DiseaseFreeStatus, 
                             as.numeric(pdac_diseasefree_y2$Cluster == "1"), tau=24)
print(diseasefree_rmst_y2)

fit.unstrat_y2 <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_y2)
cox.zph(fit.unstrat_y2)

# First 3 years
pdac_diseasefree_y3 <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 36, 0, DiseaseFreeStatus))
results_y3 <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_y3, max.num.perms = 1000, seed=42)
print(results_y3)

diseasefree_rmst_y3 <- rmst2(pdac_diseasefree_y3$DiseaseFreeMonths, pdac_diseasefree_y3$DiseaseFreeStatus, 
                             as.numeric(pdac_diseasefree_y3$Cluster == "1"), tau=36)
print(diseasefree_rmst_y3)

fit.unstrat_y3 <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_y3)
cox.zph(fit.unstrat_y3)

# First 5 years
pdac_diseasefree_y5 <- pdac_diseasefree %>%
  mutate(DiseaseFreeStatus = ifelse(DiseaseFreeMonths > 60, 0, DiseaseFreeStatus))
results_y5 <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ Cluster, 
                                            data=pdac_diseasefree_y5, max.num.perms = 1000, seed=42)
print(results_y5)

diseasefree_rmst_y5 <- rmst2(pdac_diseasefree_y5$DiseaseFreeMonths, pdac_diseasefree_y5$DiseaseFreeStatus, 
                             as.numeric(pdac_diseasefree_y5$Cluster == "1"), tau=60)
print(diseasefree_rmst_y5)

fit.unstrat_y5 <- coxph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ Cluster, data=pdac_diseasefree_y5)
cox.zph(fit.unstrat_y5)



## DFP-CNA Association ----

cna_tcga <- read_csv("data/TCGA/omics_data/raw/cancer_data_PAAD_GISTIC_Peaks-20160128.csv")

cna_tcga_ID <- cna_tcga[,c(3,13:ncol(cna_tcga))]

cna_tcga_ID <- as.data.frame(t(cna_tcga_ID))

colnames(cna_tcga_ID) <- cna_tcga_ID[1,]
cna_tcga_ID <- cna_tcga_ID[-1,]
cna_tcga_ID <- cna_tcga_ID %>% select(`21q11.2`, `17p12`, `18q21.2`, `9p21.3`)

# Retrieve patients names

patient <- rownames(cna_tcga_ID)
patient <- patient %>% str_replace(pattern = "-01A.*", replacement = "")

cna_tcga_ID <- apply(cna_tcga_ID, 2, as.numeric)

cna_tcga_ID <- as.data.frame(cna_tcga_ID)

rownames(cna_tcga_ID) <- patient

## Align survival data with CNA data

common_patients <- intersect(rownames(pdac_diseasefree), rownames(cna_tcga_ID))

cna_tcga_ID <- cna_tcga_ID[common_patients, ,drop = FALSE]
pdac_diseasefree <- pdac_diseasefree[common_patients, ,drop = FALSE]

df_results <- data.frame()
for (k in colnames(cna_tcga_ID)) {
  cna_band <- cna_tcga_ID[[k]]
  row <- list()
  row$Band <- k
  group <- as.numeric(cna_band == 0)+1
  if (length(unique(group)) == 1) {
    next
  }
  pdac_diseasefree$group <- group
  pdac_diseasefree$group_factor <- as.factor(group)
  invisible(capture.output(logrank <- logrankHeinze::logrank.heinze(DiseaseFreeMonths + DiseaseFreeStatus ~ group, data=pdac_diseasefree,
                                                                    max.num.perms = 1000, seed=42)$pvalue))
  row$logrank <- logrank
  survival_rmst <- rmst2(pdac_diseasefree[["DiseaseFreeMonths"]], pdac_diseasefree[["DiseaseFreeStatus"]], group-1)
  row$rmst_coef <- survival_rmst$unadjusted.result[1,1]
  row$rmst_pval <- survival_rmst$unadjusted.result[1,4]
  cox_model <- coxph(Surv(pdac_diseasefree[["DiseaseFreeMonths"]], pdac_diseasefree[["DiseaseFreeStatus"]]) ~ as.factor(group))
  row$cox_coef <- summary(cox_model)$coefficients[,"exp(coef)"]
  row$cox_pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
  dd <- datadist(pdac_diseasefree)
  options(datadist = 'dd')
  fit <- cph(Surv(DiseaseFreeMonths, DiseaseFreeStatus) ~ group_factor, data=pdac_diseasefree, x=TRUE, y=TRUE)
  set.seed(42)
  invisible(capture.output(v <- validate(fit, method="boot", B=1000)))
  c_index_results <- (v["Dxy", ] / 2) + 0.5
  row$c_index <- c_index_results["index.corrected"]
  df_results <- rbind(df_results, row)
}
df_results$rmst_pval_adj <- p.adjust(df_results$rmst_pval, method = "BH")
ggplot(df_results, aes(x = rmst_coef, y=reorder(Band, rmst_coef))) +
  geom_segment(aes(xend = 0, yend = Band), color = "gray80") +
  geom_point(size=6) +
  labs(x = "RMST (months)", y = "Deletion") +
  geom_text(aes(label = paste0("p =", round(rmst_pval_adj, 5))), vjust = -0.9) +
  xlim(0, 15) +
  theme_cowplot()


