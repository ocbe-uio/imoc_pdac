
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



clusters_survival <- read.csv('results/cluster_analysis/survival_diseasefree_analysis/survival_final_clusters.csv')
pdac_survival <- clusters_survival[, c("Patient.ID", "Overall.Survival.Status", 
                                       "Overall.Survival..Months.", "Cluster", "Weights")]
row.names(pdac_survival) <- pdac_survival$Patient.ID
pdac_survival <- pdac_survival[, -1]
colnames(pdac_survival) <- c("SurvivalStatus", "SurvivalMonths", "Cluster", "Weights")
pdac_survival$Cluster <- as.numeric(as.factor(pdac_survival$Cluster))


### Log-rank tests, RMST test at different timepoints and checking if PH model can be used

# Complete survival timeline
result <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival, 
                                        max.num.perms = 1000, seed=42)
print(result)

survival_rmst <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                       as.numeric(pdac_survival$Cluster == "1"), tau=max(pdac_survival$SurvivalMonths))
print(survival_rmst)

fit.unstrat <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival, weights=Weights)
cox.zph(fit.unstrat)


# First six months
pdac_survival_6m <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 6, 0, SurvivalStatus))
results_6m <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_6m, 
                                            max.num.perms = 1000, seed=42)
print(results_6m)

survival_rmst_6m <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=6)
print(survival_rmst_6m)

fit.unstrat_6m <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_6m, weights=Weights)
cox.zph(fit.unstrat_6m)


# First year
pdac_survival_y1 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 12, 0, SurvivalStatus))
results_y1 <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_y1, 
                                            max.num.perms = 1000, seed=42)
print(results_y1)

survival_rmst_y1 <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=12)
print(survival_rmst_y1)

fit.unstrat_y1 <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_y1, weights=Weights)
cox.zph(fit.unstrat_y1)


# First 18 months (1.5 years)
pdac_survival_18m <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 18, 0, SurvivalStatus))
results_18m <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_18m, 
                                            max.num.perms = 1000, seed=42)
print(results_18m)

survival_rmst_18m <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                           as.numeric(pdac_survival$Cluster == "1"), tau=18)
print(survival_rmst_18m)

fit.unstrat_18m <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_18m)
cox.zph(fit.unstrat_18m)

# First 2 years
pdac_survival_y2 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 24, 0, SurvivalStatus))
results_y2 <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_y2, 
                                            max.num.perms = 1000, seed=42)
print(results_y2)

survival_rmst_y2 <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=24)
print(survival_rmst_y2)

fit.unstrat_y2 <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_y2, weights=Weights)
cox.zph(fit.unstrat_y2)

# First 3 years
pdac_survival_y3 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 36, 0, SurvivalStatus))
results_y3 <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_y3, 
                                            max.num.perms = 1000, seed=42)
print(results_y3)

survival_rmst_y3 <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=36)
print(survival_rmst_y3)

fit.unstrat_y3 <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_y3, weights=Weights)
cox.zph(fit.unstrat_y3)

# First 5 years
pdac_survival_y5 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalMonths > 60, 0, SurvivalStatus))
results_y5 <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ Cluster, data=pdac_survival_y5, 
                                            max.num.perms = 1000, seed=42)
print(results_y5)

survival_rmst_y5 <- rmst2(pdac_survival$SurvivalMonths, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=60)
print(survival_rmst_y5)

fit.unstrat_y5 <- coxph(Surv(SurvivalMonths, SurvivalStatus) ~ Cluster, data=pdac_survival_y5, weights=Weights)
cox.zph(fit.unstrat_y5)


## OS-CNA Association ----

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

pdac_survival_cna <- pdac_survival[rownames(cna_tcga_ID),]

df_results <- data.frame()
for (k in colnames(cna_tcga_ID)) {
  cna_band <- cna_tcga_ID[[k]]
  row <- list()
  row$Band <- k
  group <- as.numeric(cna_band == 0)+1
  if (length(unique(group)) == 1) {
    next
  }
  pdac_survival_cna$group <- group
  pdac_survival_cna$group_factor <- as.factor(group)
  invisible(capture.output(logrank <- logrankHeinze::logrank.heinze(SurvivalMonths + SurvivalStatus ~ group, data=pdac_survival_cna,
                                                                    max.num.perms = 1000, seed=42)$pvalue))
  row$logrank <- logrank
  survival_rmst <- rmst2(pdac_survival_cna[["SurvivalMonths"]], pdac_survival_cna[["SurvivalStatus"]], group-1, tau = 24)
  row$rmst_coef <- survival_rmst$unadjusted.result[1,1]
  row$rmst_pval <- survival_rmst$unadjusted.result[1,4]
  cox_model <- coxph(Surv(pdac_survival_cna[["SurvivalMonths"]], pdac_survival_cna[["SurvivalStatus"]]) ~ as.factor(group))
  row$cox_coef <- summary(cox_model)$coefficients[,"exp(coef)"]
  row$cox_pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
  dd <- datadist(pdac_survival_cna)
  options(datadist = 'dd')
  fit <- cph(Surv(SurvivalMonths, SurvivalStatus) ~ group_factor, data=pdac_survival_cna, x=TRUE, y=TRUE)
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
  theme_cowplot()



