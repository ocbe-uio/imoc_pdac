
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


