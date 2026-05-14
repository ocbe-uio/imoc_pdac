
devtools::install_github("Shamir-Lab/Logrank-Inaccuracies/logrankHeinze")
install.packages("survRM2")
library("survRM2")
library("logrankHeinze")
library("dplyr")
library("survival")

# Change working directory
wd <- ".../imo_clustering/results/cluster_analysis"    # edit to set working directory
setwd(wd)

# Change file names for clinical data and RNAseq data
clinical_data <- read.csv('CPTAC_clinical_data.tsv', sep = "\t")
RNAseq_data <- read.csv('CPTAC_RNAseq.csv', sep="\t")

survival_data <- clinical_data[, c("case_id", "vital_status", "follow_up_days")]
row.names(RNAseq_data) <- RNAseq_data$X
RNAseq_data <- RNAseq_data[, -1]
RNAseq_data <- t(RNAseq_data)
RNAseq_data <- as.data.frame(RNAseq_data)

biomarkers <- c("CDKN2A", "DNAH9", "ZNF18", "MAP2K4", "CTIF", "LINC01667", 
                "TEKT4P2", "RNA5-8SN2", "MEOX2", "SIX6", "TBX2")
biomarkers_data <- RNAseq_data[, biomarkers]

for (gene in biomarkers) {
  gene_mean <- mean(biomarkers_data[[gene]], na.rm = TRUE)
  col_name <- paste0(gene, "_Group")
  biomarkers_data[[col_name]] <- ifelse(biomarkers_data[[gene]] > gene_mean, "2", "1")
}
biomarkers_data <- cbind(case_id = rownames(biomarkers_data), biomarkers_data)
biomarkers_data$case_id <- gsub("\\.", "-", biomarkers_data$case_id)
rownames(biomarkers_data) <- 1:nrow(biomarkers_data)

alldata_survival <- merge(survival_data, biomarkers_data, by = "case_id")
row.names(alldata_survival) <- alldata_survival$case_id
alldata_survival <- alldata_survival[, -1]

# RUN FROM HERE EVERY TIME TO TEST A NEW BIOMARKER
# Change gene name with whatever biomarker we are testing
gene <- "TBX2"    # change gene here
# Note: for "RNA5-8SN2", input "RNA5.8SN2"

cluster_totest <- paste(gene, "Group", sep="_")
pdac_survival <- data.frame(alldata_survival)
pdac_survival <- pdac_survival[, c("vital_status", "follow_up_days", cluster_totest)]
pdac_survival <- na.omit(pdac_survival)

colnames(pdac_survival) <- c("SurvivalStatus", "SurvivalDays", "Cluster")
pdac_survival$Cluster <- as.numeric(as.factor(pdac_survival$Cluster))
pdac_survival[pdac_survival == "Living"] <- 0
pdac_survival[pdac_survival == "Deceased"] <- 1
pdac_survival$SurvivalStatus <- as.numeric(pdac_survival$SurvivalStatus)


### Log-rank tests, RMST test at different timepoints and checking if PH model can be used

# Complete survival timeline
result <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival, 
                                        max.num.perms = 1000, seed=42)
print(result)

max_time_cluster_1 <- max(pdac_survival$SurvivalDays[pdac_survival$Cluster == "1"])
max_time_cluster_2 <- max(pdac_survival$SurvivalDays[pdac_survival$Cluster == "2"])
tau <- min(max_time_cluster_1, max_time_cluster_2)

survival_rmst <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                       as.numeric(pdac_survival$Cluster == "1"), tau=tau)
print(survival_rmst)

# Un-comment to test proportional hazards assumption (Cox model)
#fit.unstrat <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival)
# cox.zph(fit.unstrat)


# First six months
pdac_survival_6m <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 180, 0, SurvivalStatus))
results_6m <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_6m, 
                                            max.num.perms = 1000, seed=42)
print(results_6m)

survival_rmst_6m <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=180)
print(survival_rmst_6m)

#fit.unstrat_6m <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_6m)
#cox.zph(fit.unstrat_6m)


# First year
pdac_survival_y1 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 360, 0, SurvivalStatus))
results_y1 <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_y1, 
                                            max.num.perms = 1000, seed=42)
print(results_y1)

survival_rmst_y1 <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=360)
print(survival_rmst_y1)

#fit.unstrat_y1 <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_y1)
#cox.zph(fit.unstrat_y1)


# First 18 months (1.5 years)
pdac_survival_18m <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 540, 0, SurvivalStatus))
results_18m <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_18m, 
                                            max.num.perms = 1000, seed=42)
print(results_18m)

survival_rmst_18m <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                           as.numeric(pdac_survival$Cluster == "1"), tau=540)
print(survival_rmst_18m)

#fit.unstrat_18m <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_18m)
#cox.zph(fit.unstrat_18m)

# First 2 years
pdac_survival_y2 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 720, 0, SurvivalStatus))
results_y2 <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_y2, 
                                            max.num.perms = 1000, seed=42)
print(results_y2)

survival_rmst_y2 <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=720)
print(survival_rmst_y2)

#fit.unstrat_y2 <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_y2)
#cox.zph(fit.unstrat_y2)

# First 3 years
pdac_survival_y3 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 1080, 0, SurvivalStatus))
results_y3 <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_y3, 
                                            max.num.perms = 1000, seed=42)
print(results_y3)

survival_rmst_y3 <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=1080)
print(survival_rmst_y3)

#fit.unstrat_y3 <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_y3)
#cox.zph(fit.unstrat_y3)



# ADDITIONALTIME(not in CPTAC, it only gets to a bit over 3 years)
# First 5 years
pdac_survival_y5 <- pdac_survival %>%
  mutate(SurvivalStatus = ifelse(SurvivalDays > 60, 0, SurvivalStatus))
results_y5 <- logrankHeinze::logrank.heinze(SurvivalDays + SurvivalStatus ~ Cluster, data=pdac_survival_y5, 
                                            max.num.perms = 1000, seed=42)
print(results_y5)

survival_rmst_y5 <- rmst2(pdac_survival$SurvivalDays, pdac_survival$SurvivalStatus, 
                          as.numeric(pdac_survival$Cluster == "1"), tau=60)
print(survival_rmst_y5)

fit.unstrat_y5 <- coxph(Surv(SurvivalDays, SurvivalStatus) ~ Cluster, data=pdac_survival_y5)
cox.zph(fit.unstrat_y5)

