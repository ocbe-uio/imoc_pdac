pak::pkg_install("survminer")
devtools::install_github("Shamir-Lab/Logrank-Inaccuracies/logrankHeinze")
install.packages("survRM2")
install.packages("ggfortify")
install.packages("readxl")
install.packages("tidyverse")
library("survival")
library("survRM2")
library("logrankHeinze")
library("survival")
library("readxl")
library("dplyr")
library("ggfortify")
library("ggsurvfit")
library("rms")
library("cowplot")
library("tidyverse")


cna_data <- read.table("pancreas_cptac_gdc/data_cna.txt", header = 1)
colnames(cna_data) <- sub("^([^.]*\\.[^.]*)\\..*", "\\1", colnames(cna_data))
colnames(cna_data) <- gsub("\\.", "-", colnames(cna_data))
cna_data <- na.omit(cna_data)
head(cna_data)
survival_data <- read.table("pancreas_cptac_gdc/PDAC_survival.txt", header = T)
survival_data["OS_STATUS"] <- survival_data$OS_event
survival_data["OS_MONTHS"] <- survival_data$OS_days / 30
survival_data["PATIENT_ID"] <- survival_data$case_id
survival_data <- survival_data[!is.na(survival_data$OS_STATUS),]
survival_data <- survival_data[!is.na(survival_data$OS_MONTHS),]
patients <- intersect(survival_data$PATIENT_ID, colnames(cna_data))
rownames(survival_data) <- survival_data$PATIENT_ID
filtered_survival_data <- survival_data[patients,]
cna_data <- cbind(cna_data$Entrez_Gene_Id, cna_data[,patients])
filtered_survival_data <- filtered_survival_data[c("OS_STATUS", "OS_MONTHS")]
filtered_survival_data["OS_STATUS"] <- as.numeric(substr(filtered_survival_data[["OS_STATUS"]], 1, 1))
head(filtered_survival_data)
write.csv(cna_data, "pancreas_cptac_gdc/filtered_cna_data.csv", row.names = FALSE)

km_fit <- survfit2(Surv(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]]) ~ 1)
km_fit %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark(shape = 3, size = 10) +
  scale_ggsurvfit()

clusters_mean <- read.csv("clusters_mean.csv")
rownames(clusters_mean) <- clusters_mean$X
clusters_mean$X <- NULL
filtered_survival_data$cluster <- clusters_mean$X0
filtered_survival_data$cluster_factor <- as.factor(clusters_mean$X0)
row <- list()
invisible(capture.output(logrank <- logrankHeinze::logrank.heinze(OS_MONTHS + OS_STATUS ~ cluster, data=filtered_survival_data,
                                                                  max.num.perms = 1000, seed=42)$pvalue))
row$logrank <- logrank
survival_rmst <- rmst2(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]], filtered_survival_data$cluster-1)
row$rmst_coef <- survival_rmst$unadjusted.result[1,1]
row$rmst_pval <- survival_rmst$unadjusted.result[1,4]
cox_model <- coxph(Surv(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]]) ~ filtered_survival_data$cluster_factor)
row$cox_coef <- summary(cox_model)$coefficients[,"exp(coef)"]
row$cox_pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
dd <- datadist(filtered_survival_data)
options(datadist = 'dd')
fit <- cph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_factor, data=filtered_survival_data, x=TRUE, y=TRUE)
set.seed(42)
invisible(capture.output(v <- validate(fit, method="boot", B=1000)))
c_index_results <- (v["Dxy", ] / 2) + 0.5
row$c_index <- c_index_results["index.corrected"]
row
km_fit <- survfit2(Surv(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]]) ~ filtered_survival_data$cluster)
km_fit %>%
  ggsurvfit() +
  add_confidence_interval() +
  add_risktable() +
  add_censor_mark() +
  scale_ggsurvfit() +
  annotate("text", x = 50, y = 0.85, label = paste("RMST=", round(row$rmst_coef, 3))) +
  annotate("text", x = 50, y = 0.75, label = paste("p =", round(row$rmst_pval, 3)))
  


tmb <- read.delim("pancreas_cptac_gdc/data_clinical_sample.txt", comment.char = "#", sep = "\t", header = TRUE)
tmb <- distinct(tmb, PATIENT_ID, .keep_all = TRUE)
tmb <- tmb[tmb$PATIENT_ID %in% patients,]
rownames(tmb) <- tmb$PATIENT_ID
tmb <- tmb[patients,]
filtered_survival_data_tmb <- filtered_survival_data[!is.na(tmb$TMB_NONSYNONYMOUS),]
tmb <- tmb[!is.na(tmb$TMB_NONSYNONYMOUS),]
filtered_survival_data_tmb$tmb <- tmb$TMB_NONSYNONYMOUS
p <- kruskal.test(tmb ~ cluster_factor, data = filtered_survival_data_tmb)$p.value
ggplot(filtered_survival_data_tmb, aes(x=cluster_factor, y=tmb, fill=cluster_factor)) +
  geom_boxplot(notch = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Clusters", y = "Tumor Mutational Burden") +
  theme_cowplot() +
  theme(legend.position = "none") + 
  annotate("text", x = 1.5, y = 10, label = paste("p =", round(p, 3)))


X <- read.csv("X_mean.csv", header = T)
rownames(X) <- X$X
X$X <- NULL
colnames(X) <- substring(colnames(X), 2)
df_results <- data.frame()
for (cna_band_name in colnames(X)) {
  cna_band <- X[[cna_band_name]]
  row <- list()
  row$Band <- cna_band_name
  group <- as.numeric(cna_band >= 0)+1
  if (length(unique(group)) == 1) {
    next
  }
  filtered_survival_data$group <- group
  filtered_survival_data$group_factor <- as.factor(group)
  invisible(capture.output(logrank <- logrankHeinze::logrank.heinze(OS_MONTHS + OS_STATUS ~ group, data=filtered_survival_data,
                                                                    max.num.perms = 1000, seed=42)$pvalue))
  row$logrank <- logrank
  survival_rmst <- rmst2(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]], group-1)
  row$rmst_coef <- survival_rmst$unadjusted.result[1,1]
  row$rmst_pval <- survival_rmst$unadjusted.result[1,4]
  cox_model <- coxph(Surv(filtered_survival_data[["OS_MONTHS"]], filtered_survival_data[["OS_STATUS"]]) ~ as.factor(group))
  row$cox_coef <- summary(cox_model)$coefficients[,"exp(coef)"]
  row$cox_pval <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
  dd <- datadist(filtered_survival_data)
  options(datadist = 'dd')
  fit <- cph(Surv(OS_MONTHS, OS_STATUS) ~ group_factor, data=filtered_survival_data, x=TRUE, y=TRUE)
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


cna_tcga_ID <- merge(X, clusters_mean, by = 0)
cna_tcga_ID$cluster <- as.factor(cna_tcga_ID$X0)
cna_tcga_ID$X0 <- NULL
rownames(cna_tcga_ID) <- cna_tcga_ID$Row.names
cna_tcga_ID$Row.names <- NULL
cna_tcga_ID_long <- cna_tcga_ID %>% pivot_longer(cols = 1:4, names_to = "Descriptor", values_to = "CNA")
cna_tcga_ID_long$CNA[cna_tcga_ID_long$CNA >= 0] <- NA
contin_cna <- cna_tcga_ID_long %>% group_by(Descriptor, cluster) %>% summarise(Altered = sum(!is.na(CNA)))
contin_cna <- contin_cna %>% mutate(Normal = c(32,121) - Altered)
contingency_tables <- split(contin_cna, contin_cna$Descriptor) %>%
  lapply(function(sub_df) {
    table_matrix <- as.matrix(sub_df[, c("Altered", "Normal")])
    rownames(table_matrix) <- sub_df$cluster
    return(table_matrix)
  })
cna_fisher <- lapply(contingency_tables, FUN = fisher.test)
pvals <- cna_fisher %>% map_dfr(\(x){
  tibble(pval = x %>% pluck("p.value"))
})
pvals$Descriptor <- names(contingency_tables)
pvals$fisher_adjpval <- p.adjust(pvals$pval, method = "BH")
OddsRatio <- function(mat) {
  a <- mat[2] + 0.5
  b <- mat[4] + 0.5 
  c <- mat[1] + 0.5
  d <- mat[3] + 0.5
  OR <- ((a)/(b))/((c)/(d))
  SElogOR <- sqrt(1/a+1/b+1/c+1/d)
  lwr <- exp(log(OR) - 1.96*SElogOR)
  upr <- exp(log(OR) + 1.96*SElogOR)
  OR <- data.frame(OR = OR, lwr = lwr, upr = upr)
  return(OR)
}
OR <- lapply(contingency_tables, OddsRatio)
OR <- imap_dfr(OR, ~mutate(.x, Descriptor = .y))
OR <- merge(pvals, OR, by = "Descriptor")
openxlsx2::write_xlsx(cna_tcga, "results/omics_analysis/CNA/cna_tcga.xlsx")
OR %>% 
  arrange(OR) %>% 
  mutate(Descriptor = factor(Descriptor, levels = Descriptor)) %>% 
  ggplot(aes(x = OR, y = Descriptor)) +
  geom_segment(aes(xend = 0, yend = Descriptor), color = "gray80") +
  geom_point(aes(fill = -log10(fisher_adjpval)), size = 4, shape=23) +
  scale_fill_gradient(high = "#fc0352", low = "white", name = "-log10 Fisher's\ntest FDR", limits = c(0,4.5)) +
  theme_cowplot() +
  labs(title = "Odds ratio CNA", subtitle = "External validation") +
  xlab("Odds Ratio") + ylab("Deletion")

