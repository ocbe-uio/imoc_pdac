#install.packages("tidyverse")


library(tidyverse)

setwd("./data/TCGA/comparison_papers/")

con <- file("TCGA_PAAD.rds", "rb")
readBin(con, "raw", 10)
close(con)

# Upload my clusters
my_clusters <- read.table("../../patient_clusters.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
baseline_clusters <- my_clusters %>%
  pivot_longer(
    cols = -X,
    names_to = NULL,
    values_to = "Patient.ID",
    values_drop_na = TRUE
  ) %>%
  rename(Cluster = X) %>%
  filter(Patient.ID != "") %>%
  select(Patient.ID, Cluster)
baseline_clusters
baseline_clusters$Cluster[baseline_clusters$Cluster == "Cluster_0"] <- 1
baseline_clusters$Cluster[baseline_clusters$Cluster == "Cluster_1"] <- 2
baseline_clusters$Cluster <- as.numeric(baseline_clusters$Cluster)
colnames(baseline_clusters) <- c('PatientID', 'My Clusters')
baseline_clusters


data <- readRDS("TCGA_PAAD.rds")
bailey <- data.matrix(data[["sampInfo"]][["Bailey.Clusters"]])
collisson <- data.matrix(data[["sampInfo"]][["Collisson.Clusters"]])
copy_number <- data.matrix(data[["sampInfo"]][["Copy.Number.Clusters.All.150.Samples"]])
incrna <- data.matrix(data[["sampInfo"]][["lncRNA.Clusters.All.150.Samples"]])
moffitt <- data.matrix(data[["sampInfo"]][["Moffitt.Clusters"]])
moffitt_highpurity <- data.matrix(data[["sampInfo"]][["Moffitt.HighPurity.Clusters"]])
methyl_highpurity <- data.matrix(data[["sampInfo"]][["Methylation.Clusters.76.High.Purity.Samples.Only"]])
methyl <- data.matrix(data[["sampInfo"]][["Methylation.Clusters.All.150.Samples"]])
mirna <- data.matrix(data[["sampInfo"]][["miRNA.Clusters.All.150.Samples"]])
rppa_highpurity <- data.matrix(data[["sampInfo"]][["RPPA.Clusters.76.High.Purity.Samples.Only"]])
rppa <- data.matrix(data[["sampInfo"]][["RPPA.Clusters.All.150.Samples"]])
patients <- data.matrix(data[["sampInfo"]][["Tumor.Sample.ID"]])
patients <- substr(patients, 1, 12)
clusters_paper <- data.frame(cbind(patients, bailey, collisson, copy_number, incrna, moffitt, 
                        moffitt_highpurity, methyl, methyl_highpurity, mirna, rppa, rppa_highpurity))
colnames(clusters_paper) <- c('PatientID', 'bailey', 'collisson', 'copy_number', 
                              'incrna', 'moffitt', 'moffitt_highpurity', 'methyl', 
                              'methyl_highpurity', 'mirna', 'rppa', 'rppa_highpurity')

all_clusters <- merge(baseline_clusters, clusters_paper, by = "PatientID", all.x = TRUE)
all_clusters_cleaned <- all_clusters %>%
  filter(rowSums(!is.na(select(., -PatientID, -'My Clusters'))) > 0)
all_clusters_cleaned

write.csv(all_clusters_cleaned, "clusters_papers.csv", row.names = FALSE)
