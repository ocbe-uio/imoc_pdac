
#################
##  LIBRARIES  ##
#################

install.packages("Gmisc")
install.packages("UpSetR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager:: install()
BiocManager::valid()
BiocManager::install("curatedTCGAData", force = TRUE)
BiocManager::install("TCGAutils")
BiocManager::install("curatedTCGAData")
install("BiocFileCache")
install.packages("xfun")
install.packages("knitr")

BiocManager::install("foreign")
BiocManager::install("MASS")
BiocManager::install("nlme")
BiocManager::install("survival")
BiocManager::install("TCGAbiolinks")

BiocManager::install("maftools")

install.packages("devtools")
devtools::install_version("dbplyr", version = "2.3.4")

remove.packages("xfun")
install.packages("xfun", type = "source")

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

library(BiocManager)
library(curatedTCGAData)
library(TCGAutils)
library(UpSetR)
library(Gmisc)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

##############
##  SCRIPT  ##
##############

# set data directory to download_data.R location
setwd("../imoc_pdac/data")
folder_data <- "TCGA/omics_data"
folder_raw_data <- file.path(folder_data, "raw")

# list diseaseCodes
data('diseaseCodes', package = "TCGAutils")

# choose the code
TCGA_CODE <- "PAAD"
CANCER_CODE <- "PDAC"
PLATFORM_CODE <- "TCGA"

# check experiments
metadata <- curatedTCGAData(diseaseCode = TCGA_CODE, version = "2.0.1")

# filter by omics
omics <- c("RNASeq2GeneNorm*", "RPPA*", "*Methylation*", "*miRNA*", "Mutation", "GISTIC_Peaks")
curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, version = "2.0.1")

# download data
cancer_data <- curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, dry.run = FALSE, version = "2.0.1")

# filter by primary tumor
sampleTables(cancer_data)
tums <- TCGAsampleSelect(barcodes = colnames(cancer_data), sampleCodes = "01")
cancer_data <- cancer_data[,tums,]

# visualize data
upsetSamples(cancer_data)

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

# visualize filtered data
upsetSamples(cancer_data)

# Visualise type of mutations in each gene (optional)
rag <- "PAAD_Mutation-20160128"

library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)   # older version

oncoPrintTCGA(cancer_data, matchassay = rag)
# Note that hg19 was used, instead of hg38
# check oncoprint package



# save object
#filename_rds <- pathJoin(folder_raw_data, paste0("rnaseqnorm_meth_rppa_mirna_", CANCER_CODE, PLATFORM_CODE, ".rds"))
#saveRDS(cancer_data, filename_rds)
getwd()
out_dir <- file.path("TCGA", "omics_data", "raw")
exportClass(cancer_data, dir = out_dir, fmt = "csv", ext = ".csv")
# CNV data (since we have to get all the information in the matrix to analyse)
row_data_cna <- as.data.frame(rowData(cancer_data[[1]]))
assay_data_cna <- as.data.frame(assay(cancer_data[[1]]))
merged_data_cna <- cbind(row_data_cna, assay_data_cna)
write.csv(merged_data_cna, file = "TCGA/omics_data/raw/cancer_data_PAAD_CNA_GISTIC-20160128.csv", row.names = FALSE)

# Methylation data (genes associated with methylation sites)
row_data_methyl <- as.data.frame(rowData(cancer_data[[6]]))
row_data_methyl$Methylation_Site <- rownames(rowData(cancer_data[[6]]))
row_data_methyl <- row_data_methyl[, c("Methylation_Site", colnames(row_data_methyl)[1:3])]
write.csv(row_data_methyl, file = "TCGA/omics_data/raw/associated_genes_methyl.csv", row.names = FALSE)