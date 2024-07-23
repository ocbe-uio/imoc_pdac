
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

##############
##  SCRIPT  ##
##############

# data directory
folder_data <- "data_paad"
folder_raw_data <- file.path(folder_data, "raw")
# list diseaseCodes
data('diseaseCodes', package = "TCGAutils")
diseaseCodes
# choose the code
TCGA_CODE <- "PAAD"
CANCER_CODE <- "PDAC"
PLATFORM_CODE <- "TCGA"
# check experiments
metadata <- curatedTCGAData(diseaseCode = TCGA_CODE, version = "2.0.1")
metadata
# filter by omics
omics <- c("RNASeq2GeneNorm*", "RPPA*", "*Methylation*", "*miRNA*", "Mutation", "CNASNP")
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
cancer_data

# visualize filtered data
upsetSamples(cancer_data)

# Visualise type of mutations in each gene
cancer_data
rag <- "PAAD_Mutation-20160128"

library(ComplexHeatmap)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)   # older version


oncoPrintTCGA(cancer_data, matchassay = rag)
# Note that hg19 was used, instead of hg38
# check oncoprint package

###################################################

# The CNV data comes in segments rather than by genes, so we have to convert those segments into their corresponding gene(s)
# The code used to download the data and build the CNV matrices was produced by Hornung and Wright in the paper "Block Forests: random forests for blocks of clinical and omics covariate data."

# Reference:
# Hornung, R., Wright, M.N. Block Forests: random forests for blocks of clinical and omics covariate data. BMC Bioinformatics 20, 358 (2019). https://doi.org/10.1186/s12859-019-2942-y


library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

setwd("C:/Users/alepg/PycharmProjects/imo_clustering/data_paad")

# Download CNV data
query <- GDCquery(project = "TCGA-PAAD",
                  data.category = "Copy Number Variation",
                  data.type = "Copy Number Segment", 
                  sample.type = "Primary Tumor")
GDCdownload(query, files.per.chunk=20)
cnv <- GDCprepare(query)

# Gene names and information (chromosome, start, end)
load("loci_all.Rdata")
rows <- loci$ensembl_gene_id
rownames(loci) <- loci$ensembl_gene_id

# Build CNV function
build.cnv <- function(i, cnv, loci, cols) {
  chr.cnv <- cnv[which(cnv$Chromosome==loci$chromosome_name[i]), ]  # subset of cnv with only chromosome chr
  vec <- chr.cnv$Segment_Mean                   # CNV values for chromosome chr
  names(vec) <- chr.cnv$Sample
  positive.samples <- which(chr.cnv$Start<=loci$start_position[i] & chr.cnv$End>=loci$end_position[i])    # cnv values of the gene
  vec <- vec[positive.samples]
  res <- rep(NA, length(cols))
  names(res) <- cols
  ind <- which(names(res) %in% names(vec))
  res[ind] <- vec[names(res[ind])]
  return(res)
}

# Build CNV matrices
names.cnv <- cnv$Sample
cnv.out <- NULL
nn.cnv <- unique(names.cnv)
for(j in 1:length(nn.cnv)) {
  ind <- which(names.cnv==nn.cnv[j])
  tt <- table(cnv$Sample[ind])
  if(length(tt)>1) {
    cnv.out <- c(cnv.out, names(tt)[-1])
  }
}
ind <- which(!cnv$Sample %in% cnv.out)
cnv <- cnv[ind, ]

# retrieve CNV values for all ensembl IDs
cols <- unique(cnv$Sample)
mat.cnv <- apply(as.matrix(1:nrow(loci)), 1, build.cnv, cnv, loci, cols)
mat.cnv <- t(mat.cnv)
rownames(mat.cnv) <- loci$ensembl_gene_id
colnames(mat.cnv) <- cols
mat.cnv <- as.data.frame(mat.cnv)

# Merge both dataframes, keeping all information from loci (for now)
mat.cnv$ensembl_gene_id <- rownames(mat.cnv)
cnv_data <- merge(mat.cnv, loci)
dim(cnv_data)

# Remove rows containing genes in sexual chromosomes
cnv_data <- subset(cnv_data, !(chromosome_name %in% c("X", "Y")))
# Remove rows that have only NA values in patient columns
patient_columns <- grep("^TCGA-", colnames(cnv_data), value = TRUE)
cnv_data <- cnv_data[rowSums(is.na(cnv_data[, patient_columns])) < length(patient_columns), ]
head(cnv_data)
# Put ensembl IDs as index column
rownames(cnv_data) <- cnv_data$ensembl_gene_id
cnv_data$ensembl_gene_id <- NULL

# Create final dataset with gene names in index and patients as columns, getting median for information in duplicated genes
cnv_data_subset <- cnv_data[, c("external_gene_name", patient_columns)]
cnvmean_data <- cnv_data_subset %>%
  group_by(external_gene_name) %>%
  summarise(across(everything(), median, na.rm = TRUE))
cnvmean_data <- as.data.frame(cnvmean_data)
rownames(cnvmean_data) <- cnvmean_data$external_gene_name
cnvmean_data$external_gene_name <- NULL
head(cnvmean_data)
dim(cnvmean_data)


write.csv(cnvmean_data, "cancer_data_PAAD_CNA-20160128.csv")


################################################################################

# save object
#filename_rds <- pathJoin(folder_raw_data, paste0("rnaseqnorm_meth_rppa_mirna_", CANCER_CODE, PLATFORM_CODE, ".rds"))
#saveRDS(cancer_data, filename_rds)
exportClass(cancer_data, dir = folder_raw_data, fmt = "csv", ext = ".csv")

# FOR SUPERVISED LEARNING (not used in this project)

# split dataset in 60% as training set and remaining 40% as testing set
patients <- as.data.frame(patients)
dataset_type <- rep("training", nrow(patients))
dataset_type[1:as.integer(nrow(patients)*0.4)] <- "testing"
set.seed(42)
patients$dataset_type <- sample(dataset_type)
table(patients$dataset_type)
table(patients$dataset_type)/length(patients$dataset_type)
filename_samples <- pathJoin(folder_raw_data, paste0("patients_", CANCER_CODE, PLATFORM_CODE, ".csv"))
write.csv2(x = patients, file = filename_samples)
