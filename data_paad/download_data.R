
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
# filter by omics
omics <- c("RNASeq2GeneNorm*", "RPPA*", "*Methylation*", "*miRNA*", "Mutation")
curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, version = "2.0.1")
# download data
cancer_data <- curatedTCGAData(diseaseCode = TCGA_CODE, assays = omics, dry.run = FALSE, version = "2.0.1")
cancer_data
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

################################################################################

#### CNV SEGMENT MEAN

BiocManager::install("GenomicRanges", force = TRUE)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("biomaRt")

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(biomaRt)

# Extract the gene sequences (index column)
cancer_data[[1]]
CNV_mean <- assay(cancer_data[[1]], "Segment_Mean")
index_column <- rownames(CNV_mean)
print(index_column)

# Setting up a gene annotation template to use
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"), mart=mart)
genes <- genes[genes[,1]!="" & genes[,2] %in% c(1:22,"X","Y"),]
xidx <- which(genes[,2]=="X")
yidx <- which(genes[,2]=="Y")
genes[xidx, 2] <- 23
genes[yidx, 2] <- 24
genes[,2] <- sapply(genes[,2],as.integer)
genes <- genes[order(genes[,3]),]
genes <- genes[order(genes[,2]),]
colnames(genes) <- c("GeneSymbol","Chr","Start","End")

genes_GR <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

# Function to parse each segment
parse_segment <- function(segment) {
  # Split by colons and dashes
  parts <- unlist(strsplit(segment, "[:-]"))
  # Create a list with chromosome, start, end, and placeholders for extra fields
  list(
    chr = parts[1],
    start = as.numeric(parts[2]),
    end = as.numeric(parts[3])
  )
}

# Apply the function to all rows in the index column
parsed_segments <- lapply(index_column, parse_segment)

# Convert the list of parsed segments to a data frame
df <- do.call(rbind, lapply(parsed_segments, as.data.frame))

# Ensure the dataframe columns are correctly named
colnames(df) <- c("chr", "start", "end")
print(df)
dim(df)

# Store data as GenomicRanges object
df_GR <- makeGRangesFromDataFrame(df, na.rm = TRUE)   # Removed 13 NA values (can't work with them because we can't map them)
df_GR

# Overlap regions with reference dataset created
#hits <- findOverlaps(df_GR, genes_GR, type="within")
hits <- findOverlaps(genes_GR, df_GR, type="within")
df_ann <- cbind(df[subjectHits(hits),],genes[queryHits(hits), "GeneSymbol"])
df_ann <- unique(df_ann)   # remove duplicate rows
head(df_ann)
dim(df_ann)
print(df_ann)

# Check how many genes appear per sequence
# Define variables for chr, start, and end
chr_value <- 1  # Example value for chr
start_value <- 61735  # Example value for start
end_value <- 98588  # Example value for end

# Find rows where chr, start, and end match the defined values
matching_rows <- df_ann$chr == chr_value & df_ann$start == start_value & df_ann$end == end_value

# Subset df_ann to show only matching rows
matched_df_ann <- df_ann[matching_rows, ]

# Print the matched rows
print(matched_df_ann)

# Convert to a data frame if necessary (in case it's not already)
#CNV_mean_df <- as.data.frame(CNV_mean)

# Define the file path where you want to save the CSV
#file_path <- "CNV_mean.csv"

# Write the data to a CSV file
#write.csv(CNV_mean_df, file = file_path, row.names = TRUE)

# Print a message indicating the file has been saved
#cat("CNV_mean data has been saved to", file_path)

################################################################################

# save object
#filename_rds <- pathJoin(folder_raw_data, paste0("rnaseqnorm_meth_rppa_mirna_", CANCER_CODE, PLATFORM_CODE, ".rds"))
#saveRDS(cancer_data, filename_rds)
exportClass(cancer_data, dir = folder_raw_data, fmt = "csv", ext = ".csv")



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
