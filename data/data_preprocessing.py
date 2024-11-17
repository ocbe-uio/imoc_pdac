import pandas as pd
import numpy as np
from sklearn.pipeline import make_pipeline
from preprocessing_transformers import (InitialProcessing, RemoveFeaturesWithZeros, RemoveFeaturesWithNaN, RemoveFeaturesLowMAD,
                                        RemoveCorrelatedFeatures, Log2Transformation, GeneMutations, ValueImputation)

complete_sample_set = True   # True if all samples are being used, False if subset with complete information across all views
add_noise = False

# OMIC DATA TYPES

# RNAseq
RNAseq_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv").process_data()
RNAseq_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAD(percentage_to_keep=0.1, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True),
    Log2Transformation()
)
RNAseq_new = RNAseq_pipeline.fit_transform(RNAseq_data)

# Proteins (RPPA)
RPPA_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_RPPAArray-20160128.csv").process_data()
RPPA_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    ValueImputation(verbose=True)
)
RPPA_new = RPPA_pipeline.fit_transform(RPPA_data)

# miRNA
miRNA_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_miRNASeqGene-20160128.csv").process_data()
miRNA_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    Log2Transformation()
)
miRNA_new = miRNA_pipeline.fit_transform(miRNA_data)

# Methylation
methylation_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_Methylation-20160128.csv").process_data()
methylation_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAD(percentage_to_keep=0.01, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True)
)
methylation_new = methylation_pipeline.fit_transform(methylation_data)

# Mutations
mutations_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_Mutation-20160128.csv").process_data()
mutations_pipeline = make_pipeline(
    GeneMutations(verbose=True),
    RemoveFeaturesWithZeros(threshold=0.95, verbose=True)
)
mutations_new = mutations_pipeline.fit_transform(mutations_data)
if add_noise == True:
    mutations_new = mutations_new + np.random.default_rng(42).normal(scale=0.1, size=mutations_new.shape)

# Copy number
cnv_data = InitialProcessing("data/TCGA_PDAC_rawdata/cancer_data_PAAD_CNA_GISTIC-20160128.csv").process_data()
cnv_data.columns = cnv_data.loc['Descriptor']
cnv_new = cnv_data[cnv_data.index.str.contains('TCGA')]
if add_noise == True:
    cnv_new = cnv_new + np.random.default_rng(42).normal(scale=0.1, size=cnv_new.shape)


dataframes = [RNAseq_new, RPPA_new, miRNA_new, methylation_new, mutations_new, cnv_new]
if complete_sample_set == False:
    partial_samples = RPPA_new.index.intersection(miRNA_new.index).intersection(RNAseq_new.index).intersection(methylation_new.index).intersection(mutations_new.index).intersection(cnv_new.index)
    RNAseq_partial, RPPA_partial, miRNA_partial, methylation_partial, mutations_partial, cnv_partial = [df.loc[partial_samples] for df in dataframes]
elif complete_sample_set == True:
    complete_samples = pd.concat(dataframes, axis=0, join='outer').index.unique()
    RNAseq_complete, RPPA_complete, miRNA_complete, methylation_complete, mutations_complete, cnv_complete = [df.reindex(complete_samples) for df in dataframes]

RNAseq_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_RNAseq.csv', index=True)
RPPA_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_RPPA.csv', index=True)
miRNA_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_miRNA.csv', index=True)
methylation_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_Methylation.csv', index=True)
mutations_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_Mutation.csv', index=True)
cnv_complete.to_csv('data/TCGA_PDAC_all/TCGA_PDAC_all_CNA.csv', index=True)

print("Completed successfully!")
