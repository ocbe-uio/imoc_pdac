import pandas as pd
from sklearn.pipeline import make_pipeline
from preprocessing_transformers import (InitialProcessing, RemoveFeaturesWithZeros, RemoveFeaturesWithNaN, RemoveFeaturesLowMAD,
                                        RemoveCorrelatedFeatures, Log2Transformation, GeneMutations, ValueImputation, PrincipalFeatures)


# OMIC DATA TYPES

# RNA
RNAseq_data = InitialProcessing("data/raw_data/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv").process_data()
RNAseq_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAD(percentage_to_keep=0.1, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True),
    Log2Transformation()
)
RNAseq_new = RNAseq_pipeline.fit_transform(RNAseq_data)

# Proteins (RPPA)
RPPA_data = InitialProcessing("data/raw_data/cancer_data_PAAD_RPPAArray-20160128.csv").process_data()
RPPA_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    ValueImputation(verbose=True)
)
RPPA_new = RPPA_pipeline.fit_transform(RPPA_data)

# miRNA
miRNA_data = InitialProcessing("data/raw_data/cancer_data_PAAD_miRNASeqGene-20160128.csv").process_data()
miRNA_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    Log2Transformation()
)
miRNA_new = miRNA_pipeline.fit_transform(miRNA_data)

# Methylation
methylation_data = InitialProcessing("data/raw_data/cancer_data_PAAD_Methylation-20160128.csv").process_data()
methylation_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAD(percentage_to_keep=0.01, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True)
)
methylation_new = methylation_pipeline.fit_transform(methylation_data)

# Mutations
mutations_data = InitialProcessing("data/raw_data/cancer_data_PAAD_Mutation-20160128.csv").process_data()
mutations_pipeline = make_pipeline(
    GeneMutations(verbose=True),
    RemoveFeaturesWithZeros(threshold=0.05, verbose=True)
)
mutations_new = mutations_pipeline.fit_transform(mutations_data)


# Copy number
cnv_data = InitialProcessing("data/raw_data/cancer_data_PAAD_CNA-20160128.csv").process_data()
cnv_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    ValueImputation(),
    PrincipalFeatures(diff_n_features = 1000, explained_var = 0.99, verbose=True)
)
cnv_new = cnv_pipeline.fit_transform(cnv_data)

samples_complete = RPPA_new.index.intersection(miRNA_new.index).intersection(RNAseq_new.index).intersection(methylation_new.index).intersection(mutations_new.index).intersection(cnv_new.index)
RPPA_new, miRNA_new, RNAseq_new, methylation_new, mutations_new, cnv_new = RPPA_new.loc[samples_complete], miRNA_new.loc[samples_complete], RNAseq_new.loc[samples_complete], methylation_new.loc[samples_complete], mutations_new.loc[samples_complete], cnv_new.loc[samples_complete]
