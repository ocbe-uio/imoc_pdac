from sklearn.pipeline import make_pipeline
from preprocessing_transformers import (InitialProcessing, RemoveFeaturesWithZeros, RemoveFeaturesWithNaN, RemoveFeaturesLowMAE,
                                        RemoveCorrelatedFeatures, Log2Transformation, GeneMutations, ValueImputation)


# OMIC DATA TYPES

# RNA
RNAseq_data = InitialProcessing("data_paad/raw_data/cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv").process_data()
RNAseq_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAE(percentage_to_keep=0.1, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True),
    Log2Transformation()
)
RNAseq_new = RNAseq_pipeline.fit_transform(RNAseq_data)
RNAseq_new.to_csv("data_paad/preprocessed_data_missing_views/processed_cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv")

# Proteins (RPPA)
RPPA_data = InitialProcessing("data_paad/raw_data/cancer_data_PAAD_RPPAArray-20160128.csv").process_data()
RPPA_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    ValueImputation(verbose=True)
)
RPPA_new = RPPA_pipeline.fit_transform(RPPA_data)
RPPA_new.to_csv("data_paad/preprocessed_data_missing_views/processed_cancer_data_PAAD_RPPAArray-20160128.csv")

# miRNA
miRNA_data = InitialProcessing("data_paad/raw_data/cancer_data_PAAD_miRNASeqGene-20160128.csv").process_data()
miRNA_pipeline = make_pipeline(
    RemoveFeaturesWithZeros(threshold=0.2, verbose=True),
    Log2Transformation()
)
miRNA_new = miRNA_pipeline.fit_transform(miRNA_data)
miRNA_new.to_csv("data_paad/preprocessed_data_missing_views/processed_cancer_data_PAAD_miRNASeqGene-20160128.csv")

# Methylation
methylation_data = InitialProcessing("data_paad/raw_data/cancer_data_PAAD_Methylation-20160128.csv").process_data()
methylation_pipeline = make_pipeline(
    RemoveFeaturesWithNaN(threshold=0.2, verbose=True),
    RemoveFeaturesLowMAE(percentage_to_keep=0.01, verbose=True),
    RemoveCorrelatedFeatures(threshold=0.85, verbose=True),
    ValueImputation(verbose=True)
)
methylation_new = methylation_pipeline.fit_transform(methylation_data)
methylation_new.to_csv("data_paad/preprocessed_data_missing_views/processed_cancer_data_PAAD_Methylation-20160128.csv")

# Mutations
mutations_data = InitialProcessing("data_paad/raw_data/cancer_data_PAAD_Mutation-20160128.csv").process_data()
mutations_pipeline = make_pipeline(
    GeneMutations(verbose=True),
    RemoveFeaturesWithZeros(threshold=0.05, verbose=True)
)
mutations_new = mutations_pipeline.fit_transform(mutations_data)
mutations_new.to_csv("data_paad/preprocessed_data_missing_views/processed_cancer_data_PAAD_Mutation-20160128.csv")


samples_complete = RPPA_new.index.intersection(miRNA_new.index).intersection(RNAseq_new.index).intersection(methylation_new.index).intersection(mutations_new.index)
RPPA_new, miRNA_new, RNAseq_new, methylation_new, mutation_new = RPPA_new.loc[samples_complete], miRNA_new.loc[samples_complete], RNAseq_new.loc[samples_complete], methylation_new.loc[samples_complete], mutations_new.loc[samples_complete]

RNAseq_new.to_csv("data_paad/complete_samples_data/complete_cancer_data_PAAD_RNASeq2GeneNorm-20160128.csv")
RPPA_new.to_csv("data_paad/complete_samples_data/complete_cancer_data_PAAD_RPPAArray-20160128.csv")
miRNA_new.to_csv("data_paad/complete_samples_data/complete_cancer_data_PAAD_miRNASeqGene-20160128.csv")
methylation_new.to_csv("data_paad/complete_samples_data/complete_cancer_data_PAAD_Methylation-20160128.csv")
mutations_new.to_csv("data_paad/complete_samples_data/complete_cancer_data_PAAD_Mutation-20160128.csv")
