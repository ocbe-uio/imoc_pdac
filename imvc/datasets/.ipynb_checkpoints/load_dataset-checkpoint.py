import os
from os.path import dirname

import pandas as pd
from imvc.utils import DatasetUtils


class LoadDataset:

    @staticmethod
    def load_incomplete_PDAC_partial_samples(p, return_y: bool = False, shuffle: bool = True,
                                             assess_percentage: bool = True, random_state: int = None):
        r"""
        Load an incomplete multi-view version of the PDAC_partial_samples dataset, a six-view dataset from a multi-omic
        study on a subset of PDAC patients.

        Parameters
        ----------
        p: list or float
            The percentage that each view will have for missing samples. If p is float, all the views will have the
            same percentage
        return_y: bool, default False
            If True, return the label too
        shuffle: bool, default False
            If True, shuffle the dataset.
        assess_percentage: bool
            If False, each view is dropped independently.
        random_state: int, default None
            If int, random_state is the seed used by the random number generator.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.
        ys : optional list of array-likes
            Array with labels

        Notes
        -----
        This data consists of six views from a TCGA multi-omics study of 89 PDAC patients:
        - RPPA : protein expression of 192 proteins
        - miRNA: miRNA expression of 385 miRNA
        - RNAseq: gene expression of 1419 genes
        - methylation: methylation of 2185 CpG sites
        - mutations: mutations (1 = yes / 0 = no) of 71 genes
        - CNA: copy number alterations (GISTIC2.0 format) of 745 genes
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_PDAC_partial_samples(p=0)
        """
        filenames = ["PDAC_partial_samples_RPPA.csv",
                     "PDAC_partial_samples_miRNA.csv",
                     "PDAC_partial_samples_RNAseq.csv",
                     "PDAC_partial_samples_Methylation.csv",
                     "PDAC_partial_samples_Mutation.csv",
                     "PDAC_partial_samples_CNA.csv"]
        module_path = dirname(__file__)
        Xs = [pd.read_csv(os.path.join(module_path, "data", "PDAC_partial_samples", filename)) for filename in filenames]
        Xs = DatasetUtils.ampute(Xs=Xs, p=p, assess_percentage=assess_percentage, random_state=random_state)
        if shuffle:
            Xs = DatasetUtils.shuffle_imvd(Xs=Xs, random_state=random_state)
        if return_y:
            ys = None
            out = (Xs, ys)
        else:
            out = Xs
        return out


    @staticmethod
    def load_incomplete_PDAC_complete_samples(p, return_y: bool = False, shuffle: bool = True,
                                             assess_percentage: bool = True, random_state: int = None):
        r"""
        Load an incomplete multi-view version of the PDAC_complete_samples dataset, a two-view dataset from a multi-omic
         study on a complete_set of PDAC patients
        patients.

        Parameters
        ----------
        p: list or float
            The percentage that each view will have for missing samples. If p is float, all the views will have the
            same percentage
        return_y: bool, default False
            If True, return the label too
        shuffle: bool, default False
            If True, shuffle the dataset.
        assess_percentage: bool
            If False, each view is dropped independently.
        random_state: int, default None
            If int, random_state is the seed used by the random number generator.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.
        ys : optional list of array-likes
            Array with labels

        Notes
        -----
        This data consists of two views from a TCGA multi-omics study of 154 PDAC patients:
        - mutations: mutations (1 = yes / 0 = no) of 71 genes
        - CNA: copy number alterations (GISTIC2.0 format) of 745 genes
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_PDAC_complete_samples(p=0)
        """
        module_path = dirname(__file__)
        Xs = [pd.read_csv(os.path.join(module_path, "data", "PDAC_complete_samples", filename)) for filename in
              ["PDAC_complete_samples_Mutation.csv", "PDAC_complete_samples_CNA.csv"]]
        Xs = DatasetUtils.ampute(Xs=Xs, p=p, assess_percentage=assess_percentage, random_state=random_state)
        if shuffle:
            Xs = DatasetUtils.shuffle_imvd(Xs=Xs, random_state=random_state)
        if return_y:
            ys = None
            out = (Xs, ys)
        else:
            out = Xs
        return out
