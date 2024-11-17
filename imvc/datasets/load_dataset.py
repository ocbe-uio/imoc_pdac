import os
import pandas as pd
import json
from settings import DATA_FOLDER


class LoadDataset:

    def load_tcga_pdac_subset(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from a subset of pancreatic ductal adenocarcinoma (TCGA_PDAC_subset)
        patients, extracted from The Cancer Genome Atlas (TCGA) website. There are six types of omics data available:
        copy number variation (CNA), DNA methylation, gene mutations, RNAseq, proteomics and miRNA.

        Samples: 89; Views: 6; Features: [745, 2185, 71, 1419, 192, 385]

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_tcga_pdac_subset()
        """
        output = LoadDataset.load_dataset(dataset_name = "TCGA_PDAC_subset", return_y = return_y, return_metadata = return_metadata)
        return output

    def load_tcga_pdac_subset_noise(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from a subset of pancreatic ductal adenocarcinoma (TCGA_PDAC_subset)
        patients, extracted from The Cancer Genome Atlas (TCGA) website. There are six types of omics data available:
        copy number variation (CNA), DNA methylation, gene mutations, RNAseq, proteomics and miRNA. This dataset has
        CNA and gene mutation datasets with noise signal (originally they are integer values).

        Samples: 89; Views: 6; Features: [745, 2185, 71, 1419, 192, 385]

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_tcga_pdac_subset_noise()
        """
        output = LoadDataset.load_dataset(dataset_name = "TCGA_PDAC_subset_noise", return_y = return_y, return_metadata = return_metadata)
        return output


    def load_tcga_pdac_all(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from a complete set of pancreatic ductal adenocarcinoma (TCGA_PDAC_all)
        patients, extracted from The Cancer Genome Atlas (TCGA) website. There are six types of omics data available:
        copy number variation (CNA), DNA methylation, gene mutations, RNAseq, proteomics and miRNA.

        Samples: 154; Views: 6; Features: [745, 2185, 71, 1419, 192, 385]

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_tcga_pdac_all()
        """
        output = LoadDataset.load_dataset(dataset_name = "TCGA_PDAC_all", return_y = return_y, return_metadata = return_metadata)
        return output


    def load_mt_subset(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from a complete set of pancreatic ductal adenocarcinoma (TCGA_PDAC_all)
        patients, extracted from The Cancer Genome Atlas (TCGA) website. There are six types of omics data available:
        copy number variation (CNA), DNA methylation, gene mutations, RNAseq, proteomics and miRNA.

        Samples: 154; Views: 6; Features: [745, 2185, 71, 1419, 192, 385]

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_mt_subset()
        """
        output = LoadDataset.load_dataset(dataset_name = "MT_subset", return_y = return_y, return_metadata = return_metadata)
        return output

    def load_mt_all(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from a complete set of pancreatic ductal adenocarcinoma (TCGA_PDAC_all)
        patients, extracted from The Cancer Genome Atlas (TCGA) website. There are six types of omics data available:
        copy number variation (CNA), DNA methylation, gene mutations, RNAseq, proteomics and miRNA.

        Samples: 154; Views: 6; Features: [745, 2185, 71, 1419, 192, 385]

        Parameters
        ----------
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

        Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_mt_all()
        """
        output = LoadDataset.load_dataset(dataset_name = "MT_all", return_y = return_y, return_metadata = return_metadata)
        return output



    @staticmethod
    def load_dataset(dataset_name: str, return_y: bool = False, return_metadata: bool = False):
        r"""
        Load a multi-view dataset.

        Parameters
        ----------
        dataset_name: str
            Name of the dataset. It must be one of: "TCGA_PDAC_subset", "TCGA_PDAC_subset_noise", "TCGA_PDAC_all".
        return_y: bool, default=False
            If True, return the label too.
        return_metadata: bool, default=False
            If True, return the metadata.

        Returns
        -------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : optional
            Array with labels
        metadata : optional
            Dict with info about the dataset (data modality names, labels, etc.).

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name = "TCGA_PDAC_subset")
        """
        data_path = os.path.join(DATA_FOLDER, dataset_name)
        data_files = [filename for filename in os.listdir(data_path)]
        data_files = sorted(data_files)
        data_files = [os.path.join(data_path, filename) for filename in data_files if dataset_name in filename]
        Xs = [pd.read_csv(filename, index_col=0) for filename in data_files]    # check that files are in correct format
        output = (Xs,)
        if return_y:
            y = pd.read_csv(os.path.join(data_path, f"{dataset_name}_y.csv"))
            y = y.loc[Xs[0].index]
            if y.shape[1] > 1:
                y = y.squeeze()
            output = output + (y,)
        if return_metadata:
            metadata_filename = os.path.join(data_path, "metadata.json")
            if os.path.isfile(metadata_filename):
                with open(metadata_filename) as json_file:
                    metadata = json.load(json_file)
                output = output + (metadata,)
            else:
                output = output + (None,)
        if len(output) == 1:
            output = output[0]
        return output
