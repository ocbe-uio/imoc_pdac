import os
import pandas as pd
import json
from settings import DATA_FOLDER


class LoadDataset:

    def load_pdac(return_y: bool = False, return_metadata: bool = False):
        r"""
        The dataset comprises multi-omic data from pancreatic ductal adenocarcinoma (PDAC) patients, extracted from The
        Cancer Genome Atlas (TCGA) website. There are six types of omics data available: copy number variation (CNA),
        DNA methylation, gene mutations, RNAseq, proteomics and miRNA.

        Samples: 90; Views: 6; Features: [770, 2187, 71, 1565, 192, 553]

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
        >>> Xs = LoadDataset.load_pdac()
        """
        output = LoadDataset.load_dataset(dataset_name = "PDAC", return_y = return_y, return_metadata = return_metadata)
        return output

    @staticmethod
    def load_dataset(dataset_name: str, return_y: bool = False, return_metadata: bool = False):
        r"""
        Load a multi-view dataset.

        Parameters
        ----------
        dataset_name: str
            Name of the dataset. It must be one of: "bbcsport", "bdgp", "buaa", "caltech101", "digits", "metabric",
            "nuswide", "nutrimouse", "simulated_gm", "simulated_InterSIM", "simulated_netMUG", "tcga", "PDAC".
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
        >>> Xs = LoadDataset.load_dataset(dataset_name = "PDAC")
        """
        data_path = os.path.join(DATA_FOLDER, dataset_name)
        data_files = [filename for filename in os.listdir(data_path)]
        data_files = sorted(data_files)
        data_files = [os.path.join(data_path, filename) for filename in data_files if dataset_name in filename and not filename.endswith("y.csv")]
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
