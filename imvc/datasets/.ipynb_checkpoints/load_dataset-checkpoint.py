import os
from os.path import dirname

import pandas as pd
from imvc.utils import DatasetUtils


class LoadDataset:

    @staticmethod
    def load_incomplete_digits(p, return_y: bool = False, shuffle: bool = True, assess_percentage: bool = True,
                               random_state: int = None):
        r"""
        Load the UCI multiple features dataset, taken from the UCI Machine Learning Repository
        at https://archive.ics.uci.edu/ml/datasets/Multiple+Features. This data set consists of 6 views of
        handwritten digit images, with classes 0-9..

        Parameters
        ----------
        p: list or float
            The percentaje that each view will have for missing samples. If p is float, all the views will have the
            same percentaje.
        return_y: bool, default False
            If True, return the label too.
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
        y : optional
            Array with labels

        Notes
        -----
        This data consists of six views:
        - 6 Fourier coefficients of the character shapes
        - 216 profile correlations
        - 64 Karhunen-Love coefficients
        - 240 pixel averages of the images from 2x3 windows
        - 47 Zernike moments
        - 6 morphological features

        References
        ----------
        [#1Data] M. van Breukelen, et al. "Handwritten digit recognition by combined classifiers", Kybernetika, 34(4):381-386, 1998
        [#2UCI] Dheeru Dua and Casey Graff. UCI machine learning repository, 2017. URL http://archive.ics.uci.edu/ml.
        [url] Adapted function from mvlearn Package: Perry, Ronan, et al. "mvlearn: Multiview Machine Learning in
            Python." Journal of Machine Learning Research 22.109 (2021): 1-7.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_digits(p = 0.2)
        """
        filenames = ["mfeat-fou.csv", "mfeat-fac.csv", "mfeat-kar.csv",
                     "mfeat-pix.csv", "mfeat-zer.csv", "mfeat-mor.csv"]
        module_path = dirname(__file__)
        Xs = [pd.read_csv(os.path.join(module_path, "data", "digits", filename)) for filename in filenames]
        Xs = DatasetUtils.ampute(Xs=Xs, p=p, assess_percentage=assess_percentage, random_state=random_state)
        if shuffle:
            Xs = DatasetUtils.shuffle_imvd(Xs=Xs, random_state=random_state)
        if return_y:
            y = [X.iloc[:, -1] for X in Xs]
            Xs = [X.iloc[:, :-1] for X in Xs]
            out = (Xs, y)
        else:
            Xs = [X.iloc[:, :-1] for X in Xs]
            out = Xs
        return out

    @staticmethod
    def load_incomplete_nutrimouse(p, return_y: bool = False, shuffle: bool = True, assess_percentage: bool = True,
                                   random_state: int = None):
        r"""
        Load an incomplete multi-view version of the Nutrimouse dataset, a two-view dataset from a nutrition
        study on mice.

        Parameters
        ----------
        p: list or float
            The percentaje that each view will have for missing samples. If p is float, all the views will have the
            same percentaje.
        return_y: bool, default False
            If True, return the label too.
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
        This data consists of two views from a nutrition study of 40 mice:
        - gene : expressions of 120 potentially relevant genes
        - lipid : concentrations of 21 hepatic fatty acids

        References
        ----------
        [paper] P. Martin, H. Guillou, F. Lasserre, S. Déjean, A. Lan, J-M.
                Pascussi, M. San Cristobal, P. Legrand, P. Besse, T. Pineau.
                "Novel aspects of PPARalpha-mediated regulation of lipid and
                xenobiotic metabolism revealed through a nutrigenomic study."
                Hepatology, 2007.
        [paper] González I., Déjean S., Martin P.G.P and Baccini, A. (2008) CCA:
                "An R Package to Extend Canonical Correlation Analysis." Journal
                of Statistical Software, 23(12).
        [url] Adapted function from mvlearn Package: Perry, Ronan, et al. "mvlearn: Multiview Machine Learning in
            Python." Journal of Machine Learning Research 22.109 (2021): 1-7.

         Examples
        --------
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
        """
        module_path = dirname(__file__)
        Xs = [pd.read_csv(os.path.join(module_path, "data", "nutrimouse", filename)) for filename in ["gene.csv", "lipid.csv"]]
        ys = [pd.read_csv(os.path.join(module_path, "data", "nutrimouse", filename)) for filename in ["genotype.csv", "diet.csv"]]
        Xs = DatasetUtils.ampute(Xs=Xs, p=p, assess_percentage=assess_percentage, random_state=random_state)
        if shuffle:
            Xs = DatasetUtils.shuffle_imvd(Xs=Xs, random_state=random_state)
        if return_y:
            out = (Xs, ys)
        else:
            out = Xs
        return out
