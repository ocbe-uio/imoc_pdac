import os
from os.path import dirname

import numpy as np
import oct2py
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cluster import KMeans

from ..impute import get_observed_view_indicator, simple_view_imputer
from ..utils import check_Xs


class DAIMC(BaseEstimator, ClassifierMixin):
    r"""
    Doubly Aligned Incomplete Multi-view Clustering (DAIMC).

    The DAIMC algorithm integrates weighted semi-nonnegative matrix factorization (semi-NMF) to address incomplete
    multi-view clustering challenges. It leverages instance alignment information to learn a unified latent feature
    matrix across views and employs L2,1-Norm regularized regression to establish a consensus basis matrix, minimizing
    the impact of missing instances.

    It is recommended to normalize (Normalizer or NormalizerNaN in case incomplete views) the data before applying
    this algorithm.

    octave-control and octave-statistics should be installed. You can install them with
    'sudo apt install octave-control' and 'sudo apt install octave-statistics'.

    Parameters
    ----------
    n_clusters : int, default=8
        The number of clusters to generate.
    alpha : float, default=1
        nonnegative.
    beta : float, default=1
        Define the trade-off between sparsity and accuracy of regression for the i-th view.
    random_state : int, default=None
        Determines the randomness. Use an int to make the randomness deterministic.
    engine : str, default=matlab
        Engine to use for computing the model. Current options are 'matlab'. If engine == 'matlab',
        packages 'statistics' and 'control' should be installed in Octave. In linux, you can run: sudo apt-get install
        octave-statistics; sudo apt-get install octave-control.
.   verbose : bool, default=False
        Verbosity mode.

    Attributes
    ----------
    labels_ : array-like of shape (n_samples,)
        Labels of each point in training data.
    U_ : np.array
        Basis matrix.
    V_ : np.array
        Commont latent feature matrix.
    B_ : np.array
        Regression coefficient matrices.

    References
    ----------
    [paper1] Menglei Hu and Songcan Chen. 2018. Doubly aligned incomplete multi-view clustering. In Proceedings of the
            27th International Joint Conference on Artificial Intelligence (IJCAI'18). AAAI Press, 2262–2268.
    [paper2] Jie Wen, Zheng Zhang, Lunke Fei, Bob Zhang, Yong Xu, Zhao Zhang, Jinxing Li, A Survey on Incomplete
             Multi-view Clustering, IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS: SYSTEMS, 2022.
    [code]  https://github.com/DarrenZZhang/Survey_IMC

    Examples
    --------
    >>> from sklearn.pipeline import make_pipeline
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.cluster import DAIMC
    >>> from imvc.preprocessing import NormalizerNaN, MultiViewTransformer
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> normalizer = NormalizerNaN().set_output(transform="pandas")
    >>> estimator = DAIMC(n_clusters = 2)
    >>> pipeline = make_pipeline(MultiViewTransformer(normalizer), estimator)
    >>> labels = pipeline.fit_predict(Xs)
    """

    def __init__(self, n_clusters: int = 8, alpha: float = 1, beta: float = 1, random_state:int = None,
                 engine: str ="matlab", verbose = False):
        self.n_clusters = n_clusters
        self.alpha = alpha
        self.beta = beta
        self.random_state = random_state
        self.engine = engine
        self.verbose = verbose


    def fit(self, Xs, y=None):
        r"""
        Fit the transformer to the input data.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        self :  Fitted estimator.
        """
        Xs = check_Xs(Xs, force_all_finite='allow-nan')

        if self.engine=="matlab":
            matlab_folder = dirname(__file__)
            matlab_folder = os.path.join(matlab_folder, "_daimc")
            matlab_files = ["newinit.m", "litekmeans.m", "DAIMC.m", "UpdateV_DAIMC.m"]
            oc = oct2py.Oct2Py(temp_dir= matlab_folder)
            for matlab_file in matlab_files:
                with open(os.path.join(matlab_folder, matlab_file)) as f:
                    oc.eval(f.read())
            oc.eval("pkg load statistics")
            oc.eval("pkg load control")
            oc.warning("off", "Octave:possible-matlab-short-circuit-operator")

            if isinstance(Xs[0], pd.DataFrame):
                transformed_Xs = [X.values for X in Xs]
            elif isinstance(Xs[0], np.ndarray):
                transformed_Xs = Xs
            observed_view_indicator = get_observed_view_indicator(transformed_Xs)
            transformed_Xs = simple_view_imputer(transformed_Xs, value="zeros")
            transformed_Xs = [X.T for X in transformed_Xs]
            transformed_Xs = tuple(transformed_Xs)

            w = tuple([oc.diag(missing_view) for missing_view in observed_view_indicator.T])
            if self.random_state is not None:
                oc.rand('seed', self.random_state)
            u_0, v_0, b_0 = oc.newinit(transformed_Xs, w, self.n_clusters, len(transformed_Xs), nout=3)
            u, v, b, f, p, n = oc.DAIMC(transformed_Xs, w, u_0, v_0, b_0, None, self.n_clusters,
                                        len(transformed_Xs), {"afa": self.alpha, "beta": self.beta}, nout=6)
        else:
            raise ValueError("Only engine=='matlab' is currently supported.")

        model = KMeans(n_clusters= self.n_clusters, random_state= self.random_state)
        self.labels_ = model.fit_predict(X= v)
        self.U_ = u
        self.V_ = v
        self.B_ = b

        return self

    def _predict(self, Xs):
        r"""
        Return clustering results for samples.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        return self.labels_


    def fit_predict(self, Xs, y=None):
        r"""
        Fit the model and return clustering results.
        Convenience method; equivalent to calling fit(X) followed by predict(X).

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """

        labels = self.fit(Xs)._predict(Xs)
        return labels

