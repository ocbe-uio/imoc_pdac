from typing import Union

import numpy as np
import pandas as pd
import snf
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.cluster import SpectralClustering
from sklearn.manifold import spectral_embedding

from ..impute import get_observed_view_indicator
from ..utils import check_Xs


class NEMO(BaseEstimator, ClassifierMixin):
    r"""
    NEighborhood based Multi-Omics clustering (NEMO).

    NEMO is a method used for clustering data from multiple views sources. This algorithm operates
    through three main stages. Initially, it constructs a similarity matrix for each modality that represents the
    similarities between different samples. Then, it merges these individual view matrices into a unified one,
    combining the information from all views. Finally, the algorithm performs the actual clustering process on this
    integrated network, grouping similar samples together based on their multi-views data patterns.

    Parameters
    ----------
    n_clusters : int or list-of-int
        The number of clusters to generate. If it is a list, the number of clusters will be estimated by the algorithm
        with this range of number of clusters to choose between.
    num_neighbors : list or int, default=None
        The number of neighbors to use for each view. It can either be a number, a list of numbers or None. If it is a
        number, this is the number of neighbors used for all views. If this is a list, the number of neighbors are
        taken for each view from that list. If it is None, each view chooses the number of neighbors to be the number
        of samples divided by num_neighbors_ratio.
    num_neighbors_ratio : int, default=6
        The number of clusters to generate. If it is not provided, it will be estimated by the algorithm.
    metric : str or list-of-str, default="sqeuclidean"
        Distance metric to compute. Must be one of available metrics in :py:func`scipy.spatial.distance.pdist`. If
        multiple arrays a provided an equal number of metrics may be supplied.
    random_state : int, default=None
        Determines the randomness. Use an int to make the randomness deterministic.
    engine : str, default='python'
        Engine to use for computing the model. Must be one of ["python", "r"].
    verbose : bool, default=False
        Verbosity mode.

    Attributes
    ----------
    labels_ : array-like of shape (n_samples,)
        Labels of each point in training data.
    embedding_ : array-like of shape (n_samples, n_clusters)
        The final representation of the data to be used as input for the clustering step.
    n_clusters_ : int
        Final number of clusters.
    num_neighbors_ : int
        Final number of neighbors.
    affinity_matrix_ : np.array(n_samples, n_samples)
        Affinity matrix.

    References
    ----------
    .. [#nemopaper] Rappoport Nimrod, Shamir Ron. NEMO: Cancer subtyping by integration of partial multi-omic data.
                    Bioinformatics. 2019;35(18):3348–3356. doi: 10.1093/bioinformatics/btz058.
    .. [#nemocode] https://github.com/Shamir-Lab/NEMO

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.cluster import NEMO
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> estimator = NEMO(n_clusters = 2)
    >>> labels = estimator.fit_predict(Xs)
    """

    def __init__(self, n_clusters: Union[int,list] = 8, num_neighbors = None, num_neighbors_ratio: int = 6, metric='sqeuclidean',
                 random_state:int = None, engine: str = "python", verbose = False):
        self.n_clusters = n_clusters
        self.num_neighbors = num_neighbors
        self.num_neighbors_ratio = num_neighbors_ratio
        self.metric = metric
        self.random_state = random_state
        self._engines_options = ["python", "r"]
        if engine not in self._engines_options:
            raise ValueError(f"Invalid engine. Expected one of {self._engines_options}.")
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

        if self.engine == 'python':
            if not isinstance(Xs[0], pd.DataFrame):
                Xs = [pd.DataFrame(X) for X in Xs]
            observed_view_indicator = get_observed_view_indicator(Xs)
            samples = observed_view_indicator.index

            if self.num_neighbors is None:
                self.num_neighbors_ = [round(len(X)/self.num_neighbors_ratio) for X in Xs]
            elif not isinstance(self.num_neighbors, list):
                self.num_neighbors_ = [self.num_neighbors]*len(Xs)
            else:
                self.num_neighbors_ = self.num_neighbors

            affinity_matrix = pd.DataFrame(np.zeros((len(samples), len(samples))), columns = samples, index = samples)
            for X, neigh, view_idx in zip(Xs, self.num_neighbors_, range(len(Xs))):
                X = X.loc[observed_view_indicator[view_idx]]
                sim_data = pd.DataFrame(snf.make_affinity(X, metric = self.metric, K=neigh, normalize=False),
                                            index= X.index, columns= X.index)
                sim_data = sim_data.mask(sim_data.rank(axis=1, method='min', ascending=False) > neigh, 0)
                row_sum = sim_data.sum(1)

                sim_data /= row_sum
                sim_data += sim_data.T
                affinity_matrix.loc[sim_data.index, sim_data.columns] += sim_data

            affinity_matrix /= observed_view_indicator.sum(1)

            self.n_clusters_ = self.n_clusters if isinstance(self.n_clusters, int) else \
                snf.get_n_clusters(arr= affinity_matrix.values, n_clusters= self.n_clusters)[0]

            model = SpectralClustering(n_clusters= self.n_clusters_, random_state= self.random_state)
            labels = model.fit_predict(X= affinity_matrix)
            transformed_Xs = spectral_embedding(model.affinity_matrix_, n_components=self.n_clusters_,
                                                eigen_solver=model.eigen_solver, random_state=self.random_state,
                                                eigen_tol=model.eigen_tol, drop_first=False)


        elif self.engine == "R":
            from rpy2.robjects.packages import importr
            from ..utils import _convert_df_to_r_object
            nemo = importr("nemo")
            transformed_Xs = _convert_df_to_r_object(Xs)

            affinity_matrix = nemo.nemo.affinity.graph(transformed_Xs, k=self.num_neighbors)
            if (self.n_clusters is None):
                self.n_clusters = nemo.nemo.num.clusters(affinity_matrix)

            snftool = importr("SNFtool")
            labels = snftool.spectralClustering(affinity_matrix, self.n_clusters_)

        else:
            raise ValueError(f"Invalid engine. Expected one of {self._engines_options}.")

        self.labels_ = labels
        self.embedding_ = transformed_Xs
        self.affinity_matrix_ = affinity_matrix

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