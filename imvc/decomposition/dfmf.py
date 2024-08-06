from typing import Union

import pandas as pd
from sklearn.base import TransformerMixin, BaseEstimator
from sklearn.utils.validation import _generate_get_feature_names_out

from ..utils import check_Xs
from ._skfusion import fusion


class DFMF(TransformerMixin, BaseEstimator):
    r"""
    Data Fusion by Matrix Factorization (DFMF).

    DMFM is a data fusion approach with penalized matrix tri-factorization (DFMF) that simultaneously factorizes
    data matrices to reveal hidden associations. This method can deal with both view- and single-wise missing.

    Parameters
    ----------
    n_components : int, default=10
        Number of components to keep.
    max_iter : int, default=100
        Maximum number of iterations to perform.
    init_type : str or list of str, default='random_c'
        The algorithm to initialize latent matrix factors. Options are 'random', 'random_c' and 'random_vcol'. It can be
        a list, each item being for fit and transform, respectively.
    n_run: int, default=1
        Number of components to keep.
    stopping : tuple (target_matrix, eps), default=None
        Terminate iteration if the reconstruction error of target matrix improves by less than eps.
    stopping_system : float, default=None
        Terminate iteration if the reconstruction error of the fused system improves by less than eps. compute_err is
        to True to compute the error of the fused system.
    compute_err : bool, default=False
        Compute the reconstruction error of every relation matrix if True.
    callback : callable, default=None
        An optional user-supplied function to call after each iteration. Called as callback(G, S, cur_iter), where
        S and G are current latent estimates.
    fill_value : float, default=0
        Value to use to initially fill missing values.
    random_state : int, default=None
        Determines the randomness. Use an int to make the randomness deterministic.
    verbose : bool, default=False
        Verbosity mode.
    n_jobs : int, default=None
        Number of jobs to run in parallel. None means 1 unless in a joblib.parallel_backend context. -1 means
        using all processors.

    Attributes
    ----------
    fuser_ : Dfmf object
        Model.
    transformer_ : DfmfTransform object
        Object for transforming unseen data.
    t_: fusion.ObjectType
    ts_: list of fusion.ObjectType

    References
    ----------
    .. [#dfmfpaper] M. Žitnik and B. Zupan, "Data Fusion by Matrix Factorization," in IEEE Transactions on Pattern
                    Analysis and Machine Intelligence, vol. 37, no. 1, pp. 41-53, 1 Jan. 2015,
                    doi: 10.1109/TPAMI.2014.2343973.
    .. [#dfmfcode] https://github.com/mims-harvard/scikit-fusion/tree/master

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.decomposition import DFMF
    >>> from imvc.preprocessing import MultiViewTransformer
    >>> from sklearn.pipeline import make_pipeline
    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.cluster import KMeans
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> transformer = DFMF(n_components = 5).set_output(transform="pandas")
    >>> estimator = KMeans(n_clusters = 3)
    >>> pipeline = make_pipeline(MultiViewTransformer(StandardScaler().set_output(transform="pandas")), transformer, StandardScaler(), estimator)
    >>> labels = pipeline.fit_predict(Xs)
    """


    def __init__(self, n_components : int = 10, max_iter: int = 100, init_type: Union[str, list] = 'random_c', n_run: int = 1,
                 stopping=None, stopping_system=None, verbose=0, compute_err=False, callback=None,
                 random_state: int = None, n_jobs=1, fill_value=0):

        if isinstance(init_type, str):
            init_type = [init_type]
        if len(init_type) == 1:
            init_type *= 2

        self.n_components = n_components
        self.callback = callback
        self.max_iter = max_iter
        self.init_type = init_type
        self.n_run = n_run
        self.stopping = stopping
        self.stopping_system = stopping_system
        self.verbose = verbose
        self.compute_err = compute_err
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.fill_value = fill_value

        self.fuser_ = fusion.Dfmf(max_iter=max_iter, init_type=init_type[0], n_run=n_run, stopping=stopping,
                                  stopping_system=stopping_system, verbose=verbose, compute_err=compute_err,
                                  callback=callback, random_state=random_state, n_jobs=n_jobs)
        self.transformer_ = fusion.DfmfTransform(max_iter=max_iter, init_type=init_type[1], n_run=n_run,
                                                 stopping=stopping, stopping_system=stopping_system, verbose=verbose,
                                                 compute_err=compute_err, callback=callback, random_state=random_state,
                                                 n_jobs=n_jobs, fill_value=fill_value)
        self.t_ = fusion.ObjectType('Type 0', n_components)
        self.transform_ = None


    def fit(self, Xs, y = None):
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
        self :  returns and instance of self.
        """
        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        self.ts_ = [fusion.ObjectType(f'Type {i + 1}', self.n_components) for i in range(len(Xs))]
        relations = [fusion.Relation(X.values, self.t_, t) for X,t in zip(Xs, self.ts_)]
        fusion_graph = fusion.FusionGraph(relations)
        self.fuser_.fuse(fusion_graph)
        return self


    def transform(self, Xs):
        r"""
        Project data into the learned space.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs : list of array-likes, shape (n_samples, n_components)
            The projected data.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        relations = [fusion.Relation(X.values, self.t_, t) for X,t in zip(Xs, self.ts_)]
        fusion_graph = fusion.FusionGraph(relations)
        transformed_X = self.transformer_.transform(self.t_, fusion_graph, self.fuser_).factor(self.t_)
        if self.transform_ == "pandas":
            transformed_X = pd.DataFrame(transformed_X, index= Xs[0].index)
        return transformed_X


    def set_output(self, *, transform=None):
        r"""
        Set output container.

        Parameters
        ----------
        transform : str
            Only 'pandas' is currently supported.

        Returns
        -------
        self:  returns and instance of self.
        """
        self.transform_ = transform
        return self


    def get_feature_names_out(self, input_features=None):
        r"""
        Get output feature names for transformation. The feature names out will prefixed by the lowercased class name.
        For example, if the transformer outputs 3 features, then the feature names out are: ["class_name0",
        "class_name1", "class_name2"].

        Parameters
        ----------
        input_features: array-like of str or None, default=None.
            Only used to validate feature names with the names seen in fit.

        Returns
        -------
        Returns:
        feature_names_out: list of str objects
            Transformed feature names.
        """
        return _generate_get_feature_names_out(self, self.n_components, input_features=input_features)

