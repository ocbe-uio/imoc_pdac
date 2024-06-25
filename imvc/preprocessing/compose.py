import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import FunctionTransformer

from ..utils import check_Xs
from ..utils import DatasetUtils


class DropView(FunctionTransformer):
    r"""
    A transformer that drops a specified view from a multi-view dataset.

    Parameters
    ----------
    X_idx : int, default=0
        The index of the view to drop.

    Returns
    -------
    transformed_Xs : list of array-likes
        - Xs length: n_views - 1
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import DropView
    >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
    >>> transformer = DropView(X_idx = 1)
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self, X_idx: int = 0):
        self.X_idx = X_idx
        super().__init__(drop_view, kw_args={"X_idx": X_idx})


class ConcatenateViews(FunctionTransformer):
    r"""
    A transformer that concatenates all views from a multi-view dataset.

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features)

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import DropView
    >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
    >>> transformer = ConcatenateView()
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self):
        super().__init__(concatenate_views)


class SingleView(FunctionTransformer):
    r"""
    Transformer that selects a single view from multi-view data.

    Parameters
    ----------
    X_idx : int, default=0
        The index of the view to select from the input data.

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features[X_idx])

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import SingleView
    >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
    >>> transformer = SingleView(X_idx = 1)
    >>> transformer.fit_transform(Xs)

    """
    
    def __init__(self, X_idx : int = 0):
        self.X_idx = X_idx
        super().__init__(single_view, kw_args = {"X_idx": X_idx})


class AddMissingViews(TransformerMixin, BaseEstimator):
    r"""
    Transformer to add missing samples in each view, in a way that all the views will have the same samples.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
    samples : pd.Index or array-like  (n_samples,)
        list with all samples

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features)

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import AddMissingViews, MultiViewTransformer
    >>> from imvc.utils import DatasetUtils
    >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
    >>> samples = DatasetUtils.get_sample_names(Xs= Xs)
    >>> transformer = MultiViewTransformer(transformer = AddMissingViews(samples= samples))
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self, samples: pd.Index):
        self.samples = samples


    def fit(self, X, y=None):
        r"""
        Just for compatibility with other Scikit-learn preprocessing.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
        y : array-like, shape (n_samples,)
            Labels for each sample. Only used by supervised algorithms.

        Returns
        -------
        self :  returns and instance of self.
        """

        return self


    def transform(self, X):
        r"""
        Transform the input data into a non-negative matrix.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features_i)

        Returns
        -------
        transformed_X : array-like of shape (n_samples, n_features_i)
            The transformed data.
        """

        transformed_X = add_missing_views(X = X, samples= self.samples)
        return transformed_X


class SortData(FunctionTransformer):
    r"""
    Transformer that establish and assess the order of the incomplete multi-view dataset.

    Returns
    -------
    transformed_X : list of array-likes (n_samples, n_features_i)

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import SortData
    >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
    >>> transformer = SortData()
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self):
        super().__init__(sort_data)


def concatenate_views(Xs):
    r"""
    A function that merge features from a multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features)
    """

    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    transformed_X = pd.concat(Xs, axis= 1)
    return transformed_X


def drop_view(Xs, X_idx : int = 0):
    r"""
    A function that drops a specified view from a multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Returns
    -------
    transformed_Xs : array-like of shape (n_samples, n_features - n_features[X_idx])
    """

    Xs = check_Xs(Xs, allow_incomplete=True, force_all_finite='allow-nan')
    transformed_Xs = Xs[:X_idx] + Xs[X_idx+1 :]
    return transformed_Xs


def single_view(Xs, X_idx : int = 0):
    r"""
    A function that selects a specified view from a multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features[X_idx])
    """

    Xs = check_Xs(Xs, allow_incomplete=True, force_all_finite='allow-nan')
    transformed_X = Xs[X_idx]
    return transformed_X


def add_missing_views(X, samples):
    r"""
    Add missing samples in each view, in a way that all the views will have the same samples.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
    samples : array-like  (n_samples,)
        list with all samples

    Returns
    -------
    transformed_X : array-like of shape (n_samples, n_features)
    """

    transformed_X = X.T.copy()
    transformed_X[samples.difference(X.index)] = np.nan
    transformed_X = transformed_X.T
    transformed_X = transformed_X.loc[samples]
    return transformed_X


def sort_data(Xs):
    r"""
    A function that establish and assess the order of the incomplete multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Returns
    -------
    transformed_X : list of array-likes (n_samples, n_features_i)
    """

    Xs = check_Xs(Xs, allow_incomplete=True, force_all_finite='allow-nan')
    samples = DatasetUtils.get_sample_names(Xs=Xs)
    transformed_X = [X.loc[samples.intersection(X.index)] for X in Xs]
    return transformed_X

