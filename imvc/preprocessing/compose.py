import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer

from ..utils import check_Xs
from ..utils import DatasetUtils


class DropView(FunctionTransformer):
    r"""
    A transformer that drops a specified view from a multi-view dataset. Apply FunctionTransformer (from Scikit-learn)
    with drop_view as a function.

    Parameters
    ----------
    X_idx : int, default=0
        The index of the view to drop from the input data.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import DropView
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = DropView(X_idx = 1)
    >>> transformer.fit_transform(Xs)
    """

    def __init__(self, X_idx: int = 0):
        self.X_idx = X_idx
        super().__init__(drop_view, kw_args={"X_idx": X_idx})


class ConcatenateViews(FunctionTransformer):
    r"""
    A transformer that concatenates all views from a multi-view dataset. Apply FunctionTransformer (from Scikit-learn)
    with concatenate_views as a function.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import ConcatenateViews
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = ConcatenateViews()
    >>> transformer.fit_transform(Xs)
    """

    def __init__(self):
        super().__init__(concatenate_views)


class SingleView(FunctionTransformer):
    r"""
    Transformer that selects a single view from multi-view data. Apply FunctionTransformer (from Scikit-learn) with
    single_view as a function.

    Parameters
    ----------
    X_idx : int, default=0
        The index of the view to select from the input data.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import SingleView
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = SingleView(X_idx = 1)
    >>> transformer.fit_transform(Xs)
    """
    
    def __init__(self, X_idx : int = 0):
        self.X_idx = X_idx
        super().__init__(single_view, kw_args = {"X_idx": X_idx})


class AddMissingViews(FunctionTransformer):
    r"""
    Transformer to add missing samples in each view, in a way that all the views will have the same samples. Apply
    FunctionTransformer (from Scikit-learn) with add_missing_views as a function.

    This transformer is applied on individual views, so for applying in a multi-view dataset, we recommend to use it
    with MultiViewTransformer.

    Parameters
    ----------
    samples : array-like  (n_samples,)
        pd.Index with all samples

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import AddMissingViews, MultiViewTransformer
    >>> from imvc.utils import DatasetUtils
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> samples = DatasetUtils.get_sample_names(Xs= Xs)
    >>> transformer = MultiViewTransformer(transformer = AddMissingViews(samples= samples))
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self, samples: pd.Index):
        self.samples = samples
        super().__init__(add_missing_views, kw_args={"samples": samples})


class SortData(FunctionTransformer):
    r"""
    Transformer that establish and assess the order of the incomplete multi-view dataset. Apply
    FunctionTransformer (from Scikit-learn) with sort_data as a function.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import SortData
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = SortData()
    >>> transformer.fit_transform(Xs)

    """

    def __init__(self):
        super().__init__(sort_data)


def concatenate_views(Xs: list):
    r"""
    A function that concatenate all features from a multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.

    Returns
    -------
    transformed_Xs : array-like, shape (n_samples, n_features)
        The transformed dataset.
    """

    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    if isinstance(Xs[0], pd.DataFrame):
        transformed_X = pd.concat(Xs, axis= 1)
    elif isinstance(Xs[0], np.ndarray):
        transformed_X = np.concatenate(Xs, axis= 1)
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
    X_idx : int, default=0
        The index of the view to drop from the input data.

    Returns
    -------
    transformed_Xs : array-like, shape (n_samples, n_features)
        The transformed multi-view dataset.
    """
    if X_idx >= len(Xs):
        raise ValueError("X_idx out of range. Should be between 0 and n_views - 1")
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
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
    X_idx : int, default=0
        The index of the view to select from the input data.

    Returns
    -------
    transformed_Xs : array-like, shape (n_samples, n_features)
        The transformed dataset.
    """
    if X_idx >= len(Xs):
        raise ValueError("X_idx out of range. Should be between 0 and n_views - 1")
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
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
    pandas_format = isinstance(X, pd.DataFrame)
    if pandas_format:
        transformed_X = X.T.copy()
    else:
        X = pd.DataFrame(X)
        transformed_X = X.T.copy()
    transformed_X[samples.difference(X.index)] = np.nan
    transformed_X = transformed_X.T
    transformed_X = transformed_X.loc[samples]
    if not pandas_format:
        transformed_X = transformed_X.values
    return transformed_X


def sort_data(Xs: list):
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
        The transformed multi-view dataset.
    """

    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    if not isinstance(Xs[0], pd.DataFrame):
        Xs = [pd.DataFrame(X) for X in Xs]
    samples = DatasetUtils.get_sample_names(Xs=Xs)
    transformed_X = [X.loc[samples.intersection(X.index)] for X in Xs]
    return transformed_X

