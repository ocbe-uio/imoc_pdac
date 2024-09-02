import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin

from ..impute import get_observed_view_indicator
from ..utils import check_Xs


class SimpleViewImputer(BaseEstimator, TransformerMixin):
    r"""
    Fill incomplete samples of a dataset using a specified method.

    Parameters
    ----------
    value : str, optional (default='mean')
        The method to use for filling missing views. Possible values:
        - 'mean': replace missing samples with the mean of each feature in the corresponding view.
        - 'zeros': replace missing samples with zeros.

    Attributes
    ----------
    features_view_mean_list_ : array-like of shape (n_views,)
        The mean value of each feature in the corresponding view, if value='mean'
    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.impute import SimpleViewImputer
    >>> from imvc.ampute import Amputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="simulated_gm")
    >>> amp = Amputer(p=0.3, random_state=42)
    >>> Xs = amp.fit_transform(Xs)
    >>> transformer = SimpleViewImputer(value = 'mean')
    >>> transformer.fit_transform(Xs)
    """


    def __init__(self, value : str = 'mean'):

        values = ['mean', 'zeros']
        if value not in values:
            raise ValueError(f"Invalid value. Expected one of: {values}")
        self.value = value


    def fit(self, Xs : list, y=None):
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
        self :  returns an instance of self.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        if self.value == "mean":
            self.features_view_mean_list_ = [np.nanmean(X, axis=0) for X in Xs]
        elif self.value == "zeros":
            pass
        return self


    def transform(self, Xs : list):
        r"""
        Transform the input data by filling missing samples.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs : list of array-likes, shape (n_samples, n_features_i)
            The transformed data with filled missing samples.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        pandas_format = isinstance(Xs[0], pd.DataFrame)
        if pandas_format:
            samples = Xs[0].index
            features = [X.columns for X in Xs]
            dtypes = [X.dtypes.to_dict() for X in Xs]
            Xs = [X.values for X in Xs]
        observed_view_indicator = get_observed_view_indicator(Xs = Xs)
        n_samples = len(observed_view_indicator)

        transformed_Xs = []
        for X_idx, X in enumerate(Xs):
            n_features = X.shape[1]
            if self.value == "mean":
                features_view_mean = self.features_view_mean_list_[X_idx]
                transformed_X = np.tile(features_view_mean, (n_samples ,1))
            elif self.value == "zeros":
                transformed_X = np.zeros((n_samples, n_features))
            transformed_X[observed_view_indicator[:, X_idx]] = X[observed_view_indicator[:, X_idx]]
            if pandas_format:
                transformed_X = pd.DataFrame(transformed_X, index=samples, columns=features[X_idx])
                transformed_X = transformed_X.astype(dtypes[X_idx])
            transformed_Xs.append(transformed_X)
        return transformed_Xs


def simple_view_imputer(Xs : list, y = None, value : str = 'mean'):
    r"""
    Transform the input data by filling missing samples.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features)
        A list of different views.
    y : Ignored
            Not used, present here for API consistency by convention.
    value : str, optional (default='mean')
        The method to use for filling missing views. Possible values:
        - 'mean': replace missing samples with the mean of each feature in the corresponding view.
        - 'zeros': replace missing samples with zeros.

    Returns
    -------
    transformed_Xs : list of array-likes, shape (n_samples, n_features_i)
        The transformed data with filled missing samples.
    """
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    pandas_format = isinstance(Xs[0], pd.DataFrame)
    if pandas_format:
        samples = Xs[0].index
        features = [X.columns for X in Xs]
        dtypes = [X.dtypes.to_dict() for X in Xs]
        Xs = [X.values for X in Xs]
    observed_view_indicator = get_observed_view_indicator(Xs=Xs)
    n_samples = len(observed_view_indicator)

    transformed_Xs = []
    for X_idx, X in enumerate(Xs):
        n_features = X.shape[1]
        if value == "mean":
            features_view_mean = np.nanmean(X, axis=0)
            transformed_X = np.tile(features_view_mean, (n_samples, 1))
        elif value == "zeros":
            transformed_X = np.zeros((n_samples, n_features))
        else:
            raise ValueError(f"Invalid value. Expected one of: ['mean', 'zeros']")

        transformed_X[observed_view_indicator[:, X_idx]] = X[observed_view_indicator[:, X_idx]]
        if pandas_format:
            transformed_X = pd.DataFrame(transformed_X, index=samples, columns=features[X_idx])
            transformed_X = transformed_X.astype(dtypes[X_idx])
        transformed_Xs.append(transformed_X)
    return transformed_Xs
