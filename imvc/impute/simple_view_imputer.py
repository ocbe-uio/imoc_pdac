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
        - 'mean': replace missing samples with the mean of each feature in the corresponding view
        - 'zeros': replace missing samples with zeros

    Attributes
    ----------
    features_view_mean_list_ : array-like of shape (n_views,)
        The mean value of each feature in the corresponding view, if value='mean'
    Examples
    --------
    >>> from imvc.utils import DatasetUtils
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.impute import SimpleViewImputer
    >>> from imvc.ampute import Amputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="simulated_gm")
    >>> amp = Amputer(p=0.3, mechanism="EDM")
    >>> Xs = amp.fit_transform(Xs)
    >>> transformer = FillIncompleteSamples(value = 'mean')
    >>> transformer.fit_transform(Xs)
    """


    def __init__(self, value : str = 'mean'):

        values = ['mean', 'zeros']
        if value not in values:
            raise ValueError(f"Invalid value. Expected one of: {values}")
        self.value = value


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
        self :  returns and instance of self.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        if self.value == "mean":
            self.features_view_mean_list_ = [X.mean() for X in Xs]
        elif self.value == "zeros":
            pass
        return self


    def transform(self, Xs):
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
        missing_views = get_observed_view_indicator(Xs = Xs)
        n_samples = len(missing_views)
        pandas_format = isinstance(Xs[0], pd.DataFrame)

        transformed_Xs = []
        for X_idx, X in enumerate(Xs):
            n_features = X.shape[1]
            if self.value == "mean":
                features_view_mean = self.features_view_mean_list_[X_idx]
                transformed_X = np.tile(features_view_mean, (n_samples ,1))
            elif self.value == "zeros":
                transformed_X = np.zeros((n_samples, n_features))
            transformed_X = pd.DataFrame(transformed_X, index= missing_views.index, columns = X.columns)
            transformed_X[missing_views.loc[:, X_idx]] = X
            transformed_X = transformed_X.astype(X.dtypes.to_dict())
            if not pandas_format:
                transformed_X = transformed_X.values
            transformed_Xs.append(transformed_X)
        return transformed_Xs


def simple_view_imputer(Xs, y = None, value : str = 'mean'):
    r"""
    Return a binary indicator for observed views.

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
        - 'mean': replace missing samples with the mean of each feature in the corresponding view
        - 'zeros': replace missing samples with zeros

    Returns
    -------
    transformed_X : array-likes, shape (n_samples, n_views)
        The transformed data.
    """
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    missing_views = get_observed_view_indicator(Xs=Xs)
    n_samples = len(missing_views)
    pandas_format = isinstance(Xs[0], pd.DataFrame)
    if not pandas_format:
        missing_views = pd.DataFrame(missing_views)

    transformed_Xs = []
    for X_idx, X in enumerate(Xs):
        n_features = X.shape[1]
        if value == "mean":
            features_view_mean = [X.mean() for X in Xs]
            transformed_X = np.tile(features_view_mean, (n_samples, 1))
        elif value == "zeros":
            transformed_X = np.zeros((n_samples, n_features))
        else:
            raise ValueError(f"Invalid value. Expected one of: ['mean', 'zeros']")

        dataframe_X = pd.DataFrame(X)
        transformed_X = pd.DataFrame(transformed_X, index=missing_views.index, columns=dataframe_X.columns)
        transformed_X[missing_views.loc[:, X_idx]] = X[missing_views.loc[:, X_idx]]
        transformed_X = transformed_X.astype(dataframe_X.dtypes.to_dict())
        if not pandas_format:
            transformed_X = transformed_X.values
        transformed_Xs.append(transformed_X)
    return transformed_Xs
