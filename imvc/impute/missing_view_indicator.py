import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer

from ..utils import check_Xs

class MissingViewIndicator(FunctionTransformer):
    r"""
    Binary indicators for missing views. Apply FunctionTransformer (from Scikit-learn) with get_missing_view_indicator
    as a function.

    Note that this component typically should not be used in a vanilla Pipeline consisting of preprocessing and
    an estimator.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.impute import MissingViewIndicator
    >>> from imvc.ampute import Amputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> Xs = Amputer(p= 0.2, random_state=42).fit_transform(Xs)
    >>> transformer = MissingViewIndicator()
    >>> X_tr = transformer.fit_transform(Xs)
    """

    def __init__(self):
        super().__init__(get_missing_view_indicator)


def get_missing_view_indicator(Xs : list, y = None):
    r"""
    Return a binary indicator for missing views.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features)
        A list of different views.
    y : Ignored
        Not used, present here for API consistency by convention.

    Returns
    -------
    transformed_X : array-likes, shape (n_samples, n_views)
        Binary indicator for missing views.
    """
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    transformed_X = np.vstack([np.isnan(X).all(1) for X in Xs]).T
    if isinstance(Xs[0], pd.DataFrame):
        transformed_X = pd.DataFrame(transformed_X, index=Xs[0].index)
    return transformed_X

