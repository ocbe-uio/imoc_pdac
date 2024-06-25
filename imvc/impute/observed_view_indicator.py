import numpy as np
import pandas as pd
from sklearn.preprocessing import FunctionTransformer

from ..utils import check_Xs


class ObservedViewIndicator(FunctionTransformer):
    r"""
    Binary indicators for observed views.

    Note that this component typically should not be used in a vanilla Pipeline consisting of preprocessing and
    an estimator.

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.impute import ObservedViewIndicator
    >>> from imvc.ampute import Amputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> Xs = Amputer(p= 0.2, random_state=42).fit_transform(Xs)
    >>> transformer = ObservedViewIndicator()
    >>> X_tr = transformer.fit_transform(Xs)
    """


    def __init__(self):
        super().__init__(get_observed_view_indicator)


def get_observed_view_indicator(Xs, y = None):
    r"""
    Return a binary indicator for observed views.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features)
        A list of different views.

    Returns
    -------
    transformed_X : array-likes, shape (n_samples, n_views)
        The transformed data.
    """
    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    transformed_X = np.vstack([np.isnan(X).all(1) for X in Xs]).T
    transformed_X = ~transformed_X
    if isinstance(Xs[0], pd.DataFrame):
        transformed_X = pd.DataFrame(transformed_X, index=Xs[0].index)
    return transformed_X

