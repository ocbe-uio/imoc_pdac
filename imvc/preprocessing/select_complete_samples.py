from sklearn.preprocessing import FunctionTransformer

from ..impute import get_observed_view_indicator
from ..utils import check_Xs


class SelectCompleteSamples(FunctionTransformer):
    r"""
    Remove incomplete samples from a multi-view dataset.

    Parameters
    ----------
    -

    Attributes
    ----------
    -

    Examples
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import SelectCompleteSamples
    >>> from imvc.ampute import Amputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
    >>> transformer = SelectCompleteSamples()
    >>> transformer.fit_transform(Xs)
    """

    def __init__(self):
        super().__init__(select_complete_samples)


def select_complete_samples(Xs: list):
    r"""
    Remove incomplete samples from a multi-view dataset.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features)
        A list of different views.

    Returns
    -------
    transformed_Xs : list of array-likes, shape (n_samples, n_features_i)
        The transformed data.
    """

    Xs = check_Xs(Xs, force_all_finite='allow-nan')
    sample_views = get_observed_view_indicator(Xs)
    complete_samples = sample_views.all(axis= 1)
    transformed_Xs = [X.loc[complete_samples] for X in Xs]
    return transformed_Xs


