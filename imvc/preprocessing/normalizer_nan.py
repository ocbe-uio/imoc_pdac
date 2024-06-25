import numpy as np
from sklearn.preprocessing import Normalizer
from sklearn.utils import check_array


class NormalizerNaN(Normalizer):
    r"""
    Similar to sklearn.preprocessing.Normalizer but handles NaN values.

    Parameters
    ----------
    norm : {‘l1’, ‘l2’, ‘max’}, default=’l2’
        The norm to use to normalize each non zero sample. If norm=’max’ is used, values will be rescaled by
        the maximum of the absolute values.

    Attributes
    ----------
    features_view_mean_list_ : array-like of shape (n_views,)
        The mean value of each feature in the corresponding view, if value='mean'
    Examples
    --------
    >>> from imvc.utils import DatasetUtils
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.ampute import Amputer
    >>> from imvc.preprocessing import NormalizerNaN, MultiViewTransformer
    >>> Xs = LoadDataset.load_dataset(dataset_name="simulated_gm")
    >>> amp = Amputer(p=0.3, mechanism="EDM")
    >>> Xs = amp.fit_transform(Xs)
    >>> transformer = MultiViewTransformer(NormalizerNaN())
    >>> transformer.fit_transform(Xs)
    """


    def __init__(self, norm : str = 'l2'):

        super().__init__(norm)

        values = ['l1', 'l2', 'max']
        if norm not in values:
            raise ValueError(f"Invalid value. Expected one of: {values}")
        self.norm = norm


    def fit(self, X, y=None):
        r"""
        Fit the transformer to the input data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and n_features is the number of features.
        y : Ignored
                Not used, present here for API consistency by convention.
        Returns
        -------
        self :  returns and instance of self.
        """

        self._validate_data(X, **{"force_all_finite":'allow-nan'})
        return self


    def transform(self, X, y=None):
        r"""
        Scale each non zero row of X to unit norm.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and n_features is the number of features.
        y : Ignored
                Not used, present here for API consistency by convention.

        Returns
        -------
        transformed_X : array-like of shape (n_samples, n_features)
            Transformed data.
        """

        X = check_array(X, force_all_finite='allow-nan')
        if self.norm == "l1":
            norms = np.abs(X).sum(axis=1)
        elif self.norm == "l2":
            norms = (X**2).sum(axis=1)**0.5
        elif self.norm == "max":
            norms = np.abs(X).sum(axis=1).max(axis=1)
        else:
            raise ValueError(f"Invalid value. Expected one of: ['l1', 'l2', 'max']")
        norms[norms == 0] = 1.
        transformed_X = X / norms[:,None]
        return transformed_X

