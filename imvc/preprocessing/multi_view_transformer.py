from copy import deepcopy
from sklearn.base import BaseEstimator, TransformerMixin

from ..utils import check_Xs


class MultiViewTransformer(BaseEstimator, TransformerMixin):
    r"""
    A transformer that applies the same transformation to multiple views of data.

    Parameters
    ----------
    transformer : scikit-learn transformer object or list of scikit-learn transformer object
        A scikit-learn transformer object that will be used to transform each view of data. If a list is provided,
        each transformer will be applied on each view, otherwise the same transformer will be applied on each view.

    Attributes
    ----------
    transformer_list_ : list of preprocessing (n_views,)
        A list of preprocessing, one for each view of data.
    same_transformer_ : boolean
        A booleaing indicating if the same transformer will be applied on each view of data.

    Example
    --------
    >>> from imvc.datasets import LoadDataset
    >>> from imvc.preprocessing import MultiViewTransformer
    >>> from sklearn.impute import SimpleImputer
    >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
    >>> transformer = MultiViewTransformer(transformer = SimpleImputer.set_output(transform = 'pandas'))
    >>> transformer.fit_transform(Xs)
    """


    def __init__(self, transformer):
        self.same_transformer_ = False if isinstance(transformer, list) else True
        if self.same_transformer_:
            transformer_object = deepcopy(transformer)
            try:
                assert hasattr(transformer_object, "fit") and callable(getattr(transformer_object, "fit"))
            except AssertionError:
                raise ValueError("transformer must be a scikit-learn transformer like object")
        else:
            for transformer_object in transformer:
                try:
                    assert hasattr(transformer_object, "fit") and callable(getattr(transformer_object, "fit"))
                except AssertionError:
                    raise ValueError("transformer must be a scikit-learn transformer like object")

        self.transformer = transformer
        self.transformer_list_ = [] if self.same_transformer_ else transformer


    def fit(self, Xs, y = None):
        r"""
        Fit the transformer to the input data.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        y : array-like, shape (n_samples,)
            Labels for each sample. Only used by supervised algorithms.

        Returns
        -------
        self :  returns an instance of self.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        for X_idx,X in enumerate(Xs):
            if self.same_transformer_:
                self.transformer_list_.append(deepcopy(self.transformer))
            self.transformer_list_[X_idx].fit(X, y)
        return self


    def transform(self, Xs):
        r"""
        Transform the input data using the transformers.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs : list of array-likes, shape (n_samples, n_features_i)
            A list of transformed views of data, one for each input view.
        """

        Xs = check_Xs(Xs, force_all_finite='allow-nan')
        tranformed_Xs = [self.transformer_list_[X_idx].transform(X) for X_idx, X in enumerate(Xs)]
        return tranformed_Xs
