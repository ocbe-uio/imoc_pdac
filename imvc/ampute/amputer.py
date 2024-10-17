import copy
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from ..utils import DatasetUtils


class Amputer(BaseEstimator, TransformerMixin):
    r"""
    Ampute a complete multi-view datasets.

    Parameters
    ----------
    p: float, default=0.1
        Percentaje of incomplete samples.
    mechanism: str, default="edm"
        One of ["edm", 'mcar', 'mnar', 'pm'].
    weights: list, default=None
        The probabilities associated with each number of missing modalities. If not given, the sample
        assumes a uniform distribution.
    random_state: int, default=None
        If int, random_state is the seed used by the random number generator.

    Example
    --------
    >>> from imvc.ampute import Amputer
    >>> from imvc.datasets import LoadDataset
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = Amputer(p= 0.2)
    >>> transformer.fit_transform(Xs)
    """

    def __init__(self, p:float = 0.1, mechanism: str = "pm", weights: list = None, random_state: int = None):

        mechanisms_options = ["edm", "mcar", "mnar", "pm"]
        if mechanism not in mechanisms_options:
            raise ValueError(f"Invalid mechanism. Expected one of: {mechanisms_options}")

        self.mechanism = mechanism
        self.p = p
        self.weights = weights
        self.random_state = random_state


    def fit(self, Xs: list, y=None):
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
        self.n_views = len(Xs)
        return self


    def transform(self, Xs: list):
        r"""
        Ampute a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs : list of array-likes, shape (n_samples, n_features)
            The amputed multi-view dataset.
        """
        if self.p > 0:
            pandas_format = isinstance(Xs[0], pd.DataFrame)
            if pandas_format:
                rownames = Xs[0].index
                colnames = [X.columns for X in Xs]
                Xs = [X.values for X in Xs]
            sample_names = pd.Index(list(range(len(Xs[0]))))

            if self.mechanism == "edm":
                pseudo_observed_view_indicator = self._edm_mask(sample_names=sample_names)
            elif self.mechanism == "mcar":
                pseudo_observed_view_indicator = self._mcar_mask(sample_names=sample_names)
            elif self.mechanism == "pm":
                pseudo_observed_view_indicator = self._pm_mask(sample_names=sample_names)
            elif self.mechanism == "mnar":
                pseudo_observed_view_indicator = self._mnar_mask(sample_names=sample_names)

            pseudo_observed_view_indicator = pseudo_observed_view_indicator.astype(bool)
            transformed_Xs = DatasetUtils.convert_to_imvd(Xs=Xs, observed_view_indicator=pseudo_observed_view_indicator)

            if pandas_format:
                transformed_Xs = [pd.DataFrame(X, index=rownames, columns=colnames[X_idx])
                                  for X_idx, X in enumerate(transformed_Xs)]
        else:
            transformed_Xs = copy.deepcopy(Xs)

        return transformed_Xs


    def _edm_mask(self, sample_names):
        pseudo_observed_view_indicator = pd.DataFrame(np.ones((len(sample_names), self.n_views)), index=sample_names)
        common_samples = pd.Series(sample_names, index=sample_names).sample(frac=1 - self.p, replace=False,
                                                                            random_state=self.random_state).index
        sampled_names = copy.deepcopy(common_samples)
        n_missing = int(len(sample_names.difference(sampled_names)) / self.n_views)
        for X_idx in range(self.n_views):
            x_per_view = sample_names.difference(sampled_names)
            if X_idx != self.n_views - 1:
                x_per_view = pd.Series(x_per_view, index=x_per_view).sample(n=n_missing,
                                                                            replace=False,
                                                                            random_state=self.random_state).index
            sampled_names = sampled_names.append(x_per_view)
            idxs_to_remove = common_samples.append(x_per_view)
            idxs_to_remove = sample_names.difference(idxs_to_remove)
            pseudo_observed_view_indicator.loc[idxs_to_remove, X_idx] = 0
        return pseudo_observed_view_indicator


    def _mcar_mask(self, sample_names):
        pseudo_observed_view_indicator = pd.DataFrame(np.ones((len(sample_names), self.n_views)), index=sample_names)
        common_samples = pd.Series(sample_names, index=sample_names).sample(frac=1 - self.p, replace=False,
                                                                            random_state=self.random_state).index
        idxs_to_remove = sample_names.difference(common_samples)
        shape = pseudo_observed_view_indicator.loc[idxs_to_remove].shape
        mask = np.random.default_rng(self.random_state).choice(2, size=shape)
        mask = pd.DataFrame(mask, index=idxs_to_remove)
        samples_to_fix = mask.nunique(axis=1).eq(1)
        if samples_to_fix.any():
            samples_to_fix = samples_to_fix[samples_to_fix]
            views_to_fix = np.random.default_rng(self.random_state).integers(low=0, high=self.n_views,
                                                                             size=len(samples_to_fix))
            for view_idx in np.unique(views_to_fix):
                samples = views_to_fix == view_idx
                samples = samples_to_fix[samples].index
                mask.loc[samples, view_idx] = np.invert(mask.loc[samples, view_idx].astype(bool)).astype(int)

        pseudo_observed_view_indicator.loc[idxs_to_remove] = mask.astype(int)
        return pseudo_observed_view_indicator


    def _mnar_mask(self, sample_names):
        mask = pd.DataFrame(np.ones((len(sample_names), self.n_views)), index=sample_names)
        common_samples = pd.Series(sample_names, index=sample_names).sample(frac=1 - self.p, replace=False,
                                                                            random_state=self.random_state).index
        idxs_to_remove = sample_names.difference(common_samples)
        reference_var = np.random.default_rng(self.random_state).choice(range(1, self.n_views),
                                                                        p = self.weights,
                                                                        size=len(idxs_to_remove))
        reference_var = pd.Series(reference_var, index=idxs_to_remove)
        if self.random_state is None:
            random_state = np.random.choice(100000)
        else:
            random_state = self.random_state
        n_mods_to_remove = {n_mods_to_remove: np.random.default_rng(random_state + i).choice(self.n_views,
                                                                                             size=n_mods_to_remove,
                                                                                             replace=False)
                            for i,n_mods_to_remove in enumerate(np.unique(reference_var))}
        for keys,values in n_mods_to_remove.items():
            mask.loc[reference_var[reference_var == keys].index, values] = 0

        return mask


    def _pm_mask(self, sample_names):
        pseudo_observed_view_indicator = pd.DataFrame(np.ones((len(sample_names), self.n_views)), index=sample_names)
        common_samples = pd.Series(sample_names, index=sample_names).sample(frac=1 - self.p, replace=False,
                                                                            random_state=self.random_state).index
        idxs_to_remove = sample_names.difference(common_samples)
        if self.n_views == 2:
            col = np.random.default_rng(self.random_state).choice(self.n_views)
            pseudo_observed_view_indicator.loc[idxs_to_remove, col] = 0
        else:
            rand = np.random.default_rng(self.random_state)
            mask = rand.choice(2, size=(len(idxs_to_remove), self.n_views -1))
            mask = pd.DataFrame(mask, index=idxs_to_remove, columns=rand.choice(self.n_views, size=self.n_views-1,
                                                                                replace=False))
            samples_to_fix = mask.nunique(axis=1).eq(1)
            if samples_to_fix.any():
                samples_to_fix = samples_to_fix[samples_to_fix]
                views_to_fix = rand.choice(mask.columns, size=len(samples_to_fix))
                for view_idx in np.unique(views_to_fix):
                    samples = views_to_fix == view_idx
                    samples = samples_to_fix[samples].index
                    mask.loc[samples, view_idx] = np.invert(mask.loc[samples, view_idx].astype(bool)).astype(int)
            pseudo_observed_view_indicator.loc[idxs_to_remove, mask.columns] = mask.astype(int)
        return pseudo_observed_view_indicator
