import copy

import numpy as np
import pandas as pd


class DatasetUtils:
    r"""
    A utility class that provides general methods for working with incomplete multi-view datasets.
    """

    @staticmethod
    def convert_mvd_into_imvd(Xs: list, p, assess_percentage: bool = True, random_state: int = None):
        r"""
        Randomly drop samples in a multi-view dataset to convert it into an incomplete multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        p: list or float
            The percentaje that each view will have for missing samples. If p is float, all the views will have the
            same percentaje.
        assess_percentage: bool
            If False, each view is dropped independently.
        random_state: int, default None
            If int, random_state is the seed used by the random number generator.

        Returns
        -------
        imvd : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.

         Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0)
        >>> Xs = DatasetUtils.ampute(Xs = Xs, p = [0.2, 0.5])
        """
        if not isinstance(p, list):
            p = [p]
        if len(p) != len(Xs):
            p = p*len(Xs)

        if sum(p) > 0:
            if assess_percentage:
                p = [prob/len(p) for prob in p]
                sample_names = DatasetUtils.get_sample_names(Xs)
                total_len = len(sample_names)
                common_samples = pd.Series(sample_names).sample(frac= 1 - sum(p), random_state=random_state)
                sampled_names = copy.deepcopy(common_samples)
                imvd = []
                for X_idx,X in enumerate(Xs):
                    x_per_view = X.drop(sampled_names)
                    if X_idx != len(Xs)-1:
                        x_per_view = x_per_view.drop(x_per_view.sample(
                            n = int(p[X_idx]*total_len),
                            random_state = random_state + X_idx if random_state is not None else random_state).index)
                    sampled_names = pd.concat([sampled_names, pd.Series(x_per_view.index)])
                    idxs_to_keep = pd.concat([common_samples, pd.Series(x_per_view.index)])
                    imvd.append(X.loc[idxs_to_keep])
            else:
                imvd = [X.drop(X.sample(frac = p[X_idx]/len(p),
                                        random_state = random_state + X_idx if random_state is not None else random_state).index)
                        for X_idx,X in enumerate(Xs)]
        else:
            imvd = Xs
        return imvd


    @staticmethod
    def get_observed_view_indicator(Xs: list):
        r"""
        Get the missing view panel of an incomplete multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.

        Returns
        -------
        sample_view_panel: pd.DataFrame with binary values indicating missing views for each sample (1 has the view,
        0 otherwise).

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset

        >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
        >>> observed_view_indicator = ObservedViewIndicator().set_output(transform="pandas").fit_transform(Xs = Xs)
        """

        sample_view_panel = pd.concat([X.index.to_series() for X in Xs], axis = 1).sort_index()
        sample_view_panel = sample_view_panel.mask(sample_view_panel.isna(), 0).where(sample_view_panel.isna(), 1).astype(int)
        return sample_view_panel


    @staticmethod
    def get_sample_names(Xs: list):
        r"""
        Get all the samples in an incomplete multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.

        Returns
        -------
        samples: pd.Index with all samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset

        >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
        >>> samples = DatasetUtils.get_sample_names(Xs = Xs)
        """

        samples = pd.Index(set(sum([X.index.to_list() for X in Xs], [])))
        return samples


    @staticmethod
    def shuffle_imvd(Xs: list, random_state: int = None):
        r"""
        Shuffle the dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.
        random_state: int, default None
            If int, random_state is the seed used by the random number generator.

        Returns
        -------
        Xs: list of array-likes.
            Incomplete multi-view dataset with shuffled samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_incomplete_nutrimouse(p = 0.2)
        >>> Xs = DatasetUtils.shuffle_imvd(Xs = Xs)
        """

        observed_view_indicator = ObservedViewIndicator().set_output(transform="pandas").fit_transform(Xs)
        samples = observed_view_indicator.sample(frac = 1., random_state = random_state).index
        Xs = [X.loc[samples.intersection(X.index)] for X in Xs]
        return Xs







