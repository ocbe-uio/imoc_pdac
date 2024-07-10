import copy
from typing import Union

import numpy as np
import pandas as pd

from ..impute.observed_view_indicator import get_observed_view_indicator


class DatasetUtils:
    r"""
    A utility class that provides general methods for working with multi-view datasets.
    """

    @staticmethod
    def convert_to_imvd(Xs: list, observed_view_indicator) -> list:
        r"""
        Generate view missingness patterns in complete multi-view datasets.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        observed_view_indicator: pd.DataFrame
            Boolean array-like indicating observed views for each sample.

        Returns
        -------
        transformed_Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

         Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset
        >>> from imvc.ampute import Amputer
        >>> from sklearn.pipeline import make_pipeline
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> transformer = make_pipeline(Amputer(p=0.2, mechanism="MCAR", random_state=42),
                                        ObservedViewIndicator().set_output(transformed="pandas"))
        >>> Xs = transformer.fit_transform(Xs)
        >>> observed_view_indicator = get_observed_view_indicator(Xs)
        >>> DatasetUtils.convert_to_imvd(Xs = Xs, observed_view_indicator = observed_view_indicator)
        """
        transformed_Xs = []
        if not isinstance(Xs, pd.DataFrame):
            observed_view_indicator = pd.DataFrame(observed_view_indicator)
        for X_idx, X in enumerate(Xs):
            idxs_to_remove = observed_view_indicator[observed_view_indicator[X_idx] == False].index
            transformed_X = copy.deepcopy(X)
            transformed_X.loc[idxs_to_remove] = np.nan
            transformed_Xs.append(transformed_X)

        return transformed_Xs


    @staticmethod
    def get_n_views(Xs: list) -> int:
        r"""
        Get the number of views of a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        n_views: int
            Number of views.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> DatasetUtils.get_n_views(Xs = Xs)
        """

        n_views = len(Xs)
        return n_views


    @staticmethod
    def get_n_samples_by_view(Xs: list) -> int:
        r"""
        Get the number of samples in each view.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        n_samples_by_view: pd.Series
            Number of samples in each view.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> DatasetUtils.get_n_samples_by_view(Xs = Xs)
        """

        n_samples_by_view = get_observed_view_indicator(Xs)
        n_samples_by_view = n_samples_by_view.sum(axis=0)
        return n_samples_by_view


    @staticmethod
    def get_complete_sample_names(Xs: list) -> pd.Index:
        r"""
        Get complete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        samples: pd.Index
            Sample names with full data.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_complete_sample_names(Xs = Xs)
        """

        samples = get_observed_view_indicator(Xs)
        samples = samples[samples.all(1)].index
        return samples


    @staticmethod
    def get_incomplete_sample_names(Xs: list) -> pd.Index:
        r"""
        Get incomplete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        samples: pd.Index
            Sample names with incomplete data.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_incomplete_sample_names(Xs = Xs)
        """

        samples = get_observed_view_indicator(Xs)
        samples = samples[~samples.all(1)].index
        return samples


    @staticmethod
    def get_samples_by_view(Xs: list, return_as_list: bool = True) -> Union[list, dict]:
        r"""
        Get the samples for each view in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        return_as_list : bool, default=True
            If True, the function will return a list; a dict otherwise.

        Returns
        -------
        samples: list or dict of pd.Index
            If list, each element in the list is the sample names for each view. If dict, keys are the views and the
            values are the sample names.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_samples_by_view(Xs = Xs)
        """

        observed_view_indicator = get_observed_view_indicator(Xs)
        if return_as_list:
            samples = [view_profile[view_profile == 1].index for X_idx, view_profile in observed_view_indicator.items()]
        else:
            samples = {X_idx: view_profile[view_profile == 1].index for X_idx, view_profile in observed_view_indicator.items()}
        return samples


    @staticmethod
    def get_missing_samples_by_view(Xs: list, return_as_list: bool = True) -> Union[list, dict]:
        r"""
        Get the samples not present in each view in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        return_as_list : bool, default=True
            If list, each element in the list is the sample names for each view. If dict, keys are the views and the
            values are the sample names.

        Returns
        -------
        samples: dict of pd.Index
            Dictionary of missing samples for each view.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_missing_samples_by_view(Xs = Xs)
        """

        observed_view_indicator = get_observed_view_indicator(Xs)
        if isinstance(Xs[0], pd.DataFrame):
            observed_view_indicator = pd.DataFrame(observed_view_indicator, index=Xs[0].index)
        else:
            observed_view_indicator = pd.DataFrame(observed_view_indicator)
        if return_as_list:
            samples = [view_profile[view_profile == 0].index.to_list()
                       for X_idx, view_profile in observed_view_indicator.items()]
        else:
            samples = {X_idx: view_profile[view_profile == 0].index.to_list()
                       for X_idx, view_profile in observed_view_indicator.items()}
        return samples


    @staticmethod
    def get_n_complete_samples(Xs: list) -> int:
        r"""
        Get the number of complete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        int: number of complete samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_n_complete_samples(Xs = Xs)
        """

        n_samples = len(DatasetUtils.get_complete_sample_names(Xs=Xs))
        return n_samples


    @staticmethod
    def get_n_incomplete_samples(Xs: list) -> int:
        r"""
        Get the number of incomplete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        int: number of incomplete samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_n_incomplete_samples(Xs = Xs)
        """

        n_samples = len(DatasetUtils.get_incomplete_sample_names(Xs=Xs))
        return n_samples


    @staticmethod
    def get_percentage_complete_samples(Xs: list) -> float:
        r"""
        Get the percentage of complete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        float: percentage of complete samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_percentage_complete_samples(Xs = Xs)
        """

        percentage_samples = DatasetUtils.get_n_complete_samples(Xs=Xs) / len(Xs[0]) * 100
        return percentage_samples


    @staticmethod
    def get_percentage_incomplete_samples(Xs: list) -> float:
        r"""
        Get the percentage of incomplete samples in a multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        float: percentage of incomplete samples.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.get_percentage_incomplete_samples(Xs = Xs)
        """

        percentage_samples = DatasetUtils.get_n_incomplete_samples(Xs=Xs) / len(Xs[0]) * 100
        return percentage_samples


    def remove_missing_sample_from_view(Xs: list) -> list:
        r"""
        Remove missing samples from each specific views.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs: list of array-likes.
            - Xs length: n_views
            - Xs[i] shape: (n_samples_i, n_features_i)
            A list of different views.

        Examples
        --------
        >>> from imvc.utils import DatasetUtils
        >>> from imvc.ampute import Amputer
        >>> from imvc.datasets import LoadDataset
        >>> Xs = LoadDataset.load_dataset(dataset_name="nutrimouse")
        >>> Xs = Amputer(p=0.2, mechanism="MCAR", random_state=42).fit_transform(Xs)
        >>> DatasetUtils.remove_missing_sample_from_view(Xs = Xs)
        """

        transformed_Xs = [X[np.invert(np.isnan(X).all(1))] for X in Xs]
        return transformed_Xs


    def force_all_samples(Xs: list):
        r"""
        Add missing samples to each view, in a way that all the views will have all samples.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        Xs: list of array-likes.
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
        """

        samples = set(sum([X.index.to_list() for X in Xs], []))
        Xs = [pd.concat([X, pd.DataFrame(np.nan, index= pd.Index(samples).difference(X.index),
                                         columns= X.columns)]).loc[samples] for X in Xs]
        return Xs


    def convert_mvd_from_list_to_dict(Xs: list, keys: list):
        r"""
        Convert a multi-view dataset in list format to a dict format.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs: dict of array-likes.
            - Xs length: n_views
            - Xs[key] shape: (n_samples, n_features_i)
        """

        transformed_Xs = {key:X for key,X in zip(keys, Xs)}
        return transformed_Xs


    def convert_mvd_from_dict_to_list(Xs: list, keys: list):
        r"""
        Convert a multi-view dataset in list format to a dict format.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        transformed_Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.
        """

        transformed_Xs = list(Xs.values())
        return transformed_Xs


    def select_complete_samples(Xs: list):
        r"""
        Remove samples with missing views from the whole multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        Xs: list of array-likes.
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
        """

        samples = DatasetUtils.get_complete_sample_names(Xs=Xs)
        Xs = [X.loc[samples] for X in Xs]
        return Xs


    def select_incomplete_samples(Xs: list):
        r"""
        Remove samples with no missing views from the whole multi-view dataset.

        Parameters
        ----------
        Xs : list of array-likes
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
            A list of different views.

        Returns
        -------
        Xs: list of array-likes.
            - Xs length: n_views
            - Xs[i] shape: (n_samples, n_features_i)
        """

        samples = DatasetUtils.get_incomplete_sample_names(Xs=Xs)
        Xs = [X.loc[samples] for X in Xs]
        return Xs



