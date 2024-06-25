import copy
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from scipy import optimize

from ..utils import DatasetUtils


class Amputer(BaseEstimator, TransformerMixin):
    r"""
    Generate view missingness patterns in complete multi-view datasets.

    Parameters
    ----------
    Xs : list of array-likes
        - Xs length: n_views
        - Xs[i] shape: (n_samples, n_features_i)
        A list of different views.
    p: float
        Percentaje of incomplete samples.
    mechanism: str, default="EDM"
        One of ["EDM", 'MCAR', 'MAR', 'MNAR', 'PM'].
    random_state: int, default=None
        If int, random_state is the seed used by the random number generator.
    #todo fill args

    Examples
    --------
    >>> from imvc.ampute import Amputer
    >>> from imvc.datasets import LoadDataset
    >>> Xs = LoadDataset.load_dataset("nutrimouse")
    >>> transformer = Amputer(p= 0.2)
    >>> transformer.fit_transform(Xs)
    """

    def __init__(self, p:float, mechanism: str = "EDM", random_state: int = None,
                 opt:str = "logistic", p_obs:float = 0.1, q= 0.3, exclude_inputs:bool = True,
                 p_params= 0.3, cut='both', mcar:bool = False):
        possible_mechanisms = ["EDM", 'MCAR', 'MAR', 'MNAR', "PM"]
        if mechanism not in possible_mechanisms:
            raise ValueError(f"Invalid mechanism. Expected one of: {possible_mechanisms}")

        self.p = p
        self.mechanism = mechanism
        self.random_state = random_state
        self.opt = opt
        self.p_obs = p_obs
        self.q = q
        self.exclude_inputs = exclude_inputs
        self.p_params = p_params
        self.cut = cut
        self.mcar = mcar


    def fit(self, Xs: list, y=None):
        self.n_views = len(Xs)
        return self


    def transform(self, Xs: list, y=None):
        sample_names = Xs[0].index

        if self.mechanism == "EDM":
            pseudo_observed_view_indicator = self._edm_mask(sample_names=sample_names)
        elif self.mechanism == "MCAR":
            pseudo_observed_view_indicator = self._mcar_mask(sample_names=sample_names)
        elif self.mechanism == "PM":
            pseudo_observed_view_indicator = self._pm_mask(sample_names=sample_names)
        else:
            pseudo_observed_view_indicator = np.random.default_rng(self.random_state).normal(size=(len(Xs[0]), self.n_views))
            pseudo_observed_view_indicator = self._produce_missing(X= pseudo_observed_view_indicator, sample_names=sample_names)

        pseudo_observed_view_indicator = pseudo_observed_view_indicator.astype(bool)
        transformed_Xs = DatasetUtils.convert_to_imvd(Xs=Xs, observed_view_indicator=pseudo_observed_view_indicator)

        return transformed_Xs


    def _produce_missing(self, X, sample_names):
        """
        Generate missing values for specifics missing-data mechanism and proportion of missing values.

        Parameters
        ----------
        X : torch.DoubleTensor or np.ndarray, shape (n, d)
            Data for which missing values will be simulated.
            If a numpy array is provided, it will be converted to a pytorch tensor.
        p_miss : float
            Proportion of missing values to generate for variables which will have missing values.
        mecha : str,
                Indicates the missing-data mechanism to be used. "MCAR" by default, "MAR", "MNAR" or "MNARsmask"
        opt: str,
             For mecha = "MNAR", it indicates how the missing-data mechanism is generated: using a logistic regression ("logistic"), quantile censorship ("quantile") or logistic regression for generating a self-masked MNAR mechanism ("selfmasked").
        p_obs : float
                If mecha = "MAR", or mecha = "MNAR" with opt = "logistic" or "quanti", proportion of variables with *no* missing values that will be used for the logistic masking model.
        q : float
            If mecha = "MNAR" and opt = "quanti", quantile level at which the cuts should occur.

        Returns
        ----------
        A dictionnary containing:
        'X_init': the initial data matrix.
        'X_incomp': the data with the generated missing values.
        'mask': a matrix indexing the generated missing values.
        """
        missing_samples = pd.Series(sample_names, index=sample_names).sample(frac=self.p, replace=False,
                                                                             random_state=self.random_state).index
        missing_X = X[missing_samples]

        if self.mechanism == "MAR":
            missing_X = np.insert(missing_X, 0, np.random.default_rng(self.random_state).random(len(missing_X)), axis=1)
            mask = self._MAR_mask(X=missing_X, p=self.p, p_obs=self.p_obs)
            mask = mask[:, ~np.all(mask[1:] == mask[:-1], axis=0)]
            if mask.shape[1] != X.shape[1]:
                raise ValueError("p is too small for this dataset.") from None
        elif self.mechanism == "MNAR" and self.opt == "logistic":
            mask = self._MNAR_mask_logistic(X=missing_X, p=self.p, p_params=self.p_params, exclude_inputs=self.exclude_inputs)
        elif self.mechanism == "MNAR" and self.opt == "quantile":
            mask = self._MNAR_mask_quantiles(X=missing_X, p=self.p, q=self.q, p_params=self.p_params, cut=self.cut, MCAR=self.mcar)
        elif self.mechanism == "MNAR" and self.opt == "selfmasked":
            mask = self.MNAR_self_mask_logistic(X=missing_X, p=self.p)
        else:
            raise ValueError("MNAR mechanism can only be 'logistic', 'quantile' or 'selfmasked'")

        mask = pd.DataFrame(mask, index=missing_samples)
        samples_to_fix = mask.nunique(axis=1).eq(1)
        if samples_to_fix.any():
            samples_to_fix = samples_to_fix[samples_to_fix]
            views_to_fix = np.random.default_rng(self.random_state).integers(low=0, high=self.n_views,
                                                                             size=len(samples_to_fix))
            for view_idx in np.unique(views_to_fix):
                samples = views_to_fix == view_idx
                samples = samples_to_fix[samples].index
                mask.loc[samples, view_idx] = np.invert(mask.loc[samples, view_idx].astype(bool))

        X = pd.DataFrame(np.ones_like(X), index=sample_names)
        X.loc[mask.index] = mask.astype(int)
        return X

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
                mask.loc[samples, view_idx] = np.invert(mask.loc[samples, view_idx].astype(bool))

        pseudo_observed_view_indicator.loc[idxs_to_remove] = mask
        return pseudo_observed_view_indicator


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
                    mask.loc[samples, view_idx] = np.invert(mask.loc[samples, view_idx].astype(bool))
            pseudo_observed_view_indicator.loc[idxs_to_remove, mask.columns] = mask
        return pseudo_observed_view_indicator


    def _MAR_mask(self, X, p, p_obs):
        n, d = X.shape
        mask = np.zeros((n, d)).astype(bool)
        d_obs = max(int(p_obs * d), 1)  ## number of variables that will have no missing values (at least one variable)
        d_na = d - d_obs  ## number of variables that will have missing values
        ### Sample variables that will all be observed, and those with missing values:
        idxs_obs = np.random.default_rng(self.random_state).choice(d, d_obs, replace=False)
        idxs_nas = np.array([i for i in range(d) if i not in idxs_obs])
        ### Other variables will have NA proportions that depend on those observed variables, through a logistic model
        ### The parameters of this logistic model are random.
        ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
        coeffs = self._pick_coeffs(X, idxs_obs=idxs_obs, idxs_nas=idxs_nas)
        ### Pick the intercepts to have a desired amount of missing values
        intercepts = self._fit_intercepts(X[:, idxs_obs], coeffs, p)
        ps = np.dot(X[:, idxs_obs], coeffs) + intercepts
        ps = 1 / (1 + np.exp(-ps))
        ber = np.random.default_rng(self.random_state).random((n, d_na))
        mask[:, idxs_nas] = ber < ps
        return mask


    def _MNAR_mask_logistic(self, X, p, p_params=.3, exclude_inputs=True):
        n, d = X.shape
        mask = np.zeros((n, d)).astype(bool)

        d_params = max(int(p_params * d), 1) if exclude_inputs else d  ## number of variables used as inputs (at least 1)
        d_na = d - d_params if exclude_inputs else d  ## number of variables masked with the logistic model

        ### Sample variables that will be parameters for the logistic regression:
        idxs_params = np.random.default_rng(self.random_state).choice(d,
                                                                      d_params, replace=False) if exclude_inputs else np.arange(d)
        idxs_nas = np.array([i for i in range(d) if i not in idxs_params]) if exclude_inputs else np.arange(d)

        ### Other variables will have NA proportions selected by a logistic model
        ### The parameters of this logistic model are random.

        ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
        coeffs = self._pick_coeffs(X, idxs_obs=idxs_params, idxs_nas=idxs_nas)
        ### Pick the intercepts to have a desired amount of missing values
        intercepts = self._fit_intercepts(X[:, idxs_params], coeffs, p)

        ps = np.dot(X[:, idxs_params], coeffs) + intercepts
        ps = 1 / (1 + np.exp(-ps))

        ber = np.random.default_rng(self.random_state).random((n, d_na))
        mask[:, idxs_nas] = ber < ps

        ## If the inputs of the logistic model are excluded from MNAR missingness,
        ## mask some values used in the logistic model at random.
        ## This makes the missingness of other variables potentially dependent on masked values
        if exclude_inputs:
            mask[:, idxs_params] = np.random.default_rng(self.random_state).random((n, d_params)) < p

        return mask

    def MNAR_self_mask_logistic(self, X, p):
        """
        Missing not at random mechanism with a logistic self-masking model. Variables have missing values probabilities
        given by a logistic model, taking the same variable as input (hence, missingness is independent from one variable
        to another). The intercepts are selected to attain the desired missing rate.

        Parameters
        ----------
        X : torch.DoubleTensor or np.ndarray, shape (n, d)
            Data for which missing values will be simulated.
            If a numpy array is provided, it will be converted to a pytorch tensor.

        p : float
            Proportion of missing values to generate for variables which will have missing values.

        Returns
        -------
        mask : torch.BoolTensor or np.ndarray (depending on type of X)
            Mask of generated missing values (True if the value is missing).

        """

        n, d = X.shape
        ### Variables will have NA proportions that depend on those observed variables, through a logistic model
        ### The parameters of this logistic model are random.

        ### Pick coefficients so that W^Tx has unit variance (avoids shrinking)
        coeffs = self._pick_coeffs(X=X, self_mask=True)
        ### Pick the intercepts to have a desired amount of missing values
        intercepts = self._fit_intercepts(X=X, coeffs=coeffs, p=p, self_mask=True)

        ps = X*coeffs + intercepts
        ps = 1 / (1 + np.exp(-ps))

        ber = np.random.rand(n, d)
        mask = ber < ps
        return mask


    def _MNAR_mask_quantiles(self, X, p, q, p_params, cut='both', MCAR=False):
        """
        Missing not at random mechanism with quantile censorship. First, a subset of variables which will have missing
        variables is randomly selected. Then, missing values are generated on the q-quantiles at random. Since
        missingness depends on quantile information, it depends on masked values, hence this is a MNAR mechanism.

        Parameters
        ----------
        X : torch.DoubleTensor or np.ndarray, shape (n, d)
            Data for which missing values will be simulated.
            If a numpy array is provided, it will be converted to a pytorch tensor.

        p : float
            Proportion of missing values to generate for variables which will have missing values.

        q : float
            Quantile level at which the cuts should occur

        p_params : float
            Proportion of variables that will have missing values

        cut : 'both', 'upper' or 'lower', default = 'both'
            Where the cut should be applied. For instance, if q=0.25 and cut='upper', then missing values will be generated
            in the upper quartiles of selected variables.

        MCAR : bool, default = True
            If true, masks variables that were not selected for quantile censorship with a MCAR mechanism.

        Returns
        -------
        mask : torch.BoolTensor or np.ndarray (depending on type of X)
            Mask of generated missing values (True if the value is missing).

        """
        n, d = X.shape
        mask = np.zeros((n, d)).astype(bool)

        d_na = max(int(p_params * d), 1)  ## number of variables that will have NMAR values
        ### Sample variables that will have imps at the extremes
        idxs_na = np.random.choice(d, d_na, replace=False)  ### select at least one variable with missing values

        ### check if values are greater/smaller that corresponding quantiles
        if cut == 'upper':
            quants = np.partition(a= X[:, idxs_na], kth= 1 - q, axis=0)[q]
            m = X[:, idxs_na] >= quants
        elif cut == 'lower':
            quants = np.partition(a= X[:, idxs_na], kth= q, axis=0)[q]
            m = X[:, idxs_na] <= quants
        elif cut == 'both':
            u_quants = np.partition(a= X[:, idxs_na], kth= 1 - q, axis=0)[q]
            l_quants = np.partition(a= X[:, idxs_na], kth= q, axis=0)[q]
            m = (X[:, idxs_na] <= l_quants) | (X[:, idxs_na] >= u_quants)

        ### Hide some values exceeding quantiles
        ber = np.random.default_rng(self.random_state).random((n, d_na))
        mask[:, idxs_na] = (ber < p) & m

        if MCAR:
            ## Add a mcar mecanism on top
            mask = mask | (np.random.default_rng(self.random_state).random((n, d)) < p)

        return mask


    def _pick_coeffs(self, X, idxs_obs=None, idxs_nas=None, self_mask=False):
        n, d = X.shape
        if self_mask:
            coeffs = np.random.default_rng(self.random_state).normal(size=d)
            Wx = X * coeffs
            coeffs /= np.std(Wx, 0)
        else:
            d_obs = len(idxs_obs)
            d_na = len(idxs_nas)
            coeffs = np.random.default_rng(self.random_state).normal(size=(d_obs, d_na))
            Wx = np.dot(X[:, idxs_obs], coeffs)
            coeffs /= np.std(Wx, 0, keepdims=True)
        return coeffs

    def _fit_intercepts(self, X, coeffs, p, self_mask=False):
        if self_mask:
            d = len(coeffs)
            intercepts = [optimize.bisect(lambda x: (1 / (1 + np.exp(-(X * coeffs[j] + x)))).mean() - p, -50, 50) for j
                          in range(d)]
        else:
            d_obs, d_na = coeffs.shape
            intercepts = [
                optimize.bisect(lambda x: (1 / (1 + np.exp(-(np.dot(X, coeffs[:, j]) + x)))).mean() - p, -50, 50) for j
                in range(d_na)]
        return intercepts
