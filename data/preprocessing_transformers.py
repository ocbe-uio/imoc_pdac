import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import euclidean_distances
from collections import defaultdict


class InitialProcessing:
    def __init__(self, file_path):
        self.file_path = file_path
        self.df = None
        self.df_cleaned = None

    def load_and_transpose(self):
        self.df = pd.read_csv(self.file_path)
        self.df = self.df.transpose()
        new_header = self.df.iloc[0]
        self.df = self.df[1:]
        self.df.columns = new_header
        print(f"Original number of features: {len(self.df.columns)} features")

    def clean_data(self):
        self.df_cleaned = self.df.dropna(axis = 1, how = 'all')
        print(f"{self.__class__.__name__} keeping {len(self.df_cleaned.columns)} features")

    def patient_id(self):
        new_index = self.df_cleaned.index.str.split('-').str[:3].str.join('-')
        self.df_cleaned.index = new_index

    def process_data(self):
        self.load_and_transpose()
        self.clean_data()
        self.patient_id()
        return self.df_cleaned


class RemoveFeaturesWithZeros(BaseEstimator, TransformerMixin):

    def __init__(self, threshold: float, verbose: bool = False):
        self.threshold = threshold
        self.verbose = verbose

    def fit(self, X, y = None):
        self.columns_ = X.columns[((X != 0).sum(axis = 0) / len(X)) >= self.threshold]
        if self.verbose:
            print(f"{self.__class__.__name__} keeping {len(self.columns_)} features")
        return self

    def transform(self, X, y = None):
        transformed_X = X[self.columns_]
        return transformed_X


class RemoveFeaturesWithNaN(BaseEstimator, TransformerMixin):

    def __init__(self, threshold: float, verbose: bool = False):
        self.threshold = threshold
        self.verbose = verbose

    def fit(self, X, y = None):
        self.columns_ = X.columns[((X.notna()).sum(axis = 0) / len(X)) >= self.threshold]
        if self.verbose:
            print(f"{self.__class__.__name__} keeping {len(self.columns_)} features")
        return self

    def transform(self, X, y = None):
        transformed_X = X[self.columns_]
        return transformed_X


class RemoveFeaturesLowMAD(BaseEstimator, TransformerMixin):

    def __init__(self, percentage_to_keep: float, verbose: bool = False):
        self.percentage_to_keep = percentage_to_keep
        self.verbose = verbose

    def fit(self, X, y = None):
        X = X.apply(pd.to_numeric, errors='coerce')
        var = np.abs(X - np.median(X, axis = 0))
        var = np.median(var, axis = 0)
        var = pd.Series(var, index = X.columns)
        columns = var.nlargest(n = int(X.shape[1] * self.percentage_to_keep)).index
        self.columns_ = X.columns.intersection(columns)
        if self.verbose:
            print(f"{self.__class__.__name__} keeping {len(self.columns_)} features")
        return self

    def transform(self, X, y = None):
        transformed_X = X[self.columns_]
        return transformed_X


class RemoveCorrelatedFeatures(BaseEstimator, TransformerMixin):

    def __init__(self, threshold: float, verbose: bool = False):
        self.threshold = threshold
        self.verbose = verbose

    def fit(self, X, y = None):
        X = X.apply(pd.to_numeric, errors='coerce')
        corr_mat = np.abs(np.corrcoef(X.T.values))
        np.fill_diagonal(corr_mat, 0)
        upper = pd.DataFrame(corr_mat, index = X.columns, columns = X.columns)
        self.columns_to_drop_ = upper.columns[(upper > self.threshold).any()]
        if self.verbose:
            print(f"{self.__class__.__name__} keeping {len(X.columns) - len(self.columns_to_drop_)} features")
        return self

    def transform(self, X, y = None):
        transformed_X = X.drop(columns = self.columns_to_drop_)
        return transformed_X


class Log2Transformation(FunctionTransformer):

    def __init__(self):
        super().__init__(lambda x: np.log2(1 + x.astype(float)))


class GeneMutations(BaseEstimator, TransformerMixin):
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.sorted_genes_ = []

    def fit(self, X, y = None):
        stacked_genes = X.stack()
        genes = stacked_genes.unique()
        gene_counts = stacked_genes.value_counts()
        self.sorted_genes_ = gene_counts.sort_values(ascending=False).index.tolist()
        if self.verbose:
            print(f"{self.__class__.__name__} has {len(self.sorted_genes_)} genes")
        return self

    def transform(self, X, y = None):
        transformed_X = pd.DataFrame(0, index = X.index, columns = self.sorted_genes_)
        for idx, row in X.iterrows():
            for gene in row.dropna():
                if gene in transformed_X.columns:
                    transformed_X.at[idx, gene] = 1
        return transformed_X


class ValueImputation(BaseEstimator, TransformerMixin):
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.imputer = KNNImputer()
        self.scaler = StandardScaler()

    def fit(self, X, y = None):
        self.scaler.fit(X)
        X_scaled = self.scaler.transform(X)
        self.imputer.fit(X_scaled)
        return self

    def transform(self, X, y = None):
        missing_values = X.isna().sum().sum()
        X_scaled = self.scaler.transform(X)
        X_imputed = self.imputer.transform(X_scaled)
        X_inverse = self.scaler.inverse_transform(X_imputed)
        transformed_X = pd.DataFrame(X_inverse, columns = X.columns, index = X.index)
        if self.verbose:
            print(f"{self.__class__.__name__} transforming {missing_values} values")
        return transformed_X


class PrincipalFeatures(BaseEstimator, TransformerMixin):
    def __init__(self, diff_n_features, explained_var, verbose=False):
        self.verbose = verbose
        self.scaler = StandardScaler()
        self.diff_n_features = diff_n_features
        self.explained_var = explained_var

    def fit(self, X, y = None):
        X_scaled = self.scaler.fit_transform(X)
        pca = PCA().fit(X_scaled)
        explained_variance = pca.explained_variance_ratio_
        cumulative_expl_var = np.cumsum(explained_variance)
        q = np.argmax(cumulative_expl_var >= self.explained_var) + 1
        A_q = pca.components_.T[:, :q]
        clusternumber = min(q + self.diff_n_features, X.shape[1])
        kmeans = KMeans(n_clusters=clusternumber).fit(A_q)
        clusters = kmeans.predict(A_q)
        cluster_centers = kmeans.cluster_centers_
        dists = defaultdict(list)
        for i, c in enumerate(clusters):
            dist = euclidean_distances([A_q[i, :]], [cluster_centers[c, :]])[0][0]
            dists[c].append((i, dist))
        self.indices_ = [sorted(f, key=lambda x: x[1])[0][0] for f in dists.values()]
        if self.verbose:
            print(f"{self.__class__.__name__} keeping {len(self.indices_)} features")
        self.feature_names_ = X.columns if hasattr(X, 'columns') else np.arange(X.shape[1])
        self.selected_features_ = [self.feature_names_[i] for i in self.indices_]
        return self

    def transform(self, X, y = None):
        X_transformed = pd.DataFrame(X, columns=self.selected_features_, index=X.index)
        return X_transformed
