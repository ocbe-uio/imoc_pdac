import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.utils import shuffle
from validclust import dunn
from reval.utils import kuhn_munkres_algorithm

from src.utils import dbcv


class GetMetrics:


    @staticmethod
    def compute_supervised_metrics(y_true, y_pred):
        random_preds = [pd.Series(y_true).value_counts().index[0]] * len(y_true)
        perm_clust_labels = kuhn_munkres_algorithm(true_lab=y_true, pred_lab=y_pred)
        mcc, p_value = GetMetrics.compute_mcc(y_true=y_true, y_pred=perm_clust_labels)
        supervised_metrics = {
            'ACC': metrics.accuracy_score(y_true=y_true, y_pred=perm_clust_labels),
            'MCC': mcc,
            'MCC (p-value)': p_value,
            'F1': metrics.f1_score(y_true=y_true, y_pred=perm_clust_labels, average='macro'),
            'precision': metrics.precision_score(y_true=y_true, y_pred=perm_clust_labels, average='macro', zero_division=0),
            'recall': metrics.recall_score(y_true=y_true, y_pred=perm_clust_labels, average='macro', zero_division=0),
            "bal_acc": metrics.balanced_accuracy_score(y_true=y_true, y_pred=y_pred),
            "ami": metrics.adjusted_mutual_info_score(labels_true=y_true, labels_pred=y_pred),
            "ari": metrics.adjusted_rand_score(labels_true=y_true, labels_pred=y_pred),
            "completeness": metrics.completeness_score(labels_true=y_true, labels_pred=y_pred),
            "random_acc": metrics.accuracy_score(y_true=y_true, y_pred=random_preds),
            "random_f1": metrics.f1_score(y_true=y_true, y_pred=random_preds, average='macro'),
        }
        return supervised_metrics


    @staticmethod
    def compute_unsupervised_metrics(X, y_pred, random_state):
        if len(np.unique(y_pred)) == 1:
            unsupervised_metrics = {"silhouette": np.nan, "vrc": np.nan, "db": np.nan, "dunn": np.nan}
        else:
            unsupervised_metrics = {
                "silhouette": metrics.silhouette_score(X = X, labels = y_pred, random_state= random_state),
                "vrc": metrics.calinski_harabasz_score(X = X, labels = y_pred),
                "db": metrics.davies_bouldin_score(X = X, labels = y_pred),
                "dbcv": dbcv(X = X.values, labels = y_pred.values),
                "dunn": dunn(dist = metrics.pairwise_distances(X), labels = y_pred),
            }
        return unsupervised_metrics


    @staticmethod
    def compute_mcc(y_true, y_pred, n_permutations=1000):
        mcc_value = metrics.matthews_corrcoef(y_true=y_true, y_pred=y_pred)
        p_values = [metrics.matthews_corrcoef(y_true=shuffle(y_true, random_state=i),
                                              y_pred=shuffle(y_pred, random_state=i))
                    for i in range(n_permutations)]
        p_values = (pd.Series(p_values) > mcc_value).sum() / len(p_values)
        if p_values == 0:
            p_values = 1/(n_permutations + 1)
        return mcc_value, p_values
