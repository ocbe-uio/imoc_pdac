import numpy as np
import pandas as pd
from sklearn.cluster import spectral_clustering
from snf import compute
from rpy2.robjects.packages import importr

from src.utils import Utils


class Model:
    def __init__(self, alg_name, alg):
        self.alg_name = alg_name.lower() if alg_name in ["GroupPCA", "AJIVE", "NMFC", "DFMF", "MOFA", "MONET"] else "standard"
        self.method = alg_name.lower() if alg_name in ["SNF", "IntNMF", "COCA"] else "sklearn_method"
        self.alg_name = eval(f"self.{self.alg_name.lower()}")
        self.method = eval(f"self.{self.method.lower()}")
        self.alg = alg


    def sklearn_method(self, train_Xs, n_clusters, random_state, run_n):
        model, params = self.alg["alg"], self.alg["params"]
        model = self.alg_name(model=model, n_clusters=n_clusters, random_state=random_state, run_n=run_n)
        clusters = model.fit_predict(train_Xs)
        if self.alg_name in ["DAIMC", "PIMVC"]:
            transformed_Xs = model[-1].V_
        elif self.alg_name in ["EEIMVC", "LFIMVC", "MKKMIK", "OSLFIMVC"]:
            transformed_Xs = model[-1].H_
        elif self.alg_name == "IMSR":
            transformed_Xs = model[-1].Z_
        elif self.alg_name == "MSNE":
            transformed_Xs = model[-1].embeddings_
        elif self.alg_name == "OMVC":
            transformed_Xs = model[-1].U_star_loss_
        elif self.alg_name == "SIMCADC":
            transformed_Xs = model[-1].U
        elif self.alg_name in ["MVSpectralClustering", "MVCoRegSpectralClustering"]:
            transformed_Xs = model[-1].embedding_
        else:
            transformed_Xs = model[:-1].transform(train_Xs)
        return clusters, transformed_Xs


    def snf(self, train_Xs, n_clusters, random_state, run_n):
        model = self.alg["alg"]
        train_Xs = model.fit_transform(train_Xs)
        k_snf = np.ceil(len(train_Xs[0]) / 10).astype(int)
        affinities = compute.make_affinity(train_Xs, normalize=False, K=k_snf)
        fused = compute.snf(affinities, K=k_snf)
        clusters = spectral_clustering(fused, n_clusters=n_clusters, random_state=random_state + run_n)
        transformed_Xs = pd.DataFrame(fused, index=train_Xs[0].index)
        return clusters, transformed_Xs


    def intnmf(self, train_Xs, n_clusters, random_state, run_n):
        nmf = importr("IntNMF")
        model = self.alg["alg"]
        train_Xs = model.fit_transform(train_Xs)
        clusters = nmf.nmf_mnnals(dat=Utils.convert_df_to_r_object(train_Xs),
                                  k=n_clusters, seed=int(random_state + run_n))[-1]
        clusters = np.array(clusters) - 1
        return clusters, model


    def coca(self, train_Xs, n_clusters, random_state, run_n):
        base, coca = importr("base"), importr("coca")
        model = self.alg["alg"]
        train_Xs = model.fit_transform(train_Xs)
        base.set_seed(int(random_state + run_n))
        clusters = coca.buildMOC(Utils.convert_df_to_r_object(train_Xs), M=len(train_Xs), K=n_clusters)[0]
        clusters = coca.coca(clusters, K=n_clusters)[1]
        clusters = np.array(clusters) - 1
        return clusters, model


    def grouppca(self, model, n_clusters, random_state, run_n):
        model[1].set_params(n_components=n_clusters, random_state=random_state + run_n, multiview_output=False)
        model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
        return model


    def ajive(self, model, n_clusters, random_state, run_n):
        model[1].set_params(joint_rank=n_clusters, random_state=random_state + run_n)
        model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
        return model


    def nmfc(self, model, n_clusters, random_state, run_n):
        model[-1].set_params(n_components=n_clusters, random_state=random_state + run_n)
        return model


    def monet(self, model, n_clusters, random_state, run_n):
        model[-1].set_params(random_state=random_state + run_n)
        return model


    # def deepmf(self, model, n_clusters, random_state, run_n):
    #     model[1].set_params(joint_rank=n_clusters, random_state=random_state + run_n)
    #     model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
    #     return model


    def dfmf(self, model, n_clusters, random_state, run_n):
        model[1].set_params(n_components=n_clusters, random_state=random_state + run_n)
        model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
        return model


    def mofa(self, model, n_clusters, random_state, run_n):
        model[1].set_params(factors=n_clusters, random_state=random_state + run_n)
        model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
        return model


    def standard(self, model, n_clusters, random_state, run_n):
        model[-1].set_params(n_clusters=n_clusters, random_state=random_state + run_n)
        return model
