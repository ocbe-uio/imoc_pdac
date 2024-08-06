import argparse
import os
import time
from datetime import datetime
import pandas as pd
import torch
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler, FunctionTransformer
from mvlearn.decomposition import AJIVE, GroupPCA
from mvlearn.cluster import MultiviewSpectralClustering, MultiviewCoRegSpectralClustering
from imvc.cluster import OSLFIMVC, DAIMC, EEIMVC, LFIMVC, MKKMIK, MSNE, SIMCADC, PIMVC, IMSR, OMVC, OPIMC, SUMO, NEMO, \
    IMSCAGL
from imvc.cluster.monet import MONET
from imvc.decomposition import DFMF, MOFA, DeepMF, jNMF
from imvc.preprocessing import MultiViewTransformer, NormalizerNaN, ConcatenateViews

from src.models import Model
from settings import TIME_RESULTS_PATH, TIME_LOGS_PATH, TIME_ERRORS_PATH, RANDOM_STATE, DATASET_TABLE_PATH
from src.commons import CommonOperations

datasets = pd.read_csv(DATASET_TABLE_PATH)["dataset"].to_list()
if "nutrimouse" in datasets:
    pos = datasets.index('nutrimouse')
    datasets[pos:pos + 1] = ('nutrimouse_genotype', 'nutrimouse_diet')

algorithms = {
    "Concat": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  ConcatenateViews(),
                                    StandardScaler().set_output(transform='pandas'),
                                    KMeans()), "params": {}},
    "NMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  ConcatenateViews(),
                                  MinMaxScaler().set_output(transform='pandas'),
                                  NMF().set_output(transform='pandas'), StandardScaler(), KMeans()), "params": {}},
    "MVSC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform= "pandas")),
                                                  MultiviewSpectralClustering()),
                             "params": {}},
    "MVCRSC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform= "pandas")),
                                                       MultiviewCoRegSpectralClustering()),
                                  "params": {}},
    "GPCA": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler()), GroupPCA(), StandardScaler(), KMeans()),
                 "params": {}},
    "AJIVE": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler()), AJIVE(),
                                   MultiViewTransformer(FunctionTransformer(pd.DataFrame)), ConcatenateViews(),
                                   StandardScaler(), KMeans()),
              "params": {}},
    "SNF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform= "pandas"))), "params": {}},
    "DAIMC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                   MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                   DAIMC()), "params": {}},
    "EEIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    EEIMVC()), "params": {}},
    "IMSCAGL": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                   IMSCAGL()), "params": {}},
    "IMSR": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                   IMSR()), "params": {}},
    "LFIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    LFIMVC()), "params": {}},
    "MKKMIK": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    MKKMIK()), "params": {}},
    "MONET": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                   MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                   MONET()), "params": {}},
    "MSNE": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  MSNE()), "params": {}},
    "OMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  OMVC()), "params": {}},
    "OPIMC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                   MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                     OPIMC()), "params": {}},
    "OSLFIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                      MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                      OSLFIMVC()), "params": {}},
    "PIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                   MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                   PIMVC()), "params": {}},
    "SIMCADC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                     MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                     SIMCADC()), "params": {}},
    "SUMO": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  SUMO()), "params": {}},
    "DeepMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    ConcatenateViews(), StandardScaler(),
                                    FunctionTransformer(lambda x: torch.from_numpy(x).float().cuda().t()),
                                    DeepMF(), FunctionTransformer(lambda x: x.cpu().detach().numpy()),
                                    StandardScaler(), KMeans()),
                 "params": {}},
    "DFMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  DFMF().set_output(transform="pandas"),
                                  StandardScaler().set_output(transform="pandas"), KMeans()),
             "params": {}},
    "MOFA": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  MOFA().set_output(transform="pandas"),
                                  ConcatenateViews(), StandardScaler().set_output(transform="pandas"), KMeans()),
             "params": {}},
    "jNMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(MinMaxScaler().set_output(transform="pandas")),
                                  jNMF().set_output(transform="pandas"),
                                  StandardScaler().set_output(transform="pandas"), KMeans()),
             "params": {}},
}

# args = lambda: None
# args.continue_benchmarking, args.n_jobs, args.save_results = True, 2, False
parser = argparse.ArgumentParser()
parser.add_argument('-save_results', default= False, action='store_true')
args = parser.parse_args()
if not args.save_results:
    TIME_RESULTS_PATH = os.path.join("test", "time_evaluation.csv")
    TIME_LOGS_PATH = os.path.join("test", "time_logs.txt")
    TIME_ERRORS_PATH = os.path.join("test", "time_errors.txt")

if os.path.exists(TIME_RESULTS_PATH):
    results = pd.read_csv(TIME_RESULTS_PATH)
else:
    results = pd.DataFrame(algorithms.keys(), columns=["algorithm"]).merge(pd.DataFrame(datasets, columns=["dataset"]),
                                                                           how='cross')
    results["time"] = -1
    results["finished"] = False
    results["completed"] = False
    results["comments"] = "{}"

    os.remove(TIME_LOGS_PATH) if os.path.exists(TIME_LOGS_PATH) else None
    os.remove(TIME_ERRORS_PATH) if os.path.exists(TIME_ERRORS_PATH) else None
    open(TIME_LOGS_PATH, 'w').close()
    open(TIME_ERRORS_PATH, 'w').close()

for dataset_name in datasets:
    Xs, y, n_clusters = CommonOperations.load_Xs(dataset_name=dataset_name)

    for idx, row in results[(~results["finished"]) & (results["dataset"] == dataset_name)].iterrows():
        alg_name = row["algorithm"]
        if alg_name in ["IntNMF", "COCA", "jNMF", "NEMO"]:
            continue
        alg = algorithms[alg_name]
        with open(TIME_LOGS_PATH, "a") as f:
            f.write(f'\n {dataset_name} \t {alg_name} \t {datetime.now()}')

        start_time = time.perf_counter()
        try:
            clusters, _ = Model(alg_name=alg_name, alg=alg).method(train_Xs=Xs, n_clusters=n_clusters,
                                                                   random_state=RANDOM_STATE, run_n=0)
            elapsed_time = time.perf_counter() - start_time
            errors_dict = {}
            completed = True
        except Exception as exception:
            elapsed_time = time.perf_counter() - start_time
            exception_name = type(exception).__name__
            errors_dict = {f"{exception_name}": exception}
            with open(TIME_ERRORS_PATH, "a") as f:
                f.write(f'\n {dataset_name} \t {alg_name} \t {exception_name}: {exception} \t {datetime.now()}')
            completed = False

        results.loc[idx, ["time", "finished", "completed", "comments"]] = [elapsed_time, True, completed, errors_dict]
        results.to_csv(TIME_RESULTS_PATH, index=False)

print("Completed successfully!")
with open(TIME_LOGS_PATH, "a") as f:
    f.write(f'\n Completed successfully \t {datetime.now()}')

