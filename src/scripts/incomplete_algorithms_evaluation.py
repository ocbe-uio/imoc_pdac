import os.path
from pandarallel import pandarallel
from sklearn.cluster import KMeans
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from imvc.decomposition import DFMF, MOFA
from imvc.preprocessing import MultiViewTransformer, NormalizerNaN
from imvc.cluster import NEMO, PIMVC, EEIMVC, SIMCADC, LFIMVC, IMSR, MKKMIK
from settings import INCOMPLETE_RESULTS_PATH, INCOMPLETE_SUBRESULTS_PATH, INCOMPLETE_LOGS_PATH, INCOMPLETE_ERRORS_PATH, \
    TIME_RESULTS_PATH, DATASET_TABLE_PATH, amputation_mechanisms, probs, \
    imputation, runs_per_alg, INCOMPLETE_RESULTS_FILE, INCOMPLETE_LOGS_FILE, INCOMPLETE_ERRORS_FILE, \
    INCOMPLETE_SUBRESULTS_FOLDER, views, n_clusters, best_combination, run_amputation, select_datasets
from src.commons import CommonOperations

args = CommonOperations.get_args()

if not args.save_results:
    results_folder = 'test'
    INCOMPLETE_RESULTS_PATH = os.path.join(results_folder, INCOMPLETE_RESULTS_FILE)
    INCOMPLETE_LOGS_PATH = os.path.join(results_folder, INCOMPLETE_LOGS_FILE)
    INCOMPLETE_ERRORS_PATH = os.path.join(results_folder, INCOMPLETE_ERRORS_FILE)
    INCOMPLETE_SUBRESULTS_PATH = os.path.join(results_folder, INCOMPLETE_SUBRESULTS_FOLDER)

if args.n_jobs > 1:
    pandarallel.initialize(nb_workers= args.n_jobs)

algorithms = {
    "PIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                   MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                   PIMVC()), "params": {}},
    "NEMO": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  NEMO()), "params": {}},
    "DFMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  DFMF().set_output(transform="pandas"),
                                  StandardScaler().set_output(transform="pandas"), KMeans()), "params": {}},
    "MOFA": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                  MOFA().set_output(transform="pandas"),
                                  StandardScaler().set_output(transform="pandas"), KMeans()), "params": {}},
    "EEIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    EEIMVC()), "params": {}},
    "SIMCADC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                     MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                     SIMCADC()), "params": {}},
    "LFIMVC": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    LFIMVC()), "params": {}},
    "IMSR": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(NormalizerNaN().set_output(transform="pandas")),
                                  IMSR()), "params": {}},
    "MKKMIK": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                    MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    MKKMIK()), "params": {}},
}
incomplete_algorithms = True
CommonOperations.run_script(dataset_table_path=DATASET_TABLE_PATH, algorithms=algorithms, probs=probs, views=views,
                            amputation_mechanisms=amputation_mechanisms, imputation=imputation, n_clusters=n_clusters,
                            runs_per_alg=runs_per_alg, args=args, subresults_path=INCOMPLETE_SUBRESULTS_PATH,
                            logs_file=INCOMPLETE_LOGS_PATH, error_file=INCOMPLETE_ERRORS_PATH,
                            results_path=INCOMPLETE_RESULTS_PATH, time_results_path=TIME_RESULTS_PATH,
                            incomplete_algorithms=incomplete_algorithms, best_combination=best_combination,
                            run_amputation=run_amputation, select_datasets=select_datasets)
