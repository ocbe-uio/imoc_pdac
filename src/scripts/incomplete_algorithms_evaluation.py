import os.path
from pandarallel import pandarallel
from sklearn.cluster import KMeans
from sklearn.feature_selection import VarianceThreshold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from imvc.decomposition import DFMF, MOFA #, jNMF
from imvc.preprocessing import MultiViewTransformer, ConcatenateViews
from imvc.cluster import NEMO

from settings import INCOMPLETE_RESULTS_PATH, INCOMPLETE_SUBRESULTS_PATH, INCOMPLETE_LOGS_PATH, INCOMPLETE_ERRORS_PATH, \
    TIME_RESULTS_PATH, DATASET_TABLE_PATH, amputation_mechanisms, probs, \
    imputation, runs_per_alg, INCOMPLETE_RESULTS_FILE, INCOMPLETE_LOGS_FILE, INCOMPLETE_ERRORS_FILE, \
    INCOMPLETE_SUBRESULTS_FOLDER, omic_views, n_clusters, probs_zero, amputation_mechanisms_zero
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
    "NEMO": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
                                  MultiViewTransformer(StandardScaler().set_output(transform="pandas")),
                                    NEMO()), "params": {}},
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
    # "jNMF": {"alg": make_pipeline(MultiViewTransformer(VarianceThreshold().set_output(transform="pandas")),
    #                               MultiViewTransformer(MinMaxScaler().set_output(transform="pandas")),
    #                               jNMF().set_output(transform="pandas"),
    #                               StandardScaler().set_output(transform="pandas"), KMeans()),
    #          "params": {}},
}
incomplete_algorithms = True
CommonOperations.run_script(dataset_table_path=DATASET_TABLE_PATH, algorithms=algorithms, probs=probs_zero, omic_views=omic_views,
                            amputation_mechanisms=amputation_mechanisms_zero, imputation=imputation, n_clusters=n_clusters,
                            runs_per_alg=runs_per_alg, args=args, subresults_path=INCOMPLETE_SUBRESULTS_PATH,
                            logs_file=INCOMPLETE_LOGS_PATH, error_file=INCOMPLETE_ERRORS_PATH,
                            results_path=INCOMPLETE_RESULTS_PATH, time_results_path=TIME_RESULTS_PATH,
                            incomplete_algorithms=incomplete_algorithms)
