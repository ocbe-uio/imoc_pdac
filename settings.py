import os

import numpy as np

RANDOM_STATE = 42
TIME_LIMIT = 300

RESULTS_FOLDER = 'results'
AUX_DATA_FOLDER = 'aux_data'
PROFILES_FOLDER = 'profiles'
DATA_FOLDER = 'data'
PROFILES_PATH = os.path.join(AUX_DATA_FOLDER, PROFILES_FOLDER)

COMPLETE_SUBRESULTS_FOLDER = 'complete_subresults'
INCOMPLETE_SUBRESULTS_FOLDER = 'incomplete_subresults'
COMPLETE_SUBRESULTS_PATH = os.path.join(RESULTS_FOLDER, COMPLETE_SUBRESULTS_FOLDER)
INCOMPLETE_SUBRESULTS_PATH = os.path.join(RESULTS_FOLDER, INCOMPLETE_SUBRESULTS_FOLDER)

COMPLETE_RESULTS_FILE = 'complete_algorithms_evaluation.csv'
INCOMPLETE_RESULTS_FILE = 'incomplete_algorithms_evaluation.csv'
TIME_RESULTS_FILE = 'time_evaluation.csv'
DATASET_TABLE = "dataset_table.csv"
COMPLETE_RESULTS_PATH = os.path.join(RESULTS_FOLDER, COMPLETE_RESULTS_FILE)
INCOMPLETE_RESULTS_PATH = os.path.join(RESULTS_FOLDER, INCOMPLETE_RESULTS_FILE)
TIME_RESULTS_PATH = os.path.join(RESULTS_FOLDER, TIME_RESULTS_FILE)
DATASET_TABLE_PATH = os.path.join(AUX_DATA_FOLDER, DATASET_TABLE)

COMPLETE_LOGS_FILE = 'complete_logs.txt'
INCOMPLETE_LOGS_FILE = 'incomplete_logs.txt'
TIME_LOGS_FILE = 'time_logs.txt'
COMPLETE_LOGS_PATH = os.path.join(RESULTS_FOLDER, COMPLETE_LOGS_FILE)
INCOMPLETE_LOGS_PATH = os.path.join(RESULTS_FOLDER, INCOMPLETE_LOGS_FILE)
TIME_LOGS_PATH = os.path.join(RESULTS_FOLDER, TIME_LOGS_FILE)

COMPLETE_ERRORS_FILE = 'complete_errors.txt'
INCOMPLETE_ERRORS_FILE = 'incomplete_errors.txt'
TIME_ERRORS_FILE = 'time_errors.txt'
COMPLETE_ERRORS_PATH = os.path.join(RESULTS_FOLDER, COMPLETE_ERRORS_FILE)
INCOMPLETE_ERRORS_PATH = os.path.join(RESULTS_FOLDER, INCOMPLETE_ERRORS_FILE)
TIME_ERRORS_PATH = os.path.join(RESULTS_FOLDER, TIME_ERRORS_FILE)

views = ["CNA", "Methylation", "Mutations", "RNAseq", "RPPA", "miRNA"]  # All possible views in the dataset
best_combination = ["CNA", "Mutations"]       # Set to False if there is no set combination of data types to use (will run all possible combinations of data types)
run_amputation = False       # Set to False if there will be no amputation happening
select_datasets = ['PDAC_complete_samples']       # set to False if using all datasets

n_clusters = [2]
amputation_mechanisms = ["EDM", 'MCAR', 'MNAR', "PM"]
probs = [20, 40, 60, 80]
imputation = [False]
runs_per_alg = np.arange(10)
