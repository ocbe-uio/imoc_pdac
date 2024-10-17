import argparse
import itertools
import json
import os.path
import shutil
from sklearn.utils import shuffle
from imvc.ampute import Amputer
from imvc.impute import get_observed_view_indicator

from settings import (PROFILES_PATH, DATASET_TABLE_PATH, RANDOM_STATE, probs, amputation_mechanisms, runs_per_alg,
                      best_combination, views, run_amputation, select_datasets)
from src.commons import CommonOperations

parser = argparse.ArgumentParser()
parser.add_argument('-continue_indxs', default=False, action='store_true')
parser.add_argument('-save_results', default=False, action='store_true')
args = parser.parse_args()

if not args.continue_indxs:
    shutil.rmtree(PROFILES_PATH, ignore_errors=True)
    os.mkdir(PROFILES_PATH)

datasets, two_view_datasets = CommonOperations.get_list_of_datasets(DATASET_TABLE_PATH, select_datasets)

for dataset_name in datasets:

    if best_combination:
        if isinstance(best_combination, list) and len(best_combination) > 1:
            binary_combination = ['1' if view in best_combination else '0' for view in views]
            binary_combination = [''.join(binary_combination)]
        else:
            raise ValueError("best_combination must be a list with at least two views.")
    elif best_combination == False:
        binary_combination = [row for row in itertools.product([0, 1], repeat=len(views)) if sum(row) >= 2]
        binary_combination = [''.join(map(str, combination)) for combination in binary_combination]
    else:
        raise TypeError

    if run_amputation == False:
        probs = [0]
        amputation_mechanisms = ["edm"]

    for binary_combination, prob, amputation_mechanism, run_n in itertools.product(binary_combination, probs, amputation_mechanisms, runs_per_alg):
        Xs = CommonOperations.load_Xs(dataset_name=dataset_name)
        view_combinations = [bool(int(value)) for value in binary_combination]
        Xs = [view for i, view in enumerate(Xs) if view_combinations[i] == True]

        if prob == 0:
            if amputation_mechanism == "edm":
                amputation_mechanism = "No"
            else:
                continue

        path = f"{dataset_name}_{binary_combination}_{prob}_{amputation_mechanism}_{run_n}.json"
        path = os.path.join(PROFILES_PATH, path)
        if os.path.exists(path):
            continue
        try:
            random_state = RANDOM_STATE + run_n
            if (dataset_name in two_view_datasets) and (amputation_mechanism in ["mnar"]):
                continue
            *train_Xs, = shuffle(*Xs, random_state=random_state)
            p = prob/100

            if p != 0:
                amp = Amputer(p=round(p, 2), mechanism=amputation_mechanism, random_state=random_state)
                train_Xs = amp.fit_transform(train_Xs)

            observed_view_indicator = get_observed_view_indicator(train_Xs)
            assert (train_Xs[0].index == observed_view_indicator.index).all()

            dict_indxs = {
                "observed_view_indicator": observed_view_indicator.to_dict(),
                "valid": True,
            }

        except (AssertionError, ValueError) as exception:
            dict_indxs = {
                "valid" : False,
                "error": str(exception),
            }

        if args.save_results:
            with open(path, "w") as f:
                json.dump(dict_indxs, f)

print("Completed successfully!")
