import numpy as np
import pandas as pd
from tqdm import tqdm
import itertools
from sklearn.metrics.cluster import adjusted_mutual_info_score, adjusted_rand_score
from scipy.stats import kruskal, chi2_contingency
import seaborn as sns
from ptitprince import PtitPrince as pt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from statsmodels.stats.multitest import multipletests


# Function to remove clusters containing < 10% patients. Output is also list with top x outlier patients (patients most typically found in small clusters)
def remove_small_clusters(df: pd.DataFrame, n_outliers: int, verbose: bool):
    valid_results = df[df['relative_cluster_sizes'].apply(lambda d: all(value >= 0.1 for value in d.values()))]
    outlier_results = df[df['relative_cluster_sizes'].apply(lambda f: any(value < 0.1 for value in f.values()))]
    outlier_patients = {}
    for index, row in outlier_results.iterrows():
        cluster_sizes = row['relative_cluster_sizes']
        small_clusters = [cluster for cluster, size in cluster_sizes.items() if size < 0.1]
        patients = row['y_pred_idx']
        clusters = row['y_pred']
        for patient, cluster in zip(patients, clusters):
            if cluster in small_clusters:
                if patient not in outlier_patients:
                    outlier_patients[patient] = 1
                else:
                    outlier_patients[patient] += 1
    for key in outlier_patients:
        outlier_patients[key] /= len(outlier_results)
    top_outliers = sorted(outlier_patients.items(), key=lambda x: x[1], reverse=True)[:n_outliers]
    df_top_outliers = pd.DataFrame(data=top_outliers, columns=['Patient ID', 'Count'])
    top_outlier_patients = [patient for patient, count in top_outliers]
    if verbose == True:
        print(f"Top {n_outliers} outlier patients: {top_outlier_patients}")
    return valid_results, outlier_results, df_top_outliers



# Function to perform clinical enrichment of labels to each row
def clinical_enrichment(results, clinical_data):
    # Sort labels
    results["sorted_y_pred_idx"] = results["y_pred_idx"].apply(sorted)
    results["sorted_y_pred"] = results.apply(
        lambda row: [row["y_pred"][row["y_pred_idx"].index(patient_id)] for patient_id in row["sorted_y_pred_idx"]],
        axis=1)
    # Add clinical data
    clinical_data = clinical_data[['Patient ID', 'Mutation Count', 'Fraction Genome Altered', 'Diagnosis Age', 'Sex', 'Race Category', 
                                   'Adjuvant Postoperative Targeted Therapy Administered Indicator', 'Alcohol History Documented', 'Tumor resected max dimension',
                                   'American Joint Committee on Cancer Metastasis Stage Code', 'American Joint Committee on Cancer Tumor Stage Code',
                                   'Chronic Pancreatitis Personal Medical History Indicator', 'Did patient start adjuvant postoperative radiotherapy?', 
                                   'Disease Free Status', 'Family History of Cancer', 'Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code', 
                                   'Neoplasm Disease Stage American Joint Committee on Cancer Code', 'Neoplasm Histologic Grade', 'TMB (nonsynonymous)',
                                   'New Neoplasm Event Post Initial Therapy Indicator', 'Overall Survival (Months)', 'Overall Survival Status', 'Disease Free (Months)', 
                                   'Participant Personal Medical History Diabetes Mellitus Ind-3', 'Patient Primary Tumor Site', 'Prior Cancer Diagnosis Occurence', 
                                   'Surgical Margin Resection Status', 'Patient Smoking History Category', 'Person Neoplasm Status', 'Primary Therapy Outcome Success Type']]
    clinical_data.set_index('Patient ID', inplace=True)
    # Convert necessary data 
    clinical_data['Overall Survival Status'] = clinical_data['Overall Survival Status'].str.split(':').str[0].astype(int)
    clinical_data['Patient Smoking History Category'] = (clinical_data['Patient Smoking History Category']
                                                         .where(clinical_data['Patient Smoking History Category'].isna(), 
                                                                clinical_data['Patient Smoking History Category'].astype(float).astype(str)))
    clinical_data_columns = [col for col in clinical_data.columns if col != 'Patient ID']
    def get_filtered_clinical_data(patient_ids, clinical_data, column_name):
        filtered_data = clinical_data.loc[patient_ids]
        return filtered_data[column_name].values.tolist()
    for column in clinical_data_columns:
        results[column] = results.apply(lambda row: get_filtered_clinical_data(row['sorted_y_pred_idx'], clinical_data, column), axis=1)
    # Logrank test
    def calculate_logrank_pvalue(row):
        df = pd.DataFrame({
            'cluster': row['sorted_y_pred'],
            'vital_status': row['Overall Survival Status'],
            'days_to_death': row['Overall Survival (Months)']
        })
        kmf = KaplanMeierFitter()
        test_results = multivariate_logrank_test(df['days_to_death'], df['cluster'], df['vital_status'])
        return test_results.p_value
    results['pvalue_logrank'] = results.apply(calculate_logrank_pvalue, axis=1)
    # Function for p-values (Kruskal-Wallis, chi2)
    def pvalue_tests(row, clinical_data_columns):
        pvalues = []
        for variable in clinical_data_columns:
            df = pd.DataFrame({
                'cluster': row['y_pred'],
                variable: row[variable]
            })
            if pd.api.types.is_numeric_dtype(df[variable]):
                # Kruskal-Wallis test for numerical variables
                test_numerical = [df[df['cluster'] == cluster][variable].dropna().to_numpy() for cluster in df['cluster'].unique()]
                stat, p_value_kruskal = kruskal(*test_numerical)
                pvalues.append(p_value_kruskal)
            else:
                # Chi-square contingency test for categorical variables
                test_discrete = pd.crosstab(df['cluster'], df[variable])
                chi2, p_value_chi2, dof, freq = chi2_contingency(test_discrete)
                pvalues.append(p_value_chi2)
        reject, pvals_corr, asidack, abonf = multipletests(pvals=pvalues, alpha=0.05, method='fdr_bh')
        for idx, variable in enumerate(clinical_data_columns):
            row[f"pvalue_{variable}"] = pvals_corr[idx]
        return row
    clinical_data_columns = [col for col in clinical_data_columns if col not in ['Overall Survival Status', 'Overall Survival (Months)']]
    results = results.apply(lambda row: pvalue_tests(row, clinical_data_columns), axis=1)
    columns_to_check = [f'pvalue_{variable}' for variable in clinical_data_columns]
    results['n_enriched_clinical'] = (results[columns_to_check] < 0.05).sum(axis=1)
    return results



# Function to calculate stability metrics
def calculate_stability_metrics(results: pd.DataFrame, random_state=None, progress_bar=True):
    base_columns = ['dataset', 'view_combination', 'algorithm', 'n_clusters', 'missing_percentage', 'amputation_mechanism', 'imputation', 'run_n', "sorted_y_pred", 
                    "sorted_y_pred_idx", 'silhouette', 'vrc', 'db', 'dbcv', 'dunn', "dhi", "ssei", 'rsi', 'bhi']
    normalised_columns = [col for col in results.columns if 'normalised' in col]
    robustness_column = [col for col in results.columns if 'robustness' in col]
    clin_columns = [col for col in results.columns if col in ['pvalue_logrank', 'n_enriched_clinical']]
    all_columns = base_columns + normalised_columns + clin_columns + robustness_column
    alg_stability = results[all_columns]
    if alg_stability["imputation"].nunique() != 1:
        alg_stability = alg_stability.loc[
            (alg_stability["missing_percentage"] == 0) | (alg_stability["imputation"])
            ]
    alg_uns_metrics = alg_stability.drop(columns=["sorted_y_pred", "sorted_y_pred_idx",'imputation', 'run_n'])
    # Group by taking mean of metrics
    alg_uns_metrics = alg_uns_metrics.groupby(
        ["dataset", "algorithm", "missing_percentage", "amputation_mechanism", "view_combination", "n_clusters"], as_index=False).mean()
    iterator = alg_stability["dataset"].unique()
    if progress_bar:
        iterator = tqdm(iterator)
    for dataset in iterator:
        preds_dataset = alg_stability.loc[
            (alg_stability["dataset"] == dataset), ["missing_percentage", "algorithm", 'amputation_mechanism', 'n_clusters', 
                                                    'view_combination', "run_n", "sorted_y_pred", "sorted_y_pred_idx"]]
        for alg in preds_dataset["algorithm"].unique():
            pred_alg = preds_dataset[preds_dataset["algorithm"] == alg]
            for missing_percentage in pred_alg["missing_percentage"].unique():
                pred_missing_alg = pred_alg[pred_alg["missing_percentage"] == missing_percentage]
                for amputation_mechanism in pred_missing_alg["amputation_mechanism"].unique():
                    pred_missing_ampt_alg = pred_missing_alg[
                        pred_missing_alg["amputation_mechanism"] == amputation_mechanism]
                    for view in pred_missing_ampt_alg["view_combination"].unique():
                        pred_missing_ampt_alg_view = pred_missing_ampt_alg[
                            pred_missing_ampt_alg["view_combination"] == view]
                        for cluster in pred_missing_ampt_alg_view["n_clusters"].unique():
                            pred_missing_ampt_alg_view_clus = pred_missing_ampt_alg_view[
                                pred_missing_ampt_alg_view['n_clusters'] == cluster]
                            amis, aris = [], []
                            for run_1, run_2 in set(itertools.combinations(pred_missing_ampt_alg_view_clus["run_n"].unique(), 2)):
                                pred1_alg = pred_missing_ampt_alg_view_clus.loc[
                                    (pred_missing_ampt_alg_view_clus["run_n"] == run_1), "sorted_y_pred"].to_list()[0]
                                pred2_alg = pred_missing_ampt_alg_view_clus.loc[
                                    (pred_missing_ampt_alg_view_clus["run_n"] == run_2), "sorted_y_pred"].to_list()[0]
                                pred1_idx = pred_missing_ampt_alg_view_clus.loc[(
                                    pred_missing_ampt_alg_view_clus["run_n"] == run_1), "sorted_y_pred_idx"].to_list()[0]
                                pred2_idx = pred_missing_ampt_alg_view_clus.loc[(
                                    pred_missing_ampt_alg_view_clus["run_n"] == run_2), "sorted_y_pred_idx"].to_list()[0]
                                # Only select samples in common for stability metrics
                                common_samples = list(set(pred1_idx) & set(pred2_idx))
                                pred1_common = [pred1_alg[pred1_idx.index(i)] for i in common_samples]
                                pred2_common = [pred2_alg[pred2_idx.index(i)] for i in common_samples]
                                amis.append(adjusted_mutual_info_score(pred1_common, pred2_common)), aris.append(
                                    adjusted_rand_score(pred1_common, pred2_common))
                            alg_uns_metrics.loc[(alg_uns_metrics["dataset"] == dataset) &
                                                (alg_uns_metrics["missing_percentage"] == missing_percentage) &
                                                (alg_uns_metrics["amputation_mechanism"] == amputation_mechanism) &
                                                (alg_uns_metrics["algorithm"] == alg) & 
                                                (alg_uns_metrics["view_combination"] == view) & 
                                                (alg_uns_metrics["n_clusters"] == cluster),
                            ["AMI", "ARI"]] = [np.mean(amis), np.mean(aris)]
    return alg_uns_metrics



# Function to normalise metrics with respect to a variable
def add_normalised_metric(df, variable_to_normalise, metric, greater_is_better=True):
    possible_variables = ["dataset", "algorithm", "missing_percentage", "amputation_mechanism", "view_combination", "n_clusters"]
    valid_variables = [var for var in possible_variables if var != variable_to_normalise]
    def recursive_loop(subset, remaining_vars, current_filters):
        if not remaining_vars:
            scores = subset[metric].values
            relative_score = scores / scores.max()
            if not greater_is_better:
                relative_score = 1 - relative_score
            condition = True
            for key, value in current_filters.items():
                condition &= (df[key] == value)
            df.loc[condition, f'normalised_{metric}'] = relative_score
            return
        current_var = remaining_vars[0]
        for unique_value in subset[current_var].unique():
            filtered_subset = subset[subset[current_var] == unique_value]
            recursive_loop(filtered_subset, remaining_vars[1:], {**current_filters, current_var: unique_value})
    df[f'normalised_{metric}'] = float('nan')
    recursive_loop(df, valid_variables, {})
    return df



# Function to create consensus matrix
def consensus_matrix(df: pd.DataFrame):
    patients = list(set(sum(df['y_pred_idx'], [])))
    connectivity_matrix = pd.DataFrame(0, index=patients, columns=patients)
    indicator_matrix = pd.DataFrame(0, index=patients, columns=patients)
    for index, row in df.iterrows():
        y_pred_idx = row['y_pred_idx']
        y_pred = row['y_pred']
        for i in range(len(y_pred_idx)):
            patient_i = y_pred_idx[i]
            cluster_i = y_pred[i]
            for j in range(i + 1, len(y_pred_idx)):
                patient_j = y_pred_idx[j]
                cluster_j = y_pred[j]
                indicator_matrix.loc[patient_i, patient_j] += 1
                indicator_matrix.loc[patient_j, patient_i] += 1
                if cluster_i == cluster_j:
                    connectivity_matrix.loc[patient_i, patient_j] += 1
                    connectivity_matrix.loc[patient_j, patient_i] += 1
    consensus_matrix = connectivity_matrix.div(indicator_matrix).fillna(0)
    return consensus_matrix



# Function to plot boxplots for a specific variable (first benchmark)
def raincloud_plots_variables(df, column_name, metric_name, ax, ylabel):
    mean_values = df.groupby(column_name)[metric_name].mean().sort_values(ascending=False)
    sorted_categories = mean_values.index.tolist()
    metric_subsets = [df[df[column_name] == cat][metric_name].values for cat in sorted_categories]
    pt.RainCloud(x=column_name, y=metric_name, data=df, bw=0.2, palette=[colorblind_palette[0]],
                 width_viol=0.4, ax=ax, orient="v", move=0.2, order=sorted_categories, alpha=0.8)
    means = df.groupby(column_name)[metric_name].mean().loc[sorted_categories]
    sns.scatterplot(x=range(len(sorted_categories)), y=means.values, ax=ax, color=colorblind_palette[2], s=100, marker='^', zorder=10)
    ax.set_xticks(range(len(sorted_categories)))
    ax.set_xticklabels(sorted_categories)
    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.05, 1.05)
    ax.set_axisbelow(True)
    pvalue = kruskal(*metric_subsets).pvalue
    if pvalue >= 0.001:
        pvalue_text = f"Kruskal-Wallis, p = {pvalue:.3f}"
    else:
        pvalue_text = f"Kruskal-Wallis, p = {pvalue:.2e}"
    return pvalue, pvalue_text


# Function to measure robustness
def measure_robustness(results):
    # Sort labels
    results["sorted_y_pred_idx"] = results["y_pred_idx"].apply(sorted)
    results["sorted_y_pred"] = results.apply(
        lambda row: [row["y_pred"][row["y_pred_idx"].index(patient_id)] for patient_id in row["sorted_y_pred_idx"]],
        axis=1)
    # Divide into complete and missing data dataframes
    complete_df = results[results['missing_percentage'] == 0]
    missing_df = results[results['missing_percentage'] != 0]
    results['robustness'] = np.nan
    # Compare complete vs missing cluster assignments
    for _, base_row in complete_df.iterrows():
        base_idx = base_row["sorted_y_pred_idx"]
        base_pred = base_row["sorted_y_pred"]
        matching_rows = missing_df[
            (missing_df["dataset"] == base_row["dataset"]) & 
            (missing_df["algorithm"] == base_row["algorithm"]) &
            (missing_df["view_combination"] == base_row["view_combination"]) & 
            (missing_df["n_clusters"] == base_row["n_clusters"]) & 
            (missing_df["run_n"] == base_row["run_n"])]
        for i, miss_row in matching_rows.iterrows():
            miss_idx = miss_row["sorted_y_pred_idx"]
            miss_pred = miss_row["sorted_y_pred"]
            if base_idx != miss_idx:
                continue
            matches = sum(1 for b, m in zip(base_pred, miss_pred) if b == m)
            robustness_score = matches / len(base_pred)
            results.at[i, "robustness"] = robustness_score
        results.at[base_row.name, "robustness"] = 1
    return results
